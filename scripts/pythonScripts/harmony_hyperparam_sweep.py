#!/usr/bin/env python
"""
harmony_hyperparam_sweep.py

Grid-search Harmony hyperparameters on a UCE-embedded AnnData, tuned to remove
STUDY effects while keeping real species biology.

Why not a batch-mixing score: in this design every species comes from a single
study (and technology tracks Tan vs non-Tan), so "study mixing" cannot be told
apart from "species mixing" and would reward erasing species signal. Instead,
settings are judged on things a study effect distorts and over-correction
flattens:

  1. Tan check (gate)   - the Tan study profiles several species under one
                          protocol, so its cross-species, same-cell-type
                          distances (relative to between-type distances) should
                          survive correction:
                            tan_compression : mean relative distance after/before
                            tan_spearman    : rank correlation of per-pair distances
                          A setting must pass both thresholds to be recommended.
                          Blind spot: Harmony moves Tan as one batch, so this
                          cannot see non-Tan species being flattened.

  2. Batch-matched triplets (main target) - for an anchor species A and two
                          species B, C from the SAME study, the phylogenetically
                          closer of B/C should have the closer same-cell-type
                          centroid. B and C carry the same study offset relative
                          to A, so a pure batch embedding scores chance (0.5) and
                          so does an over-corrected one. A study effect large
                          enough to swamp biology also pushes this toward chance;
                          removing it restores the signal.

  3. Phylogeny score (reported, low weight) - per cell type, Spearman between
                          cross-species centroid distance and divergence time,
                          averaged over cell types. Uses every species pair, so
                          it is partly confounded with study wherever the study
                          split follows clades. `phylo_study_baseline` is the
                          score a pure study-only embedding would get on the
                          same pairs; read phylo_spearman relative to it.

  4. Cell-type structure (guard) - cLISI and silhouette on cell type.

overall = w_triplet * triplet_accuracy + w_phylo * (phylo_spearman + 1) / 2
          + w_celltype * mean(clisi_norm, (silhouette + 1) / 2)
The recommendation is the highest-overall Harmony setting that passes the Tan
check (triplet margin breaks ties).

Divergence times: a built-in tree with approximate TimeTree node ages (Mya)
covers human, macaque, mouse, rat, rabbit, dog, pig, goat. Verify these against
timetree.org before publishing, or pass --divergence-csv with columns
species1,species2,mya (species names as in --species-key). Only the ordering of
ages matters for triplets; the Spearman scores use ranks.

Outputs (in --outdir):
  sweep_results.csv          one row per setting (+ raw UCE baseline); written
                             incrementally during the sweep
  triplet_details.csv        every evaluated triplet for raw UCE and the
                             recommended setting (anchor, near, far, cell type,
                             distances, correct?)
  best_params.json           recommended setting and its scores vs raw
  sweep_tradeoff.png         triplet accuracy vs Tan compression
  sweep_overall_heatmap.png  overall score over theta x lambda per covariate set
  (with --apply-best) X_uce_harmony_best.npy + obs_names.txt for all cells

Example:
  python harmony_hyperparam_sweep.py uce_distances_combined.h5ad \
      --celltype-key cell_type_std --tan-study Tan --metric cosine \
      --covariate-sets study technology \
      --thetas 1 2 4 6 8 --lambdas 1 0.5 0.1 --subsample 60000 --apply-best
"""

import argparse
import itertools
import json
import logging
import time
from pathlib import Path

import anndata as ad
import harmonypy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

log = logging.getLogger("harmony_sweep")

# Approximate TimeTree node ages (Mya). Structure: (age, [children]); leaves are
# species names as they appear in obs. VERIFY against timetree.org.
DEFAULT_TREE = (
    94.0, [                                   # Boreoeutheria
        (90.0, [                              # Euarchontoglires
            (29.0, ["human", "macaque"]),     # Catarrhini
            (79.0, [                          # Glires
                (13.0, ["mouse", "rat"]),     # Murinae
                "rabbit",
            ]),
        ]),
        (78.0, [                              # Laurasiatheria (Carnivora vs Cetartiodactyla)
            "dog",
            (62.0, ["pig", "goat"]),          # Cetartiodactyla
        ]),
    ],
)


# -----------------------------------------------------------------------------
# Args
# -----------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("h5ad", help="Merged, UCE-embedded AnnData")
    p.add_argument("--outdir", default="harmony_sweep")
    p.add_argument("--emb-key", default="X_uce")
    p.add_argument("--study-key", default="study")
    p.add_argument("--species-key", default="species")
    p.add_argument("--celltype-key", default="cell_type_std")
    p.add_argument("--tan-study", default="Tan", help="Value of --study-key marking the Tan subset")
    p.add_argument("--divergence-csv", default=None,
                   help="Optional CSV (species1,species2,mya) replacing the built-in tree")

    # Grid
    p.add_argument("--covariate-sets", nargs="+", default=["study", "technology"],
                   help="Comma-separated obs columns per set, e.g. study technology study,technology")
    p.add_argument("--thetas", nargs="+", type=float, default=[1, 2, 4, 6, 8])
    p.add_argument("--lambdas", nargs="+", type=float, default=[1.0, 0.5, 0.1])
    p.add_argument("--sigmas", nargs="+", type=float, default=[0.1])
    p.add_argument("--nclust", type=int, default=None, help="Harmony clusters (default: harmonypy's)")
    p.add_argument("--max-iter", type=int, default=1000,
                   help="Cap on Harmony iterations. Harmony stops as soon as it converges, so this "
                        "is only a safety limit; runs that hit it are flagged and not recommended.")
    p.add_argument("--epsilon-harmony", type=float, default=1e-4,
                   help="Convergence threshold: stop when the relative change in Harmony's objective "
                        "falls below this. Pinned here so results don't depend on the harmonypy "
                        "version's default (1e-2 in harmonypy 2.x).")

    # Data size / space
    p.add_argument("--subsample", type=int, default=60000, help="Max cells for the sweep (0 = all)")
    p.add_argument("--cap-per-group", type=int, default=2000,
                   help="Max cells per study x species x cell-type group before subsampling")
    p.add_argument("--n-pcs", type=int, default=0,
                   help="Run in top-N PCs of UCE instead of native 1280-d space (0 = native)")

    # Metrics
    p.add_argument("--k", type=int, default=30, help="Neighbours for cLISI")
    p.add_argument("--min-cells", type=int, default=20, help="Min cells for a species x cell-type centroid")
    p.add_argument("--metric", default="cosine", choices=["euclidean", "cosine"],
                   help="Centroid distance metric (your core script uses cosine)")
    p.add_argument("--min-compression", type=float, default=0.8)
    p.add_argument("--min-spearman", type=float, default=0.8)
    p.add_argument("--w-triplet", type=float, default=0.6)
    p.add_argument("--w-phylo", type=float, default=0.2)
    p.add_argument("--w-celltype", type=float, default=0.2)

    p.add_argument("--apply-best", action="store_true",
                   help="Rerun the recommended setting on the full dataset and save the embedding")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
def load(args):
    log.info("Loading %s", args.h5ad)
    a = ad.read_h5ad(args.h5ad, backed="r")  # obs/obsm into memory, X stays on disk
    obs = a.obs.copy()
    X = np.asarray(a.obsm[args.emb_key], dtype=np.float64)
    names = a.obs_names.to_numpy()
    a.file.close()

    needed = {args.study_key, args.species_key, args.celltype_key}
    for s in args.covariate_sets:
        needed.update(s.split(","))
    missing = needed - set(obs.columns)
    if missing:
        raise KeyError(f"Missing obs columns: {sorted(missing)}. Available: {list(obs.columns)}")
    for c in needed:
        obs[c] = obs[c].astype(str)
    if args.tan_study not in set(obs[args.study_key]):
        raise ValueError(f"--tan-study '{args.tan_study}' not found in obs['{args.study_key}']")
    log.info("%d cells x %d dims", *X.shape)
    return X, obs, names


def stratified_subsample(obs, keys, cap, total, seed):
    rng = np.random.default_rng(seed)
    idx = []
    for g in obs.groupby(keys, observed=True).indices.values():
        idx.append(rng.choice(g, cap, replace=False) if len(g) > cap else g)
    idx = np.concatenate(idx)
    if total and len(idx) > total:
        idx = rng.choice(idx, total, replace=False)
    return np.sort(idx)


# -----------------------------------------------------------------------------
# Phylogeny
# -----------------------------------------------------------------------------
def _leaves(node):
    return [node] if isinstance(node, str) else [l for c in node[1] for l in _leaves(c)]


def tree_divergence(node, out=None):
    """Pairwise MRCA ages from a nested (age, [children]) tree."""
    out = {} if out is None else out
    if isinstance(node, str):
        return out
    age, children = node
    for a_set, b_set in itertools.combinations([_leaves(c) for c in children], 2):
        for a in a_set:
            for b in b_set:
                out[frozenset((a, b))] = float(age)
    for c in children:
        tree_divergence(c, out)
    return out


def load_divergence(args, species):
    if args.divergence_csv:
        t = pd.read_csv(args.divergence_csv)
        div = {frozenset((str(r.species1), str(r.species2))): float(r.mya) for r in t.itertuples()}
        src = args.divergence_csv
    else:
        div = tree_divergence(DEFAULT_TREE)
        src = "built-in approximate TimeTree ages"
    missing = [f"{a}-{b}" for a, b in itertools.combinations(sorted(species), 2)
               if frozenset((a, b)) not in div]
    if missing:
        raise ValueError(f"No divergence time for species pairs {missing} ({src}). "
                         "Pass --divergence-csv with columns species1,species2,mya.")
    log.info("Divergence times from %s", src)
    return div


def build_triplets(species_study, div):
    """(anchor, near, far) where near/far share the same study set and differ in
    divergence time from the anchor."""
    sp = sorted(species_study)
    out = []
    for a in sp:
        for b, c in itertools.combinations([x for x in sp if x != a], 2):
            if species_study[b] != species_study[c]:
                continue
            tb, tc = div[frozenset((a, b))], div[frozenset((a, c))]
            if tb == tc:
                continue
            near, far = (b, c) if tb < tc else (c, b)
            out.append((a, near, far))
    return out


# -----------------------------------------------------------------------------
# Harmony
# -----------------------------------------------------------------------------
def harmony_convergence(ho, max_iter, epsilon):
    """Iterations run and whether Harmony converged (vs hitting max_iter).
    objective_harmony holds the initial objective plus one value per iteration."""
    obj = [float(v) for v in ho.objective_harmony]
    rounds = list(getattr(ho, "kmeans_rounds", []))
    n_iter = len(rounds) if rounds else max(len(obj) - 1, 0)
    rel = (obj[-2] - obj[-1]) / abs(obj[-2]) if len(obj) >= 2 and obj[-2] != 0 else np.nan
    converged = n_iter < max_iter or (not np.isnan(rel) and rel < epsilon)
    return dict(harmony_n_iter=n_iter, harmony_converged=bool(converged), harmony_last_rel_change=rel)


def run_harmony(X, meta, covars, theta, lamb, sigma, nclust, max_iter, seed, epsilon=1e-4):
    """Returns (corrected embedding, convergence info)."""
    ho = harmonypy.run_harmony(
        X, meta, covars,
        theta=[theta] * len(covars),
        lamb=[lamb] * len(covars),
        sigma=sigma,
        nclust=nclust,
        max_iter_harmony=max_iter,
        epsilon_harmony=epsilon,
        verbose=False,
        random_state=seed,
    )
    info = harmony_convergence(ho, max_iter, epsilon)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] == X.shape[0]:
        return Z, info
    if Z.shape[1] == X.shape[0]:
        return Z.T, info
    raise ValueError(f"harmonypy output shape {Z.shape} does not match {X.shape[0]} cells")


# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------
def lisi(nn_idx, labels):
    """kNN inverse-Simpson diversity (uniform weights). Returns per-cell LISI."""
    codes, _ = pd.factorize(labels)
    L = codes[nn_idx]
    n, k = L.shape
    counts = np.zeros((n, codes.max() + 1))
    np.add.at(counts, (np.repeat(np.arange(n), k), L.ravel()), 1)
    p = counts / k
    return 1.0 / (p ** 2).sum(1)


def centroids(Z, species, ct, min_cells):
    groups = pd.DataFrame({"sp": species, "ct": ct}).groupby(["sp", "ct"], observed=True).indices
    keys = sorted(k for k, v in groups.items() if len(v) >= min_cells)
    return {k: Z[groups[k]].mean(0) for k in keys}


def dist(u, v, metric):
    return float(cdist(u[None], v[None], metric=metric)[0, 0])


def tan_relative_distances(cd, metric):
    """Same-type cross-species distances / mean different-type distance."""
    keys = sorted(cd)
    C = np.vstack([cd[k] for k in keys])
    D = cdist(C, C, metric=metric)
    sp = np.array([k[0] for k in keys])
    ct = np.array([k[1] for k in keys])
    i, j = np.triu_indices(len(keys), 1)
    d = D[i, j]
    same = (ct[i] == ct[j]) & (sp[i] != sp[j])
    diff = ct[i] != ct[j]
    if same.sum() == 0 or diff.sum() == 0:
        return {}
    mean_diff = d[diff].mean()
    return {(ct[a], *sorted((sp[a], sp[b]))): v / mean_diff for a, b, v in zip(i[same], j[same], d[same])}


def triplet_eval(cd, triplets, species_study, metric):
    rows = []
    cts = sorted({k[1] for k in cd})
    for ct in cts:
        for a, n, f in triplets:
            if (a, ct) in cd and (n, ct) in cd and (f, ct) in cd:
                dn = dist(cd[(a, ct)], cd[(n, ct)], metric)
                df = dist(cd[(a, ct)], cd[(f, ct)], metric)
                rows.append(dict(cell_type=ct, anchor=a, near=n, far=f, d_near=dn, d_far=df,
                                 correct=dn < df,
                                 cross_study=species_study[a] != species_study[n]))
    return pd.DataFrame(rows)


def phylo_eval(cd, div, species_study, metric, min_pairs=4):
    """Per-cell-type Spearman(distance, divergence), weighted by #pairs.
    Also the Spearman a pure study-only embedding would get on the same pairs."""
    rhos, base, w = [], [], []
    for ct in sorted({k[1] for k in cd}):
        sps = sorted(s for s, c in cd if c == ct)
        pairs = list(itertools.combinations(sps, 2))
        if len(pairs) < min_pairs:
            continue
        t = [div[frozenset(p)] for p in pairs]
        d = [dist(cd[(a, ct)], cd[(b, ct)], metric) for a, b in pairs]
        rho = spearmanr(t, d).statistic
        if np.isnan(rho):
            continue
        s = [float(species_study[a] != species_study[b]) for a, b in pairs]
        b0 = spearmanr(t, s).statistic if len(set(s)) > 1 else 0.0
        rhos.append(rho)
        base.append(b0)
        w.append(len(pairs))
    if not rhos:
        return np.nan, np.nan, 0
    return float(np.average(rhos, weights=w)), float(np.average(base, weights=w)), len(rhos)


def evaluate(Z, obs, args, ctx, tan_rel_raw):
    ct = obs[args.celltype_key].to_numpy()
    species = obs[args.species_key].to_numpy()
    n_ct = len(np.unique(ct))

    nn = NearestNeighbors(n_neighbors=args.k + 1, n_jobs=-1).fit(Z)
    nn_idx = nn.kneighbors(Z, return_distance=False)[:, 1:]
    clisi = np.median(lisi(nn_idx, ct))
    slisi = np.median(lisi(nn_idx, species))
    sil = silhouette_score(Z, ct, sample_size=min(len(Z), 10000), random_state=args.seed, metric=args.metric)

    cd = centroids(Z, species, ct, args.min_cells)
    trip = triplet_eval(cd, ctx["triplets"], ctx["species_study"], args.metric)
    phylo, phylo_base, n_phylo_ct = phylo_eval(cd, ctx["div"], ctx["species_study"], args.metric)

    tm = ctx["tan_mask"]
    tan_rel = tan_relative_distances(centroids(Z[tm], species[tm], ct[tm], args.min_cells), args.metric)

    def acc(df):
        return float(df.correct.mean()) if len(df) else np.nan

    margin = ((trip.d_far - trip.d_near) / (trip.d_far + trip.d_near)).mean() if len(trip) else np.nan
    out = dict(
        triplet_accuracy=acc(trip),
        triplet_margin=float(margin),
        triplet_acc_cross_study=acc(trip[trip.cross_study]) if len(trip) else np.nan,
        triplet_acc_within_study=acc(trip[~trip.cross_study]) if len(trip) else np.nan,
        n_triplets=len(trip),
        phylo_spearman=phylo,
        phylo_study_baseline=phylo_base,
        n_phylo_celltypes=n_phylo_ct,
        clisi_norm=(n_ct - clisi) / max(n_ct - 1, 1),
        celltype_silhouette=sil,
        ilisi_species=slisi,  # diagnostic: species mixing (rises with over-correction)
    )
    if tan_rel_raw is None:
        out.update(tan_compression=1.0, tan_spearman=1.0, n_tan_pairs=len(tan_rel))
    else:
        shared = sorted(set(tan_rel) & set(tan_rel_raw))
        raw_v = np.array([tan_rel_raw[s] for s in shared])
        new_v = np.array([tan_rel[s] for s in shared])
        out.update(
            tan_compression=new_v.mean() / raw_v.mean() if len(shared) else np.nan,
            tan_spearman=spearmanr(raw_v, new_v).statistic if len(shared) > 2 else np.nan,
            n_tan_pairs=len(shared),
        )
    return out, tan_rel, trip


def add_scores(df, args):
    df["tan_pass"] = (df.tan_compression >= args.min_compression) & (df.tan_spearman >= args.min_spearman)
    df["celltype_score"] = (df.clisi_norm + (df.celltype_silhouette + 1) / 2) / 2
    w = args.w_triplet + args.w_phylo + args.w_celltype
    df["overall"] = (args.w_triplet * df.triplet_accuracy
                     + args.w_phylo * (df.phylo_spearman.fillna(0) + 1) / 2
                     + args.w_celltype * df.celltype_score) / w
    return df


# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
def plot_tradeoff(df, args, path):
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    h = df[df.method == "harmony"]
    thetas = sorted(h.theta.unique())
    cmap = plt.get_cmap("Blues")
    shade = {t: cmap(0.35 + 0.6 * i / max(len(thetas) - 1, 1)) for i, t in enumerate(thetas)}
    markers = dict(zip(sorted(h.covariates.unique()), "osD^v<>p"))
    for _, r in h.iterrows():
        ax.scatter(r.triplet_accuracy, r.tan_compression, marker=markers[r.covariates],
                   color=shade[r.theta], s=60, edgecolor="#333333", linewidth=0.6, zorder=3)
    raw = df[df.method == "raw"].iloc[0]
    ax.scatter(raw.triplet_accuracy, raw.tan_compression, marker="X", color="#777777", s=130,
               zorder=4, label="raw UCE")
    ax.axhline(args.min_compression, ls="--", c="#999999", lw=1)
    ax.axvline(0.5, ls=":", c="#999999", lw=1)
    ax.text(0.5, ax.get_ylim()[0], " chance", color="#777777", fontsize=8, va="bottom")
    best = df[df.recommended]
    if len(best):
        ax.scatter(best.triplet_accuracy, best.tan_compression, s=320, facecolors="none",
                   edgecolors="#c0392b", linewidth=2, zorder=5, label="recommended")
    for cov, m in markers.items():
        ax.scatter([], [], marker=m, color="#888888", edgecolor="#333333", label=f"covariates: {cov}")
    for t in thetas:
        ax.scatter([], [], marker="s", color=shade[t], label=f"theta = {t:g}")
    ax.set_xlabel("Batch-matched triplet accuracy  (0.5 = chance; higher = phylogeny recovered)")
    ax.set_ylabel("Tan relative species distance (after / before)")
    ax.set_title("Harmony sweep: study-effect removal vs Tan preservation\n"
                 "(dashed line = Tan compression threshold)", fontsize=10)
    ax.grid(color="#eeeeee", zorder=0)
    ax.legend(fontsize=7.5, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_heatmaps(df, args, path):
    h = df[df.method == "harmony"]
    groups = list(h.groupby(["covariates", "sigma"]))
    fig, axes = plt.subplots(1, len(groups), figsize=(4.5 * len(groups), 4), squeeze=False)
    vmin, vmax = h.overall.min(), h.overall.max()
    im = None
    for ax, ((cov, sig), sub) in zip(axes[0], groups):
        piv = sub.pivot_table(index="theta", columns="lamb", values="overall")
        ok = sub.assign(ok=sub.tan_pass.astype(float)).pivot_table(index="theta", columns="lamb", values="ok")
        im = ax.imshow(piv.values, cmap="Blues", vmin=vmin, vmax=vmax, aspect="auto", origin="lower")
        ax.set_xticks(range(piv.shape[1]), [f"{v:g}" for v in piv.columns])
        ax.set_yticks(range(piv.shape[0]), [f"{v:g}" for v in piv.index])
        for (r, c), v in np.ndenumerate(piv.values):
            failed = ok.values[r, c] < 1
            if failed:
                ax.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1, facecolor="#dddddd",
                                           edgecolor="#999999", hatch="///", lw=0))
            dark = (v - vmin) / max(vmax - vmin, 1e-9) > 0.6 and not failed
            ax.text(c, r, f"{v:.2f}" + ("\nfails Tan" if failed else ""), ha="center", va="center",
                    fontsize=8, color="white" if dark else "#222222")
        ax.set_xlabel("lambda")
        ax.set_ylabel("theta")
        ax.set_title(f"{cov} | sigma={sig:g}", fontsize=10)
    fig.colorbar(im, ax=axes.ravel().tolist(), label="overall score")
    fig.suptitle(f"Hatched = fails Tan check (compression < {args.min_compression:g} or "
                 f"Spearman < {args.min_spearman:g}); not eligible", fontsize=9, y=1.02)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    X_full, obs_full, names_full = load(args)
    cov_sets = [s.split(",") for s in args.covariate_sets]
    for cs in cov_sets:
        if args.species_key in cs:
            log.warning("Covariate set %s includes species: this removes the signal you are "
                        "measuring, so its scores are partly circular.", cs)

    # Design-level context (from the full data, not the subsample)
    species_study = {sp: frozenset(g[args.study_key].unique())
                     for sp, g in obs_full.groupby(args.species_key, observed=True)}
    div = load_divergence(args, species_study)
    triplets = build_triplets(species_study, div)
    log.info("Species -> study: %s", {s: sorted(v) for s, v in species_study.items()})
    for c in sorted({c for s in cov_sets for c in s}):
            conf = obs_full.groupby(c, observed=True)[args.species_key].nunique()
            if (conf == 1).any():
                log.info("Note: covariate '%s' has levels containing a single species (%s); correcting "
                         "on it also removes those species' differences from the rest.",
                         c, ", ".join(conf[conf == 1].index))
    if not triplets:
        raise RuntimeError("No batch-matched triplets: no two species share a study. The triplet "
                           "score needs at least one study profiling 2+ species.")
    use = pd.Series([x for t in triplets for x in t[1:]]).value_counts() / len(triplets)
    log.info("%d batch-matched triplet templates (before cell-type expansion)", len(triplets))
    if use.iloc[0] > 0.5:
        log.warning("'%s' is the near/far species in %.0f%% of triplets; the triplet score leans "
                    "heavily on it (a species-specific quirk would swing the score).",
                    use.index[0], 100 * use.iloc[0])

    idx = stratified_subsample(obs_full, [args.study_key, args.species_key, args.celltype_key],
                               args.cap_per_group, args.subsample, args.seed)
    X, obs = X_full[idx], obs_full.iloc[idx].reset_index(drop=True)
    log.info("Sweep subsample: %d cells (%d Tan)", len(idx), (obs[args.study_key] == args.tan_study).sum())

    pca = None
    if args.n_pcs:
        pca = PCA(n_components=args.n_pcs, random_state=args.seed).fit(X)
        X = pca.transform(X)
        log.info("Working in top %d PCs (%.1f%% variance)", args.n_pcs, 100 * pca.explained_variance_ratio_.sum())

    ctx = dict(species_study=species_study, div=div, triplets=triplets,
               tan_mask=(obs[args.study_key] == args.tan_study).to_numpy())
    results_path = outdir / "sweep_results.csv"
    rows, trip_details = [], {}

    t0 = time.time()
    base, tan_rel_raw, trip_raw = evaluate(X, obs, args, ctx, None)
    rows.append(dict(method="raw", covariates="-", theta=np.nan, lamb=np.nan, sigma=np.nan,
                     runtime_s=time.time() - t0, **base))
    trip_details[0] = trip_raw
    log.info("raw UCE: triplets=%.3f (%d; cross-study %.3f, within %.3f) phylo rho=%.3f "
             "[pure-study baseline %.3f] | %d Tan pairs",
             base["triplet_accuracy"], base["n_triplets"], base["triplet_acc_cross_study"],
             base["triplet_acc_within_study"], base["phylo_spearman"], base["phylo_study_baseline"],
             base["n_tan_pairs"])
    if base["n_tan_pairs"] < 3:
        log.warning("Few Tan same-type cross-species pairs (%d); lower --min-cells or raise "
                    "--cap-per-group for a reliable Tan check.", base["n_tan_pairs"])
    if base["n_triplets"] == 0:
        raise RuntimeError("No triplet has all three species' centroids for any cell type. Lower "
                           "--min-cells or raise --cap-per-group.")

    grid = list(itertools.product(cov_sets, args.thetas, args.lambdas, args.sigmas))
    for n, (cs, theta, lamb, sigma) in enumerate(grid, 1):
        t0 = time.time()
        try:
            Z, conv = run_harmony(X, obs, cs, theta, lamb, sigma, args.nclust, args.max_iter,
                                  args.seed, args.epsilon_harmony)
            if not conv["harmony_converged"]:
                log.warning("[%d/%d] %s theta=%g lambda=%g sigma=%g did NOT converge within %d "
                            "iterations (last relative change %.2e); raise --max-iter.", n, len(grid),
                            cs, theta, lamb, sigma, args.max_iter, conv["harmony_last_rel_change"])
            m, _, trip = evaluate(Z, obs, args, ctx, tan_rel_raw)
        except Exception as e:  # keep the sweep going
            log.error("[%d/%d] %s theta=%g lambda=%g sigma=%g failed: %s", n, len(grid), cs, theta, lamb, sigma, e)
            continue
        rows.append(dict(method="harmony", covariates=",".join(cs), theta=theta, lamb=lamb,
                         sigma=sigma, runtime_s=time.time() - t0, **conv, **m))
        trip_details[len(rows) - 1] = trip
        log.info("[%d/%d] %s theta=%g lambda=%g sigma=%g | triplets=%.3f phylo=%.3f cLISI=%.3f "
                 "sil=%.3f | Tan comp=%.2f rho=%.2f | %d iter%s (%.0fs)", n, len(grid), ",".join(cs), theta,
                 lamb, sigma, m["triplet_accuracy"], m["phylo_spearman"], m["clisi_norm"],
                 m["celltype_silhouette"], m["tan_compression"], m["tan_spearman"],
                 conv["harmony_n_iter"], "" if conv["harmony_converged"] else " NOT CONVERGED",
                 time.time() - t0)
        pd.DataFrame(rows).to_csv(results_path, index=False)  # incremental

    df = add_scores(pd.DataFrame(rows), args)
    h = df[(df.method == "harmony") & df.overall.notna()]
    if h.empty:
        raise RuntimeError(f"No Harmony setting produced a complete score; see {results_path}.")
    not_conv = h[~h.harmony_converged.astype(bool)]
    if len(not_conv):
        log.warning("%d setting(s) hit --max-iter without converging and are excluded from the "
                    "recommendation (see harmony_converged in %s).", len(not_conv), results_path)
        if len(not_conv) < len(h):
            h = h[h.harmony_converged.astype(bool)]
        else:
            log.warning("No setting converged; considering all of them anyway. Raise --max-iter.")
    ok = h[h.tan_pass]
    pool = ok if len(ok) else h
    if not len(ok):
        log.warning("No setting passed the Tan check; recommending the best overall score anyway.")
    best_i = pool.sort_values(["overall", "triplet_margin"], ascending=False).index[0]
    df["recommended"] = df.index == best_i
    df.sort_values("overall", ascending=False).to_csv(results_path, index=False)

    pd.concat([trip_details[0].assign(setting="raw"),
               trip_details[best_i].assign(setting="recommended")]).to_csv(
        outdir / "triplet_details.csv", index=False)

    best, raw = df.loc[best_i], df[df.method == "raw"].iloc[0]
    keys = ["overall", "triplet_accuracy", "triplet_margin", "triplet_acc_cross_study",
            "triplet_acc_within_study", "phylo_spearman", "celltype_score", "clisi_norm",
            "celltype_silhouette", "ilisi_species", "tan_compression", "tan_spearman"]
    best_params = dict(
        covariates=best.covariates.split(","), theta=best.theta, lamb=best.lamb, sigma=best.sigma,
        nclust=args.nclust, max_iter_harmony=args.max_iter, epsilon_harmony=args.epsilon_harmony,
        n_pcs=args.n_pcs, metric=args.metric,
        harmony_n_iter_on_subsample=int(best.harmony_n_iter),
        passed_tan_check=bool(best.tan_pass),
        beats_raw_on_triplets=bool(best.triplet_accuracy > raw.triplet_accuracy),
        phylo_study_baseline=float(raw.phylo_study_baseline),
        n_triplets=int(best.n_triplets),
        recommended={k: float(best[k]) for k in keys},
        raw_uce={k: float(raw[k]) for k in keys},
    )
    (outdir / "best_params.json").write_text(json.dumps(best_params, indent=2))

    plot_tradeoff(df, args, outdir / "sweep_tradeoff.png")
    plot_heatmaps(df, args, outdir / "sweep_overall_heatmap.png")

    cols = ["covariates", "theta", "lamb", "sigma", "harmony_n_iter", "overall", "triplet_accuracy", "triplet_margin",
            "phylo_spearman", "celltype_score", "tan_compression", "tan_spearman", "tan_pass"]
    print("\nTop 10 settings:")
    print(df.sort_values("overall", ascending=False)[cols].head(10).to_string(index=False, float_format="%.3f"))
    print(f"\nPure-study baseline for phylo_spearman on these pairs: {raw.phylo_study_baseline:.3f} "
          "(a study-only embedding would score this; read phylo_spearman relative to it)")
    print(f"Recommended: {best_params['covariates']} theta={best.theta:g} lambda={best.lamb:g} "
          f"sigma={best.sigma:g}" + ("" if best.tan_pass else "  (did NOT pass Tan check)"))
    print(f"  triplet accuracy {raw.triplet_accuracy:.3f} (raw) -> {best.triplet_accuracy:.3f}; "
          f"phylo rho {raw.phylo_spearman:.3f} -> {best.phylo_spearman:.3f}")
    if not best_params["beats_raw_on_triplets"]:
        print("  NOTE: no eligible Harmony setting improves triplet accuracy over raw UCE. Raw UCE is "
              "the better choice here; the setting above is only the least harmful correction.\n"
              "  A likely cause: the batches are phylogenetically unbalanced (e.g. Tan vs non-Tan), so "
              "matching their means also removes clade-level biology.")

    if args.apply_best:
        if not best_params["beats_raw_on_triplets"]:
            log.warning("--apply-best: saving the recommended setting even though it does not beat raw UCE.")
        log.info("Applying recommended setting to all %d cells", len(X_full))
        Xa = pca.transform(X_full) if pca is not None else X_full
        Z, conv = run_harmony(Xa, obs_full, best_params["covariates"], best.theta, best.lamb,
                              best.sigma, args.nclust, args.max_iter, args.seed, args.epsilon_harmony)
        log.info("Full-data Harmony: %d iterations, converged=%s", conv["harmony_n_iter"],
                 conv["harmony_converged"])
        if not conv["harmony_converged"]:
            log.warning("Full-data Harmony did NOT converge within %d iterations; raise --max-iter.",
                        args.max_iter)
        np.save(outdir / "X_uce_harmony_best.npy", Z.astype(np.float32))
        np.savetxt(outdir / "obs_names.txt", names_full, fmt="%s")
        log.info("Saved %s  (load with adata.obsm['X_uce_harmony_best'] = np.load(...))",
                 outdir / "X_uce_harmony_best.npy")

    log.info("Done. Results in %s", outdir)


if __name__ == "__main__":
    main()
