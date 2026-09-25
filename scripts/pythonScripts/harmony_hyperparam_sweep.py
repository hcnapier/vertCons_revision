"""
harmony_hyperparam_sweep.py

Grid-search Harmony hyperparameters on a UCE-embedded AnnData and score each
setting on three axes:

  1. Batch removal     - study mixing (LISI) measured WITHIN each species x cell-type
                         group present in >= 2 studies, normalized 0-1 against the best
                         achievable mixing. Global study iLISI is also reported, but not
                         scored: when study and species are confounded it penalises
                         correctly preserved species differences.
  2. Biology retained  - cLISI + silhouette on cell type (higher = cleaner types)
  3. Tan calibration   - the Tan 6-species / one-study subset has no study effect,
                         so a well-behaved correction should leave its
                         cross-species, same-cell-type distances (relative to
                         between-type distances) roughly unchanged:
                           * tan_compression : mean relative distance after / before
                                               (1.0 = preserved, <1 = flattened)
                           * tan_spearman    : rank correlation of the per-pair
                                               distances (is species ordering kept?)

It also reports the same-vs-different-type separation ratio and species iLISI
(a diagnostic for over-correction; not scored).

The recommendation is the best overall-scoring setting among those that keep
tan_compression >= --min-compression and tan_spearman >= --min-spearman.

Outputs (in --outdir):
  sweep_results.csv          one row per setting (+ raw UCE baseline), written
                             incrementally so a long sweep can be inspected mid-run
  best_params.json           recommended setting
  sweep_tradeoff.png         batch mixing vs Tan compression
  sweep_overall_heatmap.png  overall score over theta x lambda, per covariate set
  (with --apply-best) X_uce_harmony_best.npy + obs_names.txt for the full dataset

Example:
  python harmony_hyperparam_sweep.py merged_uce.h5ad \
      --study-key study --species-key species --celltype-key cell_type \
      --tan-study Tan --covariate-sets study study,technology \
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
    p.add_argument("--celltype-key", default="cell_type")
    p.add_argument("--tan-study", default="Tan", help="Value of --study-key marking the Tan subset")

    # Grid
    p.add_argument("--covariate-sets", nargs="+", default=["study"],
                   help="Comma-separated obs columns per set, e.g. study study,technology")
    p.add_argument("--thetas", nargs="+", type=float, default=[1, 2, 4, 6, 8])
    p.add_argument("--lambdas", nargs="+", type=float, default=[1.0, 0.5, 0.1])
    p.add_argument("--sigmas", nargs="+", type=float, default=[0.1])
    p.add_argument("--nclust", type=int, default=None, help="Harmony clusters (default: harmonypy's)")
    p.add_argument("--max-iter", type=int, default=20)

    # Data size / space
    p.add_argument("--subsample", type=int, default=60000, help="Max cells for the sweep (0 = all)")
    p.add_argument("--cap-per-group", type=int, default=2000,
                   help="Max cells per study x species x cell-type group before subsampling")
    p.add_argument("--n-pcs", type=int, default=0,
                   help="Run in top-N PCs of UCE instead of native 1280-d space (0 = native)")

    # Metrics
    p.add_argument("--k", type=int, default=30, help="Neighbours for LISI")
    p.add_argument("--min-cells", type=int, default=20, help="Min cells for a species x cell-type centroid")
    p.add_argument("--metric", default="euclidean", choices=["euclidean", "cosine"],
                   help="Centroid distance metric")
    p.add_argument("--min-compression", type=float, default=0.5)
    p.add_argument("--min-spearman", type=float, default=0.5)

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
        raise KeyError(f"Missing obs columns: {sorted(missing)}")
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
# Harmony
# -----------------------------------------------------------------------------
def run_harmony(X, meta, covars, theta, lamb, sigma, nclust, max_iter, seed):
    ho = harmonypy.run_harmony(
        X, meta, covars,
        theta=[theta] * len(covars),
        lamb=[lamb] * len(covars),
        sigma=sigma,
        nclust=nclust,
        max_iter_harmony=max_iter,
        verbose=False,
        random_state=seed,
    )
    Z = np.asarray(ho.Z_corr)
    # harmonypy returns dims x cells
    if Z.shape[0] != X.shape[0] or (Z.shape[0] == Z.shape[1] and Z.shape[1] != X.shape[1]):
        Z = Z.T
    return Z


# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------
def knn_indices(Z, k):
    nn = NearestNeighbors(n_neighbors=k + 1, n_jobs=-1).fit(Z)
    return nn.kneighbors(Z, return_distance=False)[:, 1:]


def lisi(nn_idx, labels):
    """kNN inverse-Simpson diversity (uniform weights). Returns per-cell LISI."""
    codes, _ = pd.factorize(labels)
    L = codes[nn_idx]
    n, k = L.shape
    counts = np.zeros((n, codes.max() + 1))
    np.add.at(counts, (np.repeat(np.arange(n), k), L.ravel()), 1)
    p = counts / k
    return 1.0 / (p ** 2).sum(1)


def matched_study_mixing(Z, study, species, ct, k, min_cells):
    """Study mixing measured only where it is identifiable: within each
    species x cell-type group that is present in >= 2 studies. Global iLISI
    penalises preserved species differences whenever study and species are
    confounded; this does not. Returns the cell-weighted mean of normalized
    LISI (0 = studies fully separate, 1 = perfectly mixed), plus #groups."""
    df = pd.DataFrame({"st": study, "sp": species, "ct": ct})
    scores, weights = [], []
    for g in df.groupby(["sp", "ct"], observed=True).indices.values():
        st = df.st.to_numpy()[g]
        vc = pd.Series(st).value_counts()
        vc = vc[vc >= min_cells]
        if len(vc) < 2:
            continue
        keep = g[np.isin(st, vc.index)]
        kk = min(k, len(keep) - 1)
        nn = NearestNeighbors(n_neighbors=kk + 1).fit(Z[keep])
        idx = nn.kneighbors(Z[keep], return_distance=False)[:, 1:]
        l = lisi(idx, df.st.to_numpy()[keep])
        # best achievable LISI given the group's study proportions
        p = vc.to_numpy() / vc.sum()
        l_max = 1.0 / (p ** 2).sum()
        scores.append(np.clip((l.mean() - 1) / (l_max - 1), 0, 1))
        weights.append(len(keep))
    if not scores:
        return np.nan, 0
    return float(np.average(scores, weights=weights)), len(scores)


def centroids(Z, species, ct, min_cells):
    groups = pd.DataFrame({"sp": species, "ct": ct}).groupby(["sp", "ct"], observed=True).indices
    keys = sorted(k for k, v in groups.items() if len(v) >= min_cells)
    C = np.vstack([Z[groups[k]].mean(0) for k in keys])
    return keys, C


def separation(keys, C, metric):
    """Same-type cross-species vs different-type centroid distances.
    Returns summary + per-pair relative distances keyed by (cell_type, sp1, sp2)."""
    D = cdist(C, C, metric=metric)
    sp = np.array([k[0] for k in keys])
    ct = np.array([k[1] for k in keys])
    i, j = np.triu_indices(len(keys), 1)
    d = D[i, j]
    same = (ct[i] == ct[j]) & (sp[i] != sp[j])
    diff = ct[i] != ct[j]
    if same.sum() == 0 or diff.sum() == 0:
        return dict(same=np.nan, diff=np.nan, ratio=np.nan), {}
    mean_diff = d[diff].mean()
    rel = {(ct[a], *sorted((sp[a], sp[b]))): d_ / mean_diff
           for a, b, d_ in zip(i[same], j[same], d[same])}
    return dict(same=d[same].mean(), diff=mean_diff, ratio=mean_diff / d[same].mean()), rel


def evaluate(Z, obs, args, tan_mask, tan_rel_raw, seed):
    ct = obs[args.celltype_key].to_numpy()
    study = obs[args.study_key].to_numpy()
    species = obs[args.species_key].to_numpy()
    n_study, n_ct = len(np.unique(study)), len(np.unique(ct))

    nn = knn_indices(Z, args.k)
    ilisi = np.median(lisi(nn, study))
    clisi = np.median(lisi(nn, ct))
    slisi = np.median(lisi(nn, species))
    sil = silhouette_score(Z, ct, sample_size=min(len(Z), 10000), random_state=seed,
                           metric="cosine" if args.metric == "cosine" else "euclidean")

    matched, n_groups = matched_study_mixing(Z, study, species, ct, args.k, args.min_cells)

    keys, C = centroids(Z, species, ct, args.min_cells)
    sep, _ = separation(keys, C, args.metric)

    tk, tC = centroids(Z[tan_mask], species[tan_mask], ct[tan_mask], args.min_cells)
    _, tan_rel = separation(tk, tC, args.metric)

    out = dict(
        matched_study_mixing=matched,
        n_matched_groups=n_groups,
        ilisi_study=ilisi,
        ilisi_study_norm=(ilisi - 1) / max(n_study - 1, 1),
        clisi_norm=(n_ct - clisi) / max(n_ct - 1, 1),
        ilisi_species=slisi,
        celltype_silhouette=sil,
        sep_same_type=sep["same"],
        sep_diff_type=sep["diff"],
        separation_ratio=sep["ratio"],
    )
    if tan_rel_raw is None:  # baseline run
        out.update(tan_compression=1.0, tan_spearman=1.0, n_tan_pairs=len(tan_rel))
        return out, tan_rel

    shared = sorted(set(tan_rel) & set(tan_rel_raw))
    raw_v = np.array([tan_rel_raw[s] for s in shared])
    new_v = np.array([tan_rel[s] for s in shared])
    out.update(
        tan_compression=new_v.mean() / raw_v.mean() if len(shared) else np.nan,
        tan_spearman=spearmanr(raw_v, new_v).statistic if len(shared) > 2 else np.nan,
        n_tan_pairs=len(shared),
    )
    return out, tan_rel


def add_scores(df):
    pres = (1 - (1 - df.tan_compression).abs().clip(upper=1)) * 0.5 + df.tan_spearman.clip(lower=0) * 0.5
    df["tan_preservation"] = pres
    df["batch_score"] = df.matched_study_mixing
    df["bio_score"] = (df.clisi_norm + (df.celltype_silhouette + 1) / 2 + pres) / 3
    df["overall"] = 0.4 * df.batch_score + 0.6 * df.bio_score
    return df


# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
def plot_tradeoff(df, args, path):
    fig, ax = plt.subplots(figsize=(7, 5.5))
    h = df[df.method == "harmony"]
    markers = dict(zip(sorted(h.lamb.unique()), "osD^v<>p*"))
    thetas = sorted(h.theta.unique())
    cmap = plt.get_cmap("viridis", len(thetas))
    tnorm = matplotlib.colors.BoundaryNorm(np.arange(len(thetas) + 1) - 0.5, len(thetas))
    for cov, sub in h.groupby("covariates"):
        for _, r in sub.iterrows():
            ax.scatter(r.matched_study_mixing, r.tan_compression, marker=markers[r.lamb],
                       color=cmap(tnorm(thetas.index(r.theta))), s=60,
                       edgecolor="k" if cov == h.covariates.iloc[0] else "red", linewidth=0.8)
    raw = df[df.method == "raw"].iloc[0]
    ax.scatter(raw.matched_study_mixing, raw.tan_compression, marker="X", c="grey", s=120, label="raw UCE")
    ax.axhline(args.min_compression, ls="--", c="grey", lw=1)
    best = df[df.recommended]
    if len(best):
        ax.scatter(best.matched_study_mixing, best.tan_compression, s=300, facecolors="none",
                   edgecolors="crimson", linewidth=2, label="recommended")
    for lam, m in markers.items():
        ax.scatter([], [], marker=m, c="k", label=f"lambda={lam:g}")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=tnorm)
    cb = fig.colorbar(sm, ax=ax, label="theta", ticks=range(len(thetas)))
    cb.ax.set_yticklabels([f"{t:g}" for t in thetas])
    ax.set_xlabel("Study mixing within species x cell type  →  better batch removal")
    ax.set_ylabel("Tan relative species distance (after / before)")
    ax.set_title("Harmony sweep: batch removal vs species-distance preservation\n"
                 "(edge: black = first covariate set, red = others)", fontsize=10)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_heatmaps(df, args, path):
    h = df[df.method == "harmony"]
    groups = list(h.groupby(["covariates", "sigma"]))
    fig, axes = plt.subplots(1, len(groups), figsize=(4.5 * len(groups), 4), squeeze=False)
    vmin, vmax = h.overall.min(), h.overall.max()
    for ax, ((cov, sig), sub) in zip(axes[0], groups):
        piv = sub.pivot_table(index="theta", columns="lamb", values="overall")
        passes = ((sub.tan_compression >= args.min_compression) & (sub.tan_spearman >= args.min_spearman))
        ok = sub.assign(ok=passes.astype(float)).pivot_table(index="theta", columns="lamb", values="ok")
        im = ax.imshow(piv.values, cmap="magma", vmin=vmin, vmax=vmax, aspect="auto", origin="lower")
        ax.set_xticks(range(piv.shape[1]), [f"{v:g}" for v in piv.columns])
        ax.set_yticks(range(piv.shape[0]), [f"{v:g}" for v in piv.index])
        for (r, c), v in np.ndenumerate(piv.values):
            failed = ok.values[r, c] < 1
            if failed:
                ax.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1, facecolor="lightgrey",
                                           edgecolor="grey", hatch="///", lw=0))
            ax.text(c, r, f"{v:.2f}" + ("\nfails Tan" if failed else ""), ha="center", va="center",
                    fontsize=8, color="black" if failed or v >= (vmin + vmax) / 2 else "white")
        ax.set_xlabel("lambda")
        ax.set_ylabel("theta")
        ax.set_title(f"{cov} | sigma={sig:g}", fontsize=10)
    fig.colorbar(im, ax=axes.ravel().tolist(), label="overall score")
    fig.suptitle(f"Hatched = fails Tan check (compression < {args.min_compression:g} or Spearman < {args.min_spearman:g}); not eligible", fontsize=9, y=1.02)
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
                        "measuring, so its separation/Tan metrics are partly circular.", cs)

    idx = stratified_subsample(obs_full, [args.study_key, args.species_key, args.celltype_key],
                               args.cap_per_group, args.subsample, args.seed)
    X, obs = X_full[idx], obs_full.iloc[idx].reset_index(drop=True)
    log.info("Sweep subsample: %d cells (%d Tan)", len(idx), (obs[args.study_key] == args.tan_study).sum())

    pca = None
    if args.n_pcs:
        pca = PCA(n_components=args.n_pcs, random_state=args.seed).fit(X)
        X = pca.transform(X)
        log.info("Working in top %d PCs (%.1f%% variance)", args.n_pcs, 100 * pca.explained_variance_ratio_.sum())

    tan_mask = (obs[args.study_key] == args.tan_study).to_numpy()
    results_path = outdir / "sweep_results.csv"
    rows = []

    t0 = time.time()
    base, tan_rel_raw = evaluate(X, obs, args, tan_mask, None, args.seed)
    rows.append(dict(method="raw", covariates="-", theta=np.nan, lamb=np.nan, sigma=np.nan,
                     runtime_s=time.time() - t0, **base))
    log.info("raw UCE: mix=%.3f sep_ratio=%.2f (%d Tan pairs, %d matched groups)",
             base["matched_study_mixing"], base["separation_ratio"], base["n_tan_pairs"], base["n_matched_groups"])
    if base["n_tan_pairs"] < 3:
        log.warning("Few Tan same-type cross-species pairs (%d); lower --min-cells or raise "
                    "--cap-per-group for a reliable calibration.", base["n_tan_pairs"])

    grid = list(itertools.product(cov_sets, args.thetas, args.lambdas, args.sigmas))
    for n, (cs, theta, lamb, sigma) in enumerate(grid, 1):
        t0 = time.time()
        try:
            Z = run_harmony(X, obs, cs, theta, lamb, sigma, args.nclust, args.max_iter, args.seed)
            m, _ = evaluate(Z, obs, args, tan_mask, tan_rel_raw, args.seed)
        except Exception as e:  # keep the sweep going
            log.error("[%d/%d] %s theta=%g lambda=%g sigma=%g failed: %s", n, len(grid), cs, theta, lamb, sigma, e)
            continue
        rows.append(dict(method="harmony", covariates=",".join(cs), theta=theta, lamb=lamb,
                         sigma=sigma, runtime_s=time.time() - t0, **m))
        log.info("[%d/%d] %s theta=%g lambda=%g sigma=%g | mix=%.3f cLISI=%.3f sil=%.3f "
                 "sep=%.2f Tan comp=%.2f rho=%.2f (%.0fs)", n, len(grid), ",".join(cs), theta, lamb,
                 sigma, m["matched_study_mixing"], m["clisi_norm"], m["celltype_silhouette"],
                 m["separation_ratio"], m["tan_compression"], m["tan_spearman"], time.time() - t0)
        pd.DataFrame(rows).to_csv(results_path, index=False)  # incremental

    df = add_scores(pd.DataFrame(rows))
    h = df[df.method == "harmony"]
    ok = h[(h.tan_compression >= args.min_compression) & (h.tan_spearman >= args.min_spearman)]
    if len(ok):
        best_i = ok.overall.idxmax()
    else:
        log.warning("No setting met the Tan thresholds; recommending best overall score instead.")
        best_i = h.overall.idxmax()
    df["recommended"] = df.index == best_i
    df.sort_values("overall", ascending=False).to_csv(results_path, index=False)

    best = df.loc[best_i]
    best_params = dict(
        covariates=best.covariates.split(","), theta=best.theta, lamb=best.lamb, sigma=best.sigma,
        nclust=args.nclust, max_iter_harmony=args.max_iter, n_pcs=args.n_pcs,
        met_tan_thresholds=bool(len(ok)),
        metrics={k: float(best[k]) for k in ["overall", "batch_score", "bio_score", "matched_study_mixing", "ilisi_study_norm",
                                              "clisi_norm", "celltype_silhouette", "separation_ratio",
                                              "tan_compression", "tan_spearman"]},
        raw_baseline={k: float(df.loc[df.method == "raw", k].iloc[0])
                      for k in ["matched_study_mixing", "ilisi_study_norm", "clisi_norm", "separation_ratio"]},
    )
    (outdir / "best_params.json").write_text(json.dumps(best_params, indent=2))

    plot_tradeoff(df, args, outdir / "sweep_tradeoff.png")
    plot_heatmaps(df, args, outdir / "sweep_overall_heatmap.png")

    cols = ["covariates", "theta", "lamb", "sigma", "overall", "matched_study_mixing", "clisi_norm",
            "celltype_silhouette", "separation_ratio", "tan_compression", "tan_spearman"]
    print("\nTop 10 settings:")
    print(df.sort_values("overall", ascending=False)[cols].head(10).to_string(index=False, float_format="%.3f"))
    print(f"\nRecommended: {best_params['covariates']} theta={best.theta:g} lambda={best.lamb:g} "
          f"sigma={best.sigma:g}" + ("" if len(ok) else "  (did NOT meet Tan thresholds)"))

    if args.apply_best:
        log.info("Applying recommended setting to all %d cells", len(X_full))
        Xa = pca.transform(X_full) if pca is not None else X_full
        Z = run_harmony(Xa, obs_full, best_params["covariates"], best.theta, best.lamb, best.sigma,
                        args.nclust, args.max_iter, args.seed)
        np.save(outdir / "X_uce_harmony_best.npy", Z.astype(np.float32))
        np.savetxt(outdir / "obs_names.txt", names_full, fmt="%s")
        log.info("Saved %s  (load with adata.obsm['X_uce_harmony_best'] = np.load(...))",
                 outdir / "X_uce_harmony_best.npy")

    log.info("Done. Results in %s", outdir)


if __name__ == "__main__":
    main()
