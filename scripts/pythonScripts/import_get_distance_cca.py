"""
Merge the CCA-integrated embedding computed in R (integrate_cca.R) back
into the combined AnnData, then run the same distance/diagnostic analysis
on it (centroid_distance_matrix, same_vs_different_type_scores,
study_controlled_diagnostic, heatmaps) as the other scripts run on X_uce /
X_uce_harmony -- so you can directly compare species/cell-type separation
in CCA-integrated UCE space against the other representations.

python import_get_distance_cca.py
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt

# ============================== CONFIG ======================================

PATH = "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances_combined.h5ad"

SEURAT_DIR = "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/seurat_export"

# Maps the .obsm key to store each result under -> the CSV filename R wrote.
# Only files that actually exist are imported; missing ones are skipped
# with a warning.
EMBEDDINGS_TO_IMPORT = {
    "X_uce_cca": "X_uce_cca.csv",
}

# Where to save the updated AnnData. Set equal to PATH to overwrite in place.
OUT_PATH = "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances_combined_integrated.h5ad"

# --- Distance / diagnostic analysis on the imported embedding(s) ---

# Which imported .obsm key(s) to run the full distance analysis on. Usually
# just the CCA embedding, but you can add others (e.g. "X_uce",
# "X_uce_harmony") here too if they're already present in the AnnData, to
# get everything in one comparable set of outputs.
REPS_FOR_ANALYSIS = ["X_uce_cca"]

METRIC = "cosine"

SPECIES_PHYLO_ORDER = [
    "human", "macaque", "guineaPig", "rat", "mouse", "rabbit",
    "pig", "cow", "goat", "dog",
]

CELL_TYPE_ONTOLOGY_ORDER = [
    "ctb", "evt", "invasive", "stb", "s-tgc", "spt", "gc", "unc", "bnc",
    "epi", "stro", "mes", "endo", "leu", "bcell", "tcell", "nkcells",
    "mono", "mac", "dc", "neu",
]

# Study used for study_controlled_diagnostic()'s study-controlled subset.
# If None, auto-picks the study with the most distinct species.
DIAGNOSTIC_STUDY = "Tan"

OUT_PREFIX = "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances_cca"

# =============================================================================


def _make_rank_lookup(order_list, label_for_warning):
    """See uce_celltype_distance_comparison_v2.py for the full docstring."""
    ranks = {val: i for i, val in enumerate(order_list)}

    def rank(value):
        if value not in ranks:
            print(f"WARNING: '{value}' not found in your {label_for_warning} "
                  f"ordering list — placing it at the end. Add it to the "
                  f"CONFIG list if you want explicit control over its position.")
            return (len(order_list), value)
        return (ranks[value], "")

    return rank


def centroid_distance_matrix(combined, metric="cosine", use_rep="X_uce",
                              species_order=None):
    """Pairwise distances between (species, cell_type) centroids in use_rep."""
    X = combined.obsm[use_rep]
    groups = combined.obs["group"].values
    centroids = pd.DataFrame(X, index=groups).groupby(level=0).mean()

    if species_order is not None:
        species_rank = _make_rank_lookup(species_order, "SPECIES_PHYLO_ORDER")

        def sort_key(label):
            species, cell_type = label.split(" | ", 1)
            return (species_rank(species), cell_type)

        ordered_labels = sorted(centroids.index, key=sort_key)
        centroids = centroids.loc[ordered_labels]

    dist = squareform(pdist(centroids.values, metric=metric))
    return pd.DataFrame(dist, index=centroids.index, columns=centroids.index)


def same_vs_different_type_scores(combined, metric="cosine", max_cells_per_group=500,
                                   use_rep="X_uce", verbose=True):
    """Same-type-cross-dataset vs. different-type distances in use_rep.

    Returns (same_type_dists, diff_type_dists, ratio); ratio is
    diff.mean()/same.mean(), or None if there were no same-type-cross-
    dataset pairs.
    """
    X = combined.obsm[use_rep]
    obs = combined.obs.reset_index(drop=True)
    rng = np.random.default_rng(0)

    def subsample_idx(mask):
        idx = np.flatnonzero(mask.values)
        if len(idx) > max_cells_per_group:
            idx = rng.choice(idx, size=max_cells_per_group, replace=False)
        return idx

    same_type_dists = []
    for ct in obs["cell_type_std"].unique():
        ct_mask = obs["cell_type_std"] == ct
        if obs.loc[ct_mask, "dataset"].nunique() < 2:
            continue
        idx = subsample_idx(ct_mask)
        if len(idx) < 2:
            continue
        d = pairwise_distances(X[idx], metric=metric)
        iu = np.triu_indices_from(d, k=1)
        same_type_dists.append(d[iu])

    if not same_type_dists:
        if verbose:
            print("No cell type appears in 2+ datasets — can't compute "
                  "same-type-cross-dataset scores.")
        same_type_dists = np.array([])
    else:
        same_type_dists = np.concatenate(same_type_dists)

    idx_all = subsample_idx(pd.Series(True, index=obs.index))
    d_all = pairwise_distances(X[idx_all], metric=metric)
    types_all = obs["cell_type_std"].values[idx_all]
    iu = np.triu_indices_from(d_all, k=1)
    diff_mask = types_all[iu[0]] != types_all[iu[1]]
    diff_type_dists = d_all[iu][diff_mask]

    ratio = None
    if verbose:
        print("\n--- Same-type-vs-different-type separation ---")
    if len(same_type_dists) and verbose:
        print(f"Same cell type, different dataset: "
              f"mean={same_type_dists.mean():.4f}, n={len(same_type_dists)}")
    if verbose:
        print(f"Different cell type (pooled):        "
              f"mean={diff_type_dists.mean():.4f}, n={len(diff_type_dists)}")
    if len(same_type_dists):
        ratio = diff_type_dists.mean() / same_type_dists.mean()
        if verbose:
            print(f"Separation ratio (diff/same): {ratio:.2f}x")

    return same_type_dists, diff_type_dists, ratio


def study_controlled_diagnostic(combined, study_name=None, metric="cosine",
                                 max_cells_per_group=500, reps_to_compare=None):
    """Isolate species/cell-type separation from study/technology batch
    effects, using the subset of datasets sharing one multi-species study.
    See uce_celltype_distance_comparison_v2.py for the full docstring.
    """
    if reps_to_compare is None:
        reps_to_compare = [r for r in
                            ["X_uce", "X_uce_centered", "X_uce_harmony", "X_uce_cca"]
                            if r in combined.obsm]

    if study_name is None:
        species_per_study = combined.obs.groupby("study")["species"].nunique()
        candidates = species_per_study[species_per_study > 1]
        if candidates.empty:
            raise ValueError(
                "No study in your data profiles more than one species — "
                "can't build a study-controlled subset."
            )
        top_studies = candidates[candidates == candidates.max()].sort_index()
        study_name = top_studies.index[0]
        print(f"DIAGNOSTIC_STUDY not set — auto-picked study '{study_name}' "
              f"({candidates[study_name]} species).")

    if study_name not in set(combined.obs["study"]):
        raise ValueError(f"study_name '{study_name}' not found in .obs['study']. "
                          f"Available: {sorted(combined.obs['study'].unique())}")

    sub = combined[combined.obs["study"] == study_name].copy()
    n_species_sub = sub.obs["species"].nunique()
    species_list = sorted(sub.obs["species"].unique())
    print(f"\n=== Study-controlled diagnostic: study = '{study_name}' ===")
    print(f"{sub.obs['dataset'].nunique()} dataset(s), {n_species_sub} species: "
          f"{species_list}")
    if n_species_sub < 2:
        raise ValueError(
            f"Study '{study_name}' only has {n_species_sub} species — "
            f"nothing to isolate species separation from."
        )

    rows = []
    for rep in reps_to_compare:
        if rep not in combined.obsm:
            print(f"Skipping rep '{rep}': not found in combined.obsm.")
            continue
        print(f"\n--- rep = {rep} ---")
        print("Full dataset (all species, all studies):")
        _, _, full_ratio = same_vs_different_type_scores(
            combined, metric=metric, max_cells_per_group=max_cells_per_group,
            use_rep=rep, verbose=True)
        print(f"\nStudy-controlled subset ('{study_name}' only, "
              f"{n_species_sub} species):")
        _, _, sub_ratio = same_vs_different_type_scores(
            sub, metric=metric, max_cells_per_group=max_cells_per_group,
            use_rep=rep, verbose=True)
        rows.append({
            "rep": rep,
            "full_dataset_ratio": full_ratio,
            f"study_controlled_ratio ({study_name})": sub_ratio,
        })

    result = pd.DataFrame(rows)
    print("\n=== Summary: full-dataset vs. study-controlled separation ratio ===")
    print(result.to_string(index=False))
    return result


def reorder_by_celltype(dist_df, cell_type_order=None, species_order=None):
    """Reorder a (species | cell_type) distance matrix cell-type-first.
    See uce_celltype_distance_comparison_v2.py for the full docstring.
    """
    cell_type_rank = (_make_rank_lookup(cell_type_order, "CELL_TYPE_ONTOLOGY_ORDER")
                       if cell_type_order is not None else None)
    species_rank = (_make_rank_lookup(species_order, "SPECIES_PHYLO_ORDER")
                     if species_order is not None else None)

    def sort_key(label):
        species, cell_type = label.split(" | ", 1)
        ct_key = cell_type_rank(cell_type) if cell_type_rank else cell_type
        sp_key = species_rank(species) if species_rank else species
        return (ct_key, sp_key)

    new_order = sorted(dist_df.index, key=sort_key)
    relabeled = {lbl: " | ".join(reversed(lbl.split(" | ", 1))) for lbl in new_order}
    reordered = dist_df.loc[new_order, new_order].rename(index=relabeled, columns=relabeled)
    return reordered


def plot_heatmap(dist_df, out_path, title="UCE-CCA centroid distances"):
    fig, ax = plt.subplots(figsize=(0.55 * len(dist_df) + 3, 0.55 * len(dist_df) + 3))
    im = ax.imshow(dist_df.values, cmap="viridis")
    ax.set_xticks(range(len(dist_df)))
    ax.set_yticks(range(len(dist_df)))
    ax.set_xticklabels(dist_df.columns, rotation=90, fontsize=30)
    ax.set_yticklabels(dist_df.index, fontsize=30)
    cbar = fig.colorbar(im, ax=ax, label="distance")
    cbar.ax.tick_params(labelsize=30)
    cbar.set_label("distance", fontsize=30)
    ax.set_title(title, fontsize=50)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Saved heatmap to {out_path}")


def run_distance_analysis(combined, use_rep):
    """Run centroid distances, separation score, study-controlled
    diagnostic, and heatmaps for one .obsm rep, mirroring what the other
    scripts do for X_uce / X_uce_harmony."""
    print(f"\n########## Distance analysis: use_rep = '{use_rep}' ##########")

    dist_df = centroid_distance_matrix(combined, metric=METRIC, use_rep=use_rep,
                                        species_order=SPECIES_PHYLO_ORDER)
    dist_csv = f"{OUT_PREFIX}_{use_rep}_centroid_distances.csv"
    dist_df.to_csv(dist_csv)
    print(f"Saved centroid distance matrix to {dist_csv}")

    same_vs_different_type_scores(combined, metric=METRIC, use_rep=use_rep)

    diagnostic_df = study_controlled_diagnostic(
        combined, study_name=DIAGNOSTIC_STUDY, metric=METRIC,
        reps_to_compare=[use_rep])
    diag_csv = f"{OUT_PREFIX}_{use_rep}_study_controlled_diagnostic.csv"
    diagnostic_df.to_csv(diag_csv, index=False)
    print(f"Saved study-controlled diagnostic to {diag_csv}")

    plot_heatmap(dist_df, f"{OUT_PREFIX}_{use_rep}_heatmap.png",
                 title=f"{use_rep} centroid distances: species (phylogenetic order) | cell type")

    dist_df_by_celltype = reorder_by_celltype(
        dist_df, cell_type_order=CELL_TYPE_ONTOLOGY_ORDER,
        species_order=SPECIES_PHYLO_ORDER)
    plot_heatmap(dist_df_by_celltype, f"{OUT_PREFIX}_{use_rep}_heatmap_by_celltype.png",
                 title=f"{use_rep} centroid distances: cell type (ontology order) | species")

# =============================================================================


def main():
    combined = sc.read_h5ad(PATH)

    for obsm_key, fname in EMBEDDINGS_TO_IMPORT.items():
        fpath = os.path.join(SEURAT_DIR, fname)
        if not os.path.exists(fpath):
            print(f"Skipping '{obsm_key}': {fpath} not found.")
            continue

        emb = pd.read_csv(fpath, index_col=0)

        missing = set(combined.obs_names) - set(emb.index)
        extra = set(emb.index) - set(combined.obs_names)
        if missing:
            raise ValueError(
                f"{fname}: {len(missing)} cells from the AnnData are missing "
                f"in this embedding (e.g. {list(missing)[:5]}). Re-check that "
                f"export_uce_for_seurat.py and this AnnData are in sync."
            )
        if extra:
            print(f"NOTE: {fname} has {len(extra)} extra cells not in the "
                  f"AnnData -- dropping them.")

        emb = emb.loc[combined.obs_names]
        combined.obsm[obsm_key] = emb.values
        print(f"Imported '{obsm_key}': {emb.shape[0]} cells x {emb.shape[1]} dims.")

    combined.write_h5ad(OUT_PATH)
    print(f"\nSaved updated AnnData to {OUT_PATH}")

    for use_rep in REPS_FOR_ANALYSIS:
        if use_rep not in combined.obsm:
            print(f"\nSkipping distance analysis for '{use_rep}': not found "
                  f"in .obsm (check EMBEDDINGS_TO_IMPORT / REPS_FOR_ANALYSIS "
                  f"are in sync, and that the corresponding R output CSV "
                  f"existed).")
            continue
        run_distance_analysis(combined, use_rep)


if __name__ == "__main__":
    main()
