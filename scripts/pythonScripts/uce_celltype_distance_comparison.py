"""
Compare UCE embedding distances for cell types across multiple datasets.

Workflow:
  1. (Assumes you've already run UCE on each dataset, producing .h5ad files
     with UCE embeddings in .obsm["X_uce"]. See the "RUN UCE" section below
     for the commands to do that if you haven't yet.)
  2. Load and merge the embedded datasets.
  3. Compute cell-type centroid distances by species (pooling all datasets
     from the same species together).
  4. Compute a same-type-vs-different-type separation score.
  5. Save a distance heatmap and a UMAP colored by dataset / cell type.

Fill in the CONFIG section below, then run:
    python uce_celltype_distance_comparison.py
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.spatial.distance import pdist, squareform, cosine
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt

# ============================== CONFIG ======================================

# Paths to your UCE-embedded .h5ad files (output of eval_single_anndata / 
# uce-eval-single-anndata — each should already have .obsm["X_uce"]).
DATASET_PATHS = {
    "dataset1": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Bos_taurus_uce_adata.h5ad",
    "dataset2": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Canis_lupus_familiaris_uce_adata.h5ad",
    "dataset3": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Cavia_porcellus_uce_adata.h5ad",
    "dataset4": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Homo_sapiens_uce_adata.h5ad",
    "dataset5": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Macaca_fascicularis_uce_adata.h5ad",
    "dataset6": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Mus_musculus_uce_adata.h5ad",
    "dataset7": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Oryctolagus_cuniculus_uce_adata.h5ad",
    "dataset8": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Rattus_norvegicus_shrutx_uce_adata.h5ad",
    "dataset9": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Sus_scrofa_uce_adata.h5ad",
}

# Name of the .obs column holding cell type labels in EACH dataset.
# If it's the same column name everywhere, just repeat it.
CELL_TYPE_COLS = {
    "dataset1": "napierCellTypes",
    "dataset2": "napierCellTypes",
    "dataset3": "napierCellTypes",
    "dataset4": "napierCellTypes",
    "dataset5": "napierCellTypes",
    "dataset6": "napierCellTypes",
    "dataset7": "napierCellTypes",
    "dataset8": "napierCellTypes",
    "dataset9": "napierCellTypes"
}

# Species of each dataset (used only for bookkeeping/reporting here — UCE
# tokenizes genes via ESM2 protein embeddings rather than gene symbols, so
# it aligns homologous genes across species internally. No ortholog
# filtering needed; just make sure you ran UCE itself with the correct
# --species flag for each dataset when you generated the embeddings).
SPECIES = {
    "dataset1": "cow",
    "dataset2": "dog",
    "dataset3": "guineaPig",
    "dataset4": "human",
    "dataset5": "macaque",
    "dataset6": "mouse",
    "dataset7": "rabbit",
    "dataset8": "rat",
    "dataset9": "pig",
}


# Order species should appear in the species-organized heatmap — e.g. by
# phylogenetic distance from human (closest first). Use the exact species
# names you set in SPECIES above. Any species present in your data but NOT
# listed here gets appended at the end (alphabetically), with a warning —
# it won't be silently dropped, just placed at the back.
SPECIES_PHYLO_ORDER = [
    "human",
    "macaque",
    "guineaPig",
    "mouse",
    "rabbit",
    "pig", 
    "cow", 
    "dog"
]
 
# Order cell types should appear in the cell-type-organized heatmap — e.g.
# by Cell Ontology (CL) hierarchy (broad lineage groupings together, related
# types adjacent). Use the exact lowercased values that end up in
# cell_type_std (i.e. however your CELL_TYPE_COLS values look after
# .str.strip().str.lower()). Anything present in your data but not listed
# here gets appended at the end (alphabetically), with a warning.
CELL_TYPE_ONTOLOGY_ORDER = [
    "ctb",
    "evt",
    "invasive",
    "stb",
    "s-tgc",
    "gc",
    "unc", 
    "bnc",
    "epi",
    "stro",
    "mes",
    "endo",
    "leu",
    "bcell",
    "tcell",
    "nkcell",
    "mono",
    "mac",
    "DC",
    "neu"
]

# Distance metric for comparisons: "euclidean" or "cosine"
METRIC = "cosine"

# If True, mean-center each species' embeddings (subtract that species'
# overall mean X_uce vector from every one of its cells) before computing
# distances, UMAP, and the separation score. This removes a per-species
# "offset" in embedding space — useful if species identity is dominating
# over cell-type similarity in your comparisons. 
CENTER_BY_SPECIES = True

# Where to save outputs
OUT_PREFIX = "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances"

# =============================================================================

def load_and_merge(dataset_paths, cell_type_cols, species_map):
    """Load each embedded dataset, standardize cell type column name, tag
    dataset origin and species, and concatenate into one AnnData."""
    adatas = []
    for name, path in dataset_paths.items():
        a = sc.read_h5ad(path)
        if "X_uce" not in a.obsm:
            raise ValueError(f"{name} ({path}) has no .obsm['X_uce'] — "
                              f"make sure you ran UCE on it first.")
        ct_col = cell_type_cols[name]
        if ct_col not in a.obs:
            raise ValueError(f"{name}: column '{ct_col}' not found in .obs. "
                              f"Available columns: {list(a.obs.columns)}")
        a.obs["cell_type_std"] = a.obs[ct_col].astype(str).str.strip().str.lower()
        a.obs["dataset"] = name
        a.obs["species"] = species_map.get(name, "unknown")
        # keep only what we need to avoid var mismatch issues on concat
        a = ad.AnnData(X=a.X, obs=a.obs[["cell_type_std", "dataset", "species"]].copy(),
                        obsm={"X_uce": a.obsm["X_uce"]})
        adatas.append(a)
 
    combined = ad.concat(adatas, join="outer", label="batch", index_unique="-")
    combined.obs["group"] = (combined.obs["species"].astype(str) + " | " +
                              combined.obs["cell_type_std"].astype(str))
    n_species = combined.obs["species"].nunique()
    print(f"Loaded {combined.n_obs} cells across {len(dataset_paths)} datasets "
          f"({n_species} species), "
          f"{combined.obs['cell_type_std'].nunique()} unique cell type labels.")
    return combined
 
 
def center_by_species(combined):
    """
    Subtract each species' mean X_uce vector from every one of its cells.
    Stores the result in .obsm["X_uce_centered"] and leaves the original
    .obsm["X_uce"] untouched, so you can compare centered vs. uncentered
    results without re-loading anything.
    """
    X = combined.obsm["X_uce"]
    species = combined.obs["species"].values
    X_centered = np.empty_like(X)
 
    for sp in np.unique(species):
        mask = species == sp
        sp_mean = X[mask].mean(axis=0, keepdims=True)
        X_centered[mask] = X[mask] - sp_mean
 
    combined.obsm["X_uce_centered"] = X_centered
    print("Per-species centering applied: subtracted each species' mean "
          "embedding from its cells (stored in .obsm['X_uce_centered']).")
    return combined
 
 
def _make_rank_lookup(order_list, label_for_warning):
    """
    Build a dict mapping each value in order_list to its rank (0, 1, 2, ...).
    Returns a function that looks up any value's rank, falling back to
    (len(order_list), value) for anything not in order_list — i.e. unlisted
    values sort after all listed ones, alphabetically among themselves.
    """
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
    """Compute pairwise distances between (species, cell_type) centroids.
 
    Note: if a species has multiple datasets, all their cells are pooled
    together before averaging — i.e. this collapses dataset as a grouping
    dimension entirely. If you want per-dataset centroids too, use the
    "dataset" column instead of "species" when building combined.obs["group"]
    in load_and_merge().
 
    use_rep: which .obsm key to compute distances from — "X_uce" (raw) or
    "X_uce_centered" (after center_by_species()).
 
    species_order: optional list giving the desired species order (e.g.
    SPECIES_PHYLO_ORDER). Species are ordered primarily by this list,
    secondarily by cell type name. If None, falls back to alphabetical.
    """
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
    dist_df = pd.DataFrame(dist, index=centroids.index, columns=centroids.index)
    return dist_df
 
 
def same_vs_different_type_scores(combined, metric="cosine", max_cells_per_group=500,
                                   use_rep="X_uce"):
    """
    For each cell type present in 2+ datasets, compare:
      - distances between cells of the SAME type across DIFFERENT datasets
      - distances between cells of DIFFERENT types (pooled, across datasets)
    A well-aligned embedding should show same-type-cross-dataset distances
    noticeably smaller than different-type distances.
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
        print("No cell type appears in 2+ datasets — can't compute "
              "same-type-cross-dataset scores. Check your label harmonization.")
        same_type_dists = np.array([])
    else:
        same_type_dists = np.concatenate(same_type_dists)
 
    # different-type pairs: sample broadly across all cells
    idx_all = subsample_idx(pd.Series(True, index=obs.index))
    d_all = pairwise_distances(X[idx_all], metric=metric)
    types_all = obs["cell_type_std"].values[idx_all]
    iu = np.triu_indices_from(d_all, k=1)
    diff_mask = types_all[iu[0]] != types_all[iu[1]]
    diff_type_dists = d_all[iu][diff_mask]
 
    print("\n--- Same-type-vs-different-type separation ---")
    if len(same_type_dists):
        print(f"Same cell type, different dataset: "
              f"mean={same_type_dists.mean():.4f}, n={len(same_type_dists)}")
    print(f"Different cell type (pooled):        "
          f"mean={diff_type_dists.mean():.4f}, n={len(diff_type_dists)}")
    if len(same_type_dists):
        print(f"Separation ratio (diff/same): "
              f"{diff_type_dists.mean() / same_type_dists.mean():.2f}x "
              f"(>1 means the embedding separates cell types better than "
              f"it separates datasets)")
 
    return same_type_dists, diff_type_dists
 
def plot_heatmap(dist_df, out_path, title="UCE centroid distances: species | cell type"):
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
 
 
def reorder_by_celltype(dist_df, cell_type_order=None, species_order=None):
    """
    Reorder a (species | cell_type) distance matrix so cell types are
    grouped together (with species nested within each cell type), instead
    of the default species-first ordering. Groups' labels are expected in
    the "species | cell_type" format produced by centroid_distance_matrix.
 
    cell_type_order: optional list giving the desired cell type order (e.g.
    CELL_TYPE_ONTOLOGY_ORDER). If None, falls back to alphabetical.
    species_order: optional list giving the desired species order WITHIN
    each cell type block (e.g. SPECIES_PHYLO_ORDER). If None, falls back
    to alphabetical.
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
    # flip label order to "cell_type | species" for readability in this view
    relabeled = {lbl: " | ".join(reversed(lbl.split(" | ", 1))) for lbl in new_order}
    reordered = dist_df.loc[new_order, new_order].rename(index=relabeled, columns=relabeled)
    return reordered
 
 
def plot_umap(combined, out_path, use_rep="X_uce"):
    sc.pp.neighbors(combined, use_rep=use_rep)
    sc.tl.umap(combined)
 
    multi_species = combined.obs["species"].nunique() > 1
    n_panels = 3 if multi_species else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 6))
 
    sc.pl.umap(combined, color="dataset", ax=axes[0], show=False, title="Dataset")
    sc.pl.umap(combined, color="cell_type_std", ax=axes[1], show=False,
               title="Cell type", legend_fontsize=6)
    if multi_species:
        sc.pl.umap(combined, color="species", ax=axes[2], show=False, title="Species")
 
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Saved UMAP to {out_path}")
 
 
def main():
    combined = load_and_merge(DATASET_PATHS, CELL_TYPE_COLS, SPECIES)
 
    use_rep = "X_uce"
    if CENTER_BY_SPECIES:
        combined = center_by_species(combined)
        use_rep = "X_uce_centered"
 
    dist_df = centroid_distance_matrix(combined, metric=METRIC, use_rep=use_rep,
                                        species_order=SPECIES_PHYLO_ORDER)
    dist_df.to_csv(f"{OUT_PREFIX}_centroid_distances.csv")
    print(f"Saved centroid distance matrix to {OUT_PREFIX}_centroid_distances.csv")
 
    same_vs_different_type_scores(combined, metric=METRIC, use_rep=use_rep)
 
    title_suffix = " (species-centered)" if CENTER_BY_SPECIES else ""
    plot_heatmap(dist_df, f"{OUT_PREFIX}_heatmap.png",
                 title=f"UCE centroid distances: species (phylogenetic order) | cell type{title_suffix}")
 
    dist_df_by_celltype = reorder_by_celltype(dist_df,
                                               cell_type_order=CELL_TYPE_ONTOLOGY_ORDER,
                                               species_order=SPECIES_PHYLO_ORDER)
    plot_heatmap(dist_df_by_celltype, f"{OUT_PREFIX}_heatmap_by_celltype.png",
                 title=f"UCE centroid distances: cell type (ontology order) | species{title_suffix}")
 
    plot_umap(combined, f"{OUT_PREFIX}_umap.png", use_rep=use_rep)
 
    combined.write_h5ad(f"{OUT_PREFIX}_combined.h5ad")
    print(f"Saved merged AnnData to {OUT_PREFIX}_combined.h5ad")
 
 
if __name__ == "__main__":
    main()
 
