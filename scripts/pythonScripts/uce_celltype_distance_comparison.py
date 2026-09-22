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
  5. Run a study-controlled diagnostic that isolates genuine species/
     cell-type separation from study/technology batch effects, using the
     subset of datasets that share a single study (see
     study_controlled_diagnostic()).
  6. Save a distance heatmap and a UMAP colored by dataset / cell type.

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

# Paths to UCE-embedded .h5ad files (output of eval_single_anndata /
# uce-eval-single-anndata — each should already have .obsm["X_uce"]).
#
# NOTE: previously "dataset3" was accidentally used as the key for both
# Capra_hircus and Cavia_porcellus, which silently dropped the goat entry
# (a later duplicate dict key just overwrites the earlier one in Python)
# and shifted every dataset after it out of alignment with SPECIES/STUDY/
# TECHNOLOGY below. Re-keyed here as dataset1..dataset10 to match those
# dicts one-to-one — double check these paths are actually correct on disk.
DATASET_PATHS = {
    "dataset1": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Bos_taurus_uce_adata.h5ad",
    "dataset2": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Canis_lupus_familiaris_uce_adata.h5ad",
    "dataset3": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Capra_hircus_uce_adata.h5ad",
    "dataset4": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Cavia_porcellus_uce_adata.h5ad",
    "dataset5": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Homo_sapiens_uce_adata.h5ad",
    "dataset6": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Macaca_fascicularis_uce_adata.h5ad",
    "dataset7": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Mus_musculus_uce_adata.h5ad",
    "dataset8": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Oryctolagus_cuniculus_uce_adata.h5ad",
    "dataset9": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Rattus_norvegicus_shrutx_uce_adata.h5ad",
    "dataset10": "/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded/Sus_scrofa_uce_adata.h5ad",
}

# Name of the .obs column holding cell type labels in EACH dataset.
CELL_TYPE_COLS = {
    "dataset1": "napierCellTypes",
    "dataset2": "napierCellTypes",
    "dataset3": "napierCellTypes",
    "dataset4": "napierCellTypes",
    "dataset5": "napierCellTypes",
    "dataset6": "napierCellTypes",
    "dataset7": "napierCellTypes",
    "dataset8": "napierCellTypes",
    "dataset9": "napierCellTypes",
    "dataset10": "napierCellTypes",
}

# Species of each dataset (used only for bookkeeping/reporting here — UCE
# tokenizes genes via ESM2 protein embeddings rather than gene symbols, so
# it aligns homologous genes across species internally. No ortholog
# filtering needed; just make sure you ran UCE itself with the correct
# --species flag for each dataset when you generated the embeddings).
SPECIES = {
    "dataset1": "cow",
    "dataset2": "dog",
    "dataset3": "goat",
    "dataset4": "guineaPig",
    "dataset5": "human",
    "dataset6": "macaque",
    "dataset7": "mouse",
    "dataset8": "rabbit",
    "dataset9": "rat",
    "dataset10": "pig",
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
    "rat",
    "mouse",
    "rabbit",
    "pig", 
    "cow", 
    "goat",
    "dog"
]

# Study/publication each dataset came from
STUDY = {
    "dataset1": "Tan",
    "dataset2": "Tan",
    "dataset3": "Tan",
    "dataset4": "Tan", 
    "dataset5": "Tsang",
    "dataset6": "Wang",
    "dataset7": "Jiang",
    "dataset8": "Tan",
    "dataset9": "Iqbal",
    "dataset10": "Tan"
}
 
# Sequencing/profiling technology used for each dataset
TECHNOLOGY = {
    "dataset1": "BGISEQ",
    "dataset2": "BGISEQ",
    "dataset3": "BGISEQ",
    "dataset4": "BGISEQ",
    "dataset5": "Illumina",
    "dataset6": "Illumina",
    "dataset7": "Illumina",
    "dataset8": "BGISEQ",
    "dataset9": "Illumina",
    "dataset10": "BGISEQ"
}
 
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
    "spt",
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
    "nkcells",
    "mono",
    "mac",
    "dc",
    "neu"
]

# Distance metric for comparisons: "euclidean" or "cosine"
METRIC = "cosine"

# If True, run Harmony (via scanpy's harmony_integrate, falling back to
# harmonypy directly if that's unavailable) on the UCE embeddings, using
# HARMONY_BATCH_KEY as the batch variable to integrate over. This is a
# heavier, more standard batch-correction approach than CENTER_BY_SPECIES —
# it can capture non-linear/per-cell-type-specific batch effects, not just
# a global per-species shift. If both this and CENTER_BY_SPECIES are True,
# Harmony takes priority and centering is skipped (a warning is printed),
# since running both is usually redundant.
USE_HARMONY = True
 
# .obs column(s) Harmony integrates over. Can be a single string ("species")
# or a list to correct for multiple batch effects at once, e.g.
# ["species", "study", "technology"] — Harmony supports multiple batch
# variables simultaneously. Available columns after load_and_merge are:
# "dataset", "species", "study", "technology".
HARMONY_BATCH_KEY = ["species", "study", "technology"]

# Only relevant when HARMONY_BATCH_KEY is a list with 2+ entries. Controls
# HOW multiple keys get combined:
#   False (default): pass all keys to Harmony as separate covariates
#     (harmonypy's vars_use / scanpy's multi-key support). Harmony corrects
#     for each variable's effect roughly independently/additively.
#   True: combine the listed keys into a single composite batch column
#     first (e.g. species "human" + study "study1" -> "human_study1"), and
#     integrate on that one combined column instead. This treats every
#     unique combination as its own distinct batch, which captures
#     INTERACTION effects between the variables (e.g. a species x study
#     combination that behaves unusually together) that treating them as
#     separate additive covariates would miss. Generally the more thorough
#     option when you suspect the batch effects aren't independent, at the
#     cost of more, smaller batches for Harmony to work with (which can
#     hurt correction quality if any combination has very few cells).
USE_COMBINED_BATCH_KEY = False

# If True, mean-center each species' embeddings (subtract that species'
# overall mean X_uce vector from every one of its cells) before computing
# distances, UMAP, and the separation score.
CENTER_BY_SPECIES = False

# --- Study-controlled diagnostic ---
# Which .obs["study"] value to use as the "study-controlled" subset for
# study_controlled_diagnostic() — i.e. a study that profiled multiple
# species together (so study/technology are held constant and any
# remaining same-vs-different-type separation reflects genuine species/
# cell-type biology, not a study batch effect). If None, the script
# auto-picks the study with the most distinct species (ties broken
# alphabetically) and prints what it picked.
DIAGNOSTIC_STUDY = "Tan"

# Where to save outputs
OUT_PREFIX = "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances"

# =============================================================================

def load_and_merge(dataset_paths, cell_type_cols, species_map, study_map, technology_map):
    """Load each embedded dataset, standardize cell type column name, tag
    dataset origin, species, study, and technology, and concatenate into
    one AnnData."""
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
        a.obs["study"] = study_map.get(name, "unknown")
        a.obs["technology"] = technology_map.get(name, "unknown")
        # keep only what we need to avoid var mismatch issues on concat
        a = ad.AnnData(X=a.X, obs=a.obs[["cell_type_std", "dataset", "species",
                                          "study", "technology"]].copy(),
                        obsm={"X_uce": a.obsm["X_uce"]})
        adatas.append(a)
 
    combined = ad.concat(adatas, join="outer", label="batch", index_unique="-")
    combined.obs["group"] = (combined.obs["species"].astype(str) + " | " +
                              combined.obs["cell_type_std"].astype(str))
    n_species = combined.obs["species"].nunique()
    print(f"Loaded {combined.n_obs} cells across {len(dataset_paths)} datasets "
          f"({n_species} species, "
          f"{combined.obs['study'].nunique()} studies, "
          f"{combined.obs['technology'].nunique()} technologies), "
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
  
def run_harmony_integration(combined, batch_key="species", use_rep="X_uce",
                             combine_keys=False):
    """
    Run Harmony on the UCE embeddings directly (not on a PCA reduction —
    X_uce is already a compact learned representation, so Harmony is run
    on it as-is). Calls harmonypy directly rather than
    scanpy.external.pp.harmony_integrate, since that wrapper's unconditional
    output transpose is broken for current harmonypy versions (see note in
    the function body). Stores the result in .obsm["X_uce_harmony"].
 
    batch_key: a single .obs column name, or a list of them.
    combine_keys: only relevant if batch_key is a list with 2+ entries.
        If True, the listed columns are combined into one composite column
        (e.g. species + study -> "human_study1") and Harmony integrates on
        that single column, treating every unique combination as its own
        batch (captures interaction effects). If False, all keys are
        passed to Harmony as separate covariates instead.
 
    Requires the 'harmonypy' package: pip install harmonypy
    """
    batch_keys = [batch_key] if isinstance(batch_key, str) else list(batch_key)
    missing = [k for k in batch_keys if k not in combined.obs]
    if missing:
        raise ValueError(f"HARMONY_BATCH_KEY column(s) not found in .obs: "
                          f"{missing} (available: {list(combined.obs.columns)})")
 
    if combine_keys and len(batch_keys) > 1:
        combined_col = "_".join(batch_keys) + "_combined"
        combined.obs[combined_col] = combined.obs[batch_keys].astype(str).agg("_".join, axis=1)
        n_combos = combined.obs[combined_col].nunique()
        print(f"Combined {batch_keys} into '{combined_col}' "
              f"({n_combos} unique combinations). Integrating on this "
              f"single composite column.")
        small_combos = combined.obs[combined_col].value_counts()
        small_combos = small_combos[small_combos < 10]
        if len(small_combos):
            print(f"WARNING: {len(small_combos)} combination(s) have fewer "
                  f"than 10 cells — Harmony correction quality for those "
                  f"groups may be poor:\n{small_combos}")
        harmony_keys = [combined_col]
    else:
        harmony_keys = batch_keys
 
    # NOTE: we deliberately do NOT use scanpy.external.pp.harmony_integrate
    # here. That wrapper unconditionally transposes harmonypy's output
    # (Z_corr.T), which was correct for harmonypy <0.1.0 (dims x cells) but
    # is WRONG for harmonypy >=0.1.0, which already returns Z_corr as
    # (cells x dims) — see https://github.com/scverse/scanpy/issues/3940.
    # Calling harmonypy directly and checking the actual shape avoids
    # depending on which harmonypy version happens to be installed.
    try:
        import harmonypy
    except ImportError:
        raise ImportError(
            "Harmony integration requires the 'harmonypy' package. "
            "Install it with: pip install harmonypy"
        )
 
    ho = harmonypy.run_harmony(combined.obsm[use_rep], combined.obs, harmony_keys)
    Z = np.asarray(ho.Z_corr)
    n_cells = combined.n_obs
 
    if Z.shape[0] == n_cells:
        corrected = Z  # already (cells x dims) — harmonypy >=0.1.0
    elif Z.shape[1] == n_cells:
        corrected = Z.T  # (dims x cells) — older harmonypy, needs transpose
        print("Detected older harmonypy output orientation (dims x cells) — transposed.")
    else:
        raise ValueError(
            f"harmonypy output shape {Z.shape} doesn't match n_cells "
            f"({n_cells}) on either axis — can't determine correct "
            f"orientation. This may indicate a harmonypy version with a "
            f"different output format than expected; inspect ho.Z_corr "
            f"manually."
        )
    combined.obsm["X_uce_harmony"] = corrected
    print(f"Harmony integration complete via harmonypy "
          f"(batch_key(s)={harmony_keys}). Stored in .obsm['X_uce_harmony'] "
          f"with shape {corrected.shape}.")
 
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
                                   use_rep="X_uce", verbose=True):
    """
    For each cell type present in 2+ datasets, compare:
      - distances between cells of the SAME type across DIFFERENT datasets
      - distances between cells of DIFFERENT types (pooled, across datasets)
    A well-aligned embedding should show same-type-cross-dataset distances
    noticeably smaller than different-type distances.

    Returns (same_type_dists, diff_type_dists, ratio) where ratio is
    diff_type_dists.mean() / same_type_dists.mean(), or None if there
    weren't any same-type-cross-dataset pairs to compute it from.
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
 
    ratio = None
    if verbose:
        print("\n--- Same-type-vs-different-type separation ---")
    if len(same_type_dists):
        if verbose:
            print(f"Same cell type, different dataset: "
                  f"mean={same_type_dists.mean():.4f}, n={len(same_type_dists)}")
    if verbose:
        print(f"Different cell type (pooled):        "
              f"mean={diff_type_dists.mean():.4f}, n={len(diff_type_dists)}")
    if len(same_type_dists):
        ratio = diff_type_dists.mean() / same_type_dists.mean()
        if verbose:
            print(f"Separation ratio (diff/same): "
                  f"{ratio:.2f}x "
                  f"(>1 means the embedding separates cell types better than "
                  f"it separates datasets)")
 
    return same_type_dists, diff_type_dists, ratio


def study_controlled_diagnostic(combined, study_name=None, metric="cosine",
                                 max_cells_per_group=500, reps_to_compare=None):
    """
    Isolate genuine species/cell-type separation from study/technology
    batch effects.

    The same_vs_different_type_scores() ratio computed on the FULL dataset
    is ambiguous whenever species and study are confounded (i.e. whenever
    a species only appears in one study): "same cell type, different
    dataset" pairs for those species are simultaneously "different study"
    pairs, so you can't tell whether the embedding is separating species
    biology or just a study/technology batch effect.

    This function restricts the same computation to the subset of datasets
    sharing a single `study` value that itself contains multiple species
    (e.g. several species profiled together in one paper, one technology).
    Within that subset, study and technology are held constant, so "same
    cell type, different dataset" pairs are necessarily cross-species pairs
    that share study/technology — any remaining separation there reflects
    real species/cell-type biology rather than a study artifact.

    study_name: which combined.obs["study"] value defines the
        study-controlled subset. If None, auto-picks the study with the
        most distinct species among its datasets (ties broken
        alphabetically) and prints what it picked. Raises if the chosen
        study only has one species (nothing to isolate).
    reps_to_compare: list of .obsm keys to compute the diagnostic for. If
        None, uses every rep in combined.obsm that is one of
        ["X_uce", "X_uce_centered", "X_uce_harmony"] and is actually
        present, in that order.

    Returns a pandas DataFrame with one row per rep, comparing the
    full-dataset ratio to the study-controlled ratio.
    """
    if reps_to_compare is None:
        reps_to_compare = [r for r in ["X_uce", "X_uce_centered", "X_uce_harmony"]
                            if r in combined.obsm]

    if study_name is None:
        species_per_study = combined.obs.groupby("study")["species"].nunique()
        candidates = species_per_study[species_per_study > 1]
        if candidates.empty:
            raise ValueError(
                "No study in your data profiles more than one species — "
                "can't build a study-controlled subset. Set DIAGNOSTIC_STUDY "
                "explicitly if this is wrong, or skip this diagnostic."
            )
        # ties broken alphabetically by study name
        top_studies = candidates[candidates == candidates.max()].sort_index()
        study_name = top_studies.index[0]
        print(f"DIAGNOSTIC_STUDY not set — auto-picked study '{study_name}' "
              f"({candidates[study_name]} species) as the study-controlled subset.")

    if study_name not in set(combined.obs["study"]):
        raise ValueError(f"study_name '{study_name}' not found in .obs['study']. "
                          f"Available: {sorted(combined.obs['study'].unique())}")

    sub = combined[combined.obs["study"] == study_name].copy()
    n_species_sub = sub.obs["species"].nunique()
    n_datasets_sub = sub.obs["dataset"].nunique()
    species_list = sorted(sub.obs["species"].unique())

    print(f"\n=== Study-controlled diagnostic: study = '{study_name}' ===")
    print(f"{n_datasets_sub} dataset(s), {n_species_sub} species: {species_list}")
    print(f"{sub.obs['technology'].nunique()} technology value(s) in this subset: "
          f"{sorted(sub.obs['technology'].unique())}")
    if n_species_sub < 2:
        raise ValueError(
            f"Study '{study_name}' only has {n_species_sub} species in this "
            f"data — there's nothing to isolate species separation from "
            f"here. Pick a different DIAGNOSTIC_STUDY."
        )
    per_species_datasets = sub.obs.groupby("species")["dataset"].nunique()
    if (per_species_datasets > 1).any():
        print("NOTE: some species in this subset have >1 dataset, so "
              "'same cell type, different dataset' pairs within the subset "
              "aren't guaranteed to be cross-species — interpret with that "
              "in mind:\n" + per_species_datasets.to_string())
    else:
        print("Each species in this subset has exactly one dataset, so "
              "'same cell type, different dataset' pairs below are "
              "necessarily cross-species, same-study, same-technology "
              "pairs — a clean read on species/cell-type separation.")

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
    print(
        "\nInterpretation: if the study-controlled ratio is close to the "
        "full-dataset ratio, most of the apparent cell-type separation "
        "reflects genuine species/cell-type biology rather than a study "
        "batch effect. If the study-controlled ratio collapses toward 1.0 "
        "(or is much lower than the full-dataset ratio), a substantial "
        "part of the full-dataset separation was likely driven by study/"
        "technology confounds rather than biology."
    )
    return result

 
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
    combined = load_and_merge(DATASET_PATHS, CELL_TYPE_COLS, SPECIES, STUDY, TECHNOLOGY)
 
    use_rep = "X_uce"
    if USE_HARMONY:
        if CENTER_BY_SPECIES:
            print("Both USE_HARMONY and CENTER_BY_SPECIES are True — "
                  "running Harmony only and skipping centering, since "
                  "combining both is usually redundant.")
        combined = run_harmony_integration(combined, batch_key=HARMONY_BATCH_KEY,
                                            use_rep=use_rep,
                                            combine_keys=USE_COMBINED_BATCH_KEY)
        use_rep = "X_uce_harmony"
    elif CENTER_BY_SPECIES:
        combined = center_by_species(combined)
        use_rep = "X_uce_centered"
 
    dist_df = centroid_distance_matrix(combined, metric=METRIC, use_rep=use_rep,
                                        species_order=SPECIES_PHYLO_ORDER)
    dist_df.to_csv(f"{OUT_PREFIX}_centroid_distances.csv")
    print(f"Saved centroid distance matrix to {OUT_PREFIX}_centroid_distances.csv")
 
    same_vs_different_type_scores(combined, metric=METRIC, use_rep=use_rep)

    diagnostic_df = study_controlled_diagnostic(
        combined, study_name=DIAGNOSTIC_STUDY, metric=METRIC)
    diagnostic_df.to_csv(f"{OUT_PREFIX}_study_controlled_diagnostic.csv", index=False)
    print(f"Saved study-controlled diagnostic to "
          f"{OUT_PREFIX}_study_controlled_diagnostic.csv")
 
    if USE_HARMONY:
        title_suffix = " (Harmony-integrated)"
    elif CENTER_BY_SPECIES:
        title_suffix = " (species-centered)"
    else:
        title_suffix = ""
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
