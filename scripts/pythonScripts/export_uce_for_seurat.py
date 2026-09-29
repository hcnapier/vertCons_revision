"""
Export a combined, UCE-embedded AnnData for RPCA/CCA integration in Seurat.

Writes three things R will read:
  1. X_uce.csv        - the UCE embedding (cells x dims), used as Seurat's
                         `orig.reduction` for anchor finding.
  2. obs_meta.csv      - per-cell metadata (species, study, technology,
                         cell_type_std, dataset, group, and whatever
                         BATCH_KEY you integrate over).
  3. counts.mtx/.tsv   - the raw counts matrix, so Seurat has a real assay
                         to attach the embedding to. Only used to give the
                         Seurat object valid structure / support CCA's own
                         internal steps -- NOT used as the anchor space
                         (that's X_uce, via orig.reduction).

python export_uce_for_seurat.py
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

# ============================== CONFIG ======================================

# Path to the combined, embedded AnnData (must have .obsm["X_uce"] and the
# usual .obs columns: species, study, technology, cell_type_std, dataset,
# group).
PATH = "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances/center_combined.h5ad"

# Which .obsm key to export as the integration reduction.
USE_REP = "X_uce_centered"

# Which .obs column to split Seurat "layers" by (one layer per unique
# value = the unit RPCA/CCA integrates over). Use the same composite-key
# trick as USE_COMBINED_BATCH_KEY in the Harmony script if you want to
# integrate over more than one variable at once: build that column in
# Python before exporting (see build_combined_batch_column below).
BATCH_KEY = "study"

# If you want to integrate over multiple .obs columns jointly (mirroring
# USE_COMBINED_BATCH_KEY=True in the Harmony script), list them here and
# leave BATCH_KEY = None; a combined column will be built and used instead.
COMBINE_KEYS = None  # e.g. ["species", "study", "technology"]

# Where to write the exported files.
OUT_DIR = "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/seurat_export"

# =============================================================================


def build_combined_batch_column(obs, keys):
    col_name = "_".join(keys) + "_combined"
    obs[col_name] = obs[keys].astype(str).agg("_".join, axis=1)
    return obs, col_name


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    combined = sc.read_h5ad(PATH)

    if USE_REP not in combined.obsm:
        raise ValueError(f"'{USE_REP}' not found in .obsm. Available: "
                          f"{list(combined.obsm.keys())}")

    obs = combined.obs.copy()
    batch_key = BATCH_KEY
    if COMBINE_KEYS:
        obs, batch_key = build_combined_batch_column(obs, COMBINE_KEYS)
        print(f"Combined {COMBINE_KEYS} into '{batch_key}' "
              f"({obs[batch_key].nunique()} unique values).")

    if batch_key not in obs:
        raise ValueError(f"batch_key '{batch_key}' not found in .obs. "
                          f"Available: {list(obs.columns)}")

    # 1. X_uce embedding
    X = combined.obsm[USE_REP]
    dim_cols = [f"UCE_{i+1}" for i in range(X.shape[1])]
    emb_df = pd.DataFrame(X, index=combined.obs_names, columns=dim_cols)
    emb_path = os.path.join(OUT_DIR, "X_uce.csv")
    emb_df.to_csv(emb_path)
    print(f"Wrote embedding ({emb_df.shape[0]} cells x {emb_df.shape[1]} dims) "
          f"to {emb_path}")

    # 2. obs metadata
    meta_cols = [c for c in ["species", "study", "technology", "cell_type_std",
                              "dataset", "group", batch_key] if c in obs]
    meta_df = obs[meta_cols].copy()
    meta_df.index = combined.obs_names
    meta_path = os.path.join(OUT_DIR, "obs_meta.csv")
    meta_df.to_csv(meta_path)
    print(f"Wrote obs metadata to {meta_path} (batch_key column: '{batch_key}')")

    # 3. counts matrix (for Seurat object structure; NOT the anchor space)
    Xc = combined.X
    if not sp.issparse(Xc):
        Xc = sp.csr_matrix(Xc)
    counts_path = os.path.join(OUT_DIR, "counts.mtx")
    sio.mmwrite(counts_path, Xc.T)  # mtx convention: features x cells
    pd.Series(combined.obs_names).to_csv(
        os.path.join(OUT_DIR, "barcodes.tsv"), index=False, header=False)
    var_names = combined.var_names if combined.var_names is not None else \
        pd.Index([f"feature_{i}" for i in range(Xc.shape[1])])
    pd.Series(var_names).to_csv(
        os.path.join(OUT_DIR, "features.tsv"), index=False, header=False)
    print(f"Wrote counts matrix ({Xc.shape[0]} cells x {Xc.shape[1]} features) "
          f"to {counts_path} (+ barcodes.tsv, features.tsv)")

    print(f"\nBatch key for Seurat layers: '{batch_key}' "
          f"({obs[batch_key].nunique()} groups)")
    print("Done. Run integrate_rpca_cca.R next, pointing IN_DIR at:", OUT_DIR)


if __name__ == "__main__":
    main()
