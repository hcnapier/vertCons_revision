# integrate_cca.R
#
# Integrates UCE embeddings using Seurat's CCA method, with X_uce supplied
# as the precomputed `orig.reduction` (instead of Seurat computing its own
# PCA from expression). CCAIntegration computes canonical correlation
# directly on the embeddings and doesn't need gene loadings, which is why
# it's used here rather than RPCA (RPCA's reciprocal-projection step
# relies on linear PCA loadings that a neural embedding like X_uce doesn't
# have).
#
# Requires: Seurat (>=5.0), SeuratObject, Matrix
#   install.packages(c("Seurat", "SeuratObject", "Matrix"))
#
# Run: Rscript integrate_cca.R

.libPaths(c("~/R/R4.6.0_packages"))
options(repos = c(CRAN = "https://cloud.r-project.org")) 
library(Seurat)
library(SeuratObject)
library(Matrix)

# ============================== CONFIG =======================================

IN_DIR  <- "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/seurat_export"
OUT_DIR <- "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/seurat_export"

# Which .obs column (written by export_uce_for_seurat.py) to split Seurat
# layers by. Must match the batch_key printed at the end of that script.
BATCH_KEY <- "study"

# Which methods to run.
RUN_CCA <- TRUE

# ==============================================================================

message("Reading exported counts / embedding / metadata from ", IN_DIR)

counts <- Matrix::readMM(file.path(IN_DIR, "counts.mtx"))
barcodes <- readLines(file.path(IN_DIR, "barcodes.tsv"))
features <- readLines(file.path(IN_DIR, "features.tsv"))
rownames(counts) <- make.unique(features)
colnames(counts) <- barcodes

meta <- read.csv(file.path(IN_DIR, "obs_meta.csv"), row.names = 1,
                 check.names = FALSE)
meta <- meta[barcodes, , drop = FALSE]

emb <- read.csv(file.path(IN_DIR, "X_uce.csv"), row.names = 1,
                check.names = FALSE)
emb <- as.matrix(emb[barcodes, , drop = FALSE])

if (!(BATCH_KEY %in% colnames(meta))) {
  stop(sprintf("BATCH_KEY '%s' not found in obs_meta.csv columns: %s",
               BATCH_KEY, paste(colnames(meta), collapse = ", ")))
}

message(sprintf("Loaded %d cells, %d features, %d embedding dims. ",
                ncol(counts), nrow(counts), ncol(emb)),
        sprintf("Batch key '%s' has %d groups.",
                BATCH_KEY, length(unique(meta[[BATCH_KEY]]))))

# --- Build the Seurat object ---
so <- CreateSeuratObject(counts = counts, meta.data = meta)

# Split into per-batch layers (Seurat v5 pattern) so CCA treats each
# batch-key group as a separate dataset to integrate.
so[["RNA"]] <- split(so[["RNA"]], f = so[[BATCH_KEY, drop = TRUE]])

# Standard preprocessing. Even though anchor-finding will use the UCE
# reduction (not a fresh PCA), CCAIntegration still expects
# normalized/scaled data to exist per layer for its internal steps.
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)

# Attach X_uce as a DimReduc so it can be passed as orig.reduction. Seurat
# still wants a "pca" reduction to exist for some internal bookkeeping, so
# compute a real (small) PCA too -- it won't be used as the anchor space
# since we explicitly pass orig.reduction = "uce" below.
so <- RunPCA(so, npcs = min(30, ncol(emb)), verbose = FALSE)
so[["uce"]] <- CreateDimReducObject(
  embeddings = emb,
  key = "UCE_",
  assay = DefaultAssay(so)
)

results <- list()

if (RUN_CCA) {
  message("\n--- Running CCAIntegration (orig.reduction = 'uce') ---")
  so <- IntegrateLayers(
    object = so,
    method = CCAIntegration,
    orig.reduction = "uce",
    new.reduction = "integrated.cca",
    verbose = TRUE
  )
  results[["cca"]] <- Embeddings(so, "integrated.cca")
  message("CCA integration complete: ", nrow(results[["cca"]]), " cells x ",
          ncol(results[["cca"]]), " dims.")
}

# --- Write results back out for Python to re-import ---
for (name in names(results)) {
  out_path <- file.path(OUT_DIR, paste0("X_uce_", name, ".csv"))
  write.csv(results[[name]], out_path)
  message("Wrote ", out_path)
}

message("\nDone. Run import_integrated_embeddings.py next to merge these ",
        "back into your AnnData's .obsm.")