# integrate_cca.R
#
# Integrates the UCE embedding using Seurat's CCA method.
#
# --- Design note: UCE embedding AS the feature space (2026-09-24) ---
# Earlier versions of this script built the Seurat object from raw gene
# counts and bolted X_uce on as a separate `orig.reduction` via
# CreateDimReducObject(embeddings = ..., key = ...) with no loadings. That
# ran fine through anchor-finding and CCA, but failed at the final
# IntegrateEmbeddings step with:
#   Error in Loadings(object = reductions)[, dims.to.integrate] :
#     subscript out of bounds
# Seurat's final embedding-integration step pulls Loadings() off
# orig.reduction to project the correction back into the original space.
# A DimReducObject built from embeddings alone has no loadings (X_uce has
# no gene-level loadings -- it's a neural embedding, not a linear
# decomposition), so Loadings() returned an empty matrix and the dims
# subscript went out of bounds. (This corrects something said earlier in
# this project: it's not that CCA "doesn't need loadings" while RPCA does
# -- CCA's own anchor-finding doesn't, but Seurat's IntegrateEmbeddings
# step needs loadings for either method.)
#
# Fix: instead of gene counts + a bolted-on loadings-free reduction, build
# the Seurat object directly FROM the UCE embedding -- i.e. treat the 1280
# UCE dimensions as the "features" instead of genes. RunPCA on that then
# produces a real, loadings-bearing PCA (1280 x npcs), which
# IntegrateLayers/IntegrateEmbeddings can index into normally. This is the
# standard way to run Seurat's CCA/RPCA integration on any precomputed
# embedding that isn't itself a PCA. It also means raw gene counts are no
# longer needed at all for this script -- a large additional memory win on
# top of the sketch-based fix (268K cells x 1280 "features" instead of
# 268K x 14,512 genes).
#
# --- OOM fix (2026-09-24), still in effect ---
# Sketch-based integration (SketchData + IntegrateLayers on a per-layer
# subsample + ProjectIntegration back onto the full data), reference-based
# anchoring, and reduced dims/k.filter/etc. This part already ran cleanly
# through anchor-finding, CCA, and merging on the last run -- no changes
# needed here, kept as-is.
#
# Requires: Seurat (>=5.1, for SketchData/ProjectIntegration), SeuratObject
#   install.packages(c("Seurat", "SeuratObject"))
#
# Run: Rscript integrate_cca.R
.libPaths(c("~/R/R4.6.0_packages"))
options(repos = c(CRAN = "https://cloud.r-project.org")) 
library(Seurat)
library(SeuratObject)

# ============================== CONFIG =======================================

IN_DIR  <- "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/seurat_export"
OUT_DIR <- "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/seurat_export"

# Which .obs column (written by export_uce_for_seurat.py) to split Seurat
# layers by. Must match the batch_key printed at the end of that script.
BATCH_KEY <- "study"

# Which methods to run.
RUN_CCA <- TRUE

# Use Seurat's sketch-based integration (recommended for this scale, and
# already confirmed to get through anchor-finding/CCA/merging without OOM
# on the last run).
USE_SKETCH_INTEGRATION <- TRUE
SKETCH_NCELLS <- 5000

# Anchor-finding parameters.
N_DIMS       <- 20    # dims = 1:N_DIMS -- must be <= npcs below
K_FILTER     <- 100
K_ANCHOR     <- 5
K_SCORE      <- 30
MAX_FEATURES <- 100

# Number of PCs to compute on the UCE-embedding-as-features matrix. Must be
# >= N_DIMS (this is exactly what went out of bounds last time -- keep a
# buffer above N_DIMS here).
N_PCS <- N_DIMS + 10

# ==============================================================================

message("Reading exported embedding / metadata from ", IN_DIR)

# Only X_uce.csv and obs_meta.csv are needed now -- counts.mtx/features.tsv
# from export_uce_for_seurat.py are no longer used by this script (the UCE
# embedding itself is the feature space; see note at top of file).
emb <- read.csv(file.path(IN_DIR, "X_uce.csv"), row.names = 1,
                check.names = FALSE)
barcodes <- rownames(emb)

meta <- read.csv(file.path(IN_DIR, "obs_meta.csv"), row.names = 1,
                 check.names = FALSE)
meta <- meta[barcodes, , drop = FALSE]

if (!(BATCH_KEY %in% colnames(meta))) {
  stop(sprintf("BATCH_KEY '%s' not found in obs_meta.csv columns: %s",
               BATCH_KEY, paste(colnames(meta), collapse = ", ")))
}

# Transpose to features (UCE dims) x cells, as Seurat assays expect.
uce_mat <- t(as.matrix(emb))
rm(emb); gc()

message(sprintf("Loaded %d cells, %d UCE dims. Batch key '%s' has %d groups.",
                ncol(uce_mat), nrow(uce_mat), BATCH_KEY,
                length(unique(meta[[BATCH_KEY]]))))

# --- Build the Seurat object directly from the UCE embedding ---
# No NormalizeData: these are continuous embedding values (can be
# negative), not counts, so log-normalization is neither valid nor
# meaningful here. `data` is set equal to `counts` unchanged.
so <- CreateSeuratObject(counts = uce_mat, meta.data = meta, assay = "RNA")
so[["RNA"]]$data <- so[["RNA"]]$counts
rm(uce_mat); gc()

# Split into per-batch layers (Seurat v5 pattern) so CCA treats each
# batch-key group as a separate dataset to integrate.
so[["RNA"]] <- split(so[["RNA"]], f = so[[BATCH_KEY, drop = TRUE]])

layer_sizes <- table(so[[BATCH_KEY, drop = TRUE]])
reference_batch <- names(layer_sizes)[which.max(layer_sizes)]
message(sprintf("Layer sizes: %s", paste(sprintf("%s=%d", names(layer_sizes), layer_sizes), collapse = ", ")))
message(sprintf("Using largest layer as integration reference: '%s'", reference_batch))

# Reference index: IntegrateLayers' `reference` arg takes the integer
# position(s) of the reference layer(s) in the object's layer order, which
# follows the factor level order split() used (alphabetical by default).
batch_levels <- sort(unique(so[[BATCH_KEY, drop = TRUE]]))
reference_idx <- which(batch_levels == reference_batch)

results <- list()

if (RUN_CCA && USE_SKETCH_INTEGRATION) {
  
  message("\n--- Sketch-based CCA integration (UCE dims as features) ---")
  message(sprintf("Sketching up to %d cells per layer (leverage-score sampling)...",
                  SKETCH_NCELLS))
  
  so <- FindVariableFeatures(so, nfeatures = nrow(so))  # all UCE dims
  so <- SketchData(
    object = so,
    ncells = SKETCH_NCELLS,
    method = "LeverageScore",
    sketched.assay = "sketch"
  )
  DefaultAssay(so) <- "sketch"
  
  so <- FindVariableFeatures(so, nfeatures = nrow(so))
  so <- ScaleData(so, features = VariableFeatures(so))
  so <- RunPCA(so, features = VariableFeatures(so), npcs = N_PCS, verbose = FALSE)
  gc()
  
  so <- IntegrateLayers(
    object = so,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca.sketch",
    reference = reference_idx,
    dims = 1:N_DIMS,
    k.anchor = K_ANCHOR,
    k.filter = K_FILTER,
    k.score = K_SCORE,
    max.features = MAX_FEATURES,
    verbose = TRUE
  )
  gc()
  
  message("Projecting integrated sketch back onto the full dataset...")
  so <- ProjectIntegration(
    object = so,
    sketched.assay = "sketch",
    assay = "RNA",
    reduction = "integrated.cca.sketch",
    reduction.name = "integrated.cca.full"
  )
  DefaultAssay(so) <- "RNA"
  
  results[["cca"]] <- Embeddings(so, "integrated.cca.full")
  message("Sketch-based CCA integration complete: ", nrow(results[["cca"]]),
          " cells x ", ncol(results[["cca"]]), " dims.")
  
} else if (RUN_CCA) {
  
  message("\n--- Running CCAIntegration on full data (UCE dims as features) ---")
  
  so <- FindVariableFeatures(so, nfeatures = nrow(so))
  so <- ScaleData(so, features = VariableFeatures(so))
  gc()
  so <- RunPCA(so, features = VariableFeatures(so), npcs = N_PCS, verbose = FALSE)
  
  so <- IntegrateLayers(
    object = so,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    reference = reference_idx,
    dims = 1:N_DIMS,
    k.anchor = K_ANCHOR,
    k.filter = K_FILTER,
    k.score = K_SCORE,
    max.features = MAX_FEATURES,
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