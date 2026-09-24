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
# --- OOM fix (2026-09-24) ---
# A prior run (268,447 cells x 14,512 features x 5 study-layers) OOM'd at
# the "Merging objects" step, right after CCA anchor-finding, even with
# 300GB RAM. Two things were compounding:
#   1. All-pairwise anchor finding: 5 layers -> C(5,2) = 10 pairs, then a
#      full merge of all of them. Reference-based integration (anchor
#      everything to your largest layer) cuts that to 4 pairs and a much
#      cheaper merge.
#   2. ScaleData(so) with no `features=` wasn't guaranteed to stay
#      restricted to variable features once the assay was split into 5
#      layers -- if it silently scaled all 14,512 genes per layer instead
#      of ~2000 variable ones, that's a large dense matrix per layer.
#
# This version defaults to Seurat's documented fix for exactly this
# situation -- sketch-based integration (SketchData + IntegrateLayers on a
# per-layer subsample + ProjectIntegration back onto the full data) -- with
# reference-based, feature-restricted, lower-dims CCA on the FULL data
# available as a fallback (USE_SKETCH_INTEGRATION <- FALSE) if you'd rather
# not add the sketch step. Try the sketch path first; it's the actual
# documented answer for OOM on integration at this scale, not just smaller
# parameter tweaks on the same approach that OOM'd.
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
# (The run that OOM'd used "study", 5 groups -- update this if that's what
# you're integrating over now instead of "species".)
BATCH_KEY <- "study"

# Which methods to run.
RUN_CCA <- TRUE

# Use Seurat's sketch-based integration (recommended for this scale).
# TRUE:  subsample SKETCH_NCELLS representative cells per layer (via
#        leverage-score sampling), integrate on the sketch, then project
#        the full dataset onto that integrated space with
#        ProjectIntegration(). Memory scales with the sketch size, not the
#        full 268K cells.
# FALSE: run CCA directly on the full data (reference-based, feature-
#        restricted, reduced dims -- see below). Try this only if the
#        sketch path is unavailable in your Seurat version, or you've
#        confirmed you have enough memory headroom.
USE_SKETCH_INTEGRATION <- TRUE

# Cells sampled per layer when USE_SKETCH_INTEGRATION is TRUE. 5000 is
# Seurat's own vignette default; raise it if a layer is much smaller than
# this (sampling can't exceed a layer's cell count) or you want a more
# faithful sketch at the cost of more memory.
SKETCH_NCELLS <- 5000

# Anchor-finding parameters, used for both the sketch path and the
# fallback full-data path. Reduced from Seurat's defaults (dims=1:30,
# k.filter=200, max.features=200) to cut anchor-step memory; raise them
# back up if you have room and want to match Seurat's defaults more
# closely.
N_DIMS       <- 20    # dims = 1:N_DIMS
K_FILTER     <- 100
K_ANCHOR     <- 5
K_SCORE      <- 30
MAX_FEATURES <- 100

# Number of variable features to compute/scale on (full-data path only;
# the sketch path uses its own FindVariableFeatures call on the sketch).
N_VAR_FEATURES <- 2000

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
rm(counts); gc()

# Split into per-batch layers (Seurat v5 pattern) so CCA treats each
# batch-key group as a separate dataset to integrate.
so[["RNA"]] <- split(so[["RNA"]], f = so[[BATCH_KEY, drop = TRUE]])

layer_names <- names(so[["RNA"]]@layers)
layer_batches <- unique(sub("^counts\\.", "", layer_names[grepl("^counts\\.", layer_names)]))
layer_sizes <- table(so[[BATCH_KEY, drop = TRUE]])
reference_batch <- names(layer_sizes)[which.max(layer_sizes)]
message(sprintf("Layer sizes: %s", paste(sprintf("%s=%d", names(layer_sizes), layer_sizes), collapse = ", ")))
message(sprintf("Using largest layer as integration reference: '%s'", reference_batch))

# Attach X_uce as a DimReduc so it can be passed as orig.reduction on the
# full-data object (used directly in the fallback path, and re-attached to
# the sketch further down).
so[["uce"]] <- CreateDimReducObject(
  embeddings = emb,
  key = "UCE_",
  assay = DefaultAssay(so)
)

results <- list()

# Reference index: IntegrateLayers' `reference` arg takes the integer
# position(s) of the reference layer(s) in the object's layer order, which
# follows the factor level order `split()` used (alphabetical by default).
batch_levels <- sort(unique(so[[BATCH_KEY, drop = TRUE]]))
reference_idx <- which(batch_levels == reference_batch)

if (RUN_CCA && USE_SKETCH_INTEGRATION) {
  
  message("\n--- Sketch-based CCA integration ---")
  message(sprintf("Sketching up to %d cells per layer (leverage-score sampling)...",
                  SKETCH_NCELLS))
  
  so <- NormalizeData(so)
  so <- FindVariableFeatures(so)
  so <- SketchData(
    object = so,
    ncells = SKETCH_NCELLS,
    method = "LeverageScore",
    sketched.assay = "sketch"
  )
  DefaultAssay(so) <- "sketch"
  
  # Subset X_uce to the sketched cells and attach it as this assay's "uce"
  # reduction, so orig.reduction = "uce" below is computed on the actual
  # sketch, not the full-data embedding.
  sketch_cells <- Cells(so[["sketch"]])
  so[["uce.sketch"]] <- CreateDimReducObject(
    embeddings = emb[sketch_cells, , drop = FALSE],
    key = "UCESK_",
    assay = "sketch"
  )
  
  so <- FindVariableFeatures(so, nfeatures = N_VAR_FEATURES)
  so <- ScaleData(so, features = VariableFeatures(so))
  so <- RunPCA(so, npcs = min(30, N_DIMS + 5), verbose = FALSE)
  gc()
  
  so <- IntegrateLayers(
    object = so,
    method = CCAIntegration,
    orig.reduction = "uce.sketch",
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
  
  message("\n--- Running CCAIntegration on full data (reference-based, feature-restricted) ---")
  
  so <- NormalizeData(so)
  so <- FindVariableFeatures(so, nfeatures = N_VAR_FEATURES)
  so <- ScaleData(so, features = VariableFeatures(so))
  gc()
  
  # Seurat still wants a "pca" reduction to exist for some internal
  # bookkeeping; it isn't used as the anchor space below since we pass
  # orig.reduction = "uce" explicitly.
  so <- RunPCA(so, npcs = min(30, N_DIMS + 5), verbose = FALSE)
  
  so <- IntegrateLayers(
    object = so,
    method = CCAIntegration,
    orig.reduction = "uce",
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