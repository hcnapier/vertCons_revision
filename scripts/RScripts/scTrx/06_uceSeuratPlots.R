.libPaths(c("~/R/R4.6.0_packages"))
options(repos = c(CRAN = "https://cloud.r-project.org")) 
library(Seurat)
library(SeuratObject)

IN_DIR  <- "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/seurat_export"
OUT_DIR <- "/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/seurat_export"

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

so <- FindVariableFeatures(so, nfeatures = nrow(so))
so <- ScaleData(so, features = VariableFeatures(so))
so <- RunPCA(so, features = VariableFeatures(so), npcs = 30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30)
so <- FindClusters(so, resolution = resolution)
so <- RunUMAP(so, dims = 1:30)
