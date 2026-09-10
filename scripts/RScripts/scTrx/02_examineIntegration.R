require(Seurat)

# Load saved data ----
setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/processed")

speciesObj <- readRDS("clustIntSpecies.rds")
# Harmony
Idents(speciesObj) <- "clusters_integrated_harmony"
DimPlot(speciesObj, reduction = "umap_integrated_harmony", label = T)

Idents(speciesObj) <- "napierCellTypes"
DimPlot(speciesObj, reduction = "umap_integrated_harmony", label = F)

Idents(speciesObj) <- "species"
DimPlot(speciesObj, reduction = "umap_integrated_harmony", label = F)

# RPCA
Idents(speciesObj) <- "clusters_integrated_rpca"
DimPlot(speciesObj, reduction = "umap_integrated_rpca", label = T)

Idents(speciesObj) <- "napierCellTypes"
DimPlot(speciesObj, reduction = "umap_integrated_rpca", label = F)

Idents(speciesObj) <- "species"
DimPlot(speciesObj, reduction = "umap_integrated_rpca", label = F)
