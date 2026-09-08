require(Seurat)
require(harmony)

setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/processed")
speciesObj <- readRDS("mergedSpecies.rds")

Idents(speciesObj) <- "species"
DimPlot(speciesObj)
