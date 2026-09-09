# Load saved data ----
setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/processed")
speciesObj <- readRDS("integratedSpecies.rds")

Idents(speciesObj) <- "cluster"