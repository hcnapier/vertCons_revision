# 0.0 Setup ----
## Run scripts/shellScripts/getTanDatasets.sh to get datasets from Figshare

## 0.1 Load packages ----
require(Seurat)
require(dplyr)
require(ggplot2)

## 0.2 Load data ----
setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/")
speciesnames <- c("mouse", "rat", "rabbit", "guineaPig", "cow", "dog", "macaque", "goat", "pig")
speciesnames <- sort(speciesnames)
speciesList <- list()
for(currSpecies in speciesnames){
  filename <- paste(currSpecies, "rds", sep = ".")
  message(paste("reading", filename))
  speciesList[[currSpecies]] <- readRDS(filename)
  message("done")
}


# 1.0 First pass ----
## 1.1 Clustering ----
DimPlot(speciesList[["cow"]])
DimPlot(speciesList[["dog"]])
DimPlot(speciesList[["goat"]])
DimPlot(speciesList[["guineaPig"]])
DimPlot(speciesList[["macaque"]])
DimPlot(speciesList[["mouse"]])
DimPlot(speciesList[["pig"]])
DimPlot(speciesList[["rabbit"]])
DimPlot(speciesList[["rat"]])

## 1.2 Set useful celltype ident ----
# each species has a slightly different cell type name ident, so find the most useful one for each species
celltype_merge <- c("dog", "goat", "cow", "guineaPig", "macaque", "pig", "rabbit")
celltypes_merge <- c("mouse")
celltype <- c("rat")
cellNames <- c("CTB", "EVT", "STB", "UNC", "BNC", "Stro", "Endo", "Mac", "Epi", "DC", "Invasive tro", "Bcells", "Tcells", "NKcells", "Mono", "GC", "Mes", "Leu", "Neu")
for(currSpecies in speciesnames){
  DefaultAssay(speciesList[[currSpecies]]) <- "RNA"
  if(currSpecies %in% celltype_merge){
    Idents(speciesList[[currSpecies]]) <- "celltype_merge"
  }else if(currSpecies %in% celltypes_merge){
    Idents(speciesList[[currSpecies]]) <- "celltypes_merge"
  }else if(currSpecies %in% celltype){
    Idents(speciesList[[currSpecies]]) <- "celltype"
  }
  newnames <- gsub("-", "", speciesList[[currSpecies]]@active.ident)
  matchedIdx <- grepl(tolower("STGC"), tolower(newnames))
  newnames[matchedIdx] <- "giantCell"
  for(currCellType in cellNames){
    matchedIdx <- grepl(tolower(currCellType), tolower(newnames))
    newnames[matchedIdx] <- currCellType
  }
  matchedIdx <- grepl(tolower("SynT"), tolower(newnames))
  newnames[matchedIdx] <- "STB"
  matchedIdx <- grepl(tolower("NK"), tolower(newnames))
  newnames[matchedIdx] <- "NKcells"
  matchedIdx <- grepl(tolower("SpT"), tolower(newnames))
  newnames[matchedIdx] <- "SpT"
  matchedIdx <- grepl(tolower("ST"), tolower(newnames))
  newnames[matchedIdx] <- "STB"
  matchedIdx <- grepl(tolower("giantCell"), tolower(newnames))
  newnames[matchedIdx] <- "s-TGC"
  
  speciesList[[currSpecies]]$napierCellTypes <- newnames
  Idents(speciesList[[currSpecies]]) <- "napierCellTypes"
}

## 1.3 Plot un-integrated data ----
### Labeled by cell type ----
DimPlot(speciesList[["cow"]]) + ggtitle("Cow")
DimPlot(speciesList[["dog"]]) + ggtitle("Dog")
DimPlot(speciesList[["goat"]]) + ggtitle("Goat")
DimPlot(speciesList[["guineaPig"]]) + ggtitle("Guinea Pig")
DimPlot(speciesList[["macaque"]]) + ggtitle("Macaque")
DimPlot(speciesList[["mouse"]]) + ggtitle("Mouse")
DimPlot(speciesList[["pig"]]) + ggtitle("Pig")
DimPlot(speciesList[["rabbit"]]) + ggtitle("Rabbit")
DimPlot(speciesList[["rat"]]) + ggtitle("Rat")

# 2.0 Assess batch integration ----
## 2.1 Switch to integrated assay ----
for(currSpecies in speciesnames){
  if(speciesList[[currSpecies]]@assays %>% length() > 1){
    message(currSpecies, " assay set to integrated")
    DefaultAssay(speciesList[[currSpecies]]) <- "integrated"
    Idents(speciesList[[currSpecies]]) <- "sample"
  }
}

## 2.2 Plot integrated object ----
DimPlot(speciesList[["cow"]], alpha = 0.25) + ggtitle("Cow")
DimPlot(speciesList[["goat"]], alpha = 0.25) + ggtitle("Goat")
DimPlot(speciesList[["guineaPig"]], alpha = 0.25) + ggtitle("Guinea Pig")
DimPlot(speciesList[["macaque"]], alpha = 0.25) + ggtitle("Macaque")
DimPlot(speciesList[["mouse"]], alpha = 0.25) + ggtitle("Mouse")
DimPlot(speciesList[["pig"]], alpha = 0.25) + ggtitle("Pig")


# 3.0 Save processed object list ----
## THIS TAKES A LONG TIME!!
setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/processed")
filePaths <- paste(speciesnames, ".rds", sep = "")
Map(saveRDS, speciesList, filePaths)
