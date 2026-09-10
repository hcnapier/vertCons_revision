# 04_convertToAnnData
# Convert each species Seurat object into an AnnData object for UCE
# Hailey Napier
# September 10, 2026

# 0.0 Setup ----
## 0.1 Load packages ----
Sys.setenv(LD_LIBRARY_PATH = paste0(
  "/hpc/group/vertgenlab/hailey/software/miniconda3/envs/rEnv/lib:",
  Sys.getenv("LD_LIBRARY_PATH")
))
Sys.setenv(RETICULATE_CONDA = "/hpc/group/vertgenlab/hailey/software/miniconda3/bin/conda")
require(Seurat)

## 0.2 Load data ----
### Set up species Seurat objects ----
setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/")
speciesnames <- c("human", "mouse", "rat", "rabbit", "guineaPig", "cow", "dog", "macaque", "goat", "pig")
speciesnames <- sort(speciesnames)
speciesList <- list()
for(currSpecies in speciesnames){
  if(currSpecies == "human"){
    message("reading human_list.rds")
    human_list <- readRDS("human_list.rds")
    human <- CreateSeuratObject(counts = human_list$counts,
                                meta.data = human_list$metadata)
    human[["umap"]] <- CreateDimReducObject(embeddings = human_list$umap, key = "UMAP_", assay = "RNA")
    speciesList[[currSpecies]] <- human
    message("done")
  }else{
    filename <- paste(currSpecies, "rds", sep = ".")
    message(paste("reading", filename))
    speciesList[[currSpecies]] <- readRDS(filename)
    message("done")
  }
}
### Set useful celltype ident ----
celltype_merge <- c("dog", "goat", "cow", "guineaPig", "macaque", "pig", "rabbit", "human")
celltypes_merge <- c("mouse")
celltype <- c("rat")
cellNames <- c("CTB", "EVT", "STB", "UNC", "BNC", "Stro", "Endo", "Mac", "Epi", "DC", "Invasive tro", "Bcell", "Tcell", "NKcell", "Mono", "GC", "Mes", "Leu", "Neu")
for(currSpecies in speciesnames){
  speciesList[[currSpecies]]$species <- currSpecies
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
message("----- TAN DATASETS LOADED -----")


# 1.0 Convert to AnnData ----
require(reticulate)
#Sys.setenv(RETICULATE_CONDA = "/hpc/group/vertgenlab/hailey/software/miniconda3/bin/conda")
use_condaenv("/hpc/group/vertgenlab/hailey/software/miniconda3/envs/rEnv", conda = "/hpc/group/vertgenlab/hailey/software/miniconda3/bin/conda", required = TRUE)
require(scCustomize)
# note that this requires python 3.10, the newer anndata api isn't compatible with reticulate

for(currSpecies in speciesnames){
  filename <- paste(currSpecies, ".h5ad", sep = "")
  as.anndata(x = speciesList[[currSpecies]], file_path = "/work/hcn4/260630_vertCons_wd/scTrx/annDataObjs", file_name = filename)
}
