# normAndCluster.R
# Function to normalize and cluster a Seurat object
# Uses SCTransform
# Hailey Napier
# April 28, 2026

setwd("/hpc/group/vertgenlab/hailey/zebrin_evolution/251105_bannerRNAseq/scripts/humanPC_RNAseq/functions")
source("PCAbugFix.R")
require(Seurat)

normAndCluster <- function(obj, bugFix = F, resolution = 0.8, species = "human"){
  obj$seurat_clusters <- NULL
  #obj[["RNA"]] <- CreateAssayObject(counts = obj[["RNA"]]$counts) # reset assay dimensions
  if(species == "human"){
    obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
  }else if(species == "rhesus"){
    mt_genes <- c("ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
                  "COX1","COX2","COX3","CO1","CO2","CO3",
                  "ATP6","ATP8","CYTB","CYB")
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mt_genes)
  }
  options(future.globals.maxSize = 10000 * 1024^2) 
  obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = T, conserve.memory = T, vst.flavor = "v2")
  obj <- FindVariableFeatures(obj)
  if(bugFix){
    obj <- PCAbugFix(obj)
    obj <- FindNeighbors(obj, dims = 2:30)
    obj <- FindClusters(obj, resolution = resolution)
    obj <- RunUMAP(obj, dims = 2:30)
    return(obj)
  }else{
    obj <- RunPCA(obj)
    obj <- FindNeighbors(obj, dims = 1:30)
    obj <- FindClusters(obj, resolution = resolution)
    obj <- RunUMAP(obj, dims = 1:30)
    return(obj)
  }
}