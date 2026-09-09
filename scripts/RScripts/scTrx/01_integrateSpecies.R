# 0.0 Setup ----
## 0.1 Load packages ----
.libPaths(c("~/R/R4.6.0_packages"))
options(repos = c(CRAN = "https://cloud.r-project.org")) 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(remotes)
library(harmony)
library(ggplot2)
message("----- PACKAGES LOADED -----")

## 0.2 Load functions ----
setwd("/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("normAndCluster.R")
source("clusterAllInts.R")
message("----- FUNCTIONS LOADED -----")

# ## 0.3 Load data ----
# ### Set up one-to-one orthologs ----
# setwd("/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/data/orthologs")
# ortho1 <- read.delim("human_mouseGuineaPigDogCowGoatMacaque.txt")
# ortho2 <- read.delim("human_ratRabbitPig.txt")
# # Filter by one-to-one orthologs
# colnames(ortho1) <- tolower(colnames(ortho1))
# ortho1 <- ortho1 %>%
#   filter(if_all(matches("homology.type"), ~ . == "ortholog_one2one")) %>%
#   filter(if_all(matches("gene.name"), ~ . != "")) %>%
#   distinct()
# 
# colnames(ortho2) <- tolower(colnames(ortho2))
# ortho2 <- ortho2 %>%
#   filter(if_all(contains("homology.type"), ~ .x == "ortholog_one2one")) %>%
#   filter(if_all(matches("gene.name"), ~ . != "")) %>%
#   distinct()
# 
# orthoAll <- inner_join(ortho1, ortho2, by = "gene.name") %>%
#   select(matches("gene.name|homology.type")) %>%
#   distinct()
# 
# # Remove any duplicate genes
# speciesnames <- c("human", "mouse", "rat", "rabbit", "guineaPig", "cow", "dog", "macaque", "goat", "pig") %>% sort()
# orthoNames <- data.frame(speciesName = speciesnames, orthoName = c("cattle", "dog", "goat", "guinea.pig", "human", "crab.eating.macaque", "mouse", "pig", "rabbit", "norway_rat"))
# for(currSpecies in speciesnames){
#   orthoName = orthoNames$orthoName[which(orthoNames$speciesName == currSpecies)]
#   if(currSpecies == "human"){
#     next
#   }else if(currSpecies == "rat"){
#     geneListName <- tolower("Norway.rat...BN.NHsdMcwi.gene.name")
#   }else{
#     geneListName <- paste(orthoName, "gene.name", sep = ".")
#   }
#   before <- nrow(orthoAll)
#   orthoAll = orthoAll[!duplicated(orthoAll[[geneListName]]),]
#   after <- nrow(orthoAll)
#   message(currSpecies, " removed ", before-after, " duplicates")
# }
# message("----- ORTHOLOGS LOADED -----")
# 
# ### Set up species Seurat objects ----
# setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/")
# speciesnames <- c("human", "mouse", "rat", "rabbit", "guineaPig", "cow", "dog", "macaque", "goat", "pig")
# speciesnames <- sort(speciesnames)
# speciesList <- list()
# for(currSpecies in speciesnames){
#   if(currSpecies == "human"){
#     message("reading human_list.rds")
#     human_list <- readRDS("human_list.rds")
#     human <- CreateSeuratObject(counts = human_list$counts,
#                                 meta.data = human_list$metadata)
#     human[["umap"]] <- CreateDimReducObject(embeddings = human_list$umap, key = "UMAP_", assay = "RNA")
#     speciesList[[currSpecies]] <- human
#     message("done")
#   }else{
#     filename <- paste(currSpecies, "rds", sep = ".")
#     message(paste("reading", filename))
#     speciesList[[currSpecies]] <- readRDS(filename)
#     message("done")
#   }
# }
# ### Set useful celltype ident ----
# celltype_merge <- c("dog", "goat", "cow", "guineaPig", "macaque", "pig", "rabbit", "human")
# celltypes_merge <- c("mouse")
# celltype <- c("rat")
# cellNames <- c("CTB", "EVT", "STB", "UNC", "BNC", "Stro", "Endo", "Mac", "Epi", "DC", "Invasive tro", "Bcell", "Tcell", "NKcell", "Mono", "GC", "Mes", "Leu", "Neu")
# for(currSpecies in speciesnames){
#   DefaultAssay(speciesList[[currSpecies]]) <- "RNA"
#   if(currSpecies %in% celltype_merge){
#     Idents(speciesList[[currSpecies]]) <- "celltype_merge"
#   }else if(currSpecies %in% celltypes_merge){
#     Idents(speciesList[[currSpecies]]) <- "celltypes_merge"
#   }else if(currSpecies %in% celltype){
#     Idents(speciesList[[currSpecies]]) <- "celltype"
#   }
#   newnames <- gsub("-", "", speciesList[[currSpecies]]@active.ident)
#   matchedIdx <- grepl(tolower("STGC"), tolower(newnames))
#   newnames[matchedIdx] <- "giantCell"
#   for(currCellType in cellNames){
#     matchedIdx <- grepl(tolower(currCellType), tolower(newnames))
#     newnames[matchedIdx] <- currCellType
#   }
#   matchedIdx <- grepl(tolower("SynT"), tolower(newnames))
#   newnames[matchedIdx] <- "STB"
#   matchedIdx <- grepl(tolower("NK"), tolower(newnames))
#   newnames[matchedIdx] <- "NKcells"
#   matchedIdx <- grepl(tolower("SpT"), tolower(newnames))
#   newnames[matchedIdx] <- "SpT"
#   matchedIdx <- grepl(tolower("ST"), tolower(newnames))
#   newnames[matchedIdx] <- "STB"
#   matchedIdx <- grepl(tolower("giantCell"), tolower(newnames))
#   newnames[matchedIdx] <- "s-TGC"
# 
#   speciesList[[currSpecies]]$napierCellTypes <- newnames
#   Idents(speciesList[[currSpecies]]) <- "napierCellTypes"
# }
# message("----- TAN DATASETS LOADED -----")
# 
# 
# # 1.0 Create one-to-one ortholog Seurat object ----
# speciesnames <- c("human", "mouse", "rat", "rabbit", "guineaPig", "cow", "dog", "macaque", "goat", "pig") %>% sort()
# orthoNames <- data.frame(speciesName = speciesnames, orthoName = c("cattle", "dog", "goat", "guinea.pig", "human", "crab.eating.macaque", "mouse", "pig", "rabbit", "norway_rat"))
# orthoList <- list()
# for(currSpecies in speciesnames){
#   orthoName = orthoNames$orthoName[which(orthoNames$speciesName == currSpecies)]
#   if(currSpecies == "human"){
#     next
#   }else if(currSpecies == "rat"){
#     geneListName <- tolower("Norway.rat...BN.NHsdMcwi.gene.name")
#   }else{
#     geneListName <- paste(orthoName, "gene.name", sep = ".")
#   }
#   countMat <- GetAssayData(object = speciesList[[currSpecies]], assay = "RNA", layer = "counts")
#   metadata <- speciesList[[currSpecies]][[]]
#   geneNames <- as.data.frame(rownames(countMat))
#   names(geneNames) <- geneListName
#   geneNames <- left_join(geneNames, orthoAll)
#   orthoCountMat <- countMat
#   rownames(orthoCountMat) <- geneNames$gene.name
#   keep_idx <- !is.na(rownames(orthoCountMat))
#   orthoCountMat <- orthoCountMat[keep_idx, ]
#   orthoList[[currSpecies]] <- CreateSeuratObject(counts = orthoCountMat, meta.data = metadata)
#   orthoList[[currSpecies]]$species <- currSpecies
# }
# message("----- SEURAT OBJECTS FILTERED BY ONE-TO-ONE ORTHOLOGS -----")
# 
# # 2.0 Merge Seurat objects ----
# mergeNames <- speciesnames[speciesnames != "human"]
# mergeNames <- c("human", mergeNames)
# speciesObj <- merge(speciesList[["human"]], y = c(orthoList[[1]],
#                                           orthoList[[2]],
#                                           orthoList[[3]],
#                                           orthoList[[4]],
#                                           orthoList[[5]],
#                                           orthoList[[6]],
#                                           orthoList[[7]],
#                                           orthoList[[8]],
#                                           orthoList[[9]]), add.cell.ids = mergeNames)
# speciesObj <- JoinLayers(speciesObj)
# counts <- GetAssayData(speciesObj, assay = "RNA", layer = "counts")
# totals <- Matrix::colSums(counts)
# speciesObj <- subset(speciesObj, cells = colnames(counts)[totals > 0])
# # Refresh metadata
# speciesObj$nCount_RNA   <- Matrix::colSums(GetAssayData(speciesObj, assay = "RNA", layer = "counts"))
# speciesObj$nFeature_RNA <- Matrix::colSums(GetAssayData(speciesObj, assay = "RNA", layer = "counts") > 0)
# speciesObj
# speciesObj[["RNA"]] <- split(speciesObj[["RNA"]], f = speciesObj$species)
# speciesObj
# speciesObj <- normAndCluster(speciesObj)
# setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/processed")
# saveRDS(speciesObj,"mergedSpecies.rds")
# message("----- SEURAT OBJECTS MERGED -----")

# Load saved data ----
setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/processed")
speciesObj <- readRDS("integratedSpecies.rds")
message("----- SEURAT OBJECT LOADED -----")

# # 3.0 Integrate objects ----
# message("Running Harmony integration...")
# speciesObj <- IntegrateLayers(
#   object       = speciesObj,
#   method       = HarmonyIntegration,
#   orig.reduction = "pca",
#   new.reduction  = "integrated_harmony",
#   verbose      = T
# )
# setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/processed")
# saveRDS(speciesObj,"integratedSpecies.rds")
# message("----- HARMONY INTEGRATION DONE -----")

# # CCA Integration
# message("  Running CCA integration...")
# speciesObj <- IntegrateLayers(
#   object       = speciesObj,
#   method       = CCAIntegration,
#   orig.reduction = "pca",
#   new.reduction  = "integrated.cca",
#   verbose      = T, 
#   normalization.method = "SCT"
# )
# setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/processed")
# saveRDS(speciesObj,"integratedSpecies.rds")
# message("----- CCA INTEGRATION DONE -----")

# # RPCA Integration
# message("Running RPCA integration...")
# speciesObj <- IntegrateLayers(
#   object       = speciesObj,
#   method       = RPCAIntegration,
#   orig.reduction = "pca",
#   new.reduction  = "integrated.rpca",
#   verbose      = T, 
#   normalization.method = "SCT"
# )
# setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/processed")
# saveRDS(speciesObj,"integratedSpecies.rds")
# message("----- RPCA INTEGRATION DONE -----")

speciesObj <- clusterAllInts(speciesObj)
setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/processed")
saveRDS(speciesObj,"clustIntSpecies.rds")



