# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(Seurat)

## 0.2 Load data ----
### Set up one-to-one orthologs ----
setwd("/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/data/orthologs")
ortho1 <- read.delim("human_mouseGuineaPigDogCowGoatMacaque.txt")
ortho2 <- read.delim("human_ratRabbitPig.txt")
# Filter by one-to-one orthologs
colnames(ortho1) <- tolower(colnames(ortho1))
ortho1 <- ortho1 %>%
  filter(if_all(matches("homology.type"), ~ . == "ortholog_one2one")) %>%
  filter(if_all(matches("gene.name"), ~ . != "")) %>%
  distinct()

colnames(ortho2) <- tolower(colnames(ortho2))
ortho2 <- ortho2 %>%
  filter(if_all(contains("homology.type"), ~ .x == "ortholog_one2one")) %>%
  filter(if_all(matches("gene.name"), ~ . != "")) %>%
  distinct()

orthoAll <- inner_join(ortho1, ortho2, by = "gene.name") %>%
  select(matches("gene.name|homology.type")) %>%
  distinct()

# Remove any duplicate genes
speciesnames <- c("human", "mouse", "rat", "rabbit", "guineaPig", "cow", "dog", "macaque", "goat", "pig") %>% sort()
orthoNames <- data.frame(speciesName = speciesnames, orthoName = c("cattle", "dog", "goat", "guinea.pig", "human", "crab.eating.macaque", "mouse", "pig", "rabbit", "norway_rat"))
for(currSpecies in speciesnames){
  orthoName = orthoNames$orthoName[which(orthoNames$speciesName == currSpecies)]
  if(currSpecies == "human"){
    next
  }else if(currSpecies == "rat"){
    geneListName <- tolower("Norway.rat...BN.NHsdMcwi.gene.name")
  }else{
    geneListName <- paste(orthoName, "gene.name", sep = ".")
  }
  before <- nrow(orthoAll)
  orthoAll = orthoAll[!duplicated(orthoAll[[geneListName]]),]
  after <- nrow(orthoAll)
  message(currSpecies, " removed ", before-after, " duplicates")
}


# 1.0 Create one-to-one ortholog Seurat object ----
speciesnames <- c("human", "mouse", "rat", "rabbit", "guineaPig", "cow", "dog", "macaque", "goat", "pig") %>% sort()
orthoNames <- data.frame(speciesName = speciesnames, orthoName = c("cattle", "dog", "goat", "guinea.pig", "human", "crab.eating.macaque", "mouse", "pig", "rabbit", "norway_rat"))
orthoList <- list()
for(currSpecies in speciesnames){
  orthoName = orthoNames$orthoName[which(orthoNames$speciesName == currSpecies)]
  if(currSpecies == "human"){
    next
  }else if(currSpecies == "rat"){
    geneListName <- tolower("Norway.rat...BN.NHsdMcwi.gene.name")
  }else{
    geneListName <- paste(orthoName, "gene.name", sep = ".")
  }
  countMat <- GetAssayData(object = speciesList[[currSpecies]], assay = "RNA", layer = "counts")
  metadata <- speciesList[[currSpecies]][[]]
  geneNames <- as.data.frame(rownames(countMat))
  names(geneNames) <- geneListName
  geneNames <- left_join(geneNames, orthoAll)
  orthoCountMat <- countMat
  rownames(orthoCountMat) <- geneNames$gene.name
  keep_idx <- !is.na(rownames(orthoCountMat))
  orthoCountMat <- orthoCountMat[keep_idx, ]
  orthoList[[currSpecies]] <- CreateSeuratObject(counts = orthoCountMat, meta.data = metadata)
}


# 2.0 ----
