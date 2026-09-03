# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(Seurat)

## 0.2 Load data ----
### Orthologs ----
setwd("/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/data/orthologs")
ortho1 <- read.delim("human_mouseGuineaPigDogCowGoatMacaque.txt")
ortho2 <- read.delim("human_ratRabbitPig.txt")
# Filter by one 2 one orthologs
ortho1 <- ortho1 %>%
  filter(if_all(matches("homology_type"), ~ . == "ortholog_one2one")) %>%
  filter(Gene_stable_ID != "") %>%
  filter(if_all(matches("gene_name"), ~ . != "")) %>%
  distinct()

ortho2 <- ortho2 %>%
  filter(if_all(contains("homology_type"), ~ .x == "ortholog_one2one")) %>%
  filter(Gene_stable_ID != "") %>%
  filter(if_all(matches("gene_name"), ~ . != "")) %>%
  distinct()

orthoAll <- inner_join(ortho1, ortho2, by = "Gene_stable_ID") %>%
  select(matches("gene_name|Gene_stable_ID|homology_type"))

names(orthoAll) <- names(orthoAll) %>% tolower()
orthoNames <- data.frame(speciesName = speciesnames, orthoName = c("cattle", "dog", "goat", "guinea_pig", "crab.eating_macaque", "mouse", "pig", "rabbit", "norway_rat"))
for(currSpecies in speciesnames){
  orthoName = orthoNames$orthoName[which(orthoNames$speciesName == currSpecies)]
  if(currSpecies == "rat"){
    geneListName <- tolower("Norway_rat_._BN.NHsdMcwi_gene_name")
  }else{
    geneListName <- paste(orthoName, "gene_name", sep = "_")
  }
  orthoAll = orthoAll[!duplicated(orthoAll[[geneListName]]),]
}

### Seurat object list ----
setwd("/work/hcn4/260630_vertCons_wd/scTrx/rObjs/processed")
speciesnames <- c("mouse", "rat", "rabbit", "guineaPig", "cow", "dog", "macaque", "goat", "pig")
speciesnames <- sort(speciesnames)
speciesList <- list()
for(currSpecies in speciesnames){
  filename <- paste(currSpecies, "rds", sep = ".")
  message(paste("reading", filename))
  speciesList[[currSpecies]] <- readRDS(filename)
  message("done")
}


# 1.0 Pseudobulk ----
#names(orthoAll) <- names(orthoAll) %>% tolower()
#orthoNames <- data.frame(speciesName = speciesnames, orthoName = c("cattle", "dog", "goat", "guinea_pig", "crab.eating_macaque", "mouse", "pig", "rabbit", "norway_rat"))
pseudobulk <- list()
convertedToHuman <- list()
for(currSpecies in speciesnames){
  orthoName = orthoNames$orthoName[which(orthoNames$speciesName == currSpecies)]
  if(currSpecies == "rat"){
    geneListName <- tolower("Norway_rat_._BN.NHsdMcwi_gene_name")
  }else{
    geneListName <- paste(orthoName, "gene_name", sep = "_")
  }
  pseudobulk[currSpecies] <- AggregateExpression(speciesList[[currSpecies]], features = pull(orthoAll[geneListName]), return.seurat = F)
  convert = data.frame(row.names(pseudobulk[[currSpecies]]))
  names(convert) <- geneListName
  convertedToHuman[[currSpecies]] <- inner_join(convert, orthoAll) %>% distinct()
  row.names(pseudobulk[[currSpecies]]) <- convertedToHuman[[currSpecies]]$gene_stable_id
}

## 1.1 Get shared genes ----
sharedGenes <- Reduce(intersect, list(convertedToHuman[["cow"]]$gene_stable_id,
                                      convertedToHuman[["dog"]]$gene_stable_id, 
                                      convertedToHuman[["goat"]]$gene_stable_id, 
                                      convertedToHuman[["guineaPig"]]$gene_stable_id, 
                                      convertedToHuman[["macaque"]]$gene_stable_id, 
                                      convertedToHuman[["mouse"]]$gene_stable_id, 
                                      convertedToHuman[["pig"]]$gene_stable_id, 
                                      convertedToHuman[["rabbit"]]$gene_stable_id, 
                                      convertedToHuman[["rat"]]$gene_stable_id))


