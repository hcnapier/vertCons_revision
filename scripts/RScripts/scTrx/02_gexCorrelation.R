# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(Seurat)
require(reshape2)
require(ggplot2)

## 0.2 Load data ----
### Orthologs ----
setwd("/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/data/orthologs")
ortho1 <- read.delim("human_mouseGuineaPigDogCowGoatMacaque.txt")
ortho2 <- read.delim("human_ratRabbitPig.txt")
# Filter by one-to-one orthologs
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
  select(matches("gene_name|Gene_stable_ID|homology_type")) %>%
  distinct()

names(orthoAll) <- names(orthoAll) %>% tolower()
speciesnames <- c("mouse", "rat", "rabbit", "guineaPig", "cow", "dog", "macaque", "goat", "pig") %>% sort()
orthoNames <- data.frame(speciesName = speciesnames, orthoName = c("cattle", "dog", "goat", "guinea_pig", "crab.eating_macaque", "mouse", "pig", "rabbit", "norway_rat"))
for(currSpecies in speciesnames){
  orthoName = orthoNames$orthoName[which(orthoNames$speciesName == currSpecies)]
  if(currSpecies == "rat"){
    geneListName <- tolower("Norway_rat_._BN.NHsdMcwi_gene_name")
  }else{
    geneListName <- paste(orthoName, "gene_name", sep = "_")
  }
  before <- nrow(orthoAll)
  orthoAll = orthoAll[!duplicated(orthoAll[[geneListName]]),]
  after <- nrow(orthoAll)
  message(currSpecies, " removed ", before-after, " duplicates")
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
names(orthoAll) <- names(orthoAll) %>% tolower()
orthoNames <- data.frame(speciesName = speciesnames, orthoName = c("cattle", "dog", "goat", "guinea_pig", "crab.eating_macaque", "mouse", "pig", "rabbit", "norway_rat"))
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

row.names(pseudobulk[["guineaPig"]]) %>% length()
convertedToHuman[["guineaPig"]] %>% nrow()
row.names(pseudobulk[["guineaPig"]]) %in% convertedToHuman[["guineaPig"]]$gene_stable_id %>% sum()
head(row.names(pseudobulk[["guineaPig"]])) 
orthoAll[["guinea_pig_gene_name"]]


## 1.1 Get shared genes ----
sharedGenes <- Reduce(intersect, list(convertedToHuman[["cow"]]$gene_stable_id,
                                      convertedToHuman[["dog"]]$gene_stable_id, 
                                      convertedToHuman[["goat"]]$gene_stable_id, 
                                      #convertedToHuman[["guineaPig"]]$gene_stable_id, 
                                      convertedToHuman[["macaque"]]$gene_stable_id, 
                                      convertedToHuman[["mouse"]]$gene_stable_id, 
                                      convertedToHuman[["pig"]]$gene_stable_id, 
                                      convertedToHuman[["rabbit"]]$gene_stable_id, 
                                      convertedToHuman[["rat"]]$gene_stable_id))
# Subset each pseudobulked matrix to only include shared genes
speciesnames <- speciesnames[speciesnames != "guineaPig"]
pseudobulk_shared <- list()
for(currSpecies in speciesnames){
  pseudobulk_shared[[currSpecies]] <- pseudobulk[[currSpecies]][sharedGenes,]
}

# Format data ----
pseudobulk_shared <- lapply(pseudobulk_shared, as.matrix)
for(currSpecies in speciesnames){
  #colnames(pseudobulk_shared[[currSpecies]]) <- paste(colnames(pseudobulk_shared[[currSpecies]]), currSpecies, sep = "_") # order columns by cell type
  colnames(pseudobulk_shared[[currSpecies]]) <- paste(currSpecies, colnames(pseudobulk_shared[[currSpecies]]), sep = "_") # order columns by species
  colsOrdered <- sort(colnames(pseudobulk_shared[[currSpecies]]))
  pseudobulk_shared[[currSpecies]] <- pseudobulk_shared[[currSpecies]][, colsOrdered]
}

## Merge all matrices 
pseudobulkMerged <- pseudobulk_shared[["cow"]] %>% as.data.frame()
pseudobulkMerged$RowNames <- pseudobulkMerged %>% row.names()

colNames <- colnames(pseudobulkMerged)
species_ordered <- c("macaque", "rat", "mouse", "rabbit", "pig", "cow", "goat", "dog")
parts <- strsplit(colNames, "\\_")
species  <- sapply(parts, `[`, 1)
celltype <- sapply(parts, `[`, 2)
colsOrdered <- colNames[order(factor(species, levels = species_ordered), celltype)]
for(currSpecies in speciesnames[2:8]){
  temp <- pseudobulk_shared[[currSpecies]] %>% as.data.frame()
  temp$RowNames <- temp %>% row.names()
  pseudobulkMerged <- merge(pseudobulkMerged, temp)
  # order columns by cell type
  #colsOrdered <- sort(colnames(pseudobulkMerged))
  pseudobulkMerged <- pseudobulkMerged[,colsOrdered]
}

# Convert back into a matrix
row.names(pseudobulkMerged) <- pseudobulkMerged$RowNames
pseudobulkMerged$RowNames <- NULL
pseudobulkMerged <- pseudobulkMerged %>% as.matrix()


# 4.0 Pearson correlation ----
## 4.1 All pairwise comparisons ----
cormat_allPairwise <- cor(pseudobulkMerged)
melted_cormat_allPairwise <- melt(cormat_allPairwise)
allPairwise_corPlot <- ggplot(data = melted_cormat_allPairwise, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed()
allPairwise_corPlot

