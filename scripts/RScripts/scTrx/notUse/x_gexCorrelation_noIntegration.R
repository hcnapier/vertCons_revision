# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(Seurat)
require(reshape2)
require(ggplot2)
require(stringr)
require(tidyverse)


## 0.2 Load data ----
### Orthologs ----
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

# 1.0 Pseudobulk ----
names(orthoAll) <- names(orthoAll) %>% tolower()
orthoNames <- data.frame(speciesName = speciesnames, orthoName = c("cattle", "dog", "goat", "guinea.pig", "human", "crab.eating.macaque", "mouse", "pig", "rabbit", "norway.rat"))
pseudobulk <- list()
convertedToHuman <- list()
for(currSpecies in speciesnames){
  orthoName = orthoNames$orthoName[which(orthoNames$speciesName == currSpecies)]
  if(currSpecies == "human"){
    geneListName <- "gene.name"
  }else if(currSpecies == "rat"){
    geneListName <- tolower("Norway.rat...BN.NHsdMcwi.gene.name")
  }else{
    geneListName <- paste(orthoName, "gene.name", sep = ".")
  }
  pseudobulk[[currSpecies]] <- AggregateExpression(speciesList[[currSpecies]], features = pull(orthoAll[geneListName]), return.seurat = F)$RNA
  convert = data.frame(row.names(pseudobulk[[currSpecies]]))
  names(convert) <- geneListName
  convertedToHuman[[currSpecies]] <- inner_join(convert, orthoAll) %>% distinct()
  row.names(pseudobulk[[currSpecies]]) <- convertedToHuman[[currSpecies]]$gene.name
}

## 1.1 Get shared genes ----
sharedGenes <- Reduce(intersect, list(convertedToHuman[["cow"]]$gene.name,
                                      convertedToHuman[["dog"]]$gene.name, 
                                      convertedToHuman[["goat"]]$gene.name, 
                                      convertedToHuman[["guineaPig"]]$gene.name, 
                                      convertedToHuman[["macaque"]]$gene.name, 
                                      convertedToHuman[["mouse"]]$gene.name, 
                                      convertedToHuman[["pig"]]$gene.name, 
                                      convertedToHuman[["rabbit"]]$gene.name, 
                                      convertedToHuman[["rat"]]$gene.name, 
                                      row.names(pseudobulk[["human"]])))
# Subset each pseudobulked matrix to only include shared genes
pseudobulk_shared <- list()
for(currSpecies in speciesnames){
  pseudobulk_shared[[currSpecies]] <- pseudobulk[[currSpecies]][row.names(pseudobulk[[currSpecies]]) %in% sharedGenes,]
}

# Format data ----
pseudobulk_shared <- lapply(pseudobulk_shared, as.matrix)
for(currSpecies in speciesnames){
  colnames(pseudobulk_shared[[currSpecies]]) <- paste(colnames(pseudobulk_shared[[currSpecies]]), currSpecies, sep = "_") # order columns by cell type
  #colnames(pseudobulk_shared[[currSpecies]]) <- paste(currSpecies, colnames(pseudobulk_shared[[currSpecies]]), sep = "_") # order columns by species
  colsOrdered <- sort(colnames(pseudobulk_shared[[currSpecies]]))
  pseudobulk_shared[[currSpecies]] <- pseudobulk_shared[[currSpecies]][, colsOrdered]
}

## Merge all matrices 
pseudobulkMerged <- pseudobulk_shared[["cow"]] %>% as.data.frame()
pseudobulkMerged$RowNames <- pseudobulkMerged %>% row.names()
for(currSpecies in speciesnames[2:10]){
  message("merging ", currSpecies)
  temp <- pseudobulk_shared[[currSpecies]] %>% as.data.frame()
  temp$RowNames <- temp %>% row.names()
  message(head(temp$RowNames, n = 5))
  pseudobulkMerged <- merge(pseudobulkMerged, temp)
}

## order columns by cell type ##
#colsOrdered <- sort(colnames(pseudobulkMerged))
#pseudobulkMerged <- pseudobulkMerged[,colsOrdered]
####

## Order by evolutionary divergence ##
colNames <- colnames(pseudobulkMerged)
species_ordered <- c("human", "macaque", "guineaPig", "rat", "mouse", "rabbit", "pig", "cow", "goat", "dog")
parts <- strsplit(colNames, "\\_")
species  <- sapply(parts, `[`, 1)
celltype <- sapply(parts, `[`, 2)
colsOrdered <- colNames[order(factor(species, levels = species_ordered), celltype)]
pseudobulkMerged <- pseudobulkMerged[,colsOrdered]
####

# Convert back into a matrix
row.names(pseudobulkMerged) <- pseudobulkMerged$RowNames
pseudobulkMerged$RowNames <- NULL
pseudobulkMerged <- pseudobulkMerged %>% as.matrix()


# 4.0 Pearson correlation ----
## 4.1 All pairwise comparisons ----
cormat_allPairwise <- cor(pseudobulkMerged, method = "spearman")
melted_cormat_allPairwise <- melt(cormat_allPairwise)
allPairwise_corPlot <- ggplot(data = melted_cormat_allPairwise, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 5, hjust = 1)) +
  theme(axis.text.y = element_text(size = 5)) + 
  coord_fixed()
allPairwise_corPlot

# Dotplot ----
humanCompDF <- as.data.frame(cormat_allPairwise) %>%
  rownames_to_column("row_id") %>%
  pivot_longer(-row_id, names_to = "col_id", values_to = "correlation")

# --- 3. Parse species and cell type from labels ---
# Adjust the split pattern if your naming convention differs (e.g., more than one underscore)
humanCompDF <- humanCompDF %>%
  separate(row_id, into = c("row_celltype", "row_species"), sep = "_", extra = "merge") %>%
  separate(col_id, into = c("col_celltype", "col_species"), sep = "_", extra = "merge")

# --- 4. Keep only comparisons: human rows vs non-human columns ---
humanCompDF <- humanCompDF %>%
  filter(row_species == "human", col_species != "human")

species_ordered <- c("macaque", "guineaPig", "rat", "mouse", "rabbit", "pig", "cow", "goat", "dog")

humanCompDF <- humanCompDF %>%
  filter(col_species %in% species_ordered) %>%
  mutate(col_species = factor(col_species, levels = species_ordered))

celltypes_per_species <- humanCompDF %>%
  group_by(col_species) %>%
  summarise(celltypes = list(unique(col_celltype)))

shared_celltypes <- Reduce(intersect, celltypes_per_species$celltypes)

humanCompDF <- humanCompDF %>%
  filter(col_celltype %in% shared_celltypes)

# --- 5. Plot: human cell type vs other species' cell type, faceted by species ---
ggplot(humanCompDF, aes(x = col_species, y = correlation, color = row_celltype)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.8) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()) +
  labs(x = "Species", y = "Correlation",
       color = "Human cell type",
       title = "Human Cell Type Correlations Across Species (Shared Cell Types Only)")
