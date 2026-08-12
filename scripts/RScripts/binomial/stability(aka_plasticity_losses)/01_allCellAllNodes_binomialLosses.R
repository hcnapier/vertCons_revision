# 01_binomialLosses
# Instability = losses
# Christi's original calculations carry a Gaussian assumption. A Gaussian may not fit the data well, but a binomial may be more appropriate. 
# Calculate the probabilities that a given node has enrichment for each cell type option relative to the baseline enrichment (= total number of regions for that cell type over all cell type regions, calculated with the current node held out)

# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
install.packages("pheatmap")
require(pheatmap)
require(tidyr)
require(ggplot2)
require(stringr)

## 0.2 Load data ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/data")
loss_speciesNodes <- read.csv("losses_speciesPerNode.csv", header = T)
speciesLosses <- read.csv("hailey_speciesLosses_forR.csv", header = T)
speciesLosses$X.1 <- NULL
speciesLosses$X.2 <- NULL
speciesLosses$X.3 <- NULL 
speciesLosses$X <- NULL

## 0.3 Process data ----
# Add in missing species and correct mistakes 
fixNodeAssignments <- data.frame(Node = c(11,9,9,5,NA,9), Species = c("macEug2", "ailMel1", "felCat8", "calJac3", "panTro6", "musFur1"))
loss_speciesNodes <- rbind(loss_speciesNodes, fixNodeAssignments)
# Add node data to species loss counts
speciesLosses <- full_join(loss_speciesNodes, speciesLosses)
speciesLosses <- speciesLosses %>%
  filter(Species != "panTro6") %>%
  filter(!is.na(nSpeciesLosses))
# Calculate average number of losses per node and cell type pair
speciesLosses$NodeCellType <- paste(speciesLosses$Node, speciesLosses$CellType, sep = "_")
avNodeLosses <- speciesLosses %>%
  select(Node, nSpeciesLosses, CellType, NodeCellType) %>%
  group_by(NodeCellType) %>%
  summarize(avNodeLosses = mean(nSpeciesLosses))
nodeLosses <- full_join(avNodeLosses, speciesLosses)
nodeLosses <- nodeLosses %>%
  select(avNodeLosses, Node, CellType, nNodeRegions)
nodeLosses <- nodeLosses %>% distinct()
