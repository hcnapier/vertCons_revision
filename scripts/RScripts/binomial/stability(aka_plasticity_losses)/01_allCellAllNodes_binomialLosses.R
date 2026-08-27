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
MYA_losses <- c(10, 15, 20, 30, 45, 75, 85, 87, 95, 100, 160, 180, 320, 350, 415, 430)

# 0.3 Source functions ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("subsetAndBinomCellTypeList.R")

## 0.4 Process data ----
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
  summarize(nSuccess = mean(nSpeciesLosses))
avNodeLosses$nSuccess <- round(avNodeLosses$nSuccess)
nodeLosses <- full_join(avNodeLosses, speciesLosses)
nodeLosses <- nodeLosses %>%
  select(nSuccess, Node, CellType, nNodeRegions)
nodeLosses <- nodeLosses %>% distinct()
nodeLosses$propSuccess <- nodeLosses$nSuccess/nodeLosses$nNodeRegions
nodeLosses$nodeName <- paste("Node",nodeLosses$Node, sep = "")
nodeNames_losses <- nodeLosses$nodeName %>% unique()
nodeLosses$nRegions <- nodeLosses$nNodeRegions
colOrder <- c("Node1", "Node2", "Node3", "Node4", "Node5", "Node6", "Node7", "Node8", "Node9", "Node10", "Node11", "Node12", "Node13", "Node14", "Node15", "Node16")
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(nodeLosses, "nodeLosses.rds")

# Get total node regions 
totalNodeRegions_losses <- nodeLosses %>%
  select(Node, nNodeRegions) %>%
  group_by(Node) %>%
  summarize(nRegions = sum(nNodeRegions))
totalNodeRegions_losses$name <- paste("Node", totalNodeRegions_losses$Node, sep = "")
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(totalNodeRegions_losses, "totalNodeRegions_losses.rds")

# Get total cell type regions
totalCellTypeRegions_losses <- nodeLosses %>%
  select(CellType, nNodeRegions) %>%
  group_by(CellType) %>%
  summarize(nRegions = sum(nNodeRegions))
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(totalCellTypeRegions_losses, "totalCellTypeRegions_losses.rds")

# Get total losses per node
totalNodeLosses <- nodeLosses %>%
  select(Node, nSuccess) %>%
  group_by(Node) %>%
  summarize(nLosses = sum(nSuccess))


# 1.0 Binomial test, all cell types ----
## 1.1 Set up output data structures ----
nNodes <- nrow(totalNodeRegions_losses)
nCellTypes <- nrow(totalCellTypeRegions_losses)
binomPrMat_loss <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(binomPrMat_loss) <- c(nodeNames_losses)
rownames(binomPrMat_loss) <- totalCellTypeRegions_losses$CellType
nullPrMat <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(nullPrMat) <- nodeNames_losses
rownames(nullPrMat) <- totalCellTypeRegions_losses$CellType

binomPvalMat_loss <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(binomPvalMat_loss) <- nodeNames_losses
rownames(binomPvalMat_loss) <- totalCellTypeRegions_losses$CellType

## 1.2 Run binomial test for all cell types ----
nodes <- nodeLosses$Node %>% unique
for(currNode in nodes){
  currNodeName <- totalNodeRegions_losses$name[which(totalNodeRegions_losses$Node == currNode)]
  print(currNodeName)
  totalOlderNodeRegions <- totalNodeRegions_losses$nRegions[which(totalNodeRegions_losses$Node == currNode)]
  for(i in 1:nCellTypes){
    currCellType <- totalCellTypeRegions_losses$CellType[i]
    cellTypeDF <- nodeLosses %>%
      filter(CellType == currCellType)
    nTrials <- cellTypeDF$nRegions[which(cellTypeDF$Node == currNode)]
    numerator <- totalNodeLosses$nLosses[which(totalNodeLosses$Node == currNode)] - cellTypeDF$nSuccess[which(cellTypeDF$Node == currNode)] # All success regions for node of interest - success regions in current cell type
    denominator <- totalOlderNodeRegions - nTrials  # All regions - total regions in node of interest
    nSuccesses = cellTypeDF$nSuccess[which(cellTypeDF$Node == currNode)]
    nullPr <- numerator/denominator
    nullPrMat[currCellType,currNodeName] <- nullPr
    test <- binom.test(nSuccesses, nTrials, nullPr, alternative = "two.sided")
    binomPrMat_loss[currCellType, currNodeName] <- test$estimate[[1]]
    binomPvalMat_loss[currCellType, currNodeName] <- test$p.value[[1]]
  }
}
binomSigMat_loss <- binomPvalMat_loss < 0.05
sum(binomSigMat_loss)
binomSigMat_loss <- binomSigMat_loss[, colOrder]
binomPvalMat_loss <- binomPvalMat_loss[, colOrder]
binomPrMat_loss <- binomPrMat_loss[, colOrder]

setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(binomPrMat_loss, "binomPrMat_loss.rds")
saveRDS(binomPvalMat_loss, "binomPvalMat_loss.rds")
saveRDS(nullPrMat_loss, "nullPrMat_loss.rds")

# 2.0 Compute enrichment ----
# Enrichment = (successes/trials)/pr success
lossMat <- nodeLosses %>%
  select(CellType, nSuccess, Node) %>%
  pivot_wider(
    names_from = Node,   
    values_from = nSuccess 
  ) 
lossMat <- data.frame(lossMat)
rownames(lossMat) <- lossMat$CellType
lossMat$CellType <- NULL
lossMat <- as.matrix(lossMat)
colOrder <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16")
lossMat <- lossMat[, colOrder]
colnames(lossMat) <- nodeNames_losses

nMat_loss <- nodeLosses %>%
  select(CellType, nRegions, Node) %>%
  pivot_wider(
    names_from = Node,
    values_from = nRegions
  ) 
nMat_loss <- data.frame(nMat_loss)
rownames(nMat_loss) <- nMat_loss$CellType
nMat_loss$CellType <- NULL
nMat_loss <- as.matrix(nMat_loss)
colOrder <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16")
nMat_loss <- nMat_loss[, colOrder]
colnames(nMat_loss) <- nodeNames_losses
enrichMat_loss <- (lossMat/nMat_loss)/nullPrMat

# Save
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(enrichMat_loss, "enrichMat_loss.rds")


