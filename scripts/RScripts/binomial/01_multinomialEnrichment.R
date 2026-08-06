# 01_binomialEnrichment
# Christi's original calculations carry a Gaussian assumption. A Gaussian may not fit the data well, but a binomial may be more appropriate. 
# Calculate the probabilities that a given node has enrichment for each cell type option relative to the baseline enrichment (= total number of regions for that cell type over all cell type regions, calculated with the current node held out)

# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
install.packages("pheatmap")
require(pheatmap)
require(tidyr)
## 0.2 Load data ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/data")
regions <- read.csv("enrichments_forR.csv", header = T)
zscoreMat <- read.csv("zscoreEnrichMat.csv", header = T)

## 0.3 Process data ----
## Region counts
regions$X <- NULL
regions$nMapped <- regions$nRegions*regions$mappedFraction
regions$nMapped <- regions$nMapped %>% round
totalNodeRegions <- regions %>%
  select(Node, nMapped) %>%
  group_by(Node) %>%
  summarize(nRegions = sum(nMapped))
## Zscore matrix
rownames(zscoreMat) <- zscoreMat$CellType
zscoreMat$CellType <- NULL
zscoreMat <- as.matrix(zscoreMat)


# 1.0 Compute binomial probabilities ----
## 1.1 Set up output matrix ----
nNodes <- nrow(totalNodeRegions)
nCellTypes <- nrow(totalCellTypeRegions)
binomPrMat <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(binomPrMat) <- c("Human", paste("Node", seq(0,17), sep = ""))
rownames(binomPrMat) <- totalCellTypeRegions$CellType
totalNodeRegions$name <- colnames(binomPrMat)
backgroundPrMat <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(backgroundPrMat) <- c("Human", paste("Node", seq(0,17), sep = ""))
rownames(backgroundPrMat) <- totalCellTypeRegions$CellType

## 1.2 Compute probabilities
for(currNode in totalNodeRegions$Node){
  currNodeName <- totalNodeRegions$name[which(totalNodeRegions$Node == currNode)]
  currNodeRegions <- regions %>%
    filter(Node == currNode)
  numerators <- currNodeRegions$nRegions - currNodeRegions$nMapped
  denominator <- sum(totalCellTypeRegions$nRegions) - totalNodeRegions$nRegions[which(totalNodeRegions$Node == currNode)]  
  nullPrs <- numerators/denominator
  backgroundPrMat[,currNodeName] <- nullPrs
  nTrials = totalNodeRegions$nRegions[which(totalNodeRegions$Node == currNode)]
  nSuccesses = currNodeRegions$nMapped
  for(i in 1:nCellTypes){
    currCellType <- totalCellTypeRegions$CellType[i]
    binomPrMat[currCellType, currNodeName] <- dbinom(nSuccesses[i], size = nTrials, nullPrs[i])
  }
}

## 1.3 Calculate confidence intervals for binomial probabilities ----
# Using Wald method to approximate 95% confidence intervals. Interdependence of probabilities for each group makes this just an estimate. 
# Make a matrix containing the number of trials for each node for each cell type
nMat <- regions %>%
  select(CellType, nRegions, Node) %>%
  pivot_wider(
    names_from = Node,    # Column containing the new column names
    values_from = nRegions     # Column containing the values to fill cells
  ) 
nMat <- data.frame(nMat)
rownames(nMat) <- nMat$CellType
nMat$CellType <- NULL
nMat <- as.matrix(nMat)

# Calculate confidence intervals
upperBinPrCI <- sqrt((binomPrMat*(1-binomPrMat))/nMat) * 1.96 + binomPrMat
lowerBinPrCI <- sqrt((binomPrMat*(1-binomPrMat))/nMat) * 1.96 - binomPrMat

# Calculate background confidence intervals
upperBackPrCI <- sqrt((backgroundPrMat*(1-backgroundPrMat))/nMat) * 1.96 + backgroundPrMat
lowerBackPrCI <- sqrt((backgroundPrMat*(1-backgroundPrMat))/nMat) * 1.96 - backgroundPrMat

# Check to see where they overlap
overlap <- (lowerBinPrCI <= upperBackPrCI) & (lowerBackPrCI <= upperBinPrCI)
sum(!overlap)


# 2.0 Compute enrichment ----



# 3.0 Compute z-score p-values ----
# Does the zscore denote significance at the 95% CI
zscoreSigMat <- (zscoreMat >= 1.96) | (zscoreMat <= -1.96)
sum(zscoreSigMat)

# 2.0 Visualize matrix ----
## 2.1 Basic binomial probability heatmap ----
pheatmap(binomPrMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Probabilities", 
         fontsize_row = 3)

## 2.2 Binomial probability heatmap, scaled by column----
## Scaling is Z scores based on column average
pheatmap(binomPrMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Probabilities, scaled by column", 
         fontsize_row = 3, 
         scale = "column")

## 2.2 Z score heatmap ----
pheatmap(zscoreMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Z Scores", 
         fontsize_row = 3)

