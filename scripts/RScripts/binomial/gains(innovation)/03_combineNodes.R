# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(tidyr)

## 0.2 Load data ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
totalCellTypeRegions <- readRDS("totalCellTypeRegions.rds")
regions <- readRDS("regions.rds")

## 0.3 Source functions ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("subsetAndBinomCellTypeList.R")


# 1.0 Combine nodes -1-3 ----
nMapMat <- regions %>%
select(CellType, nMapped, Node) %>%
  pivot_wider(
    names_from = Node,   
    values_from = nMapped 
  ) 
toCombine <- c("-1", "0", "1", "2", "3")
combinedNodes <- nMapMat[, toCombine, drop = FALSE]
H_3 <- rowSums(combinedNodes)
keep <- nMapMat[, !(colnames(nMapMat) %in% toCombine), drop = FALSE]
keep <- as.data.frame(keep)
rownames(keep) <- keep$CellType
keep$CellType <- NULL
nMapMat_combNodes <- cbind(H_3, keep)
nMapMat_combNodes <- as.matrix(nMapMat_combNodes)
nMapLong_combNodes <- as.data.frame(as.table(nMapMat_combNodes))
names(nMapLong_combNodes) <- c("CellType", "Node", "nMapped")
nMapLong_combNodes$nRegions <- rep(totalCellTypeRegions$nRegions, 15)
# Save 
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(nMapLong_combNodes, "nMapLong_combNodes.rds")

# 2.0 All cells binomial enrichment ----
## 2.1 Set up data ----
totalNodeRegions_combNodes <- nMapLong_combNodes %>%
  select(Node, nMapped) %>%
  group_by(Node) %>%
  summarize(nRegions = sum(nMapped))
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(totalNodeRegions_combNodes, "totalNodeRegions_combNodes.rds")

## 2.2 Set up output matrices ----
nNodes <- nrow(totalNodeRegions_combNodes)
nCellTypes <- nrow(totalCellTypeRegions)
binomPrMat_combNodes <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(binomPrMat_combNodes) <- c("NodesH_3", paste("Node", seq(4,17), sep = ""))
nodeNames <- colnames(binomPrMat_combNodes)
rownames(binomPrMat_combNodes) <- totalCellTypeRegions$CellType
totalNodeRegions_combNodes$name <- colnames(binomPrMat_combNodes)
nullPrMat_combNodes <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(nullPrMat_combNodes) <- c("NodesH_3", paste("Node", seq(4,17), sep = ""))
rownames(nullPrMat_combNodes) <- totalCellTypeRegions$CellType

binomPvalMat_combNodes <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(binomPvalMat_combNodes) <- c("NodesH_3", paste("Node", seq(4,17), sep = ""))
rownames(binomPvalMat_combNodes) <- totalCellTypeRegions$CellType

testOut <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(testOut) <- c("NodesH_3", paste("Node", seq(4,17), sep = ""))
rownames(testOut) <- totalCellTypeRegions$CellType

## 2.3 Compute binomial probabilities ----
Nodes <- nMapLong_combNodes$Node %>% unique
for(currNode in Nodes){
  currNodeName <- totalNodeRegions_combNodes$name[which(totalNodeRegions_combNodes$Node == currNode)]
  print(currNodeName)
  currNodeRegions <- nMapLong_combNodes %>%
    filter(Node == currNode)
  numerators <- currNodeRegions$nRegions - currNodeRegions$nMapped
  denominator <- sum(totalCellTypeRegions$nRegions) - totalNodeRegions_combNodes$nRegions[which(totalNodeRegions_combNodes$Node == currNode)]  
  nullPrs <- numerators/denominator
  nullPrMat_combNodes[,currNodeName] <- nullPrs
  nTrials = totalNodeRegions_combNodes$nRegions[which(totalNodeRegions_combNodes$Node == currNode)]
  nSuccesses = currNodeRegions$nMapped
  for(i in 1:nCellTypes){
    currCellType <- totalCellTypeRegions$CellType[i]
    test <- binom.test(nSuccesses[i], nTrials, nullPrs[i], alternative = "two.sided")
    binomPrMat_combNodes[currCellType, currNodeName] <- test$estimate
    binomPvalMat_combNodes[currCellType, currNodeName] <- test$p.value
    testOut[currCellType, currNodeName] <- (nSuccesses[i]/nTrials)/nullPrs[i]
  }
}
binomSigMat_combNodes <- binomPvalMat_combNodes < 0.05
sum(binomSigMat_combNodes)


# 3.0 Compute enrichment ----
# Enrichment = (successes/trials)/pr success
nMat_combNodes <- nMapLong_combNodes %>%
  select(CellType, nRegions, Node) %>%
  pivot_wider(
    names_from = Node,
    values_from = nRegions
  ) 
nMat_combNodes <- data.frame(nMat_combNodes)
rownames(nMat_combNodes) <- nMat_combNodes$CellType
nMat_combNodes$CellType <- NULL
nMat_combNodes <- as.matrix(nMat_combNodes)
colnames(nMat_combNodes) <- c("H_3", seq(4,17))
enrichMat_combNodes <- (nMapMat_combNodes/nMat_combNodes)/nullPrMat_combNodes
