# 01_binomialInnovation
# Innovation = gains
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
regions <- read.csv("enrichments_forR.csv", header = T)

## OPTIONAL, add Wang 2026 placenta data ----
wang2026 <- read.csv("wang2026Placenta_gains.txt", sep = "\t")
wang2026 <- wang2026 %>% 
  filter(str_detect(Cell.Type, "400bp"))
wang2026$Cell.Type <- str_remove(wang2026$Cell.Type, "_.*")
wang2026$Node <- str_remove(wang2026$Node, "node")
wang2026$Node <- as.integer(wang2026$Node)
wang2026$Cell.Type <- paste("placental", wang2026$Cell.Type, sep = "_")
# change cell type names to match Domcke dataset
wang2026 <- wang2026 %>%
  mutate(Cell.Type = case_when(
    str_detect(Cell.Type, "Endo") ~ "endothelialPlacental",
    str_detect(Cell.Type, "HB") ~ "macrophagePlacental",
    str_detect(Cell.Type, "SCT") ~ "syncitiotrophoblastCytotrophoblast",
    str_detect(Cell.Type, "EVT") ~ "extravillousTrophoblast",
    TRUE ~ Cell.Type  # Keeps original value if no match
  ))
wang2026$Cell.Type <- paste("wang26", wang2026$Cell.Type, sep = "_")
wang2026 <- complete(wang2026, Cell.Type, Node = -1:17, fill = list(Regions = 0))
wangRegions <- wang2026 %>%
  select(Cell.Type, Regions) %>%
  group_by(Cell.Type) %>%
  summarize(nRegions = sum(Regions))
wang2026 <- full_join(wang2026, wangRegions)
names(wang2026) <- c("CellType", "Node", "Regions", "nRegions")
wang2026$mappedFraction <- wang2026$Regions/wang2026$nRegions
regions <- bind_rows(regions, wang2026)
regions$Regions <- NULL
##

## 0.4 Process data ----
## Region counts
regions$X <- NULL
regions$nMapped <- regions$nRegions*regions$mappedFraction
regions$nMapped <- regions$nMapped %>% round
totalNodeRegions <- regions %>%
  select(Node, nMapped) %>%
  group_by(Node) %>%
  summarize(nRegions = sum(nMapped))
## Cell type regions
cellTypeRegions <- regions$nRegions %>% unique()
cellTypes <- regions$CellType %>% unique()
totalCellTypeRegions <- data.frame(CellType = cellTypes, nRegions = cellTypeRegions)
## Save 
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(totalCellTypeRegions, "totalCellTypeRegions.rds")
saveRDS(totalNodeRegions, "totalNodeRegions.rds")
saveRDS(regions, "regions.rds")

# 1.0 Compute binomial probabilities ----
## 1.1 Set up output matrix ----
nNodes <- nrow(totalNodeRegions)
nCellTypes <- nrow(totalCellTypeRegions)
binomPrMat <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(binomPrMat) <- c("Human", paste("Node", seq(0,17), sep = ""))
nodeNames <- colnames(binomPrMat)
rownames(binomPrMat) <- totalCellTypeRegions$CellType
totalNodeRegions$name <- colnames(binomPrMat)
backgroundPrMat <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(backgroundPrMat) <- c("Human", paste("Node", seq(0,17), sep = ""))
rownames(backgroundPrMat) <- totalCellTypeRegions$CellType

binomPvalMat <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(binomPvalMat) <- c("Human", paste("Node", seq(0,17), sep = ""))
rownames(binomPvalMat) <- totalCellTypeRegions$CellType

## 1.2 Run binomial test for all cell types ----
nodes <- regions$Node %>% unique
for(currNode in nodes){
  currNodeName <- totalNodeRegions$name[which(totalNodeRegions$Node == currNode)]
  print(currNodeName)
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
    test <- binom.test(nSuccesses[i], nTrials, nullPrs[i], alternative = "two.sided")
    binomPrMat[currCellType, currNodeName] <- test$estimate
    binomPvalMat[currCellType, currNodeName] <- test$p.value
  }
}
binomSigMat <- binomPvalMat < 0.05
sum(binomSigMat)

# 2.0 Compute enrichment ----
# Enrichment = (successes/trials)/pr success
nMapMat <- regions %>%
  select(CellType, nMapped, Node) %>%
  pivot_wider(
    names_from = Node,   
    values_from = nMapped 
  ) 
nMapMat <- data.frame(nMapMat)
rownames(nMapMat) <- nMapMat$CellType
nMapMat$CellType <- NULL
nMapMat <- as.matrix(nMapMat)
colnames(nMapMat) <- nodeNames

nMat <- regions %>%
  select(CellType, nRegions, Node) %>%
  pivot_wider(
    names_from = Node,
    values_from = nRegions
  ) 
nMat <- data.frame(nMat)
rownames(nMat) <- nMat$CellType
nMat$CellType <- NULL
nMat <- as.matrix(nMat)
colnames(nMat) <- nodeNames
enrichMat <- (nMapMat/nMat)/backgroundPrMat

# Save output ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(enrichMat, "enrichMat.rds")
saveRDS(binomPrMat, "binomPrMat.rds")
saveRDS(binomPvalMat, "binomPvalMat.rds")
saveRDS(binomSigMat, "binomSigMat.rds")

