# 01_binomialEnrichment
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

## 0.2 Source functions
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("stouffer_cellTypes.R")

## 0.3 Load data ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/data")
regions <- read.csv("enrichments_forR.csv", header = T)
zscoreMat <- read.csv("zscoreEnrichMat.csv", header = T)

## 0.4 Process data ----
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
nodeNames <- colnames(binomPrMat)
rownames(binomPrMat) <- totalCellTypeRegions$CellType
totalNodeRegions$name <- colnames(binomPrMat)
backgroundPrMat <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(backgroundPrMat) <- c("Human", paste("Node", seq(0,17), sep = ""))
rownames(backgroundPrMat) <- totalCellTypeRegions$CellType

binomPvalMat <- matrix(ncol = nNodes, nrow = nCellTypes)
colnames(binomPvalMat) <- c("Human", paste("Node", seq(0,17), sep = ""))
rownames(binomPvalMat) <- totalCellTypeRegions$CellType

## 1.2 Run binomial test ----
  for(currNode in combMapped$Node){
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

## 1.3 Combine neurons ----
### Excitatory ----
exNeurMap <- getCombMapped(regions, "excitatoryNeuron")
exNeur_binom <- subtypeBinomTest(exNeurMap, totalNodeRegions = totalNodeRegions, totalCellTypeRegions = totalCellTypeRegions, nodeNames)
### Inhibitory ----
inhNeurMap <- getCombMapped(regions, "inhibitoryNeuron")
inhNeur_binom <- subtypeBinomTest(inhNeurMap, totalNodeRegions = totalNodeRegions, totalCellTypeRegions = totalCellTypeRegions, nodeNames)
### Combine for plotting ----
neuron_binom <- rbind(inhNeur_binom, exNeur_binom)

## 1.4 Combine smooth muscle ----
smoothMuscleMap <- getCombMapped(regions, "smoothMuscle")
smoothMuscle_binom <- subtypeBinomTest(smoothMuscleMap, totalNodeRegions = totalNodeRegions, totalCellTypeRegions = totalCellTypeRegions, nodeNames)

## 1.5 Combine cardiomyocytes ----
### Developing ----
devCardioMap <- getCombMapped(regions, "Cardiomyocyte.developing")
devCardio_binom <- subtypeBinomTest(devCardioMap, totalNodeRegions = totalNodeRegions, totalCellTypeRegions = totalCellTypeRegions, nodeNames)
### Adult ----
adCardioMap <- getCombMapped(regions, "Cardiomyocyte.adult")
adCardio_binom <- subtypeBinomTest(adCardioMap, totalNodeRegions = totalNodeRegions, totalCellTypeRegions = totalCellTypeRegions, nodeNames)
### Combine for plotting ----
cardio_binom <- rbind(devCardio_binom, adCardio_binom)


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


# 3.0 Compute z-score p-values ----
# Does the zscore denote significance at the 95% CI
zscoreSigMat <- (zscoreMat >= 1.96) | (zscoreMat <= -1.96)
sum(zscoreSigMat)


#  4.0 Heatmaps ----
## 4.1 Binomial test estimates heatmap ----
pheatmap(binomPrMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Probabilities", 
         fontsize_row = 3)

## 4.2 Binomial probability heatmap, scaled by column----
## Scaling is Z scores based on column average
pheatmap(binomPrMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Probabilities, scaled by column", 
         fontsize_row = 3, 
         scale = "column")

## 4.3 Binomial test pvalues ----
## Scaling is Z scores based on column average
pheatmap(binomPvalMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Test P-values", 
         fontsize_row = 3)
binomSigMat <- binomSigMat*1
pheatmap(binomSigMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Test Significance at 95% Confidence", 
         fontsize_row = 3)

## 4.4 Binomial enrichments ----
pheatmap(enrichMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Enrichments", 
         fontsize_row = 3)
pheatmap(enrichMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Enrichments, Normalized by Node", 
         fontsize_row = 3, 
         scale = "column")

## 4.5 Z score heatmap ----
pheatmap(zscoreMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Z Scores", 
         fontsize_row = 3)

## 4.6 Z score significant at 95% CI ----
zscoreSigMat <- zscoreSigMat*1
pheatmap(zscoreSigMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Z Score Significance at 95% Confidence", 
         fontsize_row = 3)


# 5.0 Line plots ----
## 5.1 Z score line plots ----
### Set up data ----
# Convert matrix to long dataframe 
zscoreLong <- as.data.frame(as.table(zscoreMat))
names(zscoreLong) <- c("cellType", "node", "zscore")
# Combine Z scores using Stouffer's Method 
inhNeuronComb <- stouffer_cellTypes(zscoreLong, "inhibitoryNeuron")
exNeuronComb <- stouffer_cellTypes(zscoreLong, "excitatoryNeuron")
smoothMuscleComb <- stouffer_cellTypes(zscoreLong, "smoothMuscle")
devCardioComb <- stouffer_cellTypes(zscoreLong, "Cardiomyocyte.developing")
adCardioComb <- stouffer_cellTypes(zscoreLong, "Cardiomyocyte.adult")

### Plot ----
# Neurons
neurons <- rbind(inhNeuronComb, exNeuronComb)
labels = c("excitatoryNeuron" = "Excitatory Neuron", 
"inhibitoryNeuron" = "Inhibitory Neuron")
ggplot(neurons, aes(x = node, y = combZscore, group = 1, color = cellType)) + 
  geom_line() + 
  scale_x_discrete(limits = rev) +
  theme_minimal() + 
  geom_hline(yintercept = 1.96, linetype = "dotted") +
  geom_hline(yintercept = -1.96, linetype = "dotted") + 
  geom_hline(yintercept = 1.6, linetype = "dashed") +
  geom_hline(yintercept = -1.6, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "longdash") + 
  annotate("rect", 
           xmin = -Inf, xmax = Inf, 
           ymin = -1.6, 
           ymax = 1.6, 
           fill = "gray", 
           alpha = 0.3) + 
  scale_y_continuous(limits = c(-6, 6)) + 
  scale_color_discrete(labels = labels) +
  facet_wrap(~cellType, 
             nrow = 2, 
             ncol = 1, 
             labeller = as_labeller(labels)) + 
  labs(title = "Z Score Enrichment, Neurons", 
       color = "Cell Type", 
       x = "Node", 
       y = "Z Score") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

# Smooth muscle
ggplot(smoothMuscleComb, aes(x = node, y = combZscore, group = 1)) + 
  geom_line() + 
  scale_x_discrete(limits = rev) +
  theme_minimal() + 
  geom_hline(yintercept = 1.96, linetype = "dotted") +
  geom_hline(yintercept = -1.96, linetype = "dotted") + 
  geom_hline(yintercept = 1.6, linetype = "dashed") +
  geom_hline(yintercept = -1.6, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "longdash") + 
  annotate("rect", 
           xmin = -Inf, xmax = Inf, 
           ymin = -1.6, 
           ymax = 1.6, 
           fill = "gray", 
           alpha = 0.3) + 
  scale_y_continuous(limits = c(-6, 6)) + 
  scale_color_discrete(labels = labels) +
  labs(title = "Z Score Enrichment, Smooth Muscle", 
       x = "Node", 
       y = "Z Score") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) 

# Cardiomyocytes
cardio <- rbind(devCardioComb, adCardioComb)
labels = c("Cardiomyocyte.developing" = "Developing Cardiomyocyte", 
           "Cardiomyocyte.adult" = "Adult Cardiomyocyte")
ggplot(cardio, aes(x = node, y = combZscore, group = 1, color = cellType)) + 
  geom_line() + 
  scale_x_discrete(limits = rev) +
  theme_minimal() + 
  geom_hline(yintercept = 1.96, linetype = "dotted") +
  geom_hline(yintercept = -1.96, linetype = "dotted") + 
  geom_hline(yintercept = 1.6, linetype = "dashed") +
  geom_hline(yintercept = -1.6, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "longdash") + 
  annotate("rect", 
           xmin = -Inf, xmax = Inf, 
           ymin = -1.6, 
           ymax = 1.6, 
           fill = "gray", 
           alpha = 0.3) + 
  scale_y_continuous(limits = c(-6, 6)) + 
  scale_color_discrete(labels = labels) +
  facet_wrap(~cellType, 
             nrow = 2, 
             ncol = 1, 
             labeller = as_labeller(labels)) + 
  labs(title = "Z Score Enrichment, Cardiomyocytes", 
       color = "Cell Type", 
       x = "Node", 
       y = "Z Score") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

## 5.2 Binomial enrichment line plots ----
### Neurons ----
labels = c("excitatoryNeuron" = "Excitatory Neuron", 
           "inhibitoryNeuron" = "Inhibitory Neuron")
ggplot(neuron_binom, aes(x = nodeName, y = enrich, group = 1, color = cellType)) + 
  geom_line() + 
  scale_x_discrete(limits = rev(nodeNames)) +
  theme_minimal() + 
  geom_point(
    data = filter(neuron_binom, sig),
    aes(x = nodeName[sig], y = enrich),
    shape = 8,
    size = 2,
    color = "black", 
    position = position_nudge(y = 0.3 * diff(range(neuron_binom$enrich)))) +
  geom_hline(yintercept = 1, linetype = "longdash") + 
  scale_color_discrete(labels = labels) +
  facet_wrap(~cellType, 
             nrow = 2, 
             ncol = 1, 
             labeller = as_labeller(labels)) + 
  labs(title = "Binomial Enrichment, Neurons", 
       color = "Cell Type", 
       x = "Node", 
       y = "Enrichment") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

### Cardiomyocytes ----
labels = c("Cardiomyocyte.developing" = "Developing Cardiomyocyte", 
           "Cardiomyocyte.adult" = "Adult Cardiomyocyte")
ggplot(cardio_binom, aes(x = nodeName, y = enrich, group = 1, color = cellType)) + 
  geom_line() + 
  scale_x_discrete(limits = rev(nodeNames)) +
  theme_minimal() + 
  geom_point(
    data = filter(cardio_binom, sig),
    aes(x = nodeName[sig], y = enrich),
    shape = 8,
    size = 2,
    color = "black", 
    position = position_nudge(y = 0.3 * diff(range(cardio_binom$enrich)))) +
  geom_hline(yintercept = 1, linetype = "longdash") + 
  scale_color_discrete(labels = labels) +
  facet_wrap(~cellType, 
             nrow = 2, 
             ncol = 1, 
             labeller = as_labeller(labels)) + 
  labs(title = "Binomial Enrichment, Cardiomyocytes", 
       color = "Cell Type", 
       x = "Node", 
       y = "Enrichment") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

