# 01_zScorePlots
# Remaking Christi's ZScore plots for revisions
# Hailey Napier, August 6, 2026

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
MYA <- c(20, 30, 45, 75, 85, 87, 95, 100, 160, 180, 320, 350, 415, 430, 560)
## Zscore matrix
rownames(zscoreMat) <- zscoreMat$CellType
zscoreMat$CellType <- NULL
zscoreMat <- as.matrix(zscoreMat)

# 1.0 Compute z-score p-values ----
## All Nodes ----
# Does the zscore denote significance at the 95% CI
zscoreSigMat <- (zscoreMat >= 1.96) | (zscoreMat <= -1.96)
sum(zscoreSigMat)
# Convert matrix to long dataframe 
zscoreLong <- as.data.frame(as.table(zscoreMat))
names(zscoreLong) <- c("cellType", "node", "zscore")
# Combine Z scores using Stouffer's Method 
inhNeuronComb <- stouffer_cellTypes(zscoreLong, "inhibitoryNeuron")
exNeuronComb <- stouffer_cellTypes(zscoreLong, "excitatoryNeuron")
smoothMuscleComb <- stouffer_cellTypes(zscoreLong, "smoothMuscle")
devCardioComb <- stouffer_cellTypes(zscoreLong, "Cardiomyocyte.developing")
adCardioComb <- stouffer_cellTypes(zscoreLong, "Cardiomyocyte.adult")

## Combine nodes -1-3 ----
toCombine <- c("Human", "Node0", "Node1", "Node2", "Node3")
combinedNodes <- zscoreMat[, toCombine, drop = FALSE]
H_3 <- rowSums(combinedNodes)/sqrt(length(toCombine))
zscoreCombNodesMat <- cbind(H_3,  zscoreMat[, !(colnames(zscoreMat) %in% toCombine), drop = FALSE])
# Does the zscore denote significance at the 95% CI
zscoreSigMat_combNodes <- (zscoreCombNodesMat >= 1.96) | (zscoreCombNodesMat <= -1.96)
sum(zscoreSigMat_combNodes)
# Convert matrix to long dataframe 
zscoreLong_combNodes<- as.data.frame(as.table(zscoreCombNodesMat))
names(zscoreLong_combNodes) <- c("cellType", "node", "zscore")
# Combine Z scores using Stouffer's Method 
inhNeuronComb_combNodes <- stouffer_cellTypes(zscoreLong_combNodes, "inhibitoryNeuron", MYA)
exNeuronComb_combNodes <- stouffer_cellTypes(zscoreLong_combNodes, "excitatoryNeuron", MYA)
smoothMuscleComb_combNodes <- stouffer_cellTypes(zscoreLong_combNodes, "smoothMuscle", MYA)
devCardioComb_combNodes <- stouffer_cellTypes(zscoreLong_combNodes, "Cardiomyocyte.developing", MYA)
adCardioComb_combNodes <- stouffer_cellTypes(zscoreLong_combNodes, "Cardiomyocyte.adult", MYA)


# 2.0 Heatmaps ----
## 2.1 ALL NODES ----
### Z score ----
pheatmap(zscoreMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Z Scores (all nodes and cell types)", 
         fontsize_row = 3)

### Significant at 95% CI ----
zscoreSigMat <- zscoreSigMat*1
pheatmap(zscoreSigMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Z Score Significance at 95% Confidence (all nodes and cell groups)", 
         fontsize_row = 3)

## 2.2 Combine nodes -1-3 ----
### Z score ----
pheatmap(zscoreCombNodesMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Z Scores", 
         fontsize_row = 3)

### Significant at 95% CI ----
zscoreSigMat_combNodes <- zscoreSigMat_combNodes*1
pheatmap(zscoreSigMat_combNodes, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Z Score Significance at 95% Confidence", 
         fontsize_row = 3)


# 3.0 Line plots, all nodes ----
## 3.1 Neurons ----
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

## 3.2 Smooth muscle ----
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

## 3.3 Cardiomyocytes ----
cardio <- rbind(devCardioComb_combNodes, adCardioComb_combNodes)
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

# 4.0 Line Plots, Nodes 1-3 Combined ----
## 4.1 Neurons ----
neurons_combNodes <- rbind(inhNeuronComb_combNodes, exNeuronComb_combNodes)
labels = c("excitatoryNeuron" = "Excitatory Neuron", 
           "inhibitoryNeuron" = "Inhibitory Neuron")
ggplot(neurons_combNodes, aes(x = MYA, y = combZscore, group = 1, color = cellType)) + 
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
  geom_line() + 
  theme_minimal() + 
  scale_y_continuous(limits = c(-6, 6)) + 
  scale_color_discrete(labels = labels) +
  facet_wrap(~cellType, 
             nrow = 2, 
             ncol = 1, 
             labeller = as_labeller(labels)) + 
  labs(title = "Z Score Enrichment, Neurons", 
       color = "Cell Type", 
       x = "Millions of Years Ago (MYA)", 
       y = "Z Score") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

## 4.2 Smooth muscle ----
ggplot(smoothMuscleComb_combNodes, aes(x = MYA, y = combZscore, group = 1)) + 
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
  geom_line() + 
  theme_minimal() + 
  scale_y_continuous(limits = c(-6, 6)) + 
  scale_color_discrete(labels = labels) +
  labs(title = "Z Score Enrichment, Smooth Muscle", 
       x = "Millions of Years Ago (MYA)", 
       y = "Z Score") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) 

## 4.3 Cardiomyocytes ----
cardio_combNodes <- rbind(devCardioComb_combNodes, adCardioComb_combNodes)
labels = c("Cardiomyocyte.developing" = "Developing Cardiomyocyte", 
           "Cardiomyocyte.adult" = "Adult Cardiomyocyte")
ggplot(cardio_combNodes, aes(x = MYA, y = combZscore, group = 1, color = cellType)) + 
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
  geom_line() + 
  theme_minimal() + 
  scale_y_continuous(limits = c(-6, 6)) + 
  scale_color_discrete(labels = labels) +
  facet_wrap(~cellType, 
             nrow = 2, 
             ncol = 1, 
             labeller = as_labeller(labels)) + 
  labs(title = "Z Score Enrichment, Cardiomyocytes", 
       color = "Cell Type", 
       x = "Millions of Years Ago (MYA)", 
       y = "Z Score") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

