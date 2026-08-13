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
zscoreLosses <- read.csv("nodePlasticityZScoreMat.csv", header = T)
rownames(zscoreLosses) <- zscoreLosses$CellType
colnames(zscoreLosses) <- str_to_sentence(colnames(zscoreLosses))
zscoreLosses$Celltype <- NULL
zscoreLosses <- as.matrix(zscoreLosses)
colOrder <- c("Node1", "Node2", "Node3", "Node4", "Node5", "Node6", "Node7", "Node8", "Node9", "Node10", "Node11", "Node12", "Node13", "Node14", "Node15", "Node16")
zscoreLosses <- zscoreLosses[,colOrder]
MYA_losses <- c(10, 15, 20, 30, 45, 75, 85, 87, 95, 100, 160, 180, 320, 350, 415, 430)

# 0.3 Source functions ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("stouffer_cellTypes.R")


# 1.0 Z Score calculations ----
## 1.1 Compute pvalues from zscores ----
# Does the zscore denote significance at the 95% CI
zscoreSigMat_loss <- (zscoreLosses >= 1.96) | (zscoreLosses <= -1.96)
sum(zscoreSigMat_loss)
# Convert matrix to long dataframe 
zscoreLong_loss <- as.data.frame(as.table(zscoreLosses))
names(zscoreLong_loss) <- c("cellType", "node", "zscore")

# Combine Z scores using Stouffer's Method 
z_loss_inhNeuronComb <- stouffer_cellTypes(zscoreLong_loss, "inhibitoryNeuron", MYA = MYA_losses)
z_loss_exNeuronComb <- stouffer_cellTypes(zscoreLong_loss, "excitatoryNeuron", MYA = MYA_losses)
z_loss_smoothMuscleComb <- stouffer_cellTypes(zscoreLong_loss, "smoothMuscle", MYA = MYA_losses)
z_loss_devCardioComb <- stouffer_cellTypes(zscoreLong_loss, "Cardiomyocyte.developing", MYA = MYA_losses)
z_loss_adCardioComb <- stouffer_cellTypes(zscoreLong_loss, "Cardiomyocyte.adult", MYA = MYA_losses)


# 2.0 Heatmaps ----
## 2.1 ALL NODES ----
### Z score ----
pheatmap(zscoreLosses, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Loss Z Scores (all nodes and cell types)", 
         fontsize_row = 3)

### Significant at 95% CI ----
zscoreSigMat_loss <- zscoreSigMat_loss*1
pheatmap(zscoreSigMat_loss, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Loss Z Scores Significant at 95% Confidence", 
         fontsize_row = 3)


# 3.0 Line plots, all nodes ----
## 3.1 Neurons ----
z_loss_neurons <- rbind(z_loss_inhNeuronComb, z_loss_exNeuronComb)
labels = c("excitatoryNeuron" = "Excitatory Neuron", 
           "inhibitoryNeuron" = "Inhibitory Neuron")
ggplot(z_loss_neurons, aes(x = MYA, y = combZscore, group = 1, color = cellType)) + 
  geom_line() + 
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
  scale_x_reverse() +
  scale_color_discrete(labels = labels) +
  facet_wrap(~cellType, 
             nrow = 2, 
             ncol = 1, 
             labeller = as_labeller(labels)) + 
  labs(title = "Loss Z Scores, Neurons", 
       color = "Cell Type", 
       x = "Millions of Years Ago (MYA)", 
       y = "Z Score") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

## 3.2 Smooth muscle ----
ggplot(z_loss_smoothMuscleComb, aes(x = MYA, y = combZscore, group = 1)) + 
  geom_line() + 
  scale_x_reverse() +
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
  labs(title = "Loss Z Scores, Smooth Muscle", 
       x = "Millions of Years Ago (MYA)", 
       y = "Z Score") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) 

## 3.3 Cardiomyocytes ----
cardio <- rbind(z_loss_devCardioComb, z_loss_adCardioComb)
labels = c("Cardiomyocyte.developing" = "Developing Cardiomyocyte", 
           "Cardiomyocyte.adult" = "Adult Cardiomyocyte")
ggplot(cardio, aes(x = MYA, y = combZscore, group = 1, color = cellType)) + 
  geom_line() + 
  scale_x_reverse() +
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
  labs(title = "Loss Z Scores, Cardiomyocytes", 
       color = "Cell Type", 
       x = "Millions of Years Ago (MYA)", 
       y = "Z Score") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
