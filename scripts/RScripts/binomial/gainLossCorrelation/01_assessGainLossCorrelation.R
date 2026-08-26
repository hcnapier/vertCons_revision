# 01_assessGainLossCorrelation
# Assess the correlation between enrichment for gains and enrichment for losses
# Hailey Napier
# August 25, 2026

# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(reshape2)
require(ggplot2)

## 0.2 Load data ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
enrichMat_loss <- readRDS("enrichMat_loss.rds")
enrichMat_combNodes <- readRDS("enrichMat_combNodes.rds")
binomPvalMat_combNodes <- readRDS("binomPvalMat_combNodes.rds")
binomPvalMat_loss <- readRDS("binomPvalMat_loss.rds")


# 1.0 Process data for scatterplot ----
## 1.1 Make each matrix into long format ----
lossEnr_long <- melt(enrichMat_loss, varnames = c("CellType", "nodeName"), value.name = "lossEnrich")
combNodeEnr_long <- melt(enrichMat_combNodes, varnames = c("CellType", "nodeName"), value.name = "gainEnrich")
lossPval_long <- melt(binomPvalMat_loss, varnames = c("CellType", "nodeName"), value.name = "lossPval")
combNodePval_long <- melt(binomPvalMat_combNodes, varnames = c("CellType", "nodeName"), value.name = "gainPval")
lossPval_long$nodeName <- str_remove_all(lossPval_long$nodeName, "Node")
lossEnr_long$nodeName <- str_remove_all(lossEnr_long$nodeName, "Node")
combNode_MYA <- data.frame(nodeName = combNodeEnr_long$nodeName %>% unique(), MYA = MYA)
lossNode_MYA <- data.frame(nodeName = lossEnr_long$nodeName %>% unique(), MYA = MYA_losses)
lossEnr_long <- full_join(lossEnr_long, lossNode_MYA)
lossPval_long <- full_join(lossPval_long, lossNode_MYA)
combNodeEnr_long <- full_join(combNodeEnr_long, combNode_MYA)
combNodePval_long$nodeName <- str_remove_all(combNodePval_long$nodeName, "Nodes")
combNodePval_long$nodeName <- str_remove_all(combNodePval_long$nodeName, "Node")
combNodePval_long <- full_join(combNodePval_long, combNode_MYA)

## 1.2 Merge matrices together ----
loss <- inner_join(lossPval_long, lossEnr_long)
gain <- inner_join(combNodeEnr_long, combNodePval_long)
gainLoss <- inner_join(gain, loss)

## 1.3 Normalize by node ----
gainLoss <- gainLoss %>%
  group_by(MYA) %>%
  mutate(gainZScore = scale(gainEnrich))
gainLoss <- gainLoss %>%
  group_by(MYA) %>%
  mutate(lossZScore = scale(lossEnrich))

## 1.3 Plot ----
ggplot(gainLoss, aes(x = gainZScore, y = lossZScore)) +
  geom_point(alpha = 0.5, stroke = 0, size = 3) + 
  theme_minimal()

ggplot(gainLoss, aes(x = log(gainZScore), y = lossZScore)) +
  geom_point(alpha = 0.5, stroke = 0, size = 3) + 
  theme_minimal()

ggplot(gainLoss, aes(x = gainZScore, y = log(lossZScore))) +
  geom_point(alpha = 0.5, stroke = 0, size = 3) + 
  theme_minimal()

ggplot(gainLoss, aes(x = log(gainZScore), y = log(lossZScore))) +
  geom_point(alpha = 0.5, stroke = 0, size = 3) + 
  theme_minimal()

# 2.0 Calculate correlation ----
