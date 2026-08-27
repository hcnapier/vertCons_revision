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
#enrichMat_combNodes <- readRDS("enrichMat_combNodes.rds")
enrichMat_combNodes <- testOut
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
combNodePval_long <- full_join(combNodePval_long, combNode_MYA)

## 1.2 Merge matrices together ----
loss <- inner_join(lossPval_long, lossEnr_long)
gain <- inner_join(combNodeEnr_long, combNodePval_long)
gain$nodeName <- str_remove_all(gain$nodeName, "Nodes")
gain$nodeName <- str_remove_all(gain$nodeName, "Node")
gainLoss <- inner_join(gain, loss)

## 1.3 Recent node subset ----
gainLoss_recent <- gainLoss %>%
  filter(MYA < 150)

## 1.4 Ancient node subset ----
gainLoss_ancient <- gainLoss %>%
  filter(MYA > 150)

# 2.0 Calculate correlation ----
## 2.1 ALL nodes ----
model <- lm(gainEnrich ~ lossEnrich, data = gainLoss)
summary(model)

## 2.2 Recent nodes ----
model_recent <- lm(gainEnrich ~ lossEnrich, data = gainLoss_recent)
summary(model_recent)

## 2.3 Ancient nodes ----
model_ancient <- lm(gainEnrich ~ lossEnrich, data = gainLoss_ancient)
summary(model_ancient)

# 3.0 Plot ----
## 3.1 All nodes ----
ggplot(gainLoss, aes(x = gainEnrich, y = lossEnrich)) +
  geom_point(alpha = 0.5, stroke = 0, size = 3) + 
  geom_smooth(method = "lm", color = "azure4", se = TRUE) +
  theme_minimal() 

ggplot(gainLoss, aes(x = gainEnrich, y = lossEnrich, color = MYA)) +
  geom_point(alpha = 0.7, stroke = 0, size = 3) + 
  scale_color_distiller(palette = "Spectral") +
  geom_smooth(method = "lm", color = "azure4", se = TRUE) +
  theme_minimal()

ggplot(gainLoss, aes(x = gainEnrich, y = lossEnrich, color = nodeName)) +
  geom_point(alpha = 0.5, stroke = 0, size = 3) + 
  theme_minimal()

ggplot(gainLoss, aes(x = log(gainEnrich), y = log(lossEnrich))) +
  geom_point(alpha = 0.5, stroke = 0, size = 3) + 
  theme_minimal()

ggplot(gainLoss, aes(x = log(gainEnrich), y = log(lossEnrich), color = MYA)) +
  geom_point(alpha = 0.5, stroke = 0, size = 3) + 
  theme_minimal()

ggplot(gainLoss, aes(x = log(gainEnrich), y = log(lossEnrich), color = nodeName)) +
  geom_point(alpha = 0.5, stroke = 0, size = 3) + 
  theme_minimal()

## 3.2 Recent nodes ----
ggplot(gainLoss_recent, aes(x = gainEnrich, y = lossEnrich, color = MYA)) +
  geom_point(alpha = 0.7, stroke = 0, size = 3) + 
  scale_color_distiller(palette = "Spectral") +
  geom_smooth(method = "lm", color = "azure4", se = TRUE) +
  theme_minimal()

## 3.3 Ancient nodes ----
ggplot(gainLoss_ancient, aes(x = gainEnrich, y = lossEnrich, color = MYA)) +
  geom_point(alpha = 0.7, stroke = 0, size = 3) + 
  scale_color_distiller(palette = "Spectral") +
  geom_smooth(method = "lm", color = "azure4", se = TRUE) +
  theme_minimal()
