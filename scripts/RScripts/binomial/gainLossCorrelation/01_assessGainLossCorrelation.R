# 01_assessGainLossCorrelation
# Assess the correlation between enrichment for gains and enrichment for losses
# Hailey Napier
# August 25, 2026

# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(reshape2)
require(ggplot2)
require(ggpmisc)
require(ggsignif)
require(stringr)
require(wesanderson)


## 0.2 Load data ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
enrichMat_loss <- readRDS("enrichMat_loss.rds")
enrichMat_combNodes <- readRDS("enrichMat_combNodes.rds")
binomPvalMat_combNodes <- readRDS("binomPvalMat_combNodes.rds")
binomPvalMat_loss <- readRDS("binomPvalMat_loss.rds")
MYA_losses <- c(10, 15, 20, 30, 45, 75, 85, 87, 95, 100, 160, 180, 320, 350, 415, 430)
MYA <- c(20, 30, 45, 75, 85, 87, 95, 100, 160, 180, 320, 350, 415, 430, 560)

# 1.0 Process data for scatterplot ----
## 1.1 Make each matrix into long format ----
lossEnr_long <- melt(enrichMat_loss, varnames = c("CellType", "nodeName"), value.name = "lossEnrich")
combNodeEnr_long <- melt(enrichMat_combNodes, varnames = c("CellType", "nodeName"), value.name = "gainEnrich")
lossPval_long <- melt(binomPvalMat_loss, varnames = c("CellType", "nodeName"), value.name = "lossPval")
combNodePval_long <- melt(binomPvalMat_combNodes, varnames = c("CellType", "nodeName"), value.name = "gainPval")
combNodePval_long$nodeName <- str_remove_all(combNodePval_long$nodeName, "Nodes")
combNodePval_long$nodeName <- str_remove_all(combNodePval_long$nodeName, "Node")
combNodeEnr_long$nodeName <- str_remove_all(combNodeEnr_long$nodeName, "Node")
combNodeEnr_long$nodeName <- str_remove_all(combNodeEnr_long$nodeName, "Node")
lossPval_long$nodeName <- str_remove_all(lossPval_long$nodeName, "Node")
lossEnr_long$nodeName <- str_remove_all(lossEnr_long$nodeName, "Node")
combNode_MYA <- data.frame(nodeName = combNodeEnr_long$nodeName %>% unique(), MYA = MYA)
lossNode_MYA <- data.frame(nodeName = lossEnr_long$nodeName %>% unique() %>% as.integer() %>% sort(), MYA = MYA_losses)
lossNode_MYA$nodeName <- as.character(lossNode_MYA$nodeName)
lossEnr_long <- full_join(lossEnr_long, lossNode_MYA)
lossPval_long <- full_join(lossPval_long, lossNode_MYA)
combNodeEnr_long <- full_join(combNodeEnr_long, combNode_MYA)
combNodePval_long <- full_join(combNodePval_long, combNode_MYA)

## 1.2 Merge matrices together ----
loss <- inner_join(lossPval_long, lossEnr_long)
gain <- inner_join(combNodeEnr_long, combNodePval_long)
gain$nodeName <- as.factor(gain$nodeName)
loss$nodeName <- as.factor(loss$nodeName)
gain$MYA <- as.integer(gain$MYA)
loss$MYA <- as.integer(loss$MYA)
gainLoss <- full_join(gain, loss)
gainLoss$MillionYearsAgo <- as.factor(gainLoss$MYA)

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
## Generate a 15-color continuous Zissou1 palette ----
zissou_15 <- wes_palette("Zissou1", 15, type = "continuous")
# View the hex codes
zissou_15 %>% as.character()
# Plot a quick color swatch to see them
image(1:15, 1, as.matrix(1:15), col = zissou_15, xlab = "Zissou 15-color", ylab = "", yaxt = "n")

## 3.1 All nodes ----
ggplot(gainLoss, aes(x = gainEnrich, y = lossEnrich)) +
  geom_point(alpha = 0.5, stroke = 0, size = 3) + 
  geom_smooth(method = "lm", color = "azure4", se = TRUE) +
  theme_minimal() 


allModel <- lm(gainEnrich ~ lossEnrich, data = gainLoss)
adjR2 <- summary(allModel)$adj.r.squared

allGainLoss_MYA <- ggplot(gainLoss, aes(x = gainEnrich, y = lossEnrich, color = MillionYearsAgo)) +
  geom_point(alpha = 0.7, stroke = 0, size = 3) + 
  scale_color_manual(values = rev(zissou_15[2:14])) +
  geom_smooth(method = "lm", color = "azure4", se = TRUE) +
  annotate("text", x = 2.5, y =0.75,
           label = sprintf('"Adj." ~ R^2 == "%.3f"', adjR2),
           parse = TRUE) +
  theme_minimal() + 
  labs(x = "pCRE Gain Enrichment Score", 
       y = "pCRE Loss Enrichment Score", 
       color = "MYA")

setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/figures/gainLossCorrPlots")
ggsave("allGainLossCorr.png", allGainLoss_MYA, width = 6, height = 4, bg = "transparent")


ggplot(gainLoss, aes(x = gainEnrich, y = lossEnrich, color = MillionYearsAgo)) +
  geom_point(alpha = 0.7, stroke = 0, size = 3) + 
  scale_color_manual(values = rev(zissou_15)) +
  geom_smooth(method = "lm", color = "azure2", se = TRUE) +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, label.x = "right", label.y = "top") +
  facet_wrap(~MYA) + 
  theme_dark()

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
ggplot(gainLoss_recent, aes(x = gainEnrich, y = lossEnrich, color = MillionYearsAgo)) +
  geom_point(alpha = 0.7, stroke = 0, size = 3) + 
  scale_color_manual(values = rev(zissou_15[2:14])) +
  geom_smooth(method = "lm", color = "azure4", se = TRUE) +
  theme_minimal()

## 3.3 Ancient nodes ----
ggplot(gainLoss_ancient, aes(x = gainEnrich, y = lossEnrich, color = MillionYearsAgo)) +
  geom_point(alpha = 0.7, stroke = 0, size = 3) + 
  scale_color_manual(values = rev(zissou_15[3:8])) +
  geom_smooth(method = "lm", color = "azure4", se = TRUE) +
  theme_minimal()

## 3.4 Color placenta cell types ----
placentaTypes <- c(#"placentalNeuron", 
                   #"fibroPlacental", 
                   "macrophagePlacental", 
                   "extravillousTrophoblast", 
                   "syncitiotrophoblastCytotrophoblast", 
                   "endothelialPlacental")
placentaPattern <- paste(placentaTypes, collapse = "|")
placentaTrophoblast <- c("extravillousTrophoblast", 
                         "syncitiotrophoblastCytotrophoblast")
placentaTrophoblastPattern <- paste(placentaTrophoblast, collapse = "|")
gainLoss <- gainLoss %>% 
  mutate(placenta = if_else(str_detect(CellType, placentaPattern), "TRUE", "FALSE"))
gainLoss <- gainLoss %>% 
  mutate(placentaTrophoblast = if_else(str_detect(CellType, placentaTrophoblastPattern), "TRUE", "FALSE"))


ggplot(gainLoss, aes(x = gainEnrich, y = lossEnrich, color = placenta, alpha = placenta)) + 
  geom_point(stroke = 0, size = 3) +
  #geom_smooth(method = "lm", color = "azure4", se = TRUE) + 
  scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.1)) + 
  scale_color_manual(values = c("TRUE" = "yellowgreen", "FALSE" = "azure4")) +
  ylim(c(0.5, 1.75)) +
  theme_minimal()

placentaGainLoss <- gainLoss %>%
  filter(placenta == TRUE)

placentaGainLoss$MYA <- as.factor(placentaGainLoss$MYA)
ggplot(placentaGainLoss, aes(x = gainEnrich, y = lossEnrich, color = MYA)) + 
  scale_color_manual(values = rev(zissou_15)) +
  geom_point(alpha = 0.7, stroke = 0, size = 3) +
  ylim(c(0.5, 1.75)) +
  #geom_smooth(method = "lm", color = "azure4", se = TRUE) + 
  theme_minimal() +
  labs(title = "Placenta")

ggplot(gainLoss, aes(x = gainEnrich, y = lossEnrich, color = placentaTrophoblast, alpha = placentaTrophoblast)) + 
  geom_point(stroke = 0, size = 3) +
  #geom_smooth(method = "lm", color = "azure4", se = TRUE) + 
  scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.15)) + 
  scale_color_manual(values = c("TRUE" = "darkolivegreen", "FALSE" = "darkgray")) +
  theme_minimal()

placentaTrophoblastGainLoss <- gainLoss %>%
  filter(placentaTrophoblast == TRUE)


ggplot(placentaGainLoss, aes(x = gainEnrich, y = lossEnrich, color = placentaTrophoblast, alpha = placentaTrophoblast)) + 
  geom_point(stroke = 0, size = 3) +
  #geom_smooth(method = "lm", color = "azure4", se = TRUE) + 
  scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.3)) + 
  scale_color_manual(values = c("TRUE" = "darkolivegreen", "FALSE" = "yellowgreen")) +
  ylim(c(0.5, 1.75)) +
  theme_minimal()

placentaTrophoblastGainLoss$MYA <- as.factor(placentaTrophoblastGainLoss$MYA)
ggplot(placentaTrophoblastGainLoss, aes(x = gainEnrich, y = lossEnrich, color = MYA)) + 
  scale_color_manual(values = rev(zissou_15)) +
  geom_point(alpha = 0.7, stroke = 0, size = 3) +
  ylim(c(0.5, 1.75)) +
  #geom_smooth(method = "lm", color = "azure4", se = TRUE) + 
  theme_minimal() +
  labs(title = "Placenta Trophoblasts")

## 3.5 Color placenta cell types different colors ----
# colors are based on Paul Tol's colorblind friendly muted palette
gainLoss$pointColor <- rep("NA", nrow(gainLoss))
gainLoss$legendLabel <- rep("NA", nrow(gainLoss))
for(i in 1:nrow(gainLoss)){
  if(gainLoss$placenta[i] == FALSE){
    gainLoss$pointColor[i] <- "#DDDDDD00"
  }else if(str_detect(gainLoss$CellType[i], "extravillousTrophoblast")){
    gainLoss$pointColor[i] <- "#DDAA33"
    gainLoss$legendLabel[i] <- "Trophoblast"
  }else if(str_detect(gainLoss$CellType[i], "syncitiotrophoblastCytotrophoblast")){
    gainLoss$pointColor[i] <- "#DDAA33"
    gainLoss$legendLabel[i] <- "Trophoblast"
  }else if(str_detect(gainLoss$CellType[i], "macrophagePlacental")){
    gainLoss$pointColor[i] <- "#44AA99"
    gainLoss$legendLabel[i] <- "Macrophage"
  }else if (str_detect(gainLoss$CellType[i], "endothelialPlacental")){
    gainLoss$pointColor[i] <- "#999933"
    gainLoss$legendLabel[i] <- "Endothelial Cell"
  }
}
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(gainLoss, "gainLoss.rds")

placentaPts <- filter(gainLoss, placenta == TRUE) %>%
  filter(legendLabel != "Fibroblast")
legendKey <- placentaPts %>%
  distinct(pointColor, legendLabel) %>%
  arrange(legendLabel)

placentaModel <- lm(gainEnrich ~ lossEnrich, data = placentaPts)
adjR2 <- summary(placentaModel)$adj.r.squared

allModel <- lm(gainEnrich ~ lossEnrich, data = gainLoss)
summary(allModel)

placentaGainLossCorr <- ggplot(gainLoss, aes(x = gainEnrich, y = lossEnrich)) +
  geom_point(data = filter(gainLoss, placenta == FALSE),
             size = 3, stroke = 0, color = "azure3", alpha = 0.5) +
  geom_point(data = placentaPts,
             aes(color = pointColor), size = 3) +
  geom_smooth(data = placentaPts, method = "lm", color = "azure4", fill = "azure4",se = TRUE) + 
  annotate("text", x = 2.5, y =0.75,
           label = sprintf('bold("Adj." ~ R^2 == "%.3f")', adjR2),
           parse = TRUE) +
  theme_minimal() + 
  scale_color_identity(guide = "legend",
                       breaks = legendKey$pointColor,
                       labels = legendKey$legendLabel,
                       name = "Placental Cell Type") +
  theme(axis.title.x = element_text(face = "bold", size = 12), 
        axis.text.x = element_text(face = "bold", size = 10), 
        axis.title.y = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(face = "bold", size = 10), 
        legend.title = element_text(face = "bold", size = 10), 
        legend.text = element_text(face = "bold", size = 9)) + 
  labs(x = "pCRE Gain Enrichment Score", 
        y = "pCRE Loss Enrichment Score")
  placentaGainLossCorr
setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/figures/gainLossCorrPlots")
ggsave("placentaGainLossCorr.png", placentaGainLossCorr, width = 6.7, height = 4, bg = "transparent")


# 4.0 Compute turnover magnitude score ----
## Within placental mammals for placental cell types
## 4.1 Compute turnover magnitude score ----
placentaPts$turnoverMag_mult <- placentaPts$lossEnrich * placentaPts$gainEnrich
placentaPts$turnoverMag_add <- placentaPts$lossEnrich + placentaPts$gainEnrich
placentaPts <- placentaPts %>%
  filter(!nodeName %in% c(1,2,3))
placentaPts$colMYA <- NA
for(currNode in as.numeric(placentaPts$nodeName)){
  print(currNode)
  placentaPts$colMYA[which(placentaPts$nodeName == currNode)] <- rev(zissou_15)[currNode-2]
}
placentaPts$broadCell_MYA <- paste(placentaPts$legendLabel, placentaPts$MYA, sep = "_")
placentaPts <- placentaPts %>% 
  group_by(broadCell_MYA) %>% 
  mutate(avg_turnMag_add = mean(turnoverMag_add, na.rm = TRUE)) %>% 
  ungroup()

placentaPts <- placentaPts %>% 
  group_by(broadCell_MYA) %>% 
  mutate(avg_turnMag_mult = mean(turnoverMag_mult, na.rm = TRUE)) %>% 
  ungroup()


## 4.2 Subset to include only nodes within placental mammals ----
recentPlacentaPts <- placentaPts %>%
  filter(MYA < 100) %>%
  select(c(legendLabel, avg_turnMag_add, pointColor, colMYA, MYA, avg_turnMag_mult)) %>%
  distinct()
recentPlacentaPts
placentaMYALegendKey <- distinct(recentPlacentaPts, colMYA, MYA)

placentaTurnMag_add <- ggplot(data = recentPlacentaPts, aes(y = avg_turnMag_add, x = legendLabel)) +
  geom_boxplot(aes(color = pointColor), outlier.shape = NA) + 
  geom_jitter(aes(color = colMYA), width = 0.2, size = 2.25, alpha = 0.8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_identity(guide = "legend",
                       breaks = placentaMYALegendKey$colMYA,
                       labels = placentaMYALegendKey$MYA,
                       name = "MYA") + 
  labs(x = "Placental Cell Type", 
       y = "Turnover Magnitude") + 
  geom_signif(
    comparisons = list(c("Trophoblast", "Macrophage")), 
    map_signif_level = TRUE, textsize = 6, color = "black") +
  geom_signif(
    comparisons = list(c("Endothelial Cell", "Macrophage")), 
    map_signif_level = TRUE, textsize = 4, y_position = 3.0, color = "black") +
  ylim(0, 5) + 
  theme(axis.title.x = element_text(face = "bold", size = 12), 
        axis.text.x = element_text(face = "bold", size = 10), 
        axis.title.y = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(face = "bold", size = 10), 
        legend.title = element_text(face = "bold", size = 10), 
        legend.text = element_text(face = "bold", size = 9)) 

placentaTurnMag_add
setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/figures/gainLossCorrPlots")
ggsave("placentaTurnMag_add.png", placentaTurnMag_add, width = 3.5, height = 4, bg = "transparent")



  
