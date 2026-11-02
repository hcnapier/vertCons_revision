# 1.0 Setup ----
## 0.1 Load packages ----
.libPaths(c("~/R/R4.6.0_packages"))
options(repos = c(CRAN = "https://cloud.r-project.org")) 
require(dplyr)
require(ggplot2)
require(ggsignif)
require(ape)
set.seed(42)

## 0.2 Load data ----
#setwd("/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances_pairwise")
setwd("/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances/center_pairwise")
mac <- read.csv(gzfile("pairwise_mac.csv.gz"))
endo <- read.csv(gzfile("pairwise_endothelial.csv.gz"))
troph <- read.csv(gzfile("pairwise_trophoblast.csv.gz"))
setwd("/work/hcn4/260630_vertCons_wd/scTrx/uce_distances")
subtypePairwise <- read.csv(gzfile("261005_center_pairwise.csv.gz"))

## 0.3 Load functions ----
setwd("/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("downsampleUCEDistances.R")
source("plotPairwiseUCEDist.R")

# 1.0 Data processing ----
## 1.1 Filter out same species comparisons ----
macFilt <- mac %>% 
  filter(species1 != species2)
endoFilt <- endo %>% 
  filter(species1 != species2)
trophFilt <- troph %>%
  filter(species1 != species2)

## 1.2 Downsample to 200 cells per species ----
macDown <- downsampleUCEDistances(macFilt, 100)
endoDown <- downsampleUCEDistances(endoFilt, 100)
trophDown <- downsampleUCEDistances(trophFilt, 100)

## 1.3 Merge into one dataframe ----
macDown$cellType <- "Macrophage"
endoDown$cellType <- "Endothelial Cell"
trophDown$cellType <- "Trophoblast"
macDown$color <- "#88CCEE"
endoDown$color <- "#999933"
trophDown$color <- "#DDAA33"
uceDist <- bind_rows(macDown, endoDown)
uceDist <- bind_rows(uceDist, trophDown)
legendKey <- uceDist %>%
  distinct(color, cellType) %>%
  arrange(cellType)


# 2.0 Get difference in means ----
anova_model <- aov(distance ~ cellType, data = uceDist)
summary(anova_model) # there is a difference in the means 
TukeyHSD(anova_model) # each cell type mean is different from the others


# 3.0 BW plot ----
uceDist_bwPlot <- ggplot(data = uceDist, aes(x = cellType, y = distance)) + 
  geom_violin(aes(color = color, fill = color), alpha = 0.1) + 
  geom_boxplot(aes(color = color), width = 0.5) + 
  scale_color_identity() + 
  scale_fill_identity() + 
  theme_minimal() + 
  geom_signif(
    comparisons = list(c("Trophoblast", "Macrophage")), 
    map_signif_level = TRUE, textsize = 6, color = "azure4", y_position = 1.85) +
  geom_signif(
    comparisons = list(c("Endothelial Cell", "Macrophage")), 
    map_signif_level = TRUE, textsize = 6, y_position = 1.7, color = "azure4") +
  labs(y = "UCE Embedding Distance", 
       x = "Placental Cell Type") + 
  ylim(0, 2.1) + 
  theme(text = element_text(family = "Helvetica"))
uceDist_bwPlot
setwd("/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/figures/fig4")
ggsave("uceDist_bwPlot.png", uceDist_bwPlot, width = 6, height = 4, bg = "transparent")


# 4.0 Line plot -----
uceMeans <- uceDist %>%
  group_by(color) %>%
  summarise(meanDist = mean(distance, na.rm = TRUE)) 

# Pairwise tests on mean distance between colors
pw <- pairwise.wilcox.test(uceDist$distance, uceDist$color,
                           p.adjust.method = "BH")

sigDF <- as.data.frame(as.table(pw$p.value)) %>%
  filter(!is.na(Freq)) %>%
  transmute(g1 = as.character(Var1), g2 = as.character(Var2), p = Freq) %>%
  mutate(stars = case_when(p < 0.001 ~ "***",
                           p < 0.01  ~ "**",
                           p < 0.05  ~ "*",
                           TRUE      ~ "ns")) %>%
  left_join(uceMeans, by = c("g1" = "color")) %>% rename(x1 = meanDist) %>%
  left_join(uceMeans, by = c("g2" = "color")) %>% rename(x2 = meanDist)

# Place brackets above the tallest density curve
maxDens <- uceDist %>%
  group_by(color) %>%
  summarise(peak = max(density(distance, na.rm = TRUE)$y)) %>%
  pull(peak) %>% max()
step <- 0.08 * maxDens

sigDF <- sigDF %>%
  mutate(xmin = pmin(x1, x2), xmax = pmax(x1, x2)) %>%
  arrange(xmax - xmin) %>%                 # short brackets lowest
  mutate(y = maxDens + step * row_number())

ggplot(data = uceDist, aes(x = distance)) + 
  geom_density(aes(color = color, fill = color), alpha = 0.25) + 
  geom_vline(data = uceMeans,
             aes(xintercept = meanDist, color = color),
             linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  # brackets: horizontal bar + end ticks + stars
  geom_segment(data = sigDF, aes(x = xmin, xend = xmax, y = y, yend = y),
               inherit.aes = FALSE) +
  geom_text(data = sigDF, aes(x = (xmin + xmax) / 2, y = y, label = stars),
            vjust = -0.1, size = 5, inherit.aes = FALSE) +
  scale_color_identity(guide = "legend",
                       breaks = legendKey$color,
                       labels = legendKey$cellType,
                       name = "Placental Cell Type") + 
  scale_fill_identity(guide = "legend",
                      breaks = legendKey$color,
                      labels = legendKey$cellType,
                      name = "Placental Cell Type") +
  theme_minimal() + 
  labs(y = "Density", x = "UCE Embedding Distance", title = "All Pairwise Embedding Distances") + 
  theme(axis.title.x = element_text(face = "bold", size = 12),   
        axis.title.y = element_text(face = "bold", size = 12), 
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"), 
        axis.text.x = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5)) 


## 4.1 human-macaque pairwise ----
plotPairwiseUCEDist("human", "macaque")

## 4.2 human-rat pairwise ----
plotPairwiseUCEDist("human", "rat")

## 4.3 human-mouse pairwise ----
plotPairwiseUCEDist("human", "mouse")

## 4.4 human-rabbit pairwise ----
plotPairwiseUCEDist("human", "rabbit")

## 4.4 human-rabbit pairwise ----
plotPairwiseUCEDist("human", "goat")

## 4.4 human-rabbit pairwise ----
plotPairwiseUCEDist("human", "dog")


# 5.0 Dendrogram ----
# Phylogenetic dendrogram: human, macaque, rat, mouse, rabbit, goat, dog
# Divergence times (million years ago) are approximate median estimates from TimeTree (timetree.org).
newick <- "(((Human:29,Macaque:29):61,((Mouse:12,Rat:12):70,Rabbit:82):8):6,((Goat:62,Pig:62):14,Dog:76):20);"
tree <- read.tree(text = newick)

plot_tree <- function() {
  par(mar = c(1, 1, 1, 1))
  # Bare dendrogram: no tip labels, node labels, axis, title or legend
  plot(tree, type = "phylogram", direction = "rightwards",
       show.tip.label = FALSE, edge.width = 2)
}

plot_tree()


# 6.0 Plot turnover magnitude vs embedding distance ----
## 6.1 Compute turnover magnitude score ----
setwd("/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/data/dataToPlot")
gainLoss <- readRDS("gainLoss.rds")
### Within placental mammals for placental cell types
### Compute turnover magnitude score ----
placentaPts <- filter(gainLoss, placenta == TRUE)
placentaPts$turnoverMag_mult <- placentaPts$lossEnrich * placentaPts$gainEnrich
placentaPts$turnoverMag_add <- placentaPts$lossEnrich + placentaPts$gainEnrich
placentaPts$colMYA <- NA
for(currNode in as.numeric(placentaPts$nodeName)){
  placentaPts$colMYA[which(placentaPts$nodeName == currNode)] <- rev(zissou_15)[currNode]
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

### Subset to include only nodes within placental mammals ----
recentPlacentaPts <- placentaPts %>%
  filter(MYA < 100) %>%
  select(c(legendLabel, avg_turnMag_add, pointColor, colMYA, MYA, avg_turnMag_mult)) %>%
  distinct()

turnoverMeans <- recentPlacentaPts %>%
  group_by(legendLabel) %>%
  summarise(meanTurnover = mean(avg_turnMag_add, na.rm = TRUE))
names(turnoverMeans) <- c("cellType", "meanTurnover")
turnoverMeans <- turnoverMeans %>%
  filter(cellType %in% c("Endothelial Cell", "Macrophage", "Trophoblast"))

## 6.2 Compute UCE distances ----
cellTypeMeans <- uceDist %>%
  group_by(cellType) %>%
  summarise(meanDist = mean(distance, na.rm = TRUE))

## 6.3 Merge ----
turnoverDist <- full_join(cellTypeMeans, turnoverMeans)
turnoverDist$ptColor <- c("#999933", "#88CCEE", "#DDAA33")

### LM ----
model <- lm(meanDist ~ meanTurnover, data = turnoverDist)
adjR2 <- summary(model)$adj.r.squared

## 6.4 Plot ----
legendKey <- turnoverDist %>%
  distinct(ptColor, cellType) %>%
  arrange(cellType)

ggplot(data = turnoverDist) + 
  geom_smooth(data = turnoverDist, aes(x = meanTurnover, y = meanDist), method = "lm", color = "azure3", fill = "azure3",se = TRUE) + 
  annotate("text", x = 3.5, y =0.6,
           label = sprintf('bold("Adj." ~ R^2 == "%.3f")', adjR2),
           parse = TRUE) +
  geom_point(aes(x = meanTurnover, y = meanDist, color = ptColor), size = 5) + 
  theme_minimal() + 
  labs(x = "Mean pCRE Turnover within Placental Mammals", 
       y = "Mean UCE Embedding Distance") +
  scale_color_identity(guide = "legend",
                       breaks = legendKey$ptColor,
                       labels = legendKey$cellType,
                       name = "Placental Cell Type") +
  theme(axis.title.x = element_text(face = "bold", size = 12),   
        axis.title.y = element_text(face = "bold", size = 12), 
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"), 
        axis.text.x = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5)) 
  
