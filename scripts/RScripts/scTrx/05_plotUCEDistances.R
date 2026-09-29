## 0.1 Load packages ----
.libPaths(c("~/R/R4.6.0_packages"))
options(repos = c(CRAN = "https://cloud.r-project.org")) 
require(dplyr)
require(ggplot2)
require(ggsignif)
set.seed(42)

## 0.2 Load data ----
#setwd("/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances_pairwise")
setwd("/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances/center_pairwise")
mac <- read.csv(gzfile("pairwise_mac.csv.gz"))
endo <- read.csv(gzfile("pairwise_endothelial.csv.gz"))
troph <- read.csv(gzfile("pairwise_trophoblast.csv.gz"))

## 0.3 Load functions ----
setwd("/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("downsampleUCEDistances.R")

# 1.0 Data processing ----
## 1.1 Filter out same species comparisons ----
macFilt <- mac %>% 
  filter(species1 != species2)
endoFilt <- endo %>% 
  filter(species1 != species2)
trophFilt <- troph %>%
  filter(species1 != species2)

## 1.2 Downsample to 200 cells per species ----
macDown <- downsampleUCEDistances(macFilt, 200)
endoDown <- downsampleUCEDistances(endoFilt, 200)
trophDown <- downsampleUCEDistances(trophFilt, 200)

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
ggplot(data = uceDist, aes(x = distance)) + 
  geom_density(aes(color = color, fill = color), alpha = 0.25) + 
  scale_color_identity(guide = "legend",
                       breaks = legendKey$color,
                       labels = legendKey$cellType,
                       name = "Placental Cell Type") + 
  scale_fill_identity(guide = "legend",
                      breaks = legendKey$color,
                      labels = legendKey$cellType,
                      name = "Placental Cell Type") +
  theme_minimal() + 
  labs(y = "Density", 
       x = "UCE Embedding Distance")

