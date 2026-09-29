## 0.1 Load packages ----
.libPaths(c("~/R/R4.6.0_packages"))
options(repos = c(CRAN = "https://cloud.r-project.org")) 
require(dplyr)
require(ggplot2)
set.seed(42)

## 0.2 Load data ----
setwd("/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances_pairwise")
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
uceDist <- bind_rows(macDown, endoDown)
uceDist <- bind_rows(uceDist, trophDown)


# 2.0 BW plot ----
ggplot(data = uceDist, aes(x = cellType, y = distance)) + 
  geom_violin() 


# 3.0 Line plot -----
ggplot(data = uceDist, aes(x = distance)) + 
  geom_freqpoly(aes(color = cellType))
