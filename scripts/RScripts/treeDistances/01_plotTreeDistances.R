# 01_plotTreeDistances
# Hailey Napier
# 9/18/2026
# To get tree distances, first run shellScripts/getTreeDists using the pruned 100-way alignment from UCSC.
# Note that distances will measured by substitution rate. 

# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(ggplot2)
require(stringr)
require(wesanderson)

## 0.2 Get data ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/data")
cov <- read.csv("genomeCoverage.csv")
speciesNodes <- read.csv("losses_speciesPerNode.csv")


# 1.0 Compute coverage and distance correlation ----
## 1.1 Set up data ----
# fix species left out 
tmp <- data.frame(Node = c("0", "5", "9", "9", "9", "11", "17"), Species = c("panTro6", "calJac3", "felCat8", "musFur1", "ailMel1", "macEug2", "petMar3"))
tmp$Node <- as.integer(tmp$Node)
speciesNodes <- rows_append(speciesNodes, tmp)
# remove percent signs
cov$hg38.Coverage <- as.numeric(gsub("%", "", cov$hg38.Coverage))
# remove human-human comparison
cov <- cov %>% 
  filter(Species.name != "human")
# convert species names to sentence case
cov$Species.name <- str_to_sentence(cov$Species.name)
# species name levels 
cov$Species.name <- factor(cov$Species.name, levels = cov$Species.name)
cov$Species <- cov$Assembly
cov$Assembly <- factor(cov$Assembly, levels = cov$Assembly)
# include node and MYA data 
cov <- left_join(cov, speciesNodes)
cov$Species <- NULL
MYA <- c(4, 6, 8, 20, 30, 45, 75, 85, 87, 95, 100, 160, 180, 320, 350, 415, 430, 560)

# Generate a 15-color continuous Zissou1 palette
zissou_18 <- wes_palette("Zissou1", 18, type = "continuous")
zissou_18[1]
cov$col <- NA
cov$MYA <- NA
for(currNode in cov$Node){
  cov$col[which(cov$Node == currNode)] <- zissou_18[18-currNode]
  cov$MYA[which(cov$Node == currNode)] <- MYA[currNode + 1]
}

## 1.2 Compute linear model 
model <- lm(log(hg38.Coverage) ~ Distance.from.hg38..substitution.rate., data = cov)
summary(model)
exp(coef(model)[1])
coef(model)[2] 

# 2.0 Plot -----
## 2.1 Bar plot: coverage by species ----
ggplot(cov, aes(y = hg38.Coverage, x = Assembly, fill = col)) + 
  geom_col() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(x = "Assembly", 
       y = "Hg38 Coverage") + 
  scale_fill_identity( guide = "legend",
                       name = "Node",
                       breaks = cov$col,
                       labels = cov$Node)

## 2.2 Scatterplot: coverage by distance ----
# Log transform y axis
ggplot(cov, aes(y = log(hg38.Coverage), x = Distance.from.hg38..substitution.rate., color = col)) + 
  geom_point(size = 4) + 
  scale_color_identity( guide = "legend",
                       name = "Node",
                       breaks = cov$col,
                       labels = cov$Node) + 
  theme_minimal() +
  labs(x = "Distance from hg38 (Substitution Rate)", 
       y = "Log(hg38 Coverage)")

# Raw y axis 
ggplot(cov, aes(y = hg38.Coverage, x = Distance.from.hg38..substitution.rate., color = col)) + 
  geom_point(size = 4) + 
  scale_color_identity( guide = "legend",
                        name = "Node",
                        breaks = cov$col,
                        labels = cov$Node) + 
  theme_minimal() +
  labs(x = "Distance from hg38 (Substitution Rate)", 
       y = "hg38 Coverage")
  
