# 01 binomialEnrichment
# Response to reviewer 2 comment regarding validity of Gaussian background
# Hailey Napier
# June 26, 2026

# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(pheatmap)

## 0.2 Load & process data ----
setwd("/hpc/group/vertgenlab/hailey/vertCons/data")
enrich <- read.csv("enrichments_forR.csv")
enrich$X <- NULL
### Add column for number of mapping regions ----
enrich$nMapped <- enrich$nRegions*enrich$mappedFraction
enrich$nMapped <- enrich$nMapped %>% round()


# 1.0 Take 1 ----
## 1.1 Binomial probability ----
# probability of getting exactly x mapped regions in n total accessible regions
# Compares across nodes for the same cell type
enrich$binProb <- dbinom(enrich$nMapped, enrich$nRegions, prob=1/enrich$nRegions)

## 1.2 Plot probabilities ----
enrichMat <- tapply(enrich$binProb, list(enrich$CellType, enrich$Node), identity)
pheatmap(enrichMat, 
         display_numbers = F, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_rownames = F, 
         show_colnames = T)


# 2.0 Take 2 ----
# I want to be able to compare cell types across all nodes, BUT...
## Each cell type has a different number of accessible regions, so the total number isn't comparable
## this also means the probability of each individual success is not the same


