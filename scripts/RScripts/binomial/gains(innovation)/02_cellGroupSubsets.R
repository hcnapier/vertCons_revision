# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(tidyr)

## 0.2 Load data ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
totalCellTypeRegions <- readRDS("totalCellTypeRegions.rds")
totalNodeRegions <- readRDS("totalNodeRegions.rds")
regions <- readRDS("regions.rds")
binomPrMat <- readRDS("binomPrMat.rds")

## 0.3 Source functions ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("subsetAndBinomCellTypeList.R")

# 1.0 Binomial enrichment, all nodes ----
nodeNames_full <- colnames(binomPrMat)
## 1.3 Neurons ----
neuronTypes <- c("excitatoryneuron", "inhibitoryneuron")
neuronGroupBinom <- subsetAndBinomCellTypeList(cellTypeList = neuronTypes, 
                                               fullDF = regions, 
                                               totalNodeRegions = totalNodeRegions,
                                               totalCellTypeRegions = totalCellTypeRegions,
                                               nodeNames = nodeNames_full, 
                                               groupName = "Neurons")
# Combine for plotting
neuron_binom <- do.call(rbind, neuronGroupBinom)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(neuron_binom, "neuron_bionm.rds")

## 1.4 All Muscle Subtypes ----
stageList <- c("developing", "adult")
muscleTypes <- c("smoothmuscle", "cardiomyocyte", "skeletalmyocyte")
muscleGroupBinom <- subsetAndBinomCellTypeList(cellTypeList = muscleTypes, 
                                               stageList = stageList, 
                                               fullDF = regions, 
                                               totalNodeRegions = totalNodeRegions, 
                                               totalCellTypeRegions = totalCellTypeRegions, 
                                               nodeNames = nodeNames_full, 
                                               groupName = "Muscle")
### Combine for plotting
muscle_binom <- do.call(rbind, muscleGroupBinom)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(muscle_binom, "muscle_binom.rds")

## 1.5 Cardiomyocytes ----
cardio_binom <- rbind(muscleGroupBinom[["cardiomyocyte_developing"]], muscleGroupBinom[["cardiomyocyte_adult"]])
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(cardio_binom, "cardio_binom.rds")

## 1.5 Innate immune cells ----
immuneTypes <- c("macrophage", "microglia", "mast", "naturalkillert")
immuneGroupBinom <- subsetAndBinomCellTypeList(cellTypeList = immuneTypes,
                                               stageList = stageList, 
                                               fullDF = regions, 
                                               totalNodeRegions = totalNodeRegions, 
                                               totalCellTypeRegions = totalCellTypeRegions, 
                                               nodeNames = nodeNames_full, 
                                               groupName = "InnateImmune")
### Combine for plotting
innateImmune_binom <- do.call(rbind, immuneGroupBinom)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(innateImmune_binom, "innateImmune_binom")
