# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(tidyr)

## 0.2 Load data ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
totalCellTypeRegions_losses <- readRDS("totalCellTypeRegions_losses.rds")
totalNodeRegions_losses <- readRDS("totalNodeRegions_losses.rds")
nodeLosses <- readRDS("nodeLosses.rds")
binomPrMat_loss <- readRDS("binomPrMat_loss.rds")

## 0.3 Source functions ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("subsetAndBinomCellTypeList.R")

# 1.0 Binomial enrichment, all nodes ----
nodeNames_losses <- c("Node1", "Node2", "Node3", "Node4", "Node5", "Node6", "Node7", "Node8", "Node9", "Node10", "Node11", "Node12", "Node13", "Node14", "Node15", "Node16")

## 1.1 Neurons ----
neuronTypes <- c("excitatoryneuron", "inhibitoryneuron")
neuronGroupBinom_loss <- subsetAndBinomCellTypeList(cellTypeList = neuronTypes, 
                                               fullDF = nodeLosses, 
                                               totalNodeRegions = totalNodeRegions_losses,
                                               totalCellTypeRegions = totalCellTypeRegions_losses,
                                               nodeNames = nodeNames_losses, 
                                               groupName = "Neurons", 
                                               MYA = MYA_losses)
# Combine for plotting
neuron_binom_loss <- do.call(rbind, neuronGroupBinom_loss)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(neuron_binom, "neuron_bionm_loss.rds")

## 1.4 All Muscle Subtypes ----
stageList <- c("developing", "adult")
muscleTypes <- c("smoothmuscle", "cardiomyocyte", "skeletalmyocyte")
muscleGroupBinom_loss <- subsetAndBinomCellTypeList(cellTypeList = muscleTypes, 
                                               stageList = stageList, 
                                               fullDF = nodeLosses, 
                                               totalNodeRegions = totalNodeRegions_losses, 
                                               totalCellTypeRegions = totalCellTypeRegions_losses, 
                                               nodeNames = nodeNames_losses, 
                                               groupName = "Muscle", 
                                               MYA = MYA_losses)
### Combine for plotting
muscle_binom_loss <- do.call(rbind, muscleGroupBinom_loss)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(muscle_binom, "muscle_binom_loss.rds")

## 1.5 Cardiomyocytes ----
cardio_binom_loss <- rbind(muscleGroupBinom_loss[["cardiomyocyte_developing"]], muscleGroupBinom_loss[["cardiomyocyte_adult"]])
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(cardio_binom_loss, "cardio_binom_loss.rds")

## 1.5 Innate immune cells ----
immuneTypes <- c("macrophage", "microglia", "mast", "naturalkillert")
immuneGroupBinom_loss <- subsetAndBinomCellTypeList(cellTypeList = immuneTypes,
                                               stageList = stageList, 
                                               fullDF = nodeLosses, 
                                               totalNodeRegions = totalNodeRegions_losses, 
                                               totalCellTypeRegions = totalCellTypeRegions_losses, 
                                               nodeNames = nodeNames_losses, 
                                               groupName = "InnateImmune")
### Combine for plotting
innateImmune_binom_loss <- do.call(rbind, immuneGroupBinom_loss)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(innateImmune_binom_loss, "innateImmune_binom_loss")
