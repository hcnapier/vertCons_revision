# 0.0 Setup ----
## 0.1 Load packages ----
require(dplyr)
require(tidyr)

## 0.2 Load data ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
totalCellTypeRegions <- readRDS("totalCellTypeRegions.rds")
totalNodeRegions_combNodes <- readRDS("totalNodeRegions_combNodes.rds")
nMapLong_combNodes <- readRDS("nMapLong_combNodes.rds")
MYA <- c(20, 30, 45, 75, 85, 87, 95, 100, 160, 180, 320, 350, 415, 430, 560)

## 0.3 Source functions ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("subsetAndBinomCellTypeList.R")

nodeNames <- totalNodeRegions_combNodes$name %>% unique()

## 1.3 Neurons ----
neuronTypes <- c("excitatoryneuron", "inhibitoryneuron")
### Separate subtypes ----
neuronGroupBinom_combNodes <- subsetAndBinomCellTypeList(cellTypeList = neuronTypes,
                                                         fullDF = nMapLong_combNodes, 
                                                         totalNodeRegions = totalNodeRegions_combNodes, 
                                                         totalCellTypeRegions = totalCellTypeRegions, 
                                                         nodeNames = nodeNames, 
                                                         MYA = MYA, 
                                                         groupName = "Neurons")
#### Combine for plotting
neuron_binom_combNodes <- do.call(rbind, neuronGroupBinom_combNodes)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(neuron_binom_combNodes, "neuron_binom_combNodes.rds")

### Combined subtypes ----
neuronCombDF <- getCombMapped(nMapLong_combNodes,
                              neuronTypes,
                              "Neurons", 
                              "OR")
neuronComb_binom_combNodes <- subtypeBinomTest(neuronCombDF, 
                                               totalNodeRegions_combNodes, 
                                               totalCellTypeRegions, 
                                               nodeNames = nodeNames, 
                                               MYA = MYA)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(neuronComb_binom_combNodes, "neuronComb_binom_combNodes.rds")

## 1.4 All Muscle Subtypes ----
stageList <- c("developing", "adult")
muscleTypes <- c("smoothmuscle", "cardiomyocyte", "skeletalmyocyte")
muscleGroupBinom_combNodes <- subsetAndBinomCellTypeList(cellTypeList = muscleTypes, 
                                                         stageList = stageList, 
                                                         fullDF = nMapLong_combNodes, 
                                                         totalNodeRegions = totalNodeRegions_combNodes,
                                                         totalCellTypeRegions = totalCellTypeRegions, 
                                                         nodeNames = nodeNames, 
                                                         MYA = MYA, 
                                                         groupName = "Muscle")
#### Combine for plotting
muscle_binom_combNodes <- do.call(rbind, muscleGroupBinom_combNodes)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(muscle_binom_combNodes, "muscle_binom_combNodes.rds")

### Combined subtypes ----
muscleCombDF <- getCombMapped(nMapLong_combNodes,
                              muscleTypes,
                              "Muscle", 
                              "OR")
muscleComb_binom_combNodes <- subtypeBinomTest(muscleCombDF, 
                                               totalNodeRegions_combNodes, 
                                               totalCellTypeRegions, 
                                               nodeNames = nodeNames, 
                                               MYA = MYA)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(muscleComb_binom_combNodes, "muscleComb_binom_combNodes.rds")

## 1.5 Innate immune cells ----
### Separate subtypes ----
immuneTypes <- c("macrophage", "microglia", "mast", "naturalkillert")
immuneGroupBinom_combNodes <- subsetAndBinomCellTypeList(cellTypeList = immuneTypes, 
                                                         stageList = stageList, 
                                                         fullDF = nMapLong_combNodes, 
                                                         totalNodeRegions = totalNodeRegions_combNodes, 
                                                         totalCellTypeRegions = totalCellTypeRegions, 
                                                         nodeNames = nodeNames, 
                                                         MYA = MYA, 
                                                         groupName= "InnateImmune")
#### Combine for plotting
innateImmune_binom_combNodes <- do.call(rbind, immuneGroupBinom_combNodes)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(innateImmune_binom_combNodes, "innateImmune_binom_combNodes.rds")

### Combined subtypes ----
innateImmuneCombDF <- getCombMapped(nMapLong_combNodes,
                              immuneTypes,
                              "InnateImmune", 
                              "OR")
innateImmuneComb_binom_combNodes <- subtypeBinomTest(innateImmuneCombDF, 
                                               totalNodeRegions_combNodes, 
                                               totalCellTypeRegions, 
                                               nodeNames = nodeNames, 
                                               MYA = MYA)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(innateImmuneComb_binom_combNodes, "innateImmuneComb_binom_combNodes.rds")


## 1.6 Placenta ----
placentaTypes <- c("placentalNeuron", 
                   "fibroPlacental", 
                   "macrophagePlacental", 
                   "extravillousTrophoblast", 
                   "syncitiotrophoblastCytotrophoblast", 
                   "endothelialPlacental")
### Separate subtypes ----
placentaGroupBinom_combNodes <- subsetAndBinomCellTypeList(cellTypeList = placentaTypes,
                                                         fullDF = nMapLong_combNodes, 
                                                         totalNodeRegions = totalNodeRegions_combNodes, 
                                                         totalCellTypeRegions = totalCellTypeRegions, 
                                                         nodeNames = nodeNames, 
                                                         MYA = MYA, 
                                                         groupName = "Placenta")
#### Combine for plotting
placenta_binom_combNodes <- do.call(rbind, placentaGroupBinom_combNodes)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(placenta_binom_combNodes, "placenta_binom_combNodes.rds")
