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
