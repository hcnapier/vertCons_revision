nodeNames <- totalNodeRegions_combNodes$name %>% unique()

## 1.3 Neurons ----
neuronTypes <- c("excitatoryneuron", "inhibitoryneuron")
### Separate subtypes ----
neuronGroupBinom_combNodes <- subsetAndBinomCellTypeList(cellTypeList = neuronTypes,
                                                         fullDF = nMapLong_combNodes, 
                                                         totalNodeRegions = totalNodeRegions_combNodes, 
                                                         totalCellTypeRegions = totalCellTypeRegions, 
                                                         nodeNames = nodeNames, 
                                                         MYA = MYA)
# Combine for plotting
neuron_binom_combNodes <- do.call(rbind, neuronGroupBinom_combNodes)

### Combined subtypes ----
neuronCombBinom_combNodes <- subsetAndBinomCellTypeList(cellTypeList = neuronTypes,
                                                        fullDF = nMapLong_combNodes, 
                                                        totalNodeRegions = totalNodeRegions_combNodes, 
                                                        totalCellTypeRegions = totalCellTypeRegions, 
                                                        nodeNames = nodeNames, 
                                                        MYA = MYA, 
                                                        logic = "OR", 
                                                        groupName = "Neurons")
## 1.4 All Muscle Subtypes ----
stageList <- c("developing", "adult")
muscleTypes <- c("smoothmuscle", "cardiomyocyte", "skeletalmyocyte")
muscleGroupBinom_combNodes <- subsetAndBinomCellTypeList(cellTypeList = muscleTypes, 
                                                         stageList = stageList, 
                                                         fullDF = nMapLong_combNodes, 
                                                         totalNodeRegions = totalNodeRegions_combNodes,
                                                         totalCellTypeRegions = totalCellTypeRegions, 
                                                         nodeNames = nodeNames, 
                                                         MYA = MYA)
### Combine for plotting
muscle_binom_combNodes <- do.call(rbind, muscleGroupBinom_combNodes)

## 1.5 Innate immune cells ----
immuneTypes <- c("macrophage", "microglia", "mast", "naturalkillert")
immuneGroupBinom_combNodes <- subsetAndBinomCellTypeList(cellTypeList = immuneTypes, 
                                                         stageList = stageList, 
                                                         fullDF = nMapLong_combNodes, 
                                                         totalNodeRegions = totalNodeRegions_combNodes, 
                                                         totalCellTypeRegions = totalCellTypeRegions, 
                                                         nodeNames = nodeNames, 
                                                         MYA = MYA)
### Combine for plotting
innateImmune_binom_combNodes <- do.call(rbind, immuneGroupBinom_combNodes)

