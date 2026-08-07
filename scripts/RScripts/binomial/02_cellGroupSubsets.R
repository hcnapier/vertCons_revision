nodeNames_full <- colnames(binomPrMat)
## 1.3 Neurons ----
neuronTypes <- c("excitatoryneuron", "inhibitoryneuron")
neuronGroupBinom <- subsetAndBinomCellTypeList(cellTypeList = neuronTypes, 
                                               fullDF = regions, 
                                               totalNodeRegions = totalNodeRegions,
                                               totalCellTypeRegions = totalCellTypeRegions,
                                               nodeNames = nodeNames_full)
# Combine for plotting
neuron_binom <- do.call(rbind, neuronGroupBinom)

## 1.4 All Muscle Subtypes ----
stageList <- c("developing", "adult")
muscleTypes <- c("smoothmuscle", "cardiomyocyte", "skeletalmyocyte")
muscleGroupBinom <- subsetAndBinomCellTypeList(cellTypeList = muscleTypes, 
                                               stageList = stageList, 
                                               fullDF = regions, 
                                               totalNodeRegions = totalNodeRegions, 
                                               totalCellTypeRegions = totalCellTypeRegions, 
                                               nodeNames = nodeNames_full)
### Combine for plotting
muscle_binom <- do.call(rbind, muscleGroupBinom)

## 1.5 Cardiomyocytes ----
cardio_binom <- rbind(muscleGroupBinom[["cardiomyocyte_developing"]], muscleGroupBinom[["cardiomyocyte_adult"]])

## 1.5 Innate immune cells ----
immuneTypes <- c("macrophage", "microglia", "mast", "naturalkillert")
immuneGroupBinom <- subsetAndBinomCellTypeList(cellTypeList = immuneTypes,
                                               stageList = stageList, 
                                               fullDF = regions, 
                                               totalNodeRegions = totalNodeRegions, 
                                               totalCellTypeRegions = totalCellTypeRegions, 
                                               nodeNames = nodeNames_full)
### Combine for plotting
innateImmune_binom <- do.call(rbind, immuneGroupBinom)

