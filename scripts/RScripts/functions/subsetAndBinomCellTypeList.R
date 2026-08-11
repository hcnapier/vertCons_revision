setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("getCombMapped.R")
source("subtypeBinomTest.R")
subsetAndBinomCellTypeList <- function(cellTypeList, stageList = NULL, fullDF, totalNodeRegions, totalCellTypeRegions, nodeNames, MYA = NA, groupName){
  mapsList <- list()
    for(cellType in cellTypeList){
      if(length(stageList) > 0){
        name <- paste0(c(cellType, stage), collapse = "_")
        for(stage in stageList){
          mapsList[[name]] <- getCombMapped(fullDF, c(cellType, stage), groupName = groupName)
        }
      }else{
        name <- cellType
        mapsList[[name]] <- getCombMapped(fullDF, cellType, groupName) 
      }
    }
  for(currName in names(mapsList)){
    if(all(is.na(mapsList[[currName]]))){
      mapsList[[currName]] <- NULL
      print(paste("removing ", name, sep = ""))
    }
  }
  outListBinom <- lapply(mapsList, subtypeBinomTest, totalNodeRegions = totalNodeRegions, totalCellTypeRegions = totalCellTypeRegions, nodeNames, MYA)
}