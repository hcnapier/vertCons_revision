setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/functions")
source("getCombSuccess.R")
source("subtypeBinomTest.R")
subsetAndBinomCellTypeList <- function(cellTypeList, stageList = NULL, fullDF, totalNodeRegions, totalCellTypeRegions = NULL, totalNodeLosses = NULL, nodeNames, MYA = NA, groupName, logic = "AND", type = "gain"){
  mapsList <- list()
    for(cellType in cellTypeList){
      if(length(stageList) > 0){
        for(stage in stageList){
          name <- paste0(c(cellType, stage), collapse = "_")
          mapsList[[name]] <- getCombSuccess(fullDF, c(cellType, stage), groupName = groupName, logic = logic)
        }
      }else{
        name <- cellType
        mapsList[[name]] <- getCombSuccess(fullDF, cellType, groupName, logic = logic) 
      }
    }
  for(currName in names(mapsList)){
    if(all(is.na(mapsList[[currName]]))){
      mapsList[[currName]] <- NULL
      print(paste("removing ", name, sep = ""))
    }
  }
  if(type == "gain"){
    outListBinom <- lapply(mapsList, subtypeBinomTest, totalNodeRegions = totalNodeRegions, totalCellTypeRegions = totalCellTypeRegions, nodeNames, MYA)
  }else if(type == "loss"){
    outListBinom <- lapply(mapsList, lossBinomTest, totalNodeRegions = totalNodeRegions, totalNodeLosses = totalNodeLosses, nodeNames = nodeNames, MYA = MYA)
  }
}