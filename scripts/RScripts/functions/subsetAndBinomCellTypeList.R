subsetAndBinomCellTypeList <- function(cellTypeList, stageList = NULL, fullDF, totalNodeRegions, totalCellTypeRegions, nodeNames, MYA = NA, logic = "AND", groupName = NULL){
  mapsList <- list()
    for(cellType in cellTypeList){
      if(length(stageList) > 0){
        if(logic == "AND"){
          name <- paste0(c(cellType, stage), collapse = "_")
          for(stage in stageList){
            mapsList[[name]] <- getCombMapped(fullDF, c(cellType, stage))
          }
        }else if(logic == "OR"){
          for(stage in stageList){
            mapsList[[name]] <- getCombMapped(fullDF, c(cellType, stage), logic = "OR", groupName = groupName)
          }
        }
      }else{
        name <- cellType
        if(logic == "AND"){
          mapsList[[name]] <- getCombMapped(fullDF, cellType)
        }else if(logic == "OR"){
          mapsList[[name]] <- getCombMapped(fullDF, cellType, logic = "OR", groupName = groupName)
        }
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