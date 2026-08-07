subsetAndBinomCellTypeList <- function(cellTypeList, stageList = NULL, fullDF, totalNodeRegions, totalCellTypeRegions, nodeNames, MYA = NA){
  mapsList <- list()
    for(cellType in cellTypeList){
      if(length(stageList) > 0){
        for(stage in stageList){
          name <- paste0(c(cellType, stage), collapse = "_")
          mapsList[[name]] <- getCombMapped(fullDF, c(cellType, stage))
        }
      }else{
        name <- cellType
        mapsList[[name]] <- getCombMapped(fullDF, cellType)
        print(paste(name, " found mapping regions", sep = ""))
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