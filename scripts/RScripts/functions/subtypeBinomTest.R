require(dplyr)
require(stringr)
subtypeBinomTest <- function(subtypeDF, totalNodeRegions, totalCellTypeRegions, nodeNames, MYA = NULL){
  outDF <- data.frame(Node = subtypeDF$Node)
  for(currNode in subtypeDF$Node){
    currNodeName <- totalNodeRegions$name[which(totalNodeRegions$Node == currNode)]
    print(currNodeName)
    currNodeRegions <- subtypeDF %>%
      filter(Node == currNode)
    numerator <- sum(subtypeDF$mappedSum) - subtypeDF$mappedSum[which(subtypeDF$Node == currNode)]
    denominator <- sum(totalCellTypeRegions$nRegions) - totalNodeRegions$nRegions[which(totalNodeRegions$Node == currNode)]  
    nullPr <- numerator/denominator
    nTrials = totalNodeRegions$nRegions[which(totalNodeRegions$Node == currNode)]
    nSuccesses = subtypeDF$mappedSum[which(subtypeDF$Node == currNode)]
    test <- binom.test(nSuccesses, nTrials, nullPr, alternative = "two.sided")
    outDF$estimate[which(outDF$Node == currNode)] <- test$estimate
    outDF$pval[which(outDF$Node == currNode)] <- test$p.value
    outDF$enrich[which(outDF$Node == currNode)] <- (nSuccesses/nTrials)/nullPr
    outDF$upperCI[which(outDF$Node == currNode)] <- test$conf.int[2]
    outDF$lowerCI[which(outDF$Node == currNode)] <- test$conf.int[1]
    outDF$enrUpperCI[which(outDF$Node == currNode)] <- test$conf.int[2]/test$null.value[[1]]
    outDF$enrLowerCI[which(outDF$Node == currNode)] <- test$conf.int[1]/test$null.value[[1]]
  } 
  outDF$cellType <- subtypeDF$cellType
  outDF$nodeName <- nodeNames
  outDF$sig <- outDF$pval <= 0.05
  outDF$stage <- NA
  cellType <- subtypeDF$cellType %>% unique()
  if(str_detect(tolower(cellType), "developing")){
    outDF$stage <- "developing"
  }else if(str_detect(tolower(cellType), "adult")){
    outDF$stage <- "adult"
  }
  if(nrow(outDF) == 15){
    outDF$MYA <- MYA
  }
  return(outDF)
}