require(dplyr)
require(stringr)
lossBinomTest <- function(subtypeDF, totalNodeRegions, totalNodeLosses, nodeNames, MYA = NULL){
  outDF <- data.frame(Node = subtypeDF$Node)
  for(currNode in subtypeDF$Node){
    currNodeName <- totalNodeRegions$name[which(totalNodeRegions$Node == currNode)]
    print(currNodeName)
    numerator <- totalNodeLosses$nLosses[which(totalNodeLosses$Node == currNode)] - subtypeDF$successSum[which(subtypeDF$Node == currNode)] # All success regions for node of interest - success regions in current cell type
    denominator <- totalNodeRegions$nRegions[which(totalNodeRegions$Node == currNode)] - subtypeDF$nRegions[which(subtypeDF$Node == currNode)]  # All regions - total regions in node of interest
    nullPr <- numerator/denominator # successes over trials with current cell type held out
    nTrials = subtypeDF$nRegions[which(subtypeDF$Node == currNode)]
    nSuccesses = subtypeDF$successSum[which(subtypeDF$Node == currNode)]
    test <- binom.test(nSuccesses, nTrials, nullPr, alternative = "two.sided")
    outDF$estimate[which(outDF$Node == currNode)] <- test$estimate
    outDF$pval[which(outDF$Node == currNode)] <- test$p.value
    outDF$enrich[which(outDF$Node == currNode)] <- (nSuccesses/nTrials)/nullPr
    outDF$upperCI[which(outDF$Node == currNode)] <- test$conf.int[2]
    outDF$lowerCI[which(outDF$Node == currNode)] <- test$conf.int[1]
    outDF$enrUpperCI[which(outDF$Node == currNode)] <- test$conf.int[2]/test$null.value[[1]]
    outDF$enrLowerCI[which(outDF$Node == currNode)] <- test$conf.int[1]/test$null.value[[1]]
    outDF$nullPR[which(outDF$Node == currNode)] <- nullPr
    outDF$num[which(outDF$Node == currNode)] <- numerator
    outDF$den[which(outDF$Node == currNode)] <- denominator
    outDF$nSuccess[which(outDF$Node == currNode)] <- nSuccesses
    outDF$nTrials[which(outDF$Node == currNode)] <- nTrials
  } 
  outDF$groupName <- subtypeDF$groupName
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