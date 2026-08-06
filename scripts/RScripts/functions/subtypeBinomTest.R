subtypeBinomTest <- function(subtypeDF, totalNodeRegions, totalCellTypeRegions, nodeNames){
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
  } 
  outDF$cellType <- subtypeDF$cellType
  outDF$nodeName <- nodeNames
  outDF$sig <- outDF$pval <= 0.05
  return(outDF)
}