getTotalGroupRegions <- function(groupDF){
  outDF = data.frame(nRegions = sum(groupDF$mappedSum))
  return(outDF)
}