# Combines z-scores for similar cell types usin Stouffer's method
# Hailey Napier
# 8/6/26

require(dplyr)
stouffer_cellTypes <- function(fullDF, pattern){
  tmp <- fullDF %>%
    filter(str_detect(cellType, pattern))
  nSubtypes <- length(unique(tmp$cellType))
  outDF <- tmp %>%
    select(node, zscore) %>%
    group_by(node) %>%
    summarize(Sum = sum(zscore))
  outDF$combZscore <- outDF$Sum/sqrt(nSubtypes)
  outDF$Sum <- NULL
  outDF$cellType <- pattern
  return(outDF)
}