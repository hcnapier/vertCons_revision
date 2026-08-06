require(dplyr)
getCombMapped <- function(fullDF, pattern){
  tmp <- fullDF %>%
    filter(str_detect(CellType, pattern))
  combMapped <- tmp %>%
    select(Node, nMapped) %>%
    group_by(Node) %>%
    summarize(mappedSum = sum(nMapped))
  combMapped$cellType <- pattern
  return(combMapped)
}