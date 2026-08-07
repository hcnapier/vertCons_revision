require(dplyr)
require(stringr)
getCombMapped <- function(fullDF, pattern, logic = "AND", groupName = NULL){
  if(length(pattern) == 1){
    tmp <- fullDF %>%
      filter(str_detect(tolower(CellType), tolower(pattern)))
    cellType <- pattern
  }else if(length(pattern) > 1){
    keywords <- tolower(pattern)
    if(logic == "AND"){
    pattern <- paste0("(?=.*", keywords, ")", collapse = "")
    tmp <- fullDF %>%
      filter(grepl(pattern, tolower(CellType), perl = T))
    cellType <- paste0(keywords, collapse = "_")
    }else if(logic == "OR"){
      pattern <- paste0(keywords, collapse = "|")
      tmp <- fullDF %>%
        filter(grepl(pattern, tolower(CellType)))
      cellType <- groupName
    }
  }
  combMapped <- tmp %>%
    select(Node, nMapped) %>%
    group_by(Node) %>%
    summarize(mappedSum = sum(nMapped))
  combMapped$cellType <- cellType
  return(combMapped)
}