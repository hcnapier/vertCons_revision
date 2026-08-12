require(dplyr)
require(stringr)
getCombSuccess <- function(fullDF, pattern, groupName, logic = "AND"){
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
  combSuccess <- tmp %>%
    select(Node, nSuccess) %>%
    group_by(Node) %>%
    summarize(successSum = sum(nSuccess))
  combSuccess$cellType <- cellType
  combSuccess$groupName <- groupName
  return(combSuccess)
}