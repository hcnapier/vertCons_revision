require(dplyr)
downsampleUCEDistances <- function(df, num){
set.seed(42)
# 1 Get every unique cell with its species (cells can appear in either column) ----
cells <- bind_rows(
  df %>% distinct(species = species1, cellid = cellID1),
  df %>% distinct(species = species2, cellid = cellID2)) %>%
  distinct()
  
# 2 Randomly pick 200 cells per species ----
sampled_cells <- cells %>%
  group_by(species) %>%
  slice_sample(n = num) %>%
  ungroup()

# 3 Filter dataframe to include only sampled cells ----
df_sub <- df %>%
  semi_join(sampled_cells, by = c(species1 = "species", cellID1 = "cellid")) %>%
  semi_join(sampled_cells, by = c(species2 = "species", cellID2 = "cellid"))

return(df_sub)
} 