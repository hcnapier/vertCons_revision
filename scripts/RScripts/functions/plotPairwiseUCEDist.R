require(dplyr)
require(stringr)
plotPairwiseUCEDist <- function(firstSpecies, secondSpecies){
  title <- paste(paste(firstSpecies, secondSpecies, sep = "-"), "Pairwise Embedding Distance", sep = " ") %>% str_to_title()
  
  subset_uceDist <- uceDist %>%
    filter(species1 %in% c(firstSpecies, secondSpecies) & species2 %in% c(firstSpecies, secondSpecies))
  subset_means <- subset_uceDist %>%
    group_by(color) %>%
    summarise(meanDist = mean(distance, na.rm = TRUE)) 
  
  # Pairwise tests on mean distance between colors
  pw <- pairwise.wilcox.test(subset_uceDist$distance, subset_uceDist$color,
                             p.adjust.method = "BH")
  
  sigDF <- as.data.frame(as.table(pw$p.value)) %>%
    filter(!is.na(Freq)) %>%
    transmute(g1 = as.character(Var1), g2 = as.character(Var2), p = Freq) %>%
    mutate(stars = case_when(p < 0.001 ~ "***",
                             p < 0.01  ~ "**",
                             p < 0.05  ~ "*",
                             TRUE      ~ "NS")) %>%
    left_join(subset_means, by = c("g1" = "color")) %>% rename(x1 = meanDist) %>%
    left_join(subset_means, by = c("g2" = "color")) %>% rename(x2 = meanDist)
  
  # Place brackets above the tallest density curve
  maxDens <- subset_uceDist %>%
    group_by(color) %>%
    summarise(peak = max(density(distance, na.rm = TRUE)$y)) %>%
    pull(peak) %>% max()
  step <- 0.08 * maxDens
  
  sigDF <- sigDF %>%
    mutate(xmin = pmin(x1, x2), xmax = pmax(x1, x2)) %>%
    arrange(xmax - xmin) %>%                 # short brackets lowest
    mutate(y = maxDens + step * row_number())
  
  pairwisePlot <- ggplot(data = subset_uceDist, aes(x = distance)) + 
    geom_density(aes(color = color, fill = color), alpha = 0.25) + 
    geom_vline(data = subset_means,
               aes(xintercept = meanDist, color = color),
               linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
    geom_segment(data = sigDF, aes(x = xmin, xend = xmax, y = y, yend = y),
                 inherit.aes = FALSE) +
    geom_text(data = sigDF, aes(x = (xmin + xmax) / 2, y = y, label = stars),
              vjust = -0.1, size = 5, inherit.aes = FALSE) +
    scale_color_identity(guide = "legend",
                         breaks = legendKey$color,
                         labels = legendKey$cellType,
                         name = "Placental Cell Type") + 
    scale_fill_identity(guide = "legend",
                        breaks = legendKey$color,
                        labels = legendKey$cellType,
                        name = "Placental Cell Type") +
    theme_minimal() + 
    labs(y = "Density", x = "UCE Embedding Distance", 
         title = title) + 
    theme(text = element_text(family = "Helvetica"), 
          axis.title.x = element_text(face = "bold", size = 12),   
          axis.title.y = element_text(face = "bold", size = 12), 
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"), 
          axis.text.x = element_text(face = "bold"), 
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5)) 
  
  return(pairwisePlot)
}
