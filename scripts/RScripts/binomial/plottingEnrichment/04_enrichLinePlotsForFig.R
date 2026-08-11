setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/figures/enrichmentLinePlots")
# Innate immune ----
## Combined ----
innateImmune_combined_enrichPlot <- ggplot(innateImmuneComb_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 1) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  theme(
    panel.background = element_rect(fill = "transparent", color = NA), # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),  # Transparent background
    legend.position = "none", 
    axis.text.y = element_text(size = 10, face = "bold"), 
    axis.line.y = element_line(linewidth = 0.75),
    axis.title = element_blank(), 
    axis.line.x  = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.text.x  = element_blank(), 
    ) 
innateImmune_combined_enrichPlot
ggsave("innateImmune_combined_enrichPlot.png", plot = innateImmune_combined_enrichPlot, bg = "transparent", width = 9, height = 2)      

# Neurons ----
neuron_combined_enrichPlot <- ggplot(neuronComb_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 1) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 0),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  theme(
    panel.background = element_rect(fill = "transparent", color = NA), # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),  # Transparent background
    legend.position = "none", 
    axis.text.y = element_text(size = 10, face = "bold"), 
    axis.line.y = element_line(linewidth = 0.75),
    axis.title = element_blank(), 
    axis.line.x  = element_blank(),  
    axis.ticks.x = element_blank(),  
    axis.text.x  = element_blank(), 
  ) 
neuron_combined_enrichPlot
ggsave("neuron_combined_enrichPlot.png", plot = neuron_combined_enrichPlot, bg = "transparent", width = 9, height = 2)

# Muscle ----
muscle_combined_enrichPlot <- ggplot(muscleComb_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI, fill = "CI"), alpha = 1, show.legend = TRUE) +
  scale_fill_manual(values = c("CI" = "darkgray")) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  theme(
    panel.background = element_rect(fill = "transparent", color = NA), # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),  # Transparent background
    axis.text.y = element_text(size = 10, face = "bold"), 
    axis.line.y = element_line(linewidth = 0.75),
    axis.title = element_blank(), 
    axis.line.x  = element_blank(),  
    axis.ticks.x = element_blank(),  
    axis.text.x  = element_blank(), 
  ) 
muscle_combined_enrichPlot
ggsave("muscle_combined_enrichPlot.png", plot = muscle_combined_enrichPlot, bg = "transparent", width = 9, height = 2)


