setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/figures/enrichmentLinePlots")
# -------- COMBINED CELL TYPE GROUPS --------
## Innate immune ----
innateImmune_combined_enrichPlot <- ggplot(innateImmuneComb_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 1) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  #scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  scale_y_reverse(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) +
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

## Neurons ----
neuron_combined_enrichPlot <- ggplot(neuronComb_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 1) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 0),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  #scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  scale_y_reverse(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) +
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

## Muscle ----
muscle_combined_enrichPlot <- ggplot(muscleComb_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 1, show.legend = TRUE) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  #scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  scale_y_reverse(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) +
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
muscle_combined_enrichPlot
ggsave("muscle_combined_enrichPlot.png", plot = muscle_combined_enrichPlot, bg = "transparent", width = 9, height = 2)




# -------- SEPARATE SUBTYPES --------
## Innate Immune ----
labels = c("macrophage_adult" = "Macrophage, Adult", 
           "macrophage_developing" = "Macrophage, Developing", 
           "mast_adult" = "Mast Cell, Adult", 
           "microglia_adult" = "Microglia, Adult", 
           "naturalkillert_adult" = "Natural Killer T Cell, Adult")
innateImmune_subtype_enrichPlot <- ggplot(innateImmune_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "gray", alpha = 0.5) +
  geom_line(linewidth = 1, color = "black") + 
  scale_x_reverse() +
  theme_bw() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_y_continuous(breaks = seq(0, 2, by=1), limits=c(0,2)) + 
  scale_color_discrete(labels = labels) +
  labs(title = "Binomial Enrichment, Innate Immune Cells", 
       x = "Millions of Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        legend.position = "none", 
        axis.text.y = element_text(size = 15, face = "bold"), 
        axis.text.x = element_text(size = 15, face = "bold")) + 
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 
innateImmune_subtype_enrichPlot
ggsave("innateImmune_subtype_enrichPlot.png", plot = innateImmune_subtype_enrichPlot, bg = "transparent", width = 5, height = 11)

## Muscle ----
labels = c("cardiomyocyte_adult" = "Cardiomyocyte, Adult", 
           "cardiomyocyte_developing" = "Cardiomyocyte, Developing", 
           "skeletalmyocyte_adult" = "Skeletal Myocyte, Adult", 
           "skeletalmyocyte_developing" = "Skeletal Myocyte, Developing", 
           "smoothmuscle_adult" = "Smooth Muscle, Adult")
muscle_subtype_enrichPlot <- ggplot(muscle_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "gray", alpha = 0.5) +
  geom_line(linewidth = 1, color = "black") + 
  scale_x_reverse() +
  theme_bw() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_y_continuous(breaks = seq(0, 2.5, by=1), limits=c(0,2.5)) + 
  scale_color_discrete(labels = labels) +
  labs(title = "Binomial Enrichment, Muscle", 
       x = "Millions of Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        legend.position = "none", 
        axis.text.y = element_text(size = 15, face = "bold"), 
        axis.text.x = element_text(size = 15, face = "bold")) + 
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 
muscle_subtype_enrichPlot
ggsave("muscle_subtype_enrichPlot.png", plot = muscle_subtype_enrichPlot, bg = "transparent", width = 5, height = 11)

# Neurons ----
labels = c("excitatoryneuron" = "Excitatory Neuron", 
           "inhibitoryneuron" = "Inhibitory Neuron")
neuron_subtype_enrichPlot <- ggplot(neuron_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "gray", alpha = 0.5) +
  geom_line(linewidth = 1, color = "black") + 
  scale_x_reverse() +
  theme_bw() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_y_continuous(breaks = seq(0, 2.5, by=1), limits=c(0,2.5)) + 
  scale_color_discrete(labels = labels) +
  labs(title = "Binomial Enrichment, Muscle", 
       x = "Millions of Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        panel.background = element_rect(fill = "transparent", color = NA), # Transparent panel
        plot.background = element_rect(fill = "transparent", color = NA),  # Transparent background
        legend.position = "none", 
        axis.text.y = element_text(size = 15, face = "bold"), 
        axis.text.x = element_text(size = 15, face = "bold")) + 
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 
neuron_subtype_enrichPlot
ggsave("neuron_subtype_enrichPlot.png", plot = neuron_subtype_enrichPlot, bg = "transparent", width = 5, height = 4.75)
