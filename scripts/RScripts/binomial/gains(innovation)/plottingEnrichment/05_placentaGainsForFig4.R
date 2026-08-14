setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/figures/enrichmentLinePlots")
# Separate placental subtypes for schematic ----
## Placental neurons ----
placentalNeuron_binom <- placenta_binom_combNodes %>%
  filter(cellType == "placentalNeuron")
placentalNeuron_enrichPlot <- ggplot(placentalNeuron_binom, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 1) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  #scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  scale_y_reverse(breaks = seq(0, 3.75, by=1), limits=c(0,3.75)) +
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
placentalNeuron_enrichPlot
ggsave("placentalNeuron_enrichPlot.png", plot = placentalNeuron_enrichPlot, bg = "transparent", width = 9, height = 2)  

## Placental fibroblasts ----
placentalFibro_binom <- placenta_binom_combNodes %>%
  filter(cellType == "fibroPlacental")
placentalFibro_enrichPlot <- ggplot(placentalFibro_binom, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 1) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  #scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  scale_y_reverse(breaks = seq(0, 3.75, by=1), limits=c(0,3.75)) +
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
placentalFibro_enrichPlot
ggsave("placentalFibroblast_enrichPlot.png", plot = placentalFibro_enrichPlot, bg = "transparent", width = 9, height = 2)  

## Placental macrophages ----
placentalMacro_binom <- placenta_binom_combNodes %>%
  filter(cellType == "macrophagePlacental")
placentalMacro_enrichPlot <- ggplot(placentalMacro_binom, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 1) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  #scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  scale_y_reverse(breaks = seq(0, 3.75, by=1), limits=c(0,3.75)) +
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
placentalMacro_enrichPlot
ggsave("placentalMacrophage_enrichPlot.png", plot = placentalMacro_enrichPlot, bg = "transparent", width = 9, height = 2)  

## Placental endothelial cells ----
placentalEndo_binom <- placenta_binom_combNodes %>%
  filter(cellType == "endothelialPlacental")
placentalEndo_enrichPlot <- ggplot(placentalEndo_binom, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 1) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  #scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  scale_y_reverse(breaks = seq(0, 3.75, by=1), limits=c(0,3.75)) +
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
placentalEndo_enrichPlot
ggsave("placentalEndothelial_enrichPlot.png", plot = placentalEndo_enrichPlot, bg = "transparent", width = 9, height = 2)  

## Extravillous trophoblasts ----
placentalExTro_binom <- placenta_binom_combNodes %>%
  filter(cellType == "extravillousTrophoblast")
placentalExTro_enrichPlot <- ggplot(placentalExTro_binom, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 1) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  #scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  scale_y_reverse(breaks = seq(0, 3.75, by=1), limits=c(0,3.75)) +
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
placentalExTro_enrichPlot
ggsave("placentalExtravillousTrophoblast_enrichPlot.png", plot = placentalExTro_enrichPlot, bg = "transparent", width = 9, height = 2)  

## Syncitiotrophoblasts & Cytotrophoblasts ----
placentalSynCyto_binom <- placenta_binom_combNodes %>%
  filter(cellType == "syncitiotrophoblastCytotrophoblast")
placentalSynCyto_enrichPlot <- ggplot(placentalSynCyto_binom, aes(x = MYA, y = enrich)) + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 1) +
  geom_line(linewidth = 1.2) + 
  theme_classic() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3.5, stroke = 1) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  #scale_y_continuous(breaks = seq(0.5, 1.5, by=0.5), limits=c(0.5,1.5)) + 
  scale_y_reverse(breaks = seq(0, 3.75, by=1), limits=c(0,3.75)) +
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
placentalSynCyto_enrichPlot
ggsave("placentalSynCyto_enrichPlot.png", plot = placentalSynCyto_enrichPlot, bg = "transparent", width = 9, height = 2)  


# All placental cell types in a facet plot ----
## Facet wrapped ----
labels = c("placentalNeuron" = "Placental Neuron", 
           "fibroPlacental" = "Placental Fibroblast", 
           "macrophagePlacental" = "Placental Macrophage", 
           "extravillousTrophoblast" = "Extravillous Trophoblast", 
           "syncitiotrophoblastCytotrophoblast" = "Syncitiotrophoblast and Cytotrophoblast")
ggplot(placenta_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 0.8) +
  geom_line(linewidth = 1, color = "black") + 
  theme_minimal() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Gains, Binomial Enrichment, Placenta", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  ylim(0,3.5) +
  scale_x_reverse() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  facet_wrap(~cellType, 
             ncol = 2, 
             labeller = as_labeller(labels)) 