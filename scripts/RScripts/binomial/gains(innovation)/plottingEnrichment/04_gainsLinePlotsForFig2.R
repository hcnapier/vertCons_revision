# 0.0 Setup ----
## 0.1 Load packages ----
require(ggplot2)
require(wesanderson)
require(egg)
require(grid)

## 0.2 Load data ----
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
neuron_binom_combNodes <- readRDS("neuron_binom_combNodes.rds")
muscle_binom_combNodes <- readRDS("muscle_binom_combNodes.rds")
innateImmune_binom_combNodes <- readRDS("innateImmune_binom_combNodes.rds")

# 1.0 Generate a 15-color continuous Zissou1 palette ----
zissou_15 <- wes_palette("Zissou1", 15, type = "continuous")
# View the hex codes
zissou_15 %>% as.character()
# Plot a quick color swatch to see them
image(1:15, 1, as.matrix(1:15), col = zissou_15, xlab = "Zissou 15-color", ylab = "", yaxt = "n")

# 2.0 Plots for Fig 2 ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/figures/enrichmentLinePlots")
## 2.A Muscle ----
muscle_binom_combNodes$col <- NA
for(currNode in muscle_binom_combNodes$Node){
  if(currNode == "H_3"){
    muscle_binom_combNodes$col[which(muscle_binom_combNodes$Node == currNode)] <- zissou_15[15]
  }else{
    muscle_binom_combNodes$col[which(muscle_binom_combNodes$Node == currNode)] <- zissou_15[18-as.numeric(currNode)]
  }
}
muscleSubset <- muscle_binom_combNodes %>%
  filter(cellType %in% c("cardiomyocyte_adult", "smoothmuscle_adult")) %>% 
  arrange(cellType, desc(MYA))
labels = c("cardiomyocyte_adult" = "Cardiomyocyte, Adult", 
           "smoothmuscle_adult" = "Smooth Muscle, Adult")
legend_df <- muscleSubset %>%
  distinct(MYA, col) %>%
  arrange(desc(MYA))
muscle_subtype_enrichPlot <- ggplot(muscleSubset, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI, fill = "95% CI"), alpha = 0.8) +
  scale_x_reverse() +
  theme_light() +
  geom_path(aes(group = cellType, color = col), linewidth = 2) + 
  geom_point(aes(shape = sig, color = col), size = 3, stroke = 2) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_fill_manual(name = NULL, values = c("95% CI" = "gray70")) +
  scale_y_continuous(breaks = seq(0, 2, by=1), limits=c(0,2)) + 
  scale_color_discrete(labels = labels) +
  scale_color_identity(guide = "legend",
                       name = "MYA",
                       breaks = legend_df$col,
                       labels = legend_df$MYA) +
  labs(title = "Regulatory Innovation, Muscle", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.text.y = element_text(size = 15, face = "bold"), 
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 15, face = "bold"),
        strip.text = element_text(face = "bold", color = "black", size = 15, hjust = 0.5),
        #strip.background = element_rect(fill = "gray70", color = "white", linewidth = 1),
        legend.position = "none")+ 
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 
muscle_subtype_enrichPlot

### 2.A.X Sup figure ----
labels = c("cardiomyocyte_adult" = "Cardiomyocyte, Adult", 
           "cardiomyocyte_developing" = "Cardiomyocyte, Developing", 
           "skeletalmyocyte_adult" = "Skeletal Myocyte, Adult", 
           "skeletalmyocyte_developing" = "Skeletal Myocyte, Developing", 
           "smoothmuscle_adult" = "Smooth Muscle, Adult")
muscle_suppPlot <- ggplot(muscle_binom_combNodes, aes(x = MYA, y = enrich)) + 
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
  labs(title = "Regulatory Innovation, Muscle", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        legend.position = "none", 
        axis.text.y = element_text(size = 15, face = "bold"), 
        axis.text.x = element_text(size = 15, face = "bold")) + 
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 
muscle_suppPlot
ggsave("muscle_suppPlot.png", plot = muscle_suppPlot, bg = "transparent", width = 5, height = 11)

## 2.B Neurons ----
neuron_binom_combNodes$col <- NA
for(currNode in neuron_binom_combNodes$Node){
  if(currNode == "H_3"){
    neuron_binom_combNodes$col[which(neuron_binom_combNodes$Node == currNode)] <- zissou_15[15]
  }else{
    neuron_binom_combNodes$col[which(neuron_binom_combNodes$Node == currNode)] <- zissou_15[18-as.numeric(currNode)]
  }
}
neuron_binom_combNodes <- neuron_binom_combNodes %>%
  arrange(cellType, desc(MYA))
labels = c("excitatoryneuron" = "Excitatory (Glutamatergic)", 
           "inhibitoryneuron" = "Inhibitory (GABAergic)")
legend_df <- neuron_binom_combNodes %>%
  distinct(MYA, col) %>%
  arrange(desc(MYA))
neuron_subtype_enrichPlot <- ggplot(neuron_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI, fill = "95% CI"), alpha = 0.8) +
  scale_x_reverse() +
  theme_light() +
  geom_path(aes(group = cellType, color = col), linewidth = 2) + 
  geom_point(aes(shape = sig, color = col), size = 3, stroke = 2) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_fill_manual(name = NULL, values = c("95% CI" = "gray70")) +
  scale_y_continuous(breaks = seq(0, 2, by=1), limits=c(0,2)) + 
  scale_color_discrete(labels = labels) +
  scale_color_identity(guide = "legend",
                       name = "MYA",
                       breaks = legend_df$col,
                       labels = legend_df$MYA) +
  labs(title = "Regulatory Innovation, Neurons", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.text.y = element_text(size = 15, face = "bold"), 
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 15, face = "bold"),
        strip.text = element_text(face = "bold", color = "black", size = 15, hjust = 0.5),
        #strip.background = element_rect(fill = "gray70", color = "white", linewidth = 1), 
        legend.position = "none")+ 
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 
neuron_subtype_enrichPlot

## 2.C Innate Immune ----
innateImmune_binom_combNodes$col <- NA
for(currNode in innateImmune_binom_combNodes$Node){
  if(currNode == "H_3"){
    innateImmune_binom_combNodes$col[which(innateImmune_binom_combNodes$Node == currNode)] <- zissou_15[15]
  }else{
    innateImmune_binom_combNodes$col[which(innateImmune_binom_combNodes$Node == currNode)] <- zissou_15[18-as.numeric(currNode)]
  }
}
immuneSubset <- innateImmune_binom_combNodes %>%
  filter(cellType %in% c("macrophage_adult", "mast_adult", "naturalkillert_adult")) %>% 
  arrange(cellType, desc(MYA))
labels = c("macrophage_adult" = "Macrophage, Adult", 
           "mast_adult" = "Mast Cell, Adult", 
           "naturalkillert_adult" = "Natural Killer T Cell, Adult")
legend_df <- immuneSubset %>%
  distinct(MYA, col) %>%
  arrange(desc(MYA))
innateImmune_subtype_enrichPlot <- ggplot(immuneSubset, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI, fill = "95% CI"), alpha = 0.8) +
  scale_x_reverse() +
  theme_light() +
  geom_path(aes(group = cellType, color = col), linewidth = 2) + 
  geom_point(aes(shape = sig, color = col), size = 3, stroke = 2) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_fill_manual(name = NULL, values = c("95% CI" = "gray70")) +
  scale_y_continuous(breaks = seq(0, 2, by=1), limits=c(0,2)) + 
  scale_color_discrete(labels = labels) +
  scale_color_identity(guide = "legend",
                       name = "MYA",
                       breaks = legend_df$col,
                       labels = legend_df$MYA) +
  labs(title = "Regulatory Innovation, Innate Immune Cells", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.text.y = element_text(size = 15, face = "bold"), 
        axis.text.x = element_text(size = 15, face = "bold"), 
        axis.title.x = element_text(size = 15, face = "bold"), 
        axis.title.y = element_text(size = 15, face = "bold"),
        strip.text = element_text(face = "bold", color = "black", size = 15, hjust = 0.5),
        #strip.background = element_rect(fill = "gray70", color = "white", linewidth = 1),
        legend.position = "none")+ 
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 
innateImmune_subtype_enrichPlot

# 3.0 Save panels for fig 2 ----
save_fig2_panels <- function(p, file, panel_w = 5, panel_h = 1.25) {
  g <- set_panel_size(p, width = unit(panel_w, "in"), height = unit(panel_h, "in"))
  w <- convertWidth(sum(g$widths), "in", valueOnly = TRUE)
  h <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  ggsave(file, g, width = w, height = h, bg = "transparent")
}
save_fig2_panels(innateImmune_subtype_enrichPlot, "innateImmune_subtype_enrichPlot.png")
save_fig2_panels(neuron_subtype_enrichPlot, "neuron_subtype_enrichPlot.png")
save_fig2_panels(muscle_subtype_enrichPlot, "muscle_subtype_enrichPlot.png")
