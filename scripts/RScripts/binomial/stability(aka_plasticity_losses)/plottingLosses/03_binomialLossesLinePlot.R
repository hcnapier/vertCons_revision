# 1.0 All Nodes----
## 1.1 Neurons ----
neuron_binom_loss$MYA <- MYA_losses
labels = c("excitatoryneuron" = "Excitatory Neuron", 
           "inhibitoryneuron" = "Inhibitory Neuron")
ggplot(neuron_binom_loss, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 0.8) +
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Loss, Binomial Enrichment, Neurons",
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  ylim(0,5.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 

## 1.2 Cardiomyocytes ----
### Facet wrap ----
cardio_binom_loss$MYA <- MYA_losses
labels = c("cardiomyocyte_adult" = "Cardiomyocyte, Adult", 
           "cardiomyocyte_developing" = "Cardiomyocyte, Developing")
ggplot(cardio_binom_loss, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 0.8) +
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Loss, Binomial Enrichment, Cardiomyocytes",
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  ylim(0,5.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 

## 1.3 All Muscle ----
muscle_binom_loss$MYA <- MYA_losses
labels = c("cardiomyocyte_adult" = "Cardiomyocyte, Adult", 
           "cardiomyocyte_developing" = "Cardiomyocyte, Developing", 
           "skeletalmyocyte_adult" = "Skeletal Myocyte, Adult", 
           "skeletalmyocyte_developing" = "Skeletal Myocyte, Developing", 
           "smoothmuscle_adult" = "Smooth Muscle, Adult")
ggplot(muscle_binom_loss, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 0.8) +
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Loss, Binomial Enrichment, Muscle",
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  ylim(0,5.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  facet_wrap(~cellType, 
             ncol = 2, 
             labeller = as_labeller(labels)) 




# 2.0 Combined nodes ----
nodeNames <- totalNodeRegions_combNodes$name %>% unique()
## 2.1 Neurons ----
### No facet wrap ----
labels = c("excitatoryneuron" = "Excitatory Neuron", 
           "inhibitoryneuron" = "Inhibitory Neuron")
ggplot(neuron_binom_combNodes, aes(x = MYA, y = enrich, color = cellType)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_line() + 
  theme_minimal() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Binomial Enrichment, Neurons", 
       color = "Cell Type", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  scale_color_manual(values = c("black", "darkgray"))

### Facet wrap ----
labels = c("excitatoryneuron" = "Excitatory Neuron", 
           "inhibitoryneuron" = "Inhibitory Neuron")
ggplot(neuron_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 0.8) +
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Binomial Enrichment, Neurons",
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
  ylim(0,2) + 
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 

### Combined neuron group ----
ggplot(neuronComb_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "gray", alpha = 0.5) +
  geom_line() + 
  theme_minimal() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  labs(title = "Binomial Enrichment, Neurons", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  ylim(0,2) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  scale_color_manual(values = c("black", "darkgray"))


## 2.2 Muscle ----
### No facet wrap ----
labels = c("cardiomyocyte_adult" = "Adult Cardiomyocyte", 
           "cardiomyocyte_developing" = "Developing Cardiomyocyte", 
           "skeletalmyocyte_adult" = "Adult Skeletal Myocyte", 
           "skeletalmyocyte_developing" = "Developing Skeletal Myocyte", 
           "smoothmuscle_adult" = "Adult Smooth Muscle")
ggplot(muscle_binom_combNodes, aes(x = MYA, y = enrich, color = cellType)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_line() + 
  theme_minimal() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Binomial Enrichment, Muscle", 
       color = "Cell Type", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  ylim(0,2) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

### Facet wrap ----
labels = c("cardiomyocyte_adult" = "Cardiomyocyte, Adult", 
           "cardiomyocyte_developing" = "Cardiomyocyte, Developing", 
           "skeletalmyocyte_adult" = "Skeletal Myocyte, Adult", 
           "skeletalmyocyte_developing" = "Skeletal Myocyte, Developing", 
           "smoothmuscle_adult" = "Smooth Muscle, Adult")
ggplot(muscle_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "gray", alpha = 0.5) +
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Binomial Enrichment, Muscle",
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  ylim(0,2.3) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  facet_wrap(~cellType, 
             ncol = 2, 
             labeller = as_labeller(labels)) 

### Combined Muscle group ----
ggplot(muscleComb_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "gray", alpha = 0.5) +
  geom_line() + 
  theme_minimal() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  labs(title = "Binomial Enrichment, Muscle", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  ylim(0,2) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  scale_color_manual(values = c("black", "darkgray"))


## 2.3 Innate immune ----
### No facet ----
labels = c("macrophage_adult" = "Adult Macrophage", 
           "macrophage_developing" = "Developing Macrophage", 
           "mast_adult" = "Adult Mast Cell", 
           "microglia_adult" = "Adult Microglia", 
           "naturalkillert_adult" = "Adult Natural Killer T Cells")
ggplot(innateImmune_binom_combNodes, aes(x = MYA, y = enrich, color = cellType)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_line() + 
  theme_minimal() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Binomial Enrichment, Innate Immune Cells", 
       color = "Cell Type", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  ylim(0,2) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

### Facet wrapped ----
labels = c("macrophage_adult" = "Macrophage, Adult", 
           "macrophage_developing" = "Macrophage, Developing", 
           "mast_adult" = "Mast Cell, Adult", 
           "microglia_adult" = "Microglia, Adult", 
           "naturalkillert_adult" = "Natural Killer T Cell, Adult")
ggplot(innateImmune_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "gray", alpha = 0.5) +
  geom_line(linewidth = 1, color = "black") + 
  theme_minimal() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Binomial Enrichment, Innate Immune Cells", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  ylim(0,2) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  facet_wrap(~cellType, 
             ncol = 2, 
             labeller = as_labeller(labels)) 

### Combined innate immune group ----
ggplot(innateImmuneComb_binom_combNodes, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "gray", alpha = 0.5) +
  geom_line() + 
  theme_minimal() + 
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  labs(title = "Binomial Enrichment, Innate Immune Cells", 
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  ylim(0,2) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  scale_color_manual(values = c("black", "darkgray"))

