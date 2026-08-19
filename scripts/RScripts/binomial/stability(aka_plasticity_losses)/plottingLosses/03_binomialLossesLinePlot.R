# 1.0 All Nodes----
## FACET WRAP ----
### Neurons ----
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
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 

### Cardiomyocytes ----
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
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  facet_wrap(~cellType, 
             ncol = 1, 
             labeller = as_labeller(labels)) 

### All Muscle ----
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
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  facet_wrap(~cellType, 
             ncol = 2, 
             labeller = as_labeller(labels)) 

### All Innate Immune ----
innateImmune_binom_loss$MYA <- MYA_losses
labels = c("macrophage_adult" = "Macrophage, Adult", 
           "macrophage_developing" = "Macrophage, Developing", 
           "mast_adult" = "Mast Cell, Adult", 
           "microglia_adult" = "Microglia, Adult", 
           "naturalkillert_adult" = "Natural Killer T Cell, Adult")
ggplot(innateImmune_binom_loss, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 0.8) +
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Loss, Binomial Enrichment, Innate Immune Cells",
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  facet_wrap(~cellType, 
             ncol = 2, 
             labeller = as_labeller(labels)) 

### Placenta ----
labels = c("placentalNeuron" = "Placental Neuron", 
           "fibroPlacental" = "Placental Fibroblast", 
           "macrophagePlacental" = "Placental Macrophage", 
           "extravillousTrophoblast" = "Extravillous Trophoblast", 
           "syncitiotrophoblastCytotrophoblast" = "Syncitiotrophoblast and Cytotrophoblast", 
           "endothelialPlacental" = "Placental Endothelial")
ggplot(placenta_binom_loss, aes(x = MYA, y = enrich)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI), fill = "darkgray", alpha = 0.8) +
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  scale_x_reverse() +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Loss, Binomial Enrichment, Placenta",
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  facet_wrap(~cellType, 
             ncol = 2, 
             labeller = as_labeller(labels)) 



## NO FACET ----
### Neurons ----
neuron_binom_loss$MYA <- MYA_losses
labels = c("excitatoryneuron" = "Excitatory Neuron", 
           "inhibitoryneuron" = "Inhibitory Neuron")
ggplot(neuron_binom_loss, aes(x = MYA, y = enrich, color = cellType)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI, fill = cellType), alpha = 0.25, linetype = 0) +
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  scale_x_reverse() +
  scale_color_discrete(labels = labels) +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Loss, Binomial Enrichment, Neurons",
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

### Muscle ----
labels = c("cardiomyocyte_adult" = "Cardiomyocyte, Adult", 
           "cardiomyocyte_developing" = "Cardiomyocyte, Developing", 
           "skeletalmyocyte_adult" = "Skeletal Myocyte, Adult", 
           "skeletalmyocyte_developing" = "Skeletal Myocyte, Developing", 
           "smoothmuscle_adult" = "Smooth Muscle, Adult")
ggplot(muscle_binom_loss, aes(x = MYA, y = enrich, color = cellType)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI, fill = cellType), alpha = 0.25, linetype = 0) +
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  scale_x_reverse() +
  scale_color_discrete(labels = labels) +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Loss, Binomial Enrichment, Muscle",
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

### Innate immune ----
labels = c("macrophage_adult" = "Macrophage, Adult", 
           "macrophage_developing" = "Macrophage, Developing", 
           "mast_adult" = "Mast Cell, Adult", 
           "microglia_adult" = "Microglia, Adult", 
           "naturalkillert_adult" = "Natural Killer T Cell, Adult")
ggplot(innateImmune_binom_loss, aes(x = MYA, y = enrich, color = cellType)) + 
  geom_hline(yintercept = 1, linetype = "longdash") + 
  geom_ribbon(aes(ymin = enrLowerCI, ymax = enrUpperCI, fill = cellType), alpha = 0.25, linetype = 0) +
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  scale_x_reverse() +
  scale_color_discrete(labels = labels) +
  geom_point(aes(shape = sig), size = 3) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
  scale_color_discrete(labels = labels) +
  labs(title = "Loss, Binomial Enrichment, Innate Immune Cells",
       x = "Million Years Ago (MYA)", 
       y = "Enrichment", 
       shape = paste0("p < ", 0.05)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

