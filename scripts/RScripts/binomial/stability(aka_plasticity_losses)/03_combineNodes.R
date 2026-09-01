# 1.0 Combine nodes -1-3 ----
nLossMat <- nodeLosses %>%
  select(CellType, nSuccess, Node) %>%
  pivot_wider(
    names_from = Node,   
    values_from = nSuccess 
  ) 
toCombine <- c("1", "2", "3")
combinedNodes <- nLossMat[, toCombine, drop = FALSE]
comb1_3 <- rowSums(combinedNodes)
keep <- nLossMat[, !(colnames(nLossMat) %in% toCombine), drop = FALSE]
keep <- as.data.frame(keep)
rownames(keep) <- keep$CellType
keep$CellType <- NULL
nLossMat_combNodes <- cbind(comb1_3, keep)
nLossMat_combNodes <- as.matrix(nLossMat_combNodes)
nLossLong_combNodes <- as.data.frame(as.table(nLossMat_combNodes))
names(nLossLong_combNodes) <- c("CellType", "Node", "nLosses")
nregions <- nodeLosses %>%
  filter(Node == 1) %>%
  select(nRegions)
nLossLong_combNodes$nRegions <- rep(nregions$nRegions, 14)
nLossLong_combNodes$nSuccess <- nLossLong_combNodes$nLosses
# Save 
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(nLossLong_combNodes, "nLossLong_combNodes.rds")

# Get other datasets ----
# Get total losses per node
totalNodeLosses_combNodes <- nLossLong_combNodes %>%
  select(Node, nLosses) %>%
  group_by(Node) %>%
  summarize(nLosses = sum(nLosses))

totalNodeRegions_losses_combNodes <- nLossLong_combNodes %>%
  select(Node, nRegions) %>%
  group_by(Node) %>%
  summarize(nRegions = sum(nRegions))
totalNodeRegions_losses_combNodes$name <- paste("Node", totalNodeRegions_losses_combNodes$Node, sep = "")
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(totalNodeRegions_losses_combNodes, "totalNodeRegions_losses_combNodes.rds")


# 1.0 Binomial enrichment, comb nodes ----
nodeNames_losses <- c("Nodecomb1_3", "Node4", "Node5", "Node6", "Node7", "Node8", "Node9", "Node10", "Node11", "Node12", "Node13", "Node14", "Node15", "Node16")

## 1.1 Neurons ----
neuronTypes <- c("excitatoryneuron", "inhibitoryneuron")
neuronGroupBinom_loss_combNodes <- subsetAndBinomCellTypeList(cellTypeList = neuronTypes, 
                                                    fullDF = nLossLong_combNodes, 
                                                    totalNodeRegions = totalNodeRegions_losses_combNodes,
                                                    totalNodeLosses = totalNodeLosses_combNodes,
                                                    nodeNames = nodeNames_losses, 
                                                    groupName = "Neurons", 
                                                    MYA = MYA_losses, 
                                                    type = "loss")
# Combine for plotting
neuron_binom_loss_combNodes <- do.call(rbind, neuronGroupBinom_loss)
setwd("~/Work/VertGenLab/Projects/vertCons/code/vertCons_revision/scripts/RScripts/rData")
saveRDS(neuron_binom_loss_combNodes, "neuron_bionm_loss.rds")

