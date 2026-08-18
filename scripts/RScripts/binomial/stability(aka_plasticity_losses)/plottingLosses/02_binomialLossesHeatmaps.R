# 1.0 All cells and nodes ----
## 1.1 Binomial test estimates heatmap ----
pheatmap(binomPrMat_loss, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Losses, Binomial Probabilities", 
         fontsize_row = 3)

## 1.2 Binomial probability heatmap, scaled by column----
## Scaling is Z scores based on column average
pheatmap(binomPrMat_loss, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Losses, Binomial Probabilities (Normalized by Node)", 
         fontsize_row = 3, 
         scale = "column")

## 1.3 Binomial test pvalues ----
## Scaling is Z scores based on column average
pheatmap(binomPvalMat_loss, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Test P-values", 
         fontsize_row = 3)
binomSigMat_loss <- binomSigMat_loss*1
pheatmap(binomSigMat_loss, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Test Significance at 95% Confidence", 
         fontsize_row = 3)

## 1.4 Binomial enrichments ----
pheatmap(enrichMat_loss, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Losses, Binomial Enrichments", 
         fontsize_row = 3)
pheatmap(enrichMat_loss, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Losses, Binomial Enrichments (Normalized by Node)", 
         fontsize_row = 3, 
         scale = "column")

