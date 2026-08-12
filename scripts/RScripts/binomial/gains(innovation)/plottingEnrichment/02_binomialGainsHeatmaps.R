# 1.0 All cells and nodes ----
## 1.1 Binomial test estimates heatmap ----
pheatmap(binomPrMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Probabilities", 
         fontsize_row = 3)

## 1.2 Binomial probability heatmap, scaled by column----
## Scaling is Z scores based on column average
pheatmap(binomPrMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Probabilities, scaled by column", 
         fontsize_row = 3, 
         scale = "column")

## 1.3 Binomial test pvalues ----
## Scaling is Z scores based on column average
pheatmap(binomPvalMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Test P-values", 
         fontsize_row = 3)
binomSigMat <- binomSigMat*1
pheatmap(binomSigMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Test Significance at 95% Confidence", 
         fontsize_row = 3)

## 1.4 Binomial enrichments ----
pheatmap(enrichMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Enrichments", 
         fontsize_row = 3)
pheatmap(enrichMat, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Enrichments, Normalized by Node", 
         fontsize_row = 3, 
         scale = "column")


# 2.0 Combine nodes -1 to 3 ----
## 2.1 Binomial test estimates heatmap ----
pheatmap(binomPrMat_combNodes, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Probabilities", 
         fontsize_row = 3)

## 2.2 Binomial probability heatmap, scaled by column----
## Scaling is Z scores based on column average
pheatmap(binomPrMat_combNodes, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Probabilities, scaled by column", 
         fontsize_row = 3, 
         scale = "column")

## 1.4 Binomial enrichments ----
pheatmap(enrichMat_combNodes, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Enrichments", 
         fontsize_row = 3)
pheatmap(enrichMat_combNodes, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Binomial Enrichments, Normalized by Node", 
         fontsize_row = 3, 
         scale = "column")
