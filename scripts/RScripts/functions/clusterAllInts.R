clusterAllInts <- function(
    seurat_obj,
    dims       = 1:30,
    resolution = 0.8) {
  
  # Validate input
  stopifnot("Input must be a Seurat object" = inherits(seurat_obj, "Seurat"))
  
  integration_reductions <- grep("^integrated", Reductions(seurat_obj), value = TRUE)
  
  if (length(integration_reductions) == 0) {
    stop("No reductions found with names beginning with 'integrated'.")
  }
  
  message("Found ", length(integration_reductions), " integration reduction(s): ",
          paste(integration_reductions, collapse = ", "))
  
  for (reduction in integration_reductions) {
    message("\nProcessing: ", reduction)
    
    suffix       <- gsub("\\.", "_", reduction)
    nn_name      <- paste0("nn_",       suffix)
    snn_name     <- paste0("snn_",      suffix)
    umap_name    <- paste0("umap_",     suffix)
    cluster_col  <- paste0("clusters_", suffix)
    umap_key     <- paste0(gsub("_", "", suffix), "UMAP_")
    
    seurat_obj <- FindNeighbors(
      seurat_obj,
      reduction  = reduction,
      dims       = dims,
      graph.name = c(nn_name, snn_name),
      verbose    = FALSE
    )
    
    seurat_obj <- FindClusters(
      seurat_obj,
      graph.name   = snn_name,
      resolution   = resolution,
      cluster.name = cluster_col,
      verbose      = FALSE
    )
    
    seurat_obj <- RunUMAP(
      seurat_obj,
      reduction      = reduction,
      dims           = dims,
      reduction.name = umap_name,
      reduction.key  = umap_key,
      verbose        = FALSE
    )
    
    message("  Clusters : ", cluster_col)
    message("  UMAP     : ", umap_name)
  }
  
  return(seurat_obj)
}