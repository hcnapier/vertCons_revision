#!/usr/bin/env Rscript
# 04_droplet_filter.R
#
# Remove empty droplets from a PISA raw gene expression matrix using
# DropletUtils, and write a filtered matrix in 10x MEX format.
#
# Usage:
#   Rscript 04_droplet_filter.R <raw_gene_expression_dir> <filtered_outdir> [method]
#
#   method: "inflection" (default, matches the barcodeRanks-inflection
#           approach used in the PISA tutorial) or "emptyDrops" (statistical
#           test against ambient profile; slower, gives an FDR).
#
# Requires: DropletUtils, Matrix. Optionally Yano (for ReadPISA); falls back
# to reading MEX directly with Matrix::readMM if Yano isn't installed.

suppressPackageStartupMessages({
  library(DropletUtils)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 04_droplet_filter.R <raw_gene_expression_dir> <filtered_outdir> [inflection|emptyDrops]")
}
raw_dir <- args[1]
out_dir <- args[2]
method  <- ifelse(length(args) >= 3, args[3], "inflection")

message("[", Sys.time(), "] Reading raw matrix from ", raw_dir)

read_pisa_mex <- function(dir) {
  if (requireNamespace("Yano", quietly = TRUE)) {
    return(Yano::ReadPISA(dir))
  }
  # Fallback: read standard 10x-style MEX triplet directly
  mtx  <- file.path(dir, "matrix.mtx.gz")
  bcs  <- file.path(dir, "barcodes.tsv.gz")
  feat <- file.path(dir, "features.tsv.gz")
  m <- as(Matrix::readMM(gzfile(mtx)), "CsparseMatrix")
  barcodes <- readLines(gzfile(bcs))
  features <- read.delim(gzfile(feat), header = FALSE)
  rownames(m) <- features[[1]]
  colnames(m) <- barcodes
  m
}

counts <- read_pisa_mex(raw_dir)
message("Raw matrix: ", nrow(counts), " features x ", ncol(counts), " barcodes")

sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))

if (method == "emptyDrops") {
  message("[", Sys.time(), "] Running emptyDrops()")
  set.seed(100)
  ed <- emptyDrops(counts(sce))
  is_cell <- which(ed$FDR <= 0.01)
  filter_bcs <- colnames(sce)[is_cell]
} else {
  message("[", Sys.time(), "] Running barcodeRanks() and taking inflection point")
  br.out <- barcodeRanks(sce)
  inflection <- S4Vectors::metadata(br.out)$inflection
  filter_bcs <- rownames(br.out)[br.out$total > inflection]
}

message("Cells retained: ", length(filter_bcs), " / ", ncol(counts))

filter_counts <- counts[, filter_bcs, drop = FALSE]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write10xCounts(out_dir, filter_counts, overwrite = TRUE)

message("[", Sys.time(), "] Filtered matrix written to ", out_dir)
