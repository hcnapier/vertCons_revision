#!/usr/bin/env bash
# 05_souporcell.sh
#
# Run Souporcell on the final annotated BAM + filtered barcode list to flag
# genotype-based clusters / cross-sample contamination (as used in the
# paper to address potential ambient RNA / cross-individual contamination).
#
# Souporcell needs:
#   - a coordinate-sorted, indexed BAM with CB tags (use final.bam sorted
#     from step 6 of 02_process_run.sh -- resort it since PISA corr's
#     output is not guaranteed to stay coordinate sorted)
#   - the filtered barcodes.tsv (from 04_droplet_filter.R output)
#   - the reference fasta used for alignment
#   - an estimate of the number of individuals/genotypes pooled in the run (-k)
#
# Usage:
#   ./05_souporcell.sh <run_dir> <genome_fa> <filtered_barcodes.tsv> <k_genotypes> <outdir> [threads]
#
# run_dir must contain final.bam from 02_process_run.sh.
# Requires souporcell (https://github.com/wheaton5/souporcell) available,
# typically easiest via its Singularity/Docker image.

set -euo pipefail

RUN_DIR=${1:?dir containing final.bam}
GENOME_FA=${2:?reference fasta used for STAR alignment}
FILTERED_BARCODES=${3:?filtered barcodes.tsv (unzipped, one CB per line)}
K=${4:?expected number of pooled genotypes/individuals}
OUTDIR=${5:?souporcell output dir}
THREADS=${6:-16}

mkdir -p "${OUTDIR}"

echo "[$(date)] Re-sorting final BAM by coordinate for souporcell"
samtools sort -@ "${THREADS}" -o "${RUN_DIR}/final.sorted.bam" "${RUN_DIR}/final.bam"
samtools index "${RUN_DIR}/final.sorted.bam"

echo "[$(date)] Running souporcell_pipeline.py"
souporcell_pipeline.py \
  -i "${RUN_DIR}/final.sorted.bam" \
  -b "${FILTERED_BARCODES}" \
  -f "${GENOME_FA}" \
  -t "${THREADS}" \
  -o "${OUTDIR}" \
  -k "${K}"

echo "[$(date)] Souporcell done. Cluster/genotype assignments in ${OUTDIR}/clusters.tsv"
