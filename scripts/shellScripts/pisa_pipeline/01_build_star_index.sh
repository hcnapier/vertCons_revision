#!/usr/bin/env bash
# 01_build_star_index.sh
#
# Build a STAR genome index for one species reference.
# Run once per species before processing any FASTQs for that species.
#
# Usage:
#   ./01_build_star_index.sh <species> <genome.fa> <annotation.gtf> <outdir> [threads] [sjdbOverhang]
#
# Notes:
#   - sjdbOverhang should be (read2_length - 1). Read 2 here is 100 bp -> use 99.
#   - RAM: mammalian genomes typically need ~32-64GB RAM for indexing.
#   - Download the correct FASTA/GTF for each species from ref/species_refs.tsv
#     before running this (URLs there are starting points -- confirm exact
#     assembly patch/version on the source site, since e.g. dog UU_Cfam_GSD_1.0
#     and rabbit UM_NZW_1.0 are NCBI RefSeq/GenBank assemblies whereas cow,
#     goat, pig and guinea pig are readily available on Ensembl).

set -euo pipefail

SPECIES=${1:?species label, e.g. cow}
GENOME_FA=${2:?path to genome fasta}
GTF=${3:?path to GTF annotation}
OUTDIR=${4:?output index directory}
THREADS=${5:-16}
SJDB_OVERHANG=${6:-99}

mkdir -p "${OUTDIR}"

echo "[$(date)] Building STAR index for ${SPECIES} -> ${OUTDIR}"

STAR \
  --runMode genomeGenerate \
  --runThreadN "${THREADS}" \
  --genomeDir "${OUTDIR}" \
  --genomeFastaFiles "${GENOME_FA}" \
  --sjdbGTFfile "${GTF}" \
  --sjdbOverhang "${SJDB_OVERHANG}"

echo "[$(date)] Done. Index written to ${OUTDIR}"
