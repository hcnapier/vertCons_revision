#!/usr/bin/env bash
# 03_run_all.sh
#
# Loop over samples.tsv and process every BGIseq run listed there.
#
# Usage:
#   ./03_run_all.sh samples.tsv star_index_root/ gtf_root/ outdir/ [threads] [whitelist]
#
# star_index_root/<species>/   must exist (built by 01_build_star_index.sh)
# gtf_root/<species>.gtf       must exist
#
# samples.tsv columns: run_id  species  reference_key  fastq_R1  fastq_R2

set -euo pipefail

SAMPLES=${1:?samples.tsv}
STAR_INDEX_ROOT=${2:?dir containing one STAR index subfolder per species}
GTF_ROOT=${3:?dir containing <species>.gtf files}
OUTDIR=${4:?output root dir}
THREADS=${5:-16}
WHITELIST=${6:-}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while IFS=$'\t' read -r RUN_ID SPECIES REF_KEY FQ1 FQ2; do
  [[ "${RUN_ID}" =~ ^#.*$ || -z "${RUN_ID}" ]] && continue

  STAR_INDEX="${STAR_INDEX_ROOT}/${SPECIES}"
  GTF="${GTF_ROOT}/${SPECIES}.gtf"

  if [[ ! -d "${STAR_INDEX}" ]]; then
    echo "ERROR: STAR index for species '${SPECIES}' not found at ${STAR_INDEX}" >&2
    exit 1
  fi
  if [[ ! -f "${GTF}" ]]; then
    echo "ERROR: GTF for species '${SPECIES}' not found at ${GTF}" >&2
    exit 1
  fi

  echo "##### Processing ${RUN_ID} (${SPECIES}) #####"
  "${SCRIPT_DIR}/02_process_run.sh" \
    "${RUN_ID}" "${FQ1}" "${FQ2}" "${STAR_INDEX}" "${GTF}" "${OUTDIR}" "${THREADS}" "${WHITELIST}"

done < "${SAMPLES}"

echo "All runs in ${SAMPLES} processed. Next: run 04_droplet_filter.R per run, then 05_souporcell.sh."
