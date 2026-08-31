#!/usr/bin/env bash
# 02_process_run.sh
#
# Process one BGIseq scRNA-seq run (paired FASTQ) into a raw cell x gene
# count matrix, following the method:
#
#   Read 1 (30bp): bases 1-20 = cell barcode, bases 21-30 = UMI
#   Read 2 (100bp): transcript sequence
#
#   1. PISA parse   -> FASTQ+, correcting CB against allow-list (<=1 mismatch)
#   2. STAR         -> splice-aware alignment to the species reference
#   3. PISA sam2bam -> BAM + alignment summary, restoring CB/UR tags
#   4. PISA anno    -> gene annotation (GN tag) from GTF
#   5. samtools sort -> coordinate sort (required for robust UMI correction)
#   6. PISA corr    -> UMI correction within CB+GN blocks (<=1 mismatch) -> UB tag
#   7. PISA count   -> raw cell x gene matrix (MEX format)
#
# Usage:
#   ./02_process_run.sh <run_id> <fastq_R1> <fastq_R2> <star_index_dir> <gtf> <outdir> [threads] [whitelist]
#
# whitelist: plain-text file, one cell barcode per line (20bp), used to
#            correct CB with Hamming distance <= 1. If your barcodes come
#            from a fixed DNBelab/BGI C4 allow-list, point to that file.
#            If you don't have one, omit it and PISA will still record raw
#            CB without allow-list correction (see NOTE below).

set -euo pipefail

RUN_ID=${1:?run id, e.g. SRR12345678}
FQ1=${2:?fastq R1 (30bp: CB 1-20, UMI 21-30)}
FQ2=${3:?fastq R2 (100bp transcript)}
STAR_INDEX=${4:?STAR genome index dir for this species}
GTF=${5:?GTF annotation for this species}
OUTDIR=${6:?output directory}
THREADS=${7:-16}
WHITELIST=${8:-}

mkdir -p "${OUTDIR}/${RUN_ID}"
cd "${OUTDIR}/${RUN_ID}"

echo "[$(date)] === ${RUN_ID}: 1/7 PISA parse (barcode/UMI extraction) ==="
if [[ -n "${WHITELIST}" ]]; then
  RULE="CB,R1:1-20,${WHITELIST},CB,1;UR,R1:21-30;R1,R2"
else
  echo "  NOTE: no whitelist supplied -- CB will be recorded as-is without allow-list correction."
  RULE="CB,R1:1-20;UR,R1:21-30;R1,R2"
fi

PISA parse \
  -rule "${RULE}" \
  -1 reads.fq \
  -report parse_report.csv \
  -q 5 \
  -t "${THREADS}" \
  "${FQ1}" "${FQ2}"

echo "[$(date)] === ${RUN_ID}: 2/7 STAR alignment ==="
STAR \
  --runThreadN "${THREADS}" \
  --genomeDir "${STAR_INDEX}" \
  --readFilesIn reads.fq \
  --outSAMtype SAM \
  --outFileNamePrefix star_ \
  --outSAMunmapped Within

# STAR writes star_Aligned.out.sam

echo "[$(date)] === ${RUN_ID}: 3/7 PISA sam2bam (restore CB/UR tags, mapq adjust) ==="
PISA sam2bam \
  -report alignment_report.csv \
  -@ "${THREADS}" \
  -adjust-mapq \
  -gtf "${GTF}" \
  -o aln.bam \
  star_Aligned.out.sam

echo "[$(date)] === ${RUN_ID}: 4/7 PISA anno (gene annotation -> GN tag) ==="
PISA anno \
  -gtf "${GTF}" \
  -o anno.bam \
  -t "${THREADS}" \
  -report anno_report.csv \
  aln.bam

echo "[$(date)] === ${RUN_ID}: 5/7 samtools sort ==="
samtools sort -@ "${THREADS}" -o sorted.bam anno.bam
samtools index sorted.bam

echo "[$(date)] === ${RUN_ID}: 6/7 PISA corr (UMI correction, <=1 mismatch, per CB+GN) ==="
PISA corr \
  -tag UR \
  -new-tag UB \
  -tags-block CB,GN \
  -cr \
  -o final.bam \
  -@ "${THREADS}" \
  sorted.bam

echo "[$(date)] === ${RUN_ID}: 7/7 PISA count (raw cell x gene matrix) ==="
mkdir -p raw_gene_expression
PISA count \
  -tags CB \
  -anno-tag GN \
  -umi UB \
  -outdir raw_gene_expression \
  -@ "${THREADS}" \
  final.bam

echo "[$(date)] === ${RUN_ID}: cleanup intermediate SAM ==="
rm -f star_Aligned.out.sam reads.fq

echo "[$(date)] ${RUN_ID} done. Raw matrix at ${OUTDIR}/${RUN_ID}/raw_gene_expression/"
