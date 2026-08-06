#!/bin/bash
#SBATCH --job-name=annotate_peaks
#SBATCH --output=logs/annotate_%A_%a.out
#SBATCH --error=logs/annotate_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-17

inDir="/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/data/coverageCalls/coverageCalls.hg38"
outDir="/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/data/regionAnnotations"
gencodeDir="/work/hcn4/260630_vertCons_wd"

mkdir -p $outDir logs

module load Bedtools/2.30.0 

# Build an array of input files (sorted for reproducibility)
bedList=($(ls ${inDir}/node*.bed | sort))

# Select this task's file based on array index
currInBed=${bedList[$SLURM_ARRAY_TASK_ID]}
filename=${currInBed##*/} 
currName=${filename%.bed}

#currName=$(basename $inDir .bed)
outOverlap=$outDir/${currName}_geneOverlap.bed
outIntergenic=${outDir}/${currName}_intergenic.bed

echo $currName

echo "Processing: $currInBed"

bedtools intersect -a ${currInBed} -b ${gencodeDir}/gencode_full.bed -wa -u > $outOverlap
bedtools intersect -a ${currInBed} -b ${gencodeDir}/gencode_full.bed -v > $outIntergenic

echo "done!"
