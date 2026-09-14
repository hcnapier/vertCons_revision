#!/bin/bash
#SBATCH --job-name=prot2ESM2
#SBATCH --output=logs/prot2ESM2_%A_%a.out
#SBATCH --error=logs/prot2ESM2_%A_%a.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --array=0-9
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu-common

export LD_LIBRARY_PATH=/hpc/group/vertgenlab/hailey/software/miniconda3/envs/scrna/lib/gcc/x86_64-conda-linux-gnu/15.2.0:$LD_LIBRARY_PATH

inDir="/work/hcn4/260630_vertCons_wd/scTrx/proteinFastas"
outDir="/work/hcn4/260630_vertCons_wd/scTrx/esm2Embeddings"
pythonDir="/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/pythonScripts"

# Build an array of input files (sorted for reproducibility)
fastaList=($(ls ${inDir}/*.fa | sort))

# Select this task's file based on array index
currInFasta=${fastaList[$SLURM_ARRAY_TASK_ID]}i
filename="${currInFasta##*/}"
speciesName="${filename%%.*}"

echo ${speciesName}

source  /hpc/group/vertgenlab/hailey/software/miniconda3/etc/profile.d/conda.sh
conda activate scrna

python ${pythonDir}/proteome_to_esm2_embeddings.py \
        --fasta ${currInFasta} \
        --species_name ${speciesName} \
        --out_dir ${outDir}/ \
        --gene_id_regex "gene:(\\S+)"

conda deactivate
