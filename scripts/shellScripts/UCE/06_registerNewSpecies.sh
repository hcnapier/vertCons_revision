#!/bin/bash
#SBATCH --job-name=registerNewSpecies
#SBATCH --output=logs/registerNewSpecies_%A_%a.out
#SBATCH --error=logs/resgisterNewSpecies_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --array=0-9
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G


export LD_LIBRARY_PATH=/hpc/group/vertgenlab/hailey/software/miniconda3/envs/scrna/lib/gcc/x86_64-conda-linux-gnu/15.2.0:$LD_LIBRARY_PATH

fastaDir="/work/hcn4/260630_vertCons_wd/scTrx/proteinFastas"
esm2EmbeddingDir="/work/hcn4/260630_vertCons_wd/scTrx/esm2Embeddings"
pythonDir="/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/pythonScripts"
uceDir="/hpc/group/vertgenlab/hailey/software/UCE" 

# Build an array of input files (sorted for reproducibility)
fastaList=($(ls ${fastaDir}/*.fa | sort))
embeddingList=($(ls ${esm2EmbeddingDir}/*.pt | sort))

# Select this task's file based on array index
currInFasta=${fastaList[$SLURM_ARRAY_TASK_ID]}
filename="${currInFasta##*/}"
speciesName="${filename%%.*}"
currInEmbed=${embeddingList[$SLURM_ARRAY_TASK_ID]}

source  /hpc/group/vertgenlab/hailey/software/miniconda3/etc/profile.d/conda.sh
conda activate scrna

export HF_HOME=/work/hcn4/260630_vertCons_wd/scTrx/esm2
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python register_uce_species.py \
        --species_name ${speciesName} \
        --protein_embeddings_path ${esm2EmbeddingDir}/${currInEmbed} \
        --csv_path ${uceDir}/model_files/new_species_protein_embeddings.csv

conda deactivate
