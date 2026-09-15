#!/bin/bash
#SBATCH --job-name=newSpeciesUCE
#SBATCH --output=logs/newSpeciesUCE_%A_%a.out
#SBATCH --error=logs/newSpeciesUCE_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --array=0-9
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G

export LD_LIBRARY_PATH=/hpc/group/vertgenlab/hailey/software/miniconda3/envs/scrna/lib/gcc/x86_64-conda-linux-gnu/15.2.0:$LD_LIBRARY_PATH

fastaDir="/work/hcn4/260630_vertCons_wd/scTrx/proteinFastas"
esm2EmbeddingDir="/work/hcn4/260630_vertCons_wd/scTrx/esm2Embeddings"
pythonDir="/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/pythonScripts"
outDir="/work/hcn4/260630_vertCons_wd/scTrx/uceSpeciesFiles" 

# Build an array of input files (sorted for reproducibility)
fastaList=($(ls ${fastaDir}/*.fa | sort))
embeddingList=($(ls ${esm2EmbeddingDir}/*.pt | sort))

# Build an array of NCBI taxonomy IDs
taxIDs=("9913" "9615" "9925" "10141" "9606" "9541" "10090" "9986" "10116" "9823")

# Select this task's file based on array index
currInFasta=${fastaList[$SLURM_ARRAY_TASK_ID]}
filename="${currInFasta##*/}"
speciesName="${filename%%.*}"
tmp="${filename#*.}"
assemblyName="${tmp%.pep}"
currInEmbed=${embeddingList[$SLURM_ARRAY_TASK_ID]}
currTaxID=${taxIDs[$SLURM_ARRAY_TASK_ID]}

echo ${speciesName}
echo ${assemblyName}

source  /hpc/group/vertgenlab/hailey/software/miniconda3/etc/profile.d/conda.sh
conda activate scrna

export HF_HOME=/work/hcn4/260630_vertCons_wd/scTrx/esm2
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python ${pythonDir}/create_new_species_files.py \
    --species_name ${speciesName} \
    --fasta ${currInFasta} \
    --protein_embeddings ${currInEmbed} \
    --assembly_name ${assemblyName} \
    --taxonomy_id ${currTaxID} \
    --all_tokens_path /hpc/group/vertgenlab/hailey/software/UCE/model_files/all_tokens.torch \
    --out_dir ${outDir}

conda deactivate
