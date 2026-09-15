#!/bin/bash
#SBATCH --job-name=collapseNewSpecies
#SBATCH --output=logs/collapseNewSpecies_%A_%a.out
#SBATCH --error=logs/collapseNewSpecies_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G

export LD_LIBRARY_PATH=/hpc/group/vertgenlab/hailey/software/miniconda3/envs/scrna/lib/gcc/x86_64-conda-linux-gnu/15.2.0:$LD_LIBRARY_PATH

uceDir="/hpc/group/vertgenlab/hailey/software/UCE" 
tmpDir="/work/hcn4/260630_vertCons_wd/scTrx/uce_species_temp"
pythonDir="/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/pythonScripts"

source  /hpc/group/vertgenlab/hailey/software/miniconda3/etc/profile.d/conda.sh
conda activate scrna

export HF_HOME=/work/hcn4/260630_vertCons_wd/scTrx/esm2
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python ${pythonDir}/collapse_new_species_csv.py \
        --temp_dir ${tmpDir} \
        --csv_path ${uceDir}/model_files/new_species_protein_embeddings.csv
