#!/bin/bash
#SBATCH --job-name=getUCEDistances
#SBATCH --output=logs/uceDistances_%A_%a.out
#SBATCH --error=logs/uceDistances_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

export LD_LIBRARY_PATH=/hpc/group/vertgenlab/hailey/software/miniconda3/envs/scrna/lib/gcc/x86_64-conda-linux-gnu/15.2.0:$LD_LIBRARY_PATH

source  /hpc/group/vertgenlab/hailey/software/miniconda3/etc/profile.d/conda.sh
conda activate scrna

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

pythonDir="/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/pythonScripts"

python ${pythonDir}/uce_celltype_distance_comparison.py

