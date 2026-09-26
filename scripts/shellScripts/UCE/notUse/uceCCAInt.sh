#!/bin/bash
#SBATCH --job-name=uceCCAInt
#SBATCH --output=logs/ccaInt_%A.out
#SBATCH --error=logs/ccaInt_%A.err
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

module load R/4.6.0 

pythonDir="/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/pythonScripts"
rDir="/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/RScripts/scTrx"

export LD_LIBRARY_PATH=/hpc/group/vertgenlab/hailey/software/miniconda3/envs/scrna/lib/gcc/x86_64-conda-linux-gnu/15.2.0:$LD_LIBRARY_PATH
source  /hpc/group/vertgenlab/hailey/software/miniconda3/etc/profile.d/conda.sh
conda activate scrna
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
#python ${pythonDir}/export_uce_for_seurat.py

Rscript ${rDir}/uce_cca_integration.R

python ${pythonDir}/import_get_distance_cca.py

conda deactivate