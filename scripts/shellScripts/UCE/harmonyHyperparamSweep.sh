#!/bin/bash
#SBATCH --job-name=harmonyHyperparamSweep
#SBATCH --output=logs/harmonyHyperparam_%A_%a.out
#SBATCH --error=logs/harmonyHyperparam_%A_%a.err
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=300G

export LD_LIBRARY_PATH=/hpc/group/vertgenlab/hailey/software/miniconda3/envs/scrna/lib/gcc/x86_64-conda-linux-gnu/15.2.0:$LD_LIBRARY_PATH
source  /hpc/group/vertgenlab/hailey/software/miniconda3/etc/profile.d/conda.sh
conda activate scrna
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

pythonDir="/hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/pythonScripts"
uceObj="/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/uce_distances_combined.h5ad"
outDir="/work/hcn4/260630_vertCons_wd/scTrx/uce_distances/harmonySweep"

python harmony_hyperparam_sweep.py ${uceObj} \
    --study-key study --species-key species --celltype-key cell_type \
    --tan-study Tan --covariate-sets study study,technology \
    --thetas 1 2 4 6 8 --lambdas 1 0.5 0.1 --subsample 60000 --apply-best \
    --outdir ${outDir}

conda deactivate