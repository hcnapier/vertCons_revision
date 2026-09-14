#!/bin/bash
#SBATCH --job-name=getESM2Model
#SBATCH --time=00:30:00
#SBATCH --mem=36G
#SBATCH --cpus-per-task=1

cd /work/hcn4/260630_vertCons_wd/scTrx/esm2

export LD_LIBRARY_PATH=/hpc/group/vertgenlab/hailey/software/miniconda3/envs/scrna/lib/gcc/x86_64-conda-linux-gnu/15.2.0:$LD_LIBRARY_PATH

export HF_HOME=/work/hcn4/260630_vertCons_wd/scTrx/esm2

conda activate scrna

python -c "from transformers import AutoModel; AutoModel.from_pretrained('facebook/esm2_t48_15B_UR50D')"
python -c "from transformers import AutoTokenizer; AutoTokenizer.from_pretrained('facebook/esm2_t48_15B_UR50D')"
