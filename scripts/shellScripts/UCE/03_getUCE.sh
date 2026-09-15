#!/bin/bash
#SBATCH --job-name=getTanData
#SBATCH --time=00:30:00
#SBATCH --mem=36G
#SBATCH --cpus-per-task=1

cd /hpc/group/vertgenlab/hailey/software
git clone https://github.com/snap-stanford/UCE.git

export LD_LIBRARY_PATH=/hpc/group/vertgenlab/hailey/software/miniconda3/envs/scrna/lib/gcc/x86_64-conda-linux-gnu/15.2.0:$LD_LIBRARY_PATH

source  /hpc/group/vertgenlab/hailey/software/miniconda3/etc/profile.d/conda.sh
conda activate scrna

HF_HUB_DISABLE_XET=1
export HF_HOME=/work/hcn4/260630_vertCons_wd/scTrx/esm2
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python -c "
from huggingface_hub import hf_hub_download
import shutil
path = hf_hub_download(repo_id='minwoosun/uce-misc', filename='all_tokens.torch')
shutil.copy(path, '/hpc/group/vertgenlab/hailey/software/UCE/model_files/all_tokens.torch')
"
python -c "
from huggingface_hub import hf_hub_download
import shutil
for fname in ['species_chrom.csv', 'species_offsets.pkl']:
    path = hf_hub_download(repo_id='minwoosun/uce-misc', filename=fname)
    shutil.copy(path, f'/path/to/UCE/model_files/{fname}')
"
