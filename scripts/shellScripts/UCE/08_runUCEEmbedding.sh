#!/bin/bash
#SBATCH --job-name=uceEmbedding
#SBATCH --output=logs/uceEmbedding_%A_%a.out
#SBATCH --error=logs/uceEmbedding_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --gres=gpu:5000_ada:2
#SBATCH --partition=gpu-common
#SBATCH --array=1-9

# ============================== CONFIG ======================================

# Directory containing your .h5ad datasets to embed. Filenames are expected
# to start with the species name followed by a dot, e.g.
#   bos_taurus.heart_atlas.h5ad
# matching the same naming convention used for the proteome FASTAs earlier.
inDir="/work/hcn4/260630_vertCons_wd/scTrx/annDataObjs/forUCE"

# Where the finished, UCE-embedded .h5ad files get written
outDir="/work/hcn4/260630_vertCons_wd/scTrx/uceEmbedded"

# Where UCE itself lives (the cloned repo, containing eval_single_anndata.py)
uceDir="/hpc/group/vertgenlab/hailey/software/UCE"

# Where create_new_species_files.py wrote its three output files per species
speciesFilesDir="/work/hcn4/260630_vertCons_wd/scTrx/uceSpeciesFiles"

# Path to the model weights .torch file
modelLoc="/hpc/group/vertgenlab/hailey/software/UCE/model_files/33l_8ep_1024t_1280.torch"

# =============================================================================

export LD_LIBRARY_PATH=/hpc/group/vertgenlab/hailey/software/miniconda3/envs/scrna/lib/gcc/x86_64-conda-linux-gnu/15.2.0:$LD_LIBRARY_PATH
export HF_HOME=/work/hcn4/260630_vertCons_wd/scTrx/esm2
export HF_HUB_DISABLE_XET=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

source /hpc/group/vertgenlab/hailey/software/miniconda3/etc/profile.d/conda.sh
conda activate scrna

mkdir -p "${outDir}" logs

# Build an array of input datasets (sorted for reproducibility)
datasetList=($(ls ${inDir}/*.h5ad | sort))
currInDataset=${datasetList[$SLURM_ARRAY_TASK_ID]}
filename="${currInDataset##*/}"
speciesName="${filename%%.*}"

echo "Dataset:     ${currInDataset}"
echo "Species:     ${speciesName}"

# Species-specific files produced by create_new_species_files.py
chromCsv="${speciesFilesDir}/${speciesName}_to_chrom_pos.csv"
tokenFile="${speciesFilesDir}/${speciesName}_pe_tokens.torch"
offsetPkl="${speciesFilesDir}/${speciesName}_offsets.pkl"

for f in "${chromCsv}" "${tokenFile}" "${offsetPkl}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: expected species file not found: ${f}" >&2
        echo "Did you run create_new_species_files.py for '${speciesName}' yet?" >&2
        exit 1
    fi
done

# Recompute CHROM_TOKEN_OFFSET at runtime — it's deterministic from the
# token file and chrom position csv, so there's no separate file to track:
#   CHROM_TOKEN_OFFSET = (rows in pe_tokens.torch) - (unique chromosomes)
chromTokenOffset=$(python -c "
import torch, pandas as pd
pe = torch.load('${tokenFile}')
chrom_pos = pd.read_csv('${chromCsv}')
n_uniq_chrom = chrom_pos['chromosome'].nunique()
print(pe.shape[0] - n_uniq_chrom)
")
echo "CHROM_TOKEN_OFFSET: ${chromTokenOffset}"

# Reminder: before this will work, you must have already (see
# create_new_species_files.py's printed instructions):
#   1. Added a row for ${speciesName} to
#      ${uceDir}/model_files/new_species_protein_embeddings.csv
#   2. Added ${speciesName} to the species dict around line 247 of
#      ${uceDir}/data_proc/data_utils.py

cd "${uceDir}"

accelerate launch eval_single_anndata.py \
    --adata_path "${currInDataset}" \
    --dir "${outDir}/" \
    --species "${speciesName}" \
    --model_loc "${modelLoc}" \
    --batch_size 8 \
    --CHROM_TOKEN_OFFSET "${chromTokenOffset}" \
    --spec_chrom_csv_path "${chromCsv}" \
    --token_file "${tokenFile}" \
    --offset_pkl_path "${offsetPkl}" \
    --multi_gpu True \
    --nlayers 33

conda deactivate
