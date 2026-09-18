#!/bin/bash
#SBATCH --job-name=gethg38Cov
#SBATCH --output=logs/hg38Cov_%A_%a.logs
#SBATCH --time=5:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-59%2

set -euo pipefail

cd /work/hcn4/260630_vertCons_wd/hg38Coverage/

fbPath="/hpc/group/vertgenlab/cl454/bin/x86_64/" 
archive="/hpc/group/vertgenlab/christi/vertCons/zippedAlignments/hg38.wholeGenomeAlignments.60way.tar.gz"
filelist="/work/hcn4/260630_vertCons_wd/hg38Coverage/filelist.txt"

file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$filelist")

if [[ -z "$file" ]]; then
    echo "No file found for array index $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

echo "Processing: $file"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

bed="$tmpdir/input.bed"

tar -xzOf "$archive" "$file" |
    awk 'BEGIN { OFS="\t" } NF { print $1, $2, $3, $4 }' > "$bed"

"${fbPath}/featureBits" hg38 "$bed"