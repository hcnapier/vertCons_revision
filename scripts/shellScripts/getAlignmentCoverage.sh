#!/bin/bash
#SBATCH --job-name=gethg38Cov
#SBATCH --output=logs/hg38Cov_%A_%a.logs
#SBATCH --time=5:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-59

set -euo pipefail

cd /work/hcn4/260630_vertCons_wd/hg38Coverage/

fbPath="/hpc/group/vertgenlab/cl454/bin/x86_64/" 
archive="/hpc/group/vertgenlab/christi/vertCons/zippedAlignments/hg38.wholeGenomeAlignments.60way.tar.gz"
tmpdir=$(mktemp -d)
fifo="$tmpdir/stdin.bed"
filelist="/work/hcn4/260630_vertCons_wd/hg38Coverage/filelist.txt"

mkfifo "$fifo"

# Get the file for this array task
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$filelist")

if [[ -z "$file" ]]; then
    echo "No file found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Processing: $file"

tmpdir=$(mktemp -d)
fifo="$tmpdir/stdin.bed"
mkfifo "$fifo"

# Extract, normalize whitespace to tabs, feed into FIFO
tar -xzOf "$archive" "$file" | awk '{$1=$1}1' OFS='\t' > "$fifo" &

"${fbPath}/featureBits" hg38 "$fifo"

wait
rm -rf "$tmpdir"