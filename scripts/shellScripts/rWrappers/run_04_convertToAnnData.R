#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --mail-user=hailey.napier@duke.edu
#SBATCH -c 1
#SBATCH --mem=100gb
#SBATCH -N 1
#SBATCH -n 1

module load /opt/apps/modules-bak/R/4.6.0 

cd /hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/RScripts/scTrx

Rscript 04_convertToAnnData.R