#!/bin/bash
#SBATCH --job-name=getWangData
#SBATCH --time=05:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

cd /work/hcn4/260630_vertCons_wd/thirdTriPlacenta_snATAC 
wget https://cell.ucsf.edu/snPlacenta/download/fragments.tar.gz
tar -xf fragments.tar.gz
