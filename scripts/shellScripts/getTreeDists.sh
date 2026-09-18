#!/bin/bash
#SBATCH --job-name=getTreeDists
#SBATCH --output=logs/treeDists_%A.logs
#SBATCH --time=05:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

allDists="/hpc/group/vertgenlab/cl454/src/phast/phast/bin/all_dists"
tree="/hpc/group/vertgenlab/christi/vertCons/trees/hg38.60way.nhs"
tmp="/work/hcn4/260630_vertCons_wd/60wayTree.dists.txt"
out="/work/hcn4/60630_vertCons_wd/hg38_60wayTree.dists.txt"
#set species=/work/cf189/runPairwiseAlignments/species.list

$allDists $tree > $tmp
head $out -n 59 > $out
#cat $out | cut -f1| grep -v "(total)" | uniq > $species
#tail -2 $out | head -1 | cut -f2 >> $species

echo DONE 