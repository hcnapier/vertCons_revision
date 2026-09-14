#!/bin/bash
#SBATCH --job-name=getTanData
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

module load samtools/1.21 

cd /work/hcn4/260630_vertCons_wd/scTrx/proteinFastas

#wget https://ftp.ensembl.org/pub/release-116/fasta/bos_taurus/pep/Bos_taurus.ARS-UCD2.0.pep.all.fa.gz
#gunzip Bos_taurus.ARS-UCD2.0.pep.all.fa.gz
#samtools faidx Bos_taurus.ARS-UCD2.0.pep.all.fa
#mv Bos_taurus.ARS-UCD2.0.pep.all.fa.fai Bos_taurus.ARS-UCD2.0.pep.all.fai


#wget https://ftp.ensembl.org/pub/release-116/fasta/capra_hircus/pep/Capra_hircus.ARS1.pep.all.fa.gz
#gunzip Capra_hircus.ARS1.pep.all.fa.gz
#samtools faidx Capra_hircus.ARS1.pep.all.fa
#mv Capra_hircus.ARS1.pep.all.fa.fai Capra_hircus.ARS1.pep.all.fai

#wget https://ftp.ensembl.org/pub/release-116/fasta/sus_scrofa/pep/Sus_scrofa.Sscrofa11.1.pep.all.fa.gz
#gunzip Sus_scrofa.Sscrofa11.1.pep.all.fa.gz
#samtools faidx Sus_scrofa.Sscrofa11.1.pep.all.fa
#mv Sus_scrofa.Sscrofa11.1.pep.all.fa.fai Sus_scrofa.Sscrofa11.1.pep.all.fai

#wget https://ftp.ensembl.org/pub/release-116/fasta/canis_lupus_familiaris/pep/Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa.gz
#gunzip Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa.gz
#samtools faidx Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa
#mv Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa.fai Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fai

#wget https://ftp.ensembl.org/pub/release-116/fasta/oryctolagus_cuniculus/pep/Oryctolagus_cuniculus.OryCun2.0.pep.all.fa.gz
#gunzip Oryctolagus_cuniculus.OryCun2.0.pep.all.fa.gz
#samtools faidx Oryctolagus_cuniculus.OryCun2.0.pep.all.fa
#mv Oryctolagus_cuniculus.OryCun2.0.pep.all.fa.fai Oryctolagus_cuniculus.OryCun2.0.pep.all.fai

#wget https://ftp.ensembl.org/pub/release-116/fasta/cavia_porcellus/pep/Cavia_porcellus.Cavpor3.0.pep.all.fa.gz
#gunzip Cavia_porcellus.Cavpor3.0.pep.all.fa.gz
#samtools faidx Cavia_porcellus.Cavpor3.0.pep.all.fa
#mv Cavia_porcellus.Cavpor3.0.pep.all.fa.fai Cavia_porcellus.Cavpor3.0.pep.all.fai

#wget https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
#gunzip Homo_sapiens.GRCh38.pep.all.fa.gz
#samtools faidx Homo_sapiens.GRCh38.pep.all.fa
#mv Homo_sapiens.GRCh38.pep.all.fa.fai Homo_sapiens.GRCh38.pep.all.fai

#wget https://ftp.ensembl.org/pub/release-116/fasta/macaca_fascicularis/pep/Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all.fa.gz
#gunzip Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all.fa.gz
#samtools faidx Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all.fa
#mv Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all.fa.fai Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all.fai

wget https://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz
gunzip -f Mus_musculus.GRCm39.pep.all.fa.gz
samtools faidx Mus_musculus.GRCm39.pep.all.fa
mv Mus_musculus.GRCm39.pep.all.fa.fai Mus_musculus.GRCm39.pep.all.fa

#wget https://ftp.ensembl.org/pub/release-116/fasta/rattus_norvegicus_shrutx/pep/Rattus_norvegicus_shrutx.UTH_Rnor_SHR_Utx.pep.all.fa.gz
#gunzip Rattus_norvegicus_shrutx.UTH_Rnor_SHR_Utx.pep.all.fa.gz
#samtools faidx Rattus_norvegicus_shrutx.UTH_Rnor_SHR_Utx.pep.all.fa
#mv Rattus_norvegicus_shrutx.UTH_Rnor_SHR_Utx.pep.all.fa.fai Rattus_norvegicus_shrutx.UTH_Rnor_SHR_Utx.pep.all.fai
