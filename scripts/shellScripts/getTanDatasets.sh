#!/bin/bash
#SBATCH --job-name=getTanData
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

cd /work/hcn4/260630_vertCons_wd/scTrx/rObjs

# rabbit
#wget https://api.figshare.com/v2/file/download/64585266
#mv 64585266 rabbit.rds

# mouse
#wget https://api.figshare.com/v2/file/download/65232021
#mv 65232021 mouse.rds

# rat 
#wget https://api.figshare.com/v2/file/download/65232027
#mv 65232027 rat.rds

# macaque
#wget https://api.figshare.com/v2/file/download/65242284
#mv 65242284 macaque.rds

# pig
#wget https://api.figshare.com/v2/file/download/64577394
#mv 64577394 pig.rds

# guinea pig
#wget https://api.figshare.com/v2/file/download/64577391
#mv 64577391 guineaPig.rds

# goat 
#wget https://api.figshare.com/v2/file/download/64577388
#mv 64577388 goat.rds

# dog
#wget https://api.figshare.com/v2/file/download/64577382
#mv 64577382 dog.rds

# cow
#wget https://api.figshare.com/v2/file/download/64577376
#mv 64577376 cow.rds

# human
wget https://api.figshare.com/v2/file/download/68171923
mv 68171923 human_list.rds
