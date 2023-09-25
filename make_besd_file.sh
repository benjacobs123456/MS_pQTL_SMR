#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y

# list all files 
cd /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR

./useful_files/smr_Linux --thread-num $NSLOTS --eqtl-flist filtered_flist_file.flist --make-besd --out combined_pqtl_besd 
