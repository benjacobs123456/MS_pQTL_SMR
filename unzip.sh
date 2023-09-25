#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -o /data/scratch/hmy117
#$ -t 1:2940

cd /data/scratch/hmy117/raw_ukb_pqtl_gwas/
file=$(ls | grep -E ".tar" | awk -v row=${SGE_TASK_ID} 'NR==row{print}')

cd /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/pQTLs
tar -xvf /data/scratch/hmy117/raw_ukb_pqtl_gwas/$file

wait 
echo "done"
