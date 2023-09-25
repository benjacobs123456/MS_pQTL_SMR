#! /bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1:22
#$ -o /data/scratch/hmy117

cd /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR

./useful_files/smr_Linux --bfile ../1kg_reference/hg38_EUR/chr${SGE_TASK_ID}\_hg38 \
--gwas-summary ./useful_files/ms_gwas.ma \
--chr ${SGE_TASK_ID} \
--beqtl-summary /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/combined_pqtl_besd \
--out ./outputs/ms_pqtl_smr_chr${SGE_TASK_ID} \
--thread-num $NSLOTS \
--diff-freq-prop 1
