#! /bin/bash
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1:22
#$ -o /data/scratch/hmy117

cd /data/Wolfson-UKBB-Dobson/SMR

# liftover to hg38
awk '{print "chr"$1,$4-1,$4,$2}' ../1kg_reference/filtered_chr${SGE_TASK_ID}\.bim > /data/scratch/hmy117/chr${SGE_TASK_ID}\_hg19.bed
/data/Wolfson-UKBB-Dobson/liftover/liftOver \
/data/scratch/hmy117/chr${SGE_TASK_ID}\_hg19.bed \
/data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz \
/data/scratch/hmy117/chr${SGE_TASK_ID}\_hg38.bed \
/data/scratch/hmy117/chr${SGE_TASK_ID}\_unlifted.bed

awk '{print $4,$3}' /data/scratch/hmy117/chr${SGE_TASK_ID}\_hg38.bed > /data/scratch/hmy117/chr${SGE_TASK_ID}\_hg38_positions.bed

# update positions
~/plink2 --bfile ../1kg_reference/filtered_chr${SGE_TASK_ID} \
--update-map /data/scratch/hmy117/chr${SGE_TASK_ID}\_hg38_positions.bed \
--extract /data/scratch/hmy117/chr${SGE_TASK_ID}\_hg38_positions.bed \
--make-bed \
--out /data/scratch/hmy117/chr${SGE_TASK_ID}\_hg38_positions_intermediate

# update SNP IDs
~/plink2 --bfile /data/scratch/hmy117/chr${SGE_TASK_ID}\_hg38_positions_intermediate \
--set-all-var-ids @:# \
--out /data/Wolfson-UKBB-Dobson/1kg_reference/hg38_EUR/chr${SGE_TASK_ID}\_hg38 \
--make-bed
