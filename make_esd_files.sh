#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -o /data/scratch/hmy117
#$ -t 1:2940

# list all files
cd /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/pQTLs

# get the number of rows
num_rows=$(awk 'END{print NR}' probe_names)

# for the ith job
## find the file
## do the formatting
## save to big file
echo "doing probe ${SGE_TASK_ID} of $num_rows"
# get probe name
probe_name=$(awk -v row=$((${SGE_TASK_ID})) 'NR==row{print}' probe_names)

# clean up
rm /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/esd_files/pQTL_combined_$probe_name\.esd

# get folder name
folder_name=$(echo "/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/pQTLs/"$probe_name)

# split probe name to short probe name
short_probe_name=$(echo $probe_name | awk -F_ '{print $1}')

# navigate to folder
cd $folder_name

chr=$(grep -E "$short_probe_name " /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/all_probes.epi | awk '{print $1}' | awk -F_ 'NR==1{print $1}')

# do formatting
file_to_use=$(ls | grep chr$chr)


# do formatting
# filter to MAF > 0.01 & INFO > 0.9
zcat $file_to_use |
awk 'NR==1{print "Chr","SNP","Bp","A1","A2","Freq","Beta","se","p"};NR>1{if($6 > 0.01 && $6 < 0.99 && $7 > 0.9) print $1,$1":"$2,$2,$5,$4,$6,$10,$11,10^-$13}' > \
/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/esd_files/pQTL_combined_$probe_name\.esd

module load R/3.6.1
Rscript /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/scripts/filter_out_multiallelics.R /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/esd_files/pQTL_combined_$probe_name\.esd
