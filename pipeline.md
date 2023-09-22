Plasma proteins & MS - MR analysis

# Preamble
This repo contains code used to produce the results in xxx.
Code by Shannon Healey and Ben Jacobs.

# Data availability
- pQTL summary statistics [here](https://doi.org/10.7303/syn51364943)
- MS summary statistics by request [here](https://imsgc.net/)
- SMR software [here](https://yanglab.westlake.edu.cn/software/smr/#Overview)

# Code

## Download pQTL GWAS from UKB (https://doi.org/10.7303/syn51364943)
Setup
```unix
module load python/3.10.7
python -m venv env
source ~/env/bin/activate
python3 -m pip3 install --upgrade synapseclient
```
Download
```unix
module load python/3.10.7
source ~/env/bin/activate
synapse login "xxxxxx" "xxxxxx"
cd /data/scratch/hmy117/raw_ukb_pqtl_gwas
synapse get -r syn51365303
```
Unzip
```unix
for file in *.tar
do
  tar xvf $file &
done
```
## Prepare proteomics data

Get proteomics assay data from UKB
```unix
cd /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files
wget https://biobank.ndph.ox.ac.uk/showcase/ukb/auxdata/olink_assay.dat
```
Get probe names
```unix
cd /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/pQTLs

# write list of files
find -type d | grep -v .tar | awk -F/ 'NR>1{print $2}' > probe_names
```

Make esd directory (input for SMR)
````unix
# only run one-off as it will overwrite existing
# mkdir /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/esd_files
```

Make draft epi file (input for SMR)
NB `./useful_files/protein_list` is downloaded the GENCODE v44 gene list downloaded directly from the UCSC browser.
```unix
cd /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR
awk 'NR>1{print $2,$7,9,$4,$6,$3}' ./useful_files/protein_list | awk -Fchr '{print $2}' > all_probes.epi
```

# get filelist with full names
cd /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR
find /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/pQTLs/ -type f | grep .gz > pqtl_filelist

# make the esds
qsub ./scripts/make_esd_files.sh

# navigate to wd
cd /data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/

# edit epi file and make flist file in R
Rscript ./scripts/make_flist_and_epi_files.R

# make besd file
qsub ./scripts/make_besd_file.sh

## Process 1kg data
Convert 1kg reference to hg38
```unix
qsub ./scripts/liftover_1kg_ref.sh
```

Get 1kg allele frequencies
```unix
for i in {1..22};
  do
    ~/plink2 --bfile ../1kg_reference/hg38_EUR/chr$i\_hg38 \
    --freq \
    --out ./outputs/kg_chr$i\_freq
  done
```

## Process MS GWAS data
```unix
Rscript ./scripts/make_ma_file.R
```

## Run SMR
qsub ./scripts/ms_smr.sh

# explore results
qsub ./scripts/explore_results.sh
