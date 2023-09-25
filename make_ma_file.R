# lib
library(tidyverse)

# set wd
setwd("/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/")

# read in MS GWAS
df = read_table("/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/discovery_metav3.0.meta")

# generate SE
df = df %>%
  mutate(freq = NA) %>%
  dplyr::rename("p" = P) %>%
  mutate(b = log(OR)) %>%
  mutate(z = abs(qnorm (1 - (p/2) ) )) %>%
  mutate(se = b/z) %>%
  mutate(n = NA) %>%
  dplyr::select(CHR,BP,SNP,A1,A2,freq,b,se,p,n) %>%
  filter(!is.na(p))

# remove MHC SNPs
df = df %>%
  filter(!(
    CHR == 6 & BP > 25000000 & BP < 35000000
  ))


# hg19 liftover input
hg19_input = df %>%
  dplyr::select(1,2,3) %>%
  mutate(bp1 = BP - 1) %>%
  mutate(CHR = paste0("chr",CHR)) %>%
  dplyr::select(CHR,bp1,BP,SNP)

# write to file
write_tsv(hg19_input,"/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_gwas_hg19.bed",col_names=F)

# run liftover from R
cmd = paste0("/data/Wolfson-UKBB-Dobson/liftover/liftOver ",
"/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_gwas_hg19.bed ",
"/data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz ",
"/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_gwas_hg38.bed ",
"/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_gwas_unlifted.bed ")
system(cmd)

# read in hg38 coordinates
hg38 = read_table("/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_gwas_hg38.bed",col_names=F)
hg38 = hg38 %>%
  dplyr::rename("SNP" = X4)

# merge
df = df %>% left_join(hg38,"SNP")

# generate chr:pos SNP ID using hg38 coordinates
df = df %>%
  mutate(SNP = paste0(X1,":",X3)) %>%
  dplyr::select(-X1,-X2)
df$SNP = str_remove_all(df$SNP,"chr")

# get 1kg freqs
all_freqs = list()
for(i in c(1:22)){
  all_freqs[[i]] = read_table(paste0("./outputs/kg_chr",i,"_freq.afreq"),col_types = "dcccdd")
}
all_freqs = do.call("bind_rows",all_freqs)

# check no duplicates
kg_unique_sites = all_freqs %>%
  filter(ID %in% df$SNP) %>%
  dplyr::count(ID) %>%
  filter(n==1)
all_freqs = all_freqs %>% filter(ID %in% kg_unique_sites$ID)

# match on position initially
df = df %>%
  left_join(all_freqs %>%
              dplyr::rename("SNP" = ID),
              by="SNP")

# filter to SNPs in reference
df = df %>% filter(!is.na(ALT_FREQS))

# find perfect matches & set AF
df = df %>%
  filter(A1 == REF & A2 == ALT | A1 == ALT & A2 == REF ) %>%
  mutate(freq = ifelse(A1 == ALT, ALT_FREQS, 1 - ALT_FREQS))

# just biallelics
df = df %>% filter(nchar(A1)==1 & nchar(A2)==1)

# MAF > 0.01
df = df %>%
  filter(ALT_FREQS >= 0.01 & ALT_FREQS <= 0.99)

# select cols
df = df %>%
dplyr::select(
SNP,A1,A2,freq,b,se,p,n
)

write_tsv(df,"/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_gwas.ma")


# severity
df = read_tsv("/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/imsgc_mssev_discovery.tsv")

df = df %>%
  mutate(freq = AF1) %>%
  dplyr::rename("p" = P,"se"=SE,"b"=BETA) %>%
  mutate(n = 12584) %>%
  dplyr::select(CHR,POS,SNP,A1,A2,freq,b,se,p,n) %>%
  filter(!is.na(p))


# split SNP ID - fix?
hg19_input = df %>%
  dplyr::select(CHR,POS) %>%
  mutate(SNP = paste0(CHR,":",POS)) %>%
  mutate(bp1 = POS - 1) %>%
  mutate(CHR = paste0("chr",CHR)) %>%
  dplyr::select(CHR,bp1,POS,SNP)

# write to file
write_tsv(hg19_input,"/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_sev_gwas_hg19.bed",col_names=F)

# run liftover from R
cmd = paste0("/data/Wolfson-UKBB-Dobson/liftover/liftOver ",
"/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_sev_gwas_hg19.bed ",
"/data/Wolfson-UKBB-Dobson/liftover/hg19ToHg38.over.chain.gz ",
"/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_sev_gwas_hg38.bed ",
"/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_sev_gwas_unlifted.bed ")
system(cmd)

# read in hg38 coordinates
hg38 = read_table2("/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_sev_gwas_hg38.bed",col_names=F,col_types=cols(.default="c"))
gwas = df %>%
  mutate(X4 = paste0(CHR,":",POS)) %>%
  left_join(hg38,"X4")



# just biallelics
gwas = gwas %>%
  dplyr::select(-SNP) %>%
  dplyr::rename("SNP"=X4) %>%
  filter(nchar(A1)==1 & nchar(A2)==1)


# select cols
gwas = gwas %>%
  dplyr::select(SNP,A1,A2,freq,b,se,p,n)

write_tsv(gwas,"/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/useful_files/ms_sev_gwas.ma")
