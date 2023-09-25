library(tidyverse)
args = commandArgs(trailingOnly = T)
this_file = args[1]
message(this_file)
gwas = read_table2(this_file,col_types=cols(.default="c")) %>%
  filter(nchar(A1)==1 &  nchar(A2)==1 )
multiallelics = gwas %>%
  dplyr::count(SNP) %>% filter(n>1)
message(nrow(multiallelics)," multiallelics")
gwas = gwas %>%
  filter(!SNP %in% multiallelics$SNP)
write_tsv(gwas,this_file)

