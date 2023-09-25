library(tidyverse)
setwd("/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/")

# read uniprot gene list
gene_list = read_tsv("./useful_files/uniprot_proteins_genelist.tsv")

# filter to proteins with an esd file
all_files =  list.files("./esd_files/",full.names=F)
short_probe_name = str_remove_all(str_remove_all(all_files,"pQTL_combined_"),".esd")
all_files = paste0("/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR/esd_files/",all_files)

# make codex
prot_df = data.frame(PathOfEsd = all_files,
                     short_probe_name = short_probe_name,
                     original_probe_name = short_probe_name) %>%
  separate(short_probe_name,sep="_",into=c("protein","other")) %>%
  dplyr::select(protein,original_probe_name,PathOfEsd) %>%
  filter(protein !="" & !is.na(protein))

nrow(prot_df)
nrow(prot_df %>% distinct(protein))

# prots not in gene list 
missing_proteins = prot_df %>%
  dplyr::select(1,2) %>%
  filter(!protein %in% gene_list$hgncSym)
nrow(missing_proteins)

# filter 
prot_df = prot_df %>%
  filter(protein %in% gene_list$hgncSym)
nrow(prot_df)

# filter to autosomes 
gene_list = gene_list %>%
  separate(`#chrom`, sep = "chr",into = c("extra","chr")) %>%
  separate(chr, sep = "_",into = c("#chrom","extra2")) %>%
  dplyr::select(-extra,-extra2)

non_autosomes = gene_list %>%
  filter(is.na(as.numeric(`#chrom`)))

prot_df %>%
  filter(protein %in% non_autosomes$geneName) %>%
  dplyr::select(1,2) %>%
  left_join(non_autosomes %>% 
              dplyr::rename("protein" = hgncSym),by="protein") %>%
  distinct(protein,.keep_all = T) %>%
  dplyr::count(`#chrom`)

# make epi file 
epi_file = gene_list %>%
  mutate(n = 9) %>%
  filter(hgncSym %in% prot_df$protein) %>%
  distinct(hgncSym,.keep_all = T) %>%
  mutate(`#chrom` = as.numeric(`#chrom`)) %>%
  filter(!is.na(`#chrom`))

# combine with epi file 
epi_file = epi_file %>%
  dplyr::rename("protein" = hgncSym) %>%
  filter(protein %in% prot_df$protein) %>%
  left_join(prot_df,by="protein") 

# check unique probes and proteins 
# make files 
flist_file = epi_file %>% 
  dplyr::select(1,14,n,2,7,4,15)
colnames(flist_file) = c("Chr","ProbeID","GeneticDistance","ProbeBp","Gene","Orientation","PathOfEsd")
epi_file = epi_file %>% 
  dplyr::select(1,14,n,2,7,4)
epi_file %>% distinct(protein) %>% nrow()

# save  
write_tsv(epi_file,"filtered_epi_file.epi",col_names=F)
write_tsv(flist_file,"filtered_flist_file.flist")

