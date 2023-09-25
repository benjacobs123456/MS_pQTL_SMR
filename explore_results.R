# setup
library(tidyverse)
library(ggrepel)
library(coloc)

setwd("/data/Wolfson-UKBB-Dobson/MS_pQTL_SMR")

# read in res
smr_res = list()
for(i in c(1:22)){
  smr_res[[i]] = read_table(
    paste0("./outputs/ms_pqtl_smr_chr",i,".smr"),
    col_types = "cdcdcddccdddddddddddd")
}
smr_res = do.call("bind_rows",smr_res)

nrow(smr_res %>% distinct(Gene))

# look at replication hits 
xin_hits = c("FLRT3", "TAPBL", "FCRL3", "IDUA", "MAPK3")
smr_res %>% 
  filter(Gene %in% xin_hits)

# exclude MHC
smr_res = smr_res %>%
  filter(!(topSNP_chr == 6 & topSNP_bp > 25000000 & topSNP_bp < 35000000))

# recode NAs as P value = 1 (in these cases the b_GWAS is 0)
smr_res = smr_res %>%
  mutate(p_SMR = ifelse(is.na(p_SMR) & b_GWAS == 0 & se_GWAS == 0,
                        1,
                        p_SMR))

# bonf 
bonf = 0.01 / nrow(smr_res)

# lambda
library(qqman)
lambda = median(qchisq(smr_res$p_SMR,df=1,lower.tail = F)) / qchisq(0.5,df=1,lower.tail = F)
png("./outputs/plots/qq.png",res=900,units="in",width=6,height=4)
qqman::qq(smr_res$p_SMR)
dev.off()

# MAF vs beta 
png("./outputs/plots/maf_vs_beta.png",res=900,units="in",width=4,height=4)
ggplot(smr_res %>%
         mutate(maf = ifelse(Freq<0.5,Freq,1-Freq)),
       aes(maf,abs(b_SMR)))+
  geom_point()+
  theme_bw()
dev.off()

# manhattan
## get chrom midpoints 
chr_coords = smr_res %>%
  group_by(ProbeChr) %>%
  summarise(min_bp = min(Probe_bp),
            max_bp = max(Probe_bp),
            midpoint = mean(Probe_bp)) %>%
  ungroup %>%
  mutate(cumbp = cumsum(max_bp)) %>%
  mutate(ProbeChr = ProbeChr + 1)

smr_res = smr_res %>%
  left_join(chr_coords,by="ProbeChr") %>%
  mutate(cum_pos = ifelse(is.na(cumbp),Probe_bp,Probe_bp + cumbp))

midpoints = chr_coords %>%
  mutate(midpoint_cum = cumbp - (max_bp - min_bp)/2 )

library(ggrepel)
smr_res = smr_res %>%
  separate(probeID,sep = "_",into = c("protein","extra"))
bonf = 0.01/length(smr_res$p_SMR)
smr_res = smr_res %>%
  mutate(color = case_when(
    p_SMR < bonf ~ "gwas_sig",
    ProbeChr %% 2 == 0 ~ "even",
    ProbeChr %% 2 != 0 ~ "odd",
  ))

pal = c("orange1","lightblue1","pink")
names(pal) = c("gwas_sig","even","odd")

p1=ggplot(smr_res %>% filter(p_HEIDI > 0.05),
       aes(cum_pos,-log10(p_SMR),col = color))+
  geom_point(size = 0.5)+
  scale_color_manual(values = pal) + 
  geom_hline(yintercept=-log10(bonf),col="red",linetype="dashed",alpha=0.1)+
  theme_minimal()+
  scale_x_continuous(breaks = midpoints$midpoint_cum,
                     labels = c(1:22))+
  theme(panel.grid = element_blank(),legend.position = "none")+
  labs(x="Chromosome",y="-log10(P)")+
  geom_label_repel(data = smr_res %>% filter(p_SMR < bonf & p_HEIDI > 0.05),
                   mapping = aes(label = protein),
                   max.overlaps = 10,
                   force_pull = 0,
                   nudge_y = 1,
                   size = 2,
                   min.segment.length = 0,
                   color="black",
                   segment.linetype = 2)+
  ggtitle("A")


# plot orientation of effects 
sig_hits = smr_res %>%
  filter(p_SMR < bonf & p_HEIDI> 0.05)

sig_hits = sig_hits %>% arrange(desc(b_SMR))
sig_hits$protein = factor(sig_hits$protein,levels = sig_hits$protein,ordered=T)
p2=ggplot(sig_hits,
       aes(b_SMR,protein))+
  geom_point()+
  geom_vline(xintercept = 0,linetype = "dashed")+
  geom_errorbarh(mapping = aes(y=protein,xmin = b_SMR - 1.96*se_SMR,xmax = b_SMR + 1.96*se_SMR),height=0.01)+
  scale_color_manual(values = c("blue","red"))+
  labs(x = "SMR estimate \n(Effect on MS risk of increased protein level)",y = "Protein")+
  theme_bw()+
  scale_x_continuous(limits = c(-3,3))+
  annotate("text",label = "Decreases MS risk",x=-1.5,y=4.5,color="blue")+
  annotate("text",label = "Increases MS risk",x=1.5,y=4.5,color="red")+
  ggtitle("B")

png("./outputs/plots/effect_directions.png",res=900,unit="in",height=4,width=6)
p2
dev.off()

png("./outputs/plots/manhattan.png",res=900,unit="in",height=4,width=8)
p1
dev.off()

# adjust P vals
smr_res$fdr_SMR = p.adjust(smr_res$p_SMR,method="fdr")
sig_hits_for_table = smr_res %>%
  filter(fdr_SMR < 0.05 & p_HEIDI > 0.05) %>%
  dplyr::select(protein,ProbeChr,Probe_bp,topSNP,A1,A2,b_GWAS,p_GWAS,p_eQTL,b_SMR,p_SMR,p_HEIDI) %>%
  mutate(MS_risk_allele = ifelse(b_GWAS > 0, A1,A2)) %>%
  dplyr::rename("CHR" = ProbeChr,
                "BP (hg38)" = Probe_bp,
                "SNP" = topSNP)
write_csv(sig_hits_for_table,"./outputs/smr_sig_hits.csv")


# MUTLI-SMR
# read in res
multi_smr_res = list()
for(i in c(1:22)){
  multi_smr_res[[i]] = read_table(
    paste0("./outputs/ms_pqtl_smr_multi_chr",i,".msmr"),
    col_types = "cdcdcddccdddddddddddd")
}
multi_smr_res = do.call("bind_rows",multi_smr_res)

# exclude MHC
multi_smr_res = multi_smr_res %>%
  filter(!(topSNP_chr == 6 & topSNP_bp > 25000000 & topSNP_bp < 35000000))

# exclue NAs 
multi_smr_res = multi_smr_res %>%
  filter(!is.na(p_SMR))


# get bonf
multi_smr_res = multi_smr_res %>%
  mutate(bonf_SMR = p.adjust(p_SMR_multi ,method="bonf"))

# get fdr
multi_smr_res = multi_smr_res %>%
  mutate(fdr_SMR = p.adjust(p_SMR_multi ,method="fdr"))

# filter
multi_smr_res_sig_hits = multi_smr_res %>%
  filter(bonf_SMR < 0.01 & p_HEIDI > 0.05)
write_csv(multi_smr_res_sig_hits,"./outputs/smr_sig_hits_multi.csv")


# locus plot

# read in gwas
ms_gwas  = read_table("useful_files/ms_gwas.ma",col_types = "cccddddd") %>% 
  mutate(trait = "MS") %>% 
  dplyr::select(SNP,p,trait,A1,A2)

# define coloc function
do_coloc = function(locus_file,locus,input_gwas=ms_gwas_for_coloc){
  # read in pqtls
  pqtl = read_table(locus_file,col_types = "dcdccdddd")  %>% 
    mutate(trait = locus) 
  
  # read in sig hit data 
  topsnp = read_csv("./outputs/smr_sig_hits.csv",col_types = cols(.default = "c")) %>%
    filter(protein == locus)
  
  
  # filter to everything within 1MB of top snp
  pqtl = pqtl %>% filter(Chr == as.numeric(topsnp$CHR) & Bp > (as.numeric(topsnp$`BP (hg38)`) - 5e5) & 
                           Bp < (as.numeric(topsnp$`BP (hg38)`) + 5e5))
  
  
  # filter the MS gwas
  ms_gwas = input_gwas %>% filter(snp %in% pqtl$SNP)
  
  # filter pQTL gwas
  pqtl = pqtl %>% filter(SNP %in% ms_gwas$snp)
  
  # check alleles compatible 
  
  snps_to_keep = pqtl %>%
    dplyr::rename("snp" = SNP) %>%
    left_join(ms_gwas,by="snp") %>%
    filter(
      (A1.x == A1.y & A2.x == A2.y) | 
        (A1.x == A2.y & A2.x == A1.y)
    )
  
  # find SNPs that need flipping of effects 
  snps_to_flip = pqtl %>%
    dplyr::rename("snp" = SNP) %>%
    left_join(ms_gwas,by="snp") %>%
    filter(
      (A1.x == A2.y & A2.x == A1.y)
    )
  
  # flip betas for MS GWAS to align effect alleles
  ms_gwas = ms_gwas %>%
    mutate(beta = ifelse(snp %in% snps_to_flip$snp,
                         beta * -1,
                         beta))
  # filter pQTLs 
  pqtl = pqtl %>%
    filter(SNP %in% snps_to_keep$snp) %>%
    dplyr::rename("snp" = SNP) %>%
    mutate(maf = ifelse(Freq < 0.5,Freq,1-Freq)) %>%
    mutate(n = 33533) %>% 
    dplyr::select(snp,Bp,Beta,se,p,n,maf) %>%
    mutate(varbeta = se^2,
           type = "quant") %>%
    rename("beta" = "Beta") %>% 
    rename("position" = "Bp")
  
  # final inputs need to be in List format:
  ms_gwas = ms_gwas %>%
    mutate(position_temp = snp) %>%
    separate(position_temp,sep=":",into = c("chr","position"))
  
  # ms gwas 
  n_case = 14802
  n_cont = 26703 
  prop_case = n_case / (n_case + n_cont)
  
  ms_gwas_coloc_input =  list(
    beta = ms_gwas$beta, 
    varbeta = ms_gwas$varbeta,
    snp = ms_gwas$snp,
    type=ms_gwas$type[1],
    MAF = ms_gwas$maf,
    N = n_case+n_cont,
    s = prop_case,
    position = as.numeric(ms_gwas$position)
  )
  pqtl_gwas_coloc_input =  list(
    beta = pqtl$beta, 
    varbeta = pqtl$varbeta,
    snp = pqtl$snp,
    type="quant",
    pvalues=pqtl$p,
    MAF = pqtl$maf,
    N = pqtl$n[1],
    position = pqtl$position
  )
  
  # use check function to ensure ready to input
  check_dataset(ms_gwas_coloc_input,warn.minp = 0.05)
  check_dataset(pqtl_gwas_coloc_input,warn.minp = 0.05)
  
  # run coloc 
  coloc = coloc.abf(ms_gwas_coloc_input,pqtl_gwas_coloc_input)
  coloc_res = tibble(coloc$results)
  
  p=ggplot(coloc_res,aes(position,`SNP.PP.H4`))+
    geom_point()+
    geom_text_repel(data = coloc_res %>% 
                      slice_max(`SNP.PP.H4`,
                                with_ties = F),
                    mapping = aes(label = snp),
                    show.legend = F,min.segment.length = 0)+
    theme_bw()+
    labs(x=paste0("Genomic position (hg38)\non chromosome ",ms_gwas$chr[1]))+
    theme(legend.position = "none")+
    ggtitle("Colocalisation")
  return(p)
}  

# define locus plot function
make_locus_plot = function(locus_file,locus){
# read in pqtls
pqtl = read_table(locus_file,col_types = "dcdccdddd")  %>% 
  mutate(trait = locus) %>% 
  dplyr::select(SNP,Bp,Chr,p,trait,A1,A2)

topsnp = read_csv("./outputs/smr_sig_hits.csv",col_types = cols(.default = "c")) %>%
  filter(protein == locus)


# filter to everything within 1MB of top snp
message("Filtering pQTLs")
pqtl = pqtl %>% filter(Chr == as.numeric(topsnp$CHR) & Bp > (as.numeric(topsnp$`BP (hg38)`) - 5e5) & 
                         Bp < (as.numeric(topsnp$`BP (hg38)`) + 5e5))

# filter the MS gwas
message("Filtering MS GWAS")

ms_gwas = ms_gwas %>% filter(SNP %in% pqtl$SNP)

# filter pQTL gwas
pqtl = pqtl %>% filter(SNP %in% ms_gwas$SNP)

# check alleles compatible 

snps_to_keep = pqtl %>%
  left_join(ms_gwas,by="SNP") %>%
  filter(
    (A1.x == A1.y & A2.x == A2.y) | 
      (A1.x == A2.y & A2.x == A1.y)
  )

# filter to compatible SNPs
combo = ms_gwas %>% 
  bind_rows(pqtl) %>%
  filter(SNP %in% snps_to_keep$SNP)

# separate snp name
combo = combo %>% separate(SNP,sep = ":",into=c("CHR","BP"))
combo = combo %>% mutate(SNP = paste0(CHR,":",BP))

topsnps_both_traits = combo %>%
  group_by(trait) %>%
  slice_min(p,with_ties = F) %>%
  ungroup() %>%
  distinct(SNP)
snps_to_keep = topsnps_both_traits$SNP


# split plots 
p1=ggplot(combo %>% filter(trait=="MS"),
          aes(as.numeric(BP),-log10(p)))+
  geom_point(color="red")+
  geom_text_repel(data = combo %>% 
                    filter(trait=="MS") %>%
                    filter(SNP %in% snps_to_keep),
                  mapping = aes(label = SNP),
                  show.legend = F,min.segment.length = 0)+
  theme_bw()+
  labs(x=paste0("Genomic position (hg38)\non chromosome ",combo$CHR[1]))+
  theme(legend.position = "none")+
    ggtitle("MS")
p2=ggplot(combo %>% filter(trait!="MS"),
          aes(as.numeric(BP),-log10(p)))+
  geom_point(color="blue")+
  geom_text_repel(data = combo %>% 
                    filter(trait==locus) %>%
                    filter(SNP %in% snps_to_keep),
                  mapping = aes(label = SNP),
                  show.legend = F,min.segment.length = 0)+
  theme_bw()+
  labs(x=paste0("Genomic position (hg38)\non chromosome ",combo$CHR[1]))+
  theme(legend.position = "none")+
  ggtitle(locus)

message("Made locus plots, now doing coloc")
# do coloc 
# read in the gwas data & format for coloc
ms_gwas_for_coloc  = read_table("useful_files/ms_gwas.ma",col_types = "cccddddd") %>%
  mutate(trait = "MS") %>% 
  mutate(maf = ifelse(freq < 0.5,freq,1-freq)) %>%
  dplyr::select(SNP,b,se,p,n,maf,A1,A2) %>% 
  mutate(varbeta = se^2) %>%
  mutate(type = "cc") %>%
  rename("snp" = "SNP") %>% 
  rename("beta" = "b") %>%
  filter(!is.na(varbeta)) %>%
  filter(varbeta != 0)

coloc_plot = do_coloc(locus_file = locus_file, locus = locus, input_gwas = ms_gwas_for_coloc)

png(paste0("./outputs/plots/locus_plot_split_",locus,".png"),
    res=900,units="in",height=6,width=4)
print(cowplot::plot_grid(p1, p2,coloc_plot, align = "v", nrow = 3, rel_heights = rep(1/3,3)))

dev.off()
}

make_locus_plot(
  locus_file = "./esd_files/pQTL_combined_TNFRSF1A_P19438_OID21155_v1_Neurology.esd",
  locus = "TNFRSF1A"
)

make_locus_plot(
  locus_file = "./esd_files/pQTL_combined_CD40_P25942_OID20724_v1_Inflammation.esd",
  locus = "CD40"
)

make_locus_plot(
  locus_file = "./esd_files/pQTL_combined_FCRL3_Q96P31_OID20443_v1_Inflammation.esd",
  locus = "FCRL3"
)

make_locus_plot(
  locus_file = "./esd_files/pQTL_combined_CD58_P19256_OID20716_v1_Inflammation.esd",
  locus = "CD58"
)




