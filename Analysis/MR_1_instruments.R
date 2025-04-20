library(data.table)
library(dplyr)
library(stringr)
library(readxl)
library(MendelianRandomization)
library(TwoSampleMR)
library(MRPRESSO)

rm(list=ls())

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

##-----------------------------------------
## Instruments data
##-----------------------------------------

##--------------------------
## Chen NatGenet 2023
##--------------------------

### Instruments selected by the original study

chen <- read_excel("data/MR/instruments/Chen_NatGenet_2023/Chen_NatGenet_2023_instruments_selected_by_original_study.xlsx") %>%
  subset(select=c("Genetic variants","Metabolites or metabolite ratios",
                  "Beta","Se","Effect Allele","Non-Effect Allele","Effect Allele Frequency","P Value","Effector gene(s)","N",
                  "CHR","POS_hg38","POS_hg37"))
colnames(chen)
colnames(chen) <- c("SNP","exposure",
                    "beta.exposure","se.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","pval.exposure","gene.exposure","samplesize.exposure",
                    "chr.exposure","position","position_hg37")

unique(chen$exposure)
sum(grepl("^rs",chen$SNP))
chen <- chen %>% 
  mutate(SNP = case_when(grepl("^rs",SNP) ~ SNP,
                         !grepl("^rs",SNP) ~ paste0(chr.exposure,":",position_hg37)))
chen_pheno_list <- unique(chen$exposure)
chen_n <- as.data.frame(table(chen$exposure)); colnames(chen_n) <- c("pheno","n_instruments")

instru_chen_paper_list <- list()
for (pheno in chen_pheno_list){
  instru_chen_paper_list[[pheno]] <- chen %>% filter(exposure == pheno) %>% mutate(id.exposure = pheno)
}

### All significant associations for metabolites

chen_supp_met <- read_excel("data/MR/instruments/Chen_NatGenet_2023/Chen_NatGenet_2023_all_identified_associations_metabolites.xlsx") %>%
  subset(select=c("SNP","Metabolite",
                  "BETA","SE","Effect Allele","Non-Effect Allele","FREQ","P*","Effector genes","N",
                  "CHR","Position(hg38)","liftover_output_b37"))
colnames(chen_supp_met)
colnames(chen_supp_met) <- c("SNP","exposure",
                            "beta.exposure","se.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","pval.exposure","gene.exposure","samplesize.exposure",
                            "chr.exposure","position","chr_position_hg37")

chen_supp_met <- chen_supp_met %>% filter(grepl("^rs",SNP))

chen_supp_met_pheno_list <- unique(chen_supp_met$exposure)
chen_supp_met_n <- as.data.frame(table(chen_supp_met$exposure)); colnames(chen_supp_met_n) <- c("pheno","n_instruments")

instru_chen_supp_met_list <- list()
for (pheno in chen_supp_met_pheno_list){
  instru_chen_supp_met_list[[pheno]] <- chen_supp_met %>% filter(exposure == pheno) %>% mutate(id.exposure = pheno)
}

### All significant associations for metabolite ratios

chen_supp_met_ratio <- read_excel("data/MR/instruments/Chen_NatGenet_2023/Chen_NatGenet_2023_all_identified_associations_metabolite_ratio.xlsx") %>%
  subset(select=c("SNP","Metabolite ratios",
                  "BETA","SE","Effect Allele","Non-Effect Allele","FREQ","P*","Effector genes","N",
                  "CHR","POS","liftover_output_b37"))
colnames(chen_supp_met_ratio)
colnames(chen_supp_met_ratio) <- c("SNP","exposure",
                                  "beta.exposure","se.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","pval.exposure","gene.exposure","samplesize.exposure",
                                  "chr.exposure","position","chr_position_hg37")

chen_supp_met_ratio <- chen_supp_met_ratio %>% filter(grepl("^rs",SNP))

chen_supp_met_ratio_pheno_list <- unique(chen_supp_met_ratio$exposure)
chen_supp_met_ratio_n <- as.data.frame(table(chen_supp_met_ratio$exposure)); colnames(chen_supp_met_ratio_n) <- c("pheno","n_instruments")

instru_chen_supp_met_ratio_list <- list()
for (pheno in chen_supp_met_ratio_pheno_list){
  instru_chen_supp_met_ratio_list[[pheno]] <- chen_supp_met_ratio %>% filter(exposure == pheno) %>% mutate(id.exposure = pheno)
}

##-----------------------------------------
## Save data
##-----------------------------------------

gdata::keep(chen, chen_pheno_list, chen_n, instru_chen_paper_list,
            chen_supp_met, chen_supp_met_pheno_list, chen_supp_met_n, instru_chen_supp_met_list,
            chen_supp_met_ratio, chen_supp_met_ratio_pheno_list, chen_supp_met_ratio_n, instru_chen_supp_met_ratio_list,
            sure=T)

save.image("data/MR/instruments.RData")
