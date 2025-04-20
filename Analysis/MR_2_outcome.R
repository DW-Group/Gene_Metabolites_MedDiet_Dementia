library(data.table)
library(dplyr)
library(stringr)
library(readxl)
library(TwoSampleMR)

rm(list=ls())

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

##-----------------------------------------
## Instruments data
##-----------------------------------------

load("data/MR/instruments.RData")

##-----------------------------------------
## Outcome data (local)
##-----------------------------------------

##--------------------------
## Formatting 
##--------------------------

out_wightman <- read_outcome_data(
  filename = "data/MR/outcome/wightman_instruments.txt",
  phenotype_col = "phenotype",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  samplesize_col = "N",
  min_pval = 0,
  chr_col = "chromosome",
  pos_col = "base_pair_location"
)

##--------------------------
## Create outcome list 
##--------------------------

out_local_list <- list()
out_local_list[["wightman_ad"]] <- out_wightman

out_local <- Reduce(full_join, out_local_list)

### Chen_NatGenet_2023

# Instruments selected by the original study
outcome_local_chen_paper_list <- list()
for (pheno in names(instru_chen_paper_list)){
  outcome_local_chen_paper_list[[pheno]] <- out_local %>% filter(SNP %in% instru_chen_paper_list[[pheno]][["SNP"]])
}

# All significant associations for metabolites
outcome_local_chen_supp_met_list <- list()
for (pheno in names(instru_chen_supp_met_list)){
  outcome_local_chen_supp_met_list[[pheno]] <- out_local %>% filter(SNP %in% instru_chen_supp_met_list[[pheno]][["SNP"]])
}

# All significant associations for metabolite ratios
outcome_local_chen_supp_met_ratio_list <- list()
for (pheno in names(instru_chen_supp_met_ratio_list)){
  outcome_local_chen_supp_met_ratio_list[[pheno]] <- out_local %>% filter(SNP %in% instru_chen_supp_met_ratio_list[[pheno]][["SNP"]])
}

##-----------------------------------------
## Outcome data (server)
##-----------------------------------------

outcome_server_list <- c("finn-b-KRA_PSY_DEMENTIA_EXMORE",
                          "finn-b-F5_VASCDEM",
                          "ebi-a-GCST006572")

##--------------------------
## Chen_NatGenet_2023
##--------------------------

### Instruments selected by the original study

chen_paper <- as.data.frame(do.call("rbind", instru_chen_paper_list))
length(unique(chen_paper$SNP))
outcome_chen_paper <- extract_outcome_data_yl(snps=unique(chen_paper$SNP), rsq=0.8, outcomes=outcome_server_list)

outcome_server_chen_paper_list <- list()
for (pheno in names(instru_chen_paper_list)){
  outcome_server_chen_paper_list[[pheno]] <- outcome_chen_paper %>% filter(SNP %in% instru_chen_paper_list[[pheno]][["SNP"]])
}

### All significant associations for metabolites

chen_supp_met <- as.data.frame(do.call("rbind", instru_chen_supp_met_list))
length(unique(chen_supp_met$SNP))
outcome_chen_supp_met <- extract_outcome_data_yl(snps=unique(chen_supp_met$SNP), rsq=0.8, outcomes=outcome_server_list)

outcome_server_chen_supp_met_list <- list()
for (pheno in names(instru_chen_supp_met_list)){
  outcome_server_chen_supp_met_list[[pheno]] <- outcome_chen_supp_met %>% filter(SNP %in% instru_chen_supp_met_list[[pheno]][["SNP"]])
}

### All significant associations for metabolite ratios

chen_supp_met_ratio <- as.data.frame(do.call("rbind", instru_chen_supp_met_ratio_list))
length(unique(chen_supp_met_ratio$SNP))
outcome_chen_supp_met_ratio <- extract_outcome_data_yl(snps=unique(chen_supp_met_ratio$SNP), rsq=0.8, outcomes=outcome_server_list)

outcome_server_chen_supp_met_ratio_list <- list()
for (pheno in names(instru_chen_supp_met_ratio_list)){
  outcome_server_chen_supp_met_ratio_list[[pheno]] <- outcome_chen_supp_met_ratio %>% filter(SNP %in% instru_chen_supp_met_ratio_list[[pheno]][["SNP"]])
}

##-----------------------------------------
## Save data
##-----------------------------------------

gdata::keep(out_local_list, out_local, outcome_server_list,
            outcome_local_chen_paper_list, outcome_local_chen_supp_met_list, outcome_local_chen_supp_met_ratio_list, 
            outcome_chen_paper, outcome_server_chen_paper_list,
            outcome_chen_supp_met, outcome_server_chen_supp_met_list,
            outcome_chen_supp_met_ratio, outcome_server_chen_supp_met_ratio_list,
            sure=T)

save.image("data/MR/outcome.RData")
