library(data.table)
library(dplyr)
library(stringr)
library(readxl)
library(TwoSampleMR)

rm(list=ls())

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

##-----------------------------------------
## Instruments and outcome data
##-----------------------------------------

load("data/MR/instruments.RData")
load("data/MR/outcome.RData")

##-----------------------------------------
## Check numbers
##-----------------------------------------

### chen_paper

check_local_chen_paper_list <- c(); i <- 1
for(met in names(outcome_local_chen_paper_list)){
  tmp <- outcome_local_chen_paper_list[[met]] %>% 
    filter(outcome=="wightman_ad") 
  if (dim(tmp)[1]!=0){tmp <- tmp %>% mutate(exposure=met)}
  check_local_chen_paper_list[[i]] <- tmp; i <- i + 1
}
check_local_chen_paper <- as.data.frame(do.call("rbind",check_local_chen_paper_list)) %>%
  subset(select=c(outcome,exposure,SNP))

check_server_chen_paper_list <- c(); i <- 1
for(met in names(outcome_server_chen_paper_list)){
  tmp <- outcome_server_chen_paper_list[[met]] %>% 
    filter(outcome %in% c("Any dementia (more controls excluded) || id:finn-b-KRA_PSY_DEMENTIA_EXMORE",
                          "Cognitive performance || id:ebi-a-GCST006572",
                          "Vascular dementia (F5_VASCDEM) || id:finn-b-F5_VASCDEM")) 
  if (dim(tmp)[1]!=0){tmp <- tmp %>% mutate(exposure=met)}
  check_server_chen_paper_list[[i]] <- tmp; i <- i + 1
}
check_server_chen_paper <- as.data.frame(do.call("rbind",check_server_chen_paper_list)) %>%
  subset(select=c(outcome,exposure,SNP))

check_chen_paper <- as.data.frame(rbind(check_local_chen_paper,check_server_chen_paper))

chen_supp_met <- read_excel("data/MR/instruments/Chen_NatGenet_2023/Chen_NatGenet_2023_all_identified_associations_metabolites.xlsx") %>%
  subset(select=c("SNP","Metabolite",
                  "BETA","SE","Effect Allele","Non-Effect Allele","FREQ","P*","Effector genes","N",
                  "CHR","Position(hg38)","liftover_output_b37"))
chen_supp_met_ratio <- read_excel("data/MR/instruments/Chen_NatGenet_2023/Chen_NatGenet_2023_all_identified_associations_metabolite_ratio.xlsx") %>%
  subset(select=c("SNP","Metabolite ratios",
                  "BETA","SE","Effect Allele","Non-Effect Allele","FREQ","P*","Effector genes","N",
                  "CHR","POS","liftover_output_b37"))
check_chen_paper_met <- check_chen_paper %>% filter(exposure %in% chen_supp_met$Metabolite)
check_chen_paper_met_ratio <- check_chen_paper %>% filter(exposure %in% chen_supp_met_ratio$`Metabolite ratios`)

### chen_supp_met

check_local_chen_supp_met_list <- c(); i <- 1
for(met in names(outcome_local_chen_supp_met_list)){
  tmp <- outcome_local_chen_supp_met_list[[met]] %>% 
    filter(outcome=="wightman_ad") 
  if (dim(tmp)[1]!=0){tmp <- tmp %>% mutate(exposure=met)}
  check_local_chen_supp_met_list[[i]] <- tmp; i <- i + 1
}
check_local_chen_supp_met <- as.data.frame(do.call("rbind",check_local_chen_supp_met_list)) %>%
  subset(select=c(outcome,exposure,SNP))

check_server_chen_supp_met_list <- c(); i <- 1
for(met in names(outcome_server_chen_supp_met_list)){
  tmp <- outcome_server_chen_supp_met_list[[met]] %>% 
    filter(outcome %in% c("Any dementia (more controls excluded) || id:finn-b-KRA_PSY_DEMENTIA_EXMORE",
                          "Cognitive performance || id:ebi-a-GCST006572",
                          "Vascular dementia (F5_VASCDEM) || id:finn-b-F5_VASCDEM")) 
  if (dim(tmp)[1]!=0){tmp <- tmp %>% mutate(exposure=met)}
  check_server_chen_supp_met_list[[i]] <- tmp; i <- i + 1
}
check_server_chen_supp_met <- as.data.frame(do.call("rbind",check_server_chen_supp_met_list)) %>%
  subset(select=c(outcome,exposure,SNP))

check_chen_supp_met <- as.data.frame(rbind(check_local_chen_supp_met,check_server_chen_supp_met)) %>%
  filter(!(exposure %in% names(instru_chen_paper_list)))

### chen_supp_met_ratio

check_local_chen_supp_met_ratio_list <- c(); i <- 1
for(met in names(outcome_local_chen_supp_met_ratio_list)){
  tmp <- outcome_local_chen_supp_met_ratio_list[[met]] %>% 
    filter(outcome=="wightman_ad")
  if (dim(tmp)[1]!=0){tmp <- tmp %>% mutate(exposure=met)}
  check_local_chen_supp_met_ratio_list[[i]] <- tmp; i <- i + 1
}
check_local_chen_supp_met_ratio <- as.data.frame(do.call("rbind",check_local_chen_supp_met_ratio_list)) %>%
  subset(select=c(outcome,exposure,SNP))

check_server_chen_supp_met_ratio_list <- c(); i <- 1
for(met in names(outcome_server_chen_supp_met_ratio_list)){
  tmp <- outcome_server_chen_supp_met_ratio_list[[met]] %>% 
    filter(outcome %in% c("Any dementia (more controls excluded) || id:finn-b-KRA_PSY_DEMENTIA_EXMORE",
                          "Cognitive performance || id:ebi-a-GCST006572",
                          "Vascular dementia (F5_VASCDEM) || id:finn-b-F5_VASCDEM")) 
  if (dim(tmp)[1]!=0){tmp <- tmp %>% mutate(exposure=met)}
  check_server_chen_supp_met_ratio_list[[i]] <- tmp; i <- i + 1
}
check_server_chen_supp_met_ratio <- as.data.frame(do.call("rbind",check_server_chen_supp_met_ratio_list)) %>%
  subset(select=c(outcome,exposure,SNP))

check_chen_supp_met_ratio <- as.data.frame(rbind(check_local_chen_supp_met_ratio,check_server_chen_supp_met_ratio)) %>%
  filter(!(exposure %in% names(instru_chen_paper_list)))

##-----------------------------------------
## Harmonize data
##-----------------------------------------

##--------------------------
## Chen_NatGenet_2023
##--------------------------

### Instruments selected by the original study

harmonized_outcome_server_chen_paper_list <- list()
for (pheno in names(outcome_server_chen_paper_list)){
  harmonized_outcome_server_chen_paper_list[[pheno]] <- harmonise_data(exposure_dat=instru_chen_paper_list[[pheno]], outcome_dat=outcome_server_chen_paper_list[[pheno]])
}

harmonized_outcome_local_chen_paper_list <- list()
for (pheno in names(outcome_local_chen_paper_list)){
  harmonized_outcome_local_chen_paper_list[[pheno]] <- harmonise_data(exposure_dat=instru_chen_paper_list[[pheno]], outcome_dat=outcome_local_chen_paper_list[[pheno]])
}

### All significant associations for metabolites

harmonized_outcome_server_chen_supp_met_list <- list()
for (pheno in names(outcome_server_chen_supp_met_list)){
  harmonized_outcome_server_chen_supp_met_list[[pheno]] <- harmonise_data(exposure_dat=instru_chen_supp_met_list[[pheno]], outcome_dat=outcome_server_chen_supp_met_list[[pheno]])
}

harmonized_outcome_local_chen_supp_met_list <- list()
for (pheno in names(outcome_local_chen_supp_met_list)){
  harmonized_outcome_local_chen_supp_met_list[[pheno]] <- harmonise_data(exposure_dat=instru_chen_supp_met_list[[pheno]], outcome_dat=outcome_local_chen_supp_met_list[[pheno]])
}

### All significant associations for metabolite ratios

harmonized_outcome_server_chen_supp_met_ratio_list <- list()
for (pheno in names(outcome_server_chen_supp_met_ratio_list)){
  harmonized_outcome_server_chen_supp_met_ratio_list[[pheno]] <- harmonise_data(exposure_dat=instru_chen_supp_met_ratio_list[[pheno]], outcome_dat=outcome_server_chen_supp_met_ratio_list[[pheno]])
}

harmonized_outcome_local_chen_supp_met_ratio_list <- list()
for (pheno in names(outcome_local_chen_supp_met_ratio_list)){
  harmonized_outcome_local_chen_supp_met_ratio_list[[pheno]] <- harmonise_data(exposure_dat=instru_chen_supp_met_ratio_list[[pheno]], outcome_dat=outcome_local_chen_supp_met_ratio_list[[pheno]])
}

##-----------------------------------------
## Save data
##-----------------------------------------

gdata::keep(harmonized_outcome_local_chen_paper_list, harmonized_outcome_server_chen_paper_list,
            harmonized_outcome_local_chen_supp_met_list, harmonized_outcome_server_chen_supp_met_list,
            harmonized_outcome_local_chen_supp_met_ratio_list, harmonized_outcome_server_chen_supp_met_ratio_list,
            sure=T)

save.image("data/MR/harmonized.RData")
