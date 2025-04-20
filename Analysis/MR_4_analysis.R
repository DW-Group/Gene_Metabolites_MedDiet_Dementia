library(data.table)
library(dplyr)
library(stringr)
library(readxl)
library(MendelianRandomization)
library(TwoSampleMR)
library(MRInstruments)
library(MRPRESSO)

rm(list=ls())

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

##-----------------------------------------
## Harmonized data
##-----------------------------------------

load("data/MR/harmonized.RData")

##-----------------------------------------
## Check numbers
##-----------------------------------------

### chen_paper

check_local_chen_paper_list <- c(); i <- 1
for(met in names(harmonized_outcome_local_chen_paper_list)){
  tmp <- harmonized_outcome_local_chen_paper_list[[met]] %>% 
    filter(outcome=="wightman_ad") 
  if (dim(tmp)[1]!=0){tmp <- tmp %>% mutate(exposure=met)}
  check_local_chen_paper_list[[i]] <- tmp; i <- i + 1
}
check_local_chen_paper <- as.data.frame(do.call("rbind",check_local_chen_paper_list)) %>%
  subset(select=c(outcome,exposure,SNP))

check_server_chen_paper_list <- c(); i <- 1
for(met in names(harmonized_outcome_server_chen_paper_list)){
  tmp <- harmonized_outcome_server_chen_paper_list[[met]] %>% 
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
for(met in names(harmonized_outcome_local_chen_supp_met_list)){
  tmp <- harmonized_outcome_local_chen_supp_met_list[[met]] %>% 
    filter(outcome=="wightman_ad") 
  if (dim(tmp)[1]!=0){tmp <- tmp %>% mutate(exposure=met)}
  check_local_chen_supp_met_list[[i]] <- tmp; i <- i + 1
}
check_local_chen_supp_met <- as.data.frame(do.call("rbind",check_local_chen_supp_met_list)) %>%
  subset(select=c(outcome,exposure,SNP))

check_server_chen_supp_met_list <- c(); i <- 1
for(met in names(harmonized_outcome_server_chen_supp_met_list)){
  tmp <- harmonized_outcome_server_chen_supp_met_list[[met]] %>% 
    filter(outcome %in% c("Any dementia (more controls excluded) || id:finn-b-KRA_PSY_DEMENTIA_EXMORE",
                          "Cognitive performance || id:ebi-a-GCST006572",
                          "Vascular dementia (F5_VASCDEM) || id:finn-b-F5_VASCDEM")) 
  if (dim(tmp)[1]!=0){tmp <- tmp %>% mutate(exposure=met)}
  check_server_chen_supp_met_list[[i]] <- tmp; i <- i + 1
}
check_server_chen_supp_met <- as.data.frame(do.call("rbind",check_server_chen_supp_met_list)) %>%
  subset(select=c(outcome,exposure,SNP))

check_chen_supp_met <- as.data.frame(rbind(check_local_chen_supp_met,check_server_chen_supp_met)) %>%
  filter(!(exposure %in% names(harmonized_outcome_local_chen_paper_list)))

### chen_supp_met_ratio

check_local_chen_supp_met_ratio_list <- c(); i <- 1
for(met in names(harmonized_outcome_local_chen_supp_met_ratio_list)){
  tmp <- harmonized_outcome_local_chen_supp_met_ratio_list[[met]] %>% 
    filter(outcome=="wightman_ad")
  if (dim(tmp)[1]!=0){tmp <- tmp %>% mutate(exposure=met)}
  check_local_chen_supp_met_ratio_list[[i]] <- tmp; i <- i + 1
}
check_local_chen_supp_met_ratio <- as.data.frame(do.call("rbind",check_local_chen_supp_met_ratio_list)) %>%
  subset(select=c(outcome,exposure,SNP))

check_server_chen_supp_met_ratio_list <- c(); i <- 1
for(met in names(harmonized_outcome_server_chen_supp_met_ratio_list)){
  tmp <- harmonized_outcome_server_chen_supp_met_ratio_list[[met]] %>% 
    filter(outcome %in% c("Any dementia (more controls excluded) || id:finn-b-KRA_PSY_DEMENTIA_EXMORE",
                          "Cognitive performance || id:ebi-a-GCST006572",
                          "Vascular dementia (F5_VASCDEM) || id:finn-b-F5_VASCDEM")) 
  if (dim(tmp)[1]!=0){tmp <- tmp %>% mutate(exposure=met)}
  check_server_chen_supp_met_ratio_list[[i]] <- tmp; i <- i + 1
}
check_server_chen_supp_met_ratio <- as.data.frame(do.call("rbind",check_server_chen_supp_met_ratio_list)) %>%
  subset(select=c(outcome,exposure,SNP))

check_chen_supp_met_ratio <- as.data.frame(rbind(check_local_chen_supp_met_ratio,check_server_chen_supp_met_ratio)) %>%
  filter(!(exposure %in% names(harmonized_outcome_local_chen_paper_list)))

##-----------------------------------------
## MR analysis
##-----------------------------------------

method_list <- c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median")

##--------------------------
## Chen_NatGenet_2023
##--------------------------

### Instruments selected by the original study

mr_res_outcome_server_chen_paper_list <- list()
for (pheno in names(harmonized_outcome_server_chen_paper_list)){
  mr_res_outcome_server_chen_paper_list[[pheno]] <- mr(harmonized_outcome_server_chen_paper_list[[pheno]], method_list=method_list)
}

mr_res_outcome_local_chen_paper_list <- list()
for (pheno in names(harmonized_outcome_local_chen_paper_list)){
  mr_res_outcome_local_chen_paper_list[[pheno]] <- mr(harmonized_outcome_local_chen_paper_list[[pheno]], method_list=method_list)
}

### All significant associations for metabolites

mr_res_outcome_server_chen_supp_met_list <- list()
for (pheno in names(harmonized_outcome_server_chen_supp_met_list)){
  mr_res_outcome_server_chen_supp_met_list[[pheno]] <- mr(harmonized_outcome_server_chen_supp_met_list[[pheno]], method_list=method_list)
}

mr_res_outcome_local_chen_supp_met_list <- list()
for (pheno in names(harmonized_outcome_local_chen_supp_met_list)){
  mr_res_outcome_local_chen_supp_met_list[[pheno]] <- mr(harmonized_outcome_local_chen_supp_met_list[[pheno]], method_list=method_list)
}

### All significant associations for metabolite ratios

mr_res_outcome_server_chen_supp_met_ratio_list <- list()
for (pheno in names(harmonized_outcome_server_chen_supp_met_ratio_list)){
  mr_res_outcome_server_chen_supp_met_ratio_list[[pheno]] <- mr(harmonized_outcome_server_chen_supp_met_ratio_list[[pheno]], method_list=method_list)
}

mr_res_outcome_local_chen_supp_met_ratio_list <- list()
for (pheno in names(harmonized_outcome_local_chen_supp_met_ratio_list)){
  mr_res_outcome_local_chen_supp_met_ratio_list[[pheno]] <- mr(harmonized_outcome_local_chen_supp_met_ratio_list[[pheno]], method_list=method_list)
}

##-----------------------------------------
## Combine results
##-----------------------------------------

mr_res_outcome_server_chen_paper <- as.data.frame(do.call("rbind",mr_res_outcome_server_chen_paper_list)) %>%
  mutate(exposure_source = "chen_paper", outcome_source = "server")
mr_res_outcome_server_chen_supp_met <- as.data.frame(do.call("rbind",mr_res_outcome_server_chen_supp_met_list)) %>%
  mutate(exposure_source = "chen_supp_met", outcome_source = "server")
mr_res_outcome_server_chen_supp_met_ratio <- as.data.frame(do.call("rbind",mr_res_outcome_server_chen_supp_met_ratio_list)) %>%
  mutate(exposure_source = "chen_supp_met_ratio", outcome_source = "server")

mr_res_outcome_local_chen_paper <- as.data.frame(do.call("rbind",mr_res_outcome_local_chen_paper_list)) %>%
  mutate(exposure_source = "chen_paper", outcome_source = "local")
mr_res_outcome_local_chen_supp_met <- as.data.frame(do.call("rbind",mr_res_outcome_local_chen_supp_met_list)) %>%
  mutate(exposure_source = "chen_supp_met", outcome_source = "local")
mr_res_outcome_local_chen_supp_met_ratio <- as.data.frame(do.call("rbind",mr_res_outcome_local_chen_supp_met_ratio_list)) %>%
  mutate(exposure_source = "chen_supp_met_ratio", outcome_source = "local")

mr_res_outcome_server <- as.data.frame(do.call("rbind",list(mr_res_outcome_server_chen_paper,
                                                            mr_res_outcome_server_chen_supp_met,
                                                            mr_res_outcome_server_chen_supp_met_ratio)))

mr_res_outcome_local <- as.data.frame(do.call("rbind",list(mr_res_outcome_local_chen_paper,
                                                            mr_res_outcome_local_chen_supp_met,
                                                            mr_res_outcome_local_chen_supp_met_ratio)))

mr_res_all <- as.data.frame(rbind(mr_res_outcome_server,mr_res_outcome_local))

##-----------------------------------------
## Save data
##-----------------------------------------

gdata::keep(mr_res_outcome_server, mr_res_outcome_local, mr_res_all,
            sure=T)

save.image("results/MR/mr_results.RData")
