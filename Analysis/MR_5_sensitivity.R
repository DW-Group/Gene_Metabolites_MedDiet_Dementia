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
source("scripts/source_outcome_sample_size.R")

##-----------------------------------------
## Sensitivity analysis
##-----------------------------------------

method_heter_list <- c("mr_ivw","mr_egger_regression")

##--------------------------
## Chen_NatGenet_2023
##--------------------------

### Instruments selected by the original study

# Heterogeneity 

mr_heter_outcome_server_chen_paper_list <- list()
for (pheno in names(harmonized_outcome_server_chen_paper_list)){
  mr_heter_outcome_server_chen_paper_list[[pheno]] <- mr_heterogeneity(harmonized_outcome_server_chen_paper_list[[pheno]], method_list=method_heter_list)
}

mr_heter_outcome_local_chen_paper_list <- list()
for (pheno in names(harmonized_outcome_local_chen_paper_list)){
  mr_heter_outcome_local_chen_paper_list[[pheno]] <- mr_heterogeneity(harmonized_outcome_local_chen_paper_list[[pheno]], method_list=method_heter_list)
}

# Horizontal pleiotropy

mr_pleio_outcome_server_chen_paper_list <- list()
for (pheno in names(harmonized_outcome_server_chen_paper_list)){
  mr_pleio_outcome_server_chen_paper_list[[pheno]] <- mr_pleiotropy_test(harmonized_outcome_server_chen_paper_list[[pheno]])
}

mr_pleio_outcome_local_chen_paper_list <- list()
for (pheno in names(harmonized_outcome_local_chen_paper_list)){
  mr_pleio_outcome_local_chen_paper_list[[pheno]] <- mr_pleiotropy_test(harmonized_outcome_local_chen_paper_list[[pheno]])
}

# MR Steiger directionality test

mr_steiger_outcome_server_chen_paper_list <- list()
for (pheno in names(harmonized_outcome_server_chen_paper_list)){
  tmp_dat <- merge(harmonized_outcome_server_chen_paper_list[[pheno]], n_missing, by="id.outcome", suffixes=c("","_new"), all.x=T) %>%
    mutate(samplesize.outcome = ifelse(is.na(samplesize.outcome), samplesize.outcome_new, samplesize.outcome),
           pval.exposure = ifelse(pval.exposure==0, 1e-200, pval.exposure),
           pval.outcome = ifelse(pval.outcome==0, 1e-200, pval.outcome))
  mr_steiger_outcome_server_chen_paper_list[[pheno]] <- directionality_test(tmp_dat)
}

mr_steiger_outcome_local_chen_paper_list <- list()
for (pheno in names(harmonized_outcome_local_chen_paper_list)){
  tmp_dat <- merge(harmonized_outcome_local_chen_paper_list[[pheno]], n_missing, by.x="outcome",by.y="id.outcome", suffixes=c("","_new"), all.x=T) %>%
    mutate(samplesize.outcome = ifelse(is.na(samplesize.outcome), samplesize.outcome_new, samplesize.outcome),
           pval.exposure = ifelse(pval.exposure==0, 1e-200, pval.exposure),
           pval.outcome = ifelse(pval.outcome==0, 1e-200, pval.outcome))
  mr_steiger_outcome_local_chen_paper_list[[pheno]] <- directionality_test(tmp_dat)
}

### All significant associations for metabolites

# Heterogeneity 

mr_heter_outcome_server_chen_supp_met_list <- list()
for (pheno in names(harmonized_outcome_server_chen_supp_met_list)){
  mr_heter_outcome_server_chen_supp_met_list[[pheno]] <- mr_heterogeneity(harmonized_outcome_server_chen_supp_met_list[[pheno]], method_list=method_heter_list)
}

mr_heter_outcome_local_chen_supp_met_list <- list()
for (pheno in names(harmonized_outcome_local_chen_supp_met_list)){
  mr_heter_outcome_local_chen_supp_met_list[[pheno]] <- mr_heterogeneity(harmonized_outcome_local_chen_supp_met_list[[pheno]], method_list=method_heter_list)
}

# Horizontal pleiotropy

mr_pleio_outcome_server_chen_supp_met_list <- list()
for (pheno in names(harmonized_outcome_server_chen_supp_met_list)){
  mr_pleio_outcome_server_chen_supp_met_list[[pheno]] <- mr_pleiotropy_test(harmonized_outcome_server_chen_supp_met_list[[pheno]])
}

mr_pleio_outcome_local_chen_supp_met_list <- list()
for (pheno in names(harmonized_outcome_local_chen_supp_met_list)){
  mr_pleio_outcome_local_chen_supp_met_list[[pheno]] <- mr_pleiotropy_test(harmonized_outcome_local_chen_supp_met_list[[pheno]])
}

# MR Steiger directionality test

mr_steiger_outcome_server_chen_supp_met_list <- list()
for (pheno in names(harmonized_outcome_server_chen_supp_met_list)){
  tmp_dat <- merge(harmonized_outcome_server_chen_supp_met_list[[pheno]], n_missing, by="id.outcome", suffixes=c("","_new"), all.x=T) %>%
    mutate(samplesize.outcome = ifelse(is.na(samplesize.outcome), samplesize.outcome_new, samplesize.outcome),
           pval.exposure = ifelse(pval.exposure==0, 1e-200, pval.exposure),
           pval.outcome = ifelse(pval.outcome==0, 1e-200, pval.outcome))
  mr_steiger_outcome_server_chen_supp_met_list[[pheno]] <- directionality_test(tmp_dat)
}

mr_steiger_outcome_local_chen_supp_met_list <- list()
for (pheno in names(harmonized_outcome_local_chen_supp_met_list)){
  tmp_dat <- merge(harmonized_outcome_local_chen_supp_met_list[[pheno]], n_missing, by.x="outcome",by.y="id.outcome", suffixes=c("","_new"), all.x=T) %>%
    mutate(samplesize.outcome = ifelse(is.na(samplesize.outcome), samplesize.outcome_new, samplesize.outcome),
           pval.exposure = ifelse(pval.exposure==0, 1e-200, pval.exposure),
           pval.outcome = ifelse(pval.outcome==0, 1e-200, pval.outcome))
  mr_steiger_outcome_local_chen_supp_met_list[[pheno]] <- directionality_test(tmp_dat)
}

### All significant associations for metabolite ratios

# Heterogeneity 

mr_heter_outcome_server_chen_supp_met_ratio_list <- list()
for (pheno in names(harmonized_outcome_server_chen_supp_met_ratio_list)){
  mr_heter_outcome_server_chen_supp_met_ratio_list[[pheno]] <- mr_heterogeneity(harmonized_outcome_server_chen_supp_met_ratio_list[[pheno]], method_list=method_heter_list)
}

mr_heter_outcome_local_chen_supp_met_ratio_list <- list()
for (pheno in names(harmonized_outcome_local_chen_supp_met_ratio_list)){
  mr_heter_outcome_local_chen_supp_met_ratio_list[[pheno]] <- mr_heterogeneity(harmonized_outcome_local_chen_supp_met_ratio_list[[pheno]], method_list=method_heter_list)
}

# Horizontal pleiotropy

mr_pleio_outcome_server_chen_supp_met_ratio_list <- list()
for (pheno in names(harmonized_outcome_server_chen_supp_met_ratio_list)){
  mr_pleio_outcome_server_chen_supp_met_ratio_list[[pheno]] <- mr_pleiotropy_test(harmonized_outcome_server_chen_supp_met_ratio_list[[pheno]])
}

mr_pleio_outcome_local_chen_supp_met_ratio_list <- list()
for (pheno in names(harmonized_outcome_local_chen_supp_met_ratio_list)){
  mr_pleio_outcome_local_chen_supp_met_ratio_list[[pheno]] <- mr_pleiotropy_test(harmonized_outcome_local_chen_supp_met_ratio_list[[pheno]])
}

# MR Steiger directionality test

mr_steiger_outcome_server_chen_supp_met_ratio_list <- list()
for (pheno in names(harmonized_outcome_server_chen_supp_met_ratio_list)){
  tmp_dat <- merge(harmonized_outcome_server_chen_supp_met_ratio_list[[pheno]], n_missing, by="id.outcome", suffixes=c("","_new"), all.x=T) %>%
    mutate(samplesize.outcome = ifelse(is.na(samplesize.outcome), samplesize.outcome_new, samplesize.outcome),
           pval.exposure = ifelse(pval.exposure==0, 1e-200, pval.exposure),
           pval.outcome = ifelse(pval.outcome==0, 1e-200, pval.outcome))
  mr_steiger_outcome_server_chen_supp_met_ratio_list[[pheno]] <- directionality_test(tmp_dat)
}

mr_steiger_outcome_local_chen_supp_met_ratio_list <- list()
for (pheno in names(harmonized_outcome_local_chen_supp_met_ratio_list)){
  tmp_dat <- merge(harmonized_outcome_local_chen_supp_met_ratio_list[[pheno]], n_missing, by.x="outcome",by.y="id.outcome", suffixes=c("","_new"), all.x=T) %>%
    mutate(samplesize.outcome = ifelse(is.na(samplesize.outcome), samplesize.outcome_new, samplesize.outcome),
           pval.exposure = ifelse(pval.exposure==0, 1e-200, pval.exposure),
           pval.outcome = ifelse(pval.outcome==0, 1e-200, pval.outcome))
  mr_steiger_outcome_local_chen_supp_met_ratio_list[[pheno]] <- directionality_test(tmp_dat)
}

##-----------------------------------------
## Combine results
##-----------------------------------------

mr_heter_outcome_server_chen_paper <- as.data.frame(do.call("rbind",mr_heter_outcome_server_chen_paper_list)) %>%
  mutate(exposure_source = "chen_paper", outcome_source = "server")
mr_heter_outcome_server_chen_supp_met <- as.data.frame(do.call("rbind",mr_heter_outcome_server_chen_supp_met_list)) %>%
  mutate(exposure_source = "chen_supp_met", outcome_source = "server")
mr_heter_outcome_server_chen_supp_met_ratio <- as.data.frame(do.call("rbind",mr_heter_outcome_server_chen_supp_met_ratio_list)) %>%
  mutate(exposure_source = "chen_supp_met_ratio", outcome_source = "server")

mr_heter_outcome_local_chen_paper <- as.data.frame(do.call("rbind",mr_heter_outcome_local_chen_paper_list)) %>%
  mutate(exposure_source = "chen_paper", outcome_source = "local")
mr_heter_outcome_local_chen_supp_met <- as.data.frame(do.call("rbind",mr_heter_outcome_local_chen_supp_met_list)) %>%
  mutate(exposure_source = "chen_supp_met", outcome_source = "local")
mr_heter_outcome_local_chen_supp_met_ratio <- as.data.frame(do.call("rbind",mr_heter_outcome_local_chen_supp_met_ratio_list)) %>%
  mutate(exposure_source = "chen_supp_met_ratio", outcome_source = "local")

mr_heter_outcome_server <- as.data.frame(do.call("rbind",list(mr_heter_outcome_server_chen_paper,
                                                            mr_heter_outcome_server_chen_supp_met,
                                                            mr_heter_outcome_server_chen_supp_met_ratio)))

mr_heter_outcome_local <- as.data.frame(do.call("rbind",list(mr_heter_outcome_local_chen_paper,
                                                           mr_heter_outcome_local_chen_supp_met,
                                                           mr_heter_outcome_local_chen_supp_met_ratio)))

mr_heter_all <- as.data.frame(rbind(mr_heter_outcome_server,mr_heter_outcome_local))

mr_pleio_outcome_server_chen_paper <- as.data.frame(do.call("rbind",mr_pleio_outcome_server_chen_paper_list)) %>%
  mutate(exposure_source = "chen_paper", outcome_source = "server")
mr_pleio_outcome_server_chen_supp_met <- as.data.frame(do.call("rbind",mr_pleio_outcome_server_chen_supp_met_list)) %>%
  mutate(exposure_source = "chen_supp_met", outcome_source = "server")
mr_pleio_outcome_server_chen_supp_met_ratio <- as.data.frame(do.call("rbind",mr_pleio_outcome_server_chen_supp_met_ratio_list)) %>%
  mutate(exposure_source = "chen_supp_met_ratio", outcome_source = "server")

mr_pleio_outcome_local_chen_paper <- as.data.frame(do.call("rbind",mr_pleio_outcome_local_chen_paper_list)) %>%
  mutate(exposure_source = "chen_paper", outcome_source = "local")
mr_pleio_outcome_local_chen_supp_met <- as.data.frame(do.call("rbind",mr_pleio_outcome_local_chen_supp_met_list)) %>%
  mutate(exposure_source = "chen_supp_met", outcome_source = "local")
mr_pleio_outcome_local_chen_supp_met_ratio <- as.data.frame(do.call("rbind",mr_pleio_outcome_local_chen_supp_met_ratio_list)) %>%
  mutate(exposure_source = "chen_supp_met_ratio", outcome_source = "local")

mr_pleio_outcome_server <- as.data.frame(do.call("rbind",list(mr_pleio_outcome_server_chen_paper,
                                                              mr_pleio_outcome_server_chen_supp_met,
                                                              mr_pleio_outcome_server_chen_supp_met_ratio)))

mr_pleio_outcome_local <- as.data.frame(do.call("rbind",list(mr_pleio_outcome_local_chen_paper,
                                                             mr_pleio_outcome_local_chen_supp_met,
                                                             mr_pleio_outcome_local_chen_supp_met_ratio)))

mr_pleio_all <- as.data.frame(rbind(mr_pleio_outcome_server,mr_pleio_outcome_local))

mr_steiger_outcome_server_chen_paper <- as.data.frame(do.call("rbind",mr_steiger_outcome_server_chen_paper_list)) %>%
  mutate(exposure_source = "chen_paper", outcome_source = "server")
mr_steiger_outcome_server_chen_supp_met <- as.data.frame(do.call("rbind",mr_steiger_outcome_server_chen_supp_met_list)) %>%
  mutate(exposure_source = "chen_supp_met", outcome_source = "server")
mr_steiger_outcome_server_chen_supp_met_ratio <- as.data.frame(do.call("rbind",mr_steiger_outcome_server_chen_supp_met_ratio_list)) %>%
  mutate(exposure_source = "chen_supp_met_ratio", outcome_source = "server")

mr_steiger_outcome_local_chen_paper <- as.data.frame(do.call("rbind",mr_steiger_outcome_local_chen_paper_list)) %>%
  mutate(exposure_source = "chen_paper", outcome_source = "local")
mr_steiger_outcome_local_chen_supp_met <- as.data.frame(do.call("rbind",mr_steiger_outcome_local_chen_supp_met_list)) %>%
  mutate(exposure_source = "chen_supp_met", outcome_source = "local")
mr_steiger_outcome_local_chen_supp_met_ratio <- as.data.frame(do.call("rbind",mr_steiger_outcome_local_chen_supp_met_ratio_list)) %>%
  mutate(exposure_source = "chen_supp_met_ratio", outcome_source = "local")

mr_steiger_outcome_server <- as.data.frame(do.call("rbind",list(mr_steiger_outcome_server_chen_paper,
                                                              mr_steiger_outcome_server_chen_supp_met,
                                                              mr_steiger_outcome_server_chen_supp_met_ratio)))

mr_steiger_outcome_local <- as.data.frame(do.call("rbind",list(mr_steiger_outcome_local_chen_paper,
                                                             mr_steiger_outcome_local_chen_supp_met,
                                                             mr_steiger_outcome_local_chen_supp_met_ratio)))

mr_steiger_all <- as.data.frame(rbind(mr_steiger_outcome_server,mr_steiger_outcome_local))

##-----------------------------------------
## Save data
##-----------------------------------------

gdata::keep(mr_heter_outcome_server, mr_heter_outcome_local, mr_heter_all,
            mr_pleio_outcome_server, mr_pleio_outcome_local, mr_pleio_all,
            mr_steiger_outcome_server, mr_steiger_outcome_local, mr_steiger_all,
            sure=T)

save.image("results/MR/mr_results_sensitivity.RData")
