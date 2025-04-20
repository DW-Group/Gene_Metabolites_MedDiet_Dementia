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
## Local LD matrix
##-----------------------------------------

load("data/MR/instruments/all_instru_ld_mat.RData") 

class(all_instru_ld_mat)

all_instru_ld_mat_dict <- as.data.frame(str_split(colnames(all_instru_ld_mat), "_", simplify=T)) %>%
  mutate(rsid = V1,
         A1 = V2,
         A2 = V3,
         rsid_with_alleles = colnames(all_instru_ld_mat)) %>%
  subset(select=c(rsid_with_alleles,rsid,A1,A2))

##-----------------------------------------
## IVW MR analysis accounting for LD
##-----------------------------------------

method_list <- c("mr_ivw")

##--------------------------
## Chen_NatGenet_2023
##--------------------------

### Instruments selected by the original study

mr_ivw_ld_res_outcome_server_chen_paper_list <- list(); i <- 1
for (pheno in names(harmonized_outcome_server_chen_paper_list)){
  tmp_harmonized <- harmonized_outcome_server_chen_paper_list[[pheno]]
  if (length(unique(tmp_harmonized$SNP)) > 1){
    tmp_MRInput <- dat_to_MRInput_yl(tmp_harmonized, 
                                     get_correlations=T, local_correlations=T, local_ld_mat=all_instru_ld_mat, local_ld_mat_dict=all_instru_ld_mat_dict)
    for (out in names(tmp_MRInput)){
      tmp_ivw_res <- tryCatch(MendelianRandomization::mr_ivw(tmp_MRInput[[out]], correl=TRUE),error=function(e) e,warning=function(w) w)
      mr_ivw_ld_res_outcome_server_chen_paper_list[[i]] <- c("server", "ivw_random_ld", "chen_paper", pheno, gsub(paste0(pheno,"."),"",out, fixed=T), if (length(tmp_ivw_res) == 1) c(tmp_ivw_res$Estimate, tmp_ivw_res$StdError, tmp_ivw_res$Pvalue, tmp_ivw_res$SNPs) else c(rep(NA,2),tmp_ivw_res$message,NA)); i <- i + 1
    }
  }
}

mr_ivw_ld_res_outcome_local_chen_paper_list <- list(); i <- 1
for (pheno in names(harmonized_outcome_local_chen_paper_list)){
  tmp_harmonized <- harmonized_outcome_local_chen_paper_list[[pheno]]
  if (length(unique(tmp_harmonized$SNP)) > 1){
    tmp_MRInput <- dat_to_MRInput_yl(tmp_harmonized, 
                                     get_correlations=T, local_correlations=T, local_ld_mat=all_instru_ld_mat, local_ld_mat_dict=all_instru_ld_mat_dict)
    for (out in names(tmp_MRInput)){
      tmp_ivw_res <- tryCatch(MendelianRandomization::mr_ivw(tmp_MRInput[[out]], correl=TRUE),error=function(e) e,warning=function(w) w)
      mr_ivw_ld_res_outcome_local_chen_paper_list[[i]] <- c("local", "ivw_random_ld", "chen_paper", pheno, gsub(paste0(pheno,"."),"",out, fixed=T), if (length(tmp_ivw_res) == 1) c(tmp_ivw_res$Estimate, tmp_ivw_res$StdError, tmp_ivw_res$Pvalue, tmp_ivw_res$SNPs) else c(rep(NA,2),tmp_ivw_res$message,NA)); i <- i + 1
    }
  }
}

### All significant associations for metabolites

mr_ivw_ld_res_outcome_server_chen_supp_met_list <- list(); i <- 1
for (pheno in names(harmonized_outcome_server_chen_supp_met_list)){
  tmp_harmonized <- harmonized_outcome_server_chen_supp_met_list[[pheno]]
  if (length(unique(tmp_harmonized$SNP)) > 1){
    tmp_MRInput <- dat_to_MRInput_yl(tmp_harmonized, 
                                     get_correlations=T, local_correlations=T, local_ld_mat=all_instru_ld_mat, local_ld_mat_dict=all_instru_ld_mat_dict)
    for (out in names(tmp_MRInput)){
      tmp_ivw_res <- tryCatch(MendelianRandomization::mr_ivw(tmp_MRInput[[out]], correl=TRUE),error=function(e) e,warning=function(w) w)
      mr_ivw_ld_res_outcome_server_chen_supp_met_list[[i]] <- c("server", "ivw_random_ld", "chen_supp_met", pheno, gsub(paste0(pheno,"."),"",out, fixed=T), if (length(tmp_ivw_res) == 1) c(tmp_ivw_res$Estimate, tmp_ivw_res$StdError, tmp_ivw_res$Pvalue, tmp_ivw_res$SNPs) else c(rep(NA,2),tmp_ivw_res$message,NA)); i <- i + 1
    }
  }
}

mr_ivw_ld_res_outcome_local_chen_supp_met_list <- list(); i <- 1
for (pheno in names(harmonized_outcome_local_chen_supp_met_list)){
  tmp_harmonized <- harmonized_outcome_local_chen_supp_met_list[[pheno]]
  if (length(unique(tmp_harmonized$SNP)) > 1){
    tmp_MRInput <- dat_to_MRInput_yl(tmp_harmonized, 
                                     get_correlations=T, local_correlations=T, local_ld_mat=all_instru_ld_mat, local_ld_mat_dict=all_instru_ld_mat_dict)
    for (out in names(tmp_MRInput)){
      tmp_ivw_res <- tryCatch(MendelianRandomization::mr_ivw(tmp_MRInput[[out]], correl=TRUE),error=function(e) e,warning=function(w) w)
      mr_ivw_ld_res_outcome_local_chen_supp_met_list[[i]] <- c("local", "ivw_random_ld", "chen_supp_met", pheno, gsub(paste0(pheno,"."),"",out, fixed=T), if (length(tmp_ivw_res) == 1) c(tmp_ivw_res$Estimate, tmp_ivw_res$StdError, tmp_ivw_res$Pvalue, tmp_ivw_res$SNPs) else c(rep(NA,2),tmp_ivw_res$message,NA)); i <- i + 1
    }
  }
}

### All significant associations for metabolite ratios

mr_ivw_ld_res_outcome_server_chen_supp_met_ratio_list <- list(); i <- 1
for (pheno in names(harmonized_outcome_server_chen_supp_met_ratio_list)){
  tmp_harmonized <- harmonized_outcome_server_chen_supp_met_ratio_list[[pheno]]
  if (length(unique(tmp_harmonized$SNP)) > 1){
    tmp_MRInput <- dat_to_MRInput_yl(tmp_harmonized, 
                                     get_correlations=T, local_correlations=T, local_ld_mat=all_instru_ld_mat, local_ld_mat_dict=all_instru_ld_mat_dict)
    for (out in names(tmp_MRInput)){
      tmp_ivw_res <- tryCatch(MendelianRandomization::mr_ivw(tmp_MRInput[[out]], correl=TRUE),error=function(e) e,warning=function(w) w)
      mr_ivw_ld_res_outcome_server_chen_supp_met_ratio_list[[i]] <- c("server", "ivw_random_ld", "chen_supp_met_ratio", pheno, gsub(paste0(pheno,"."),"",out, fixed=T), if (length(tmp_ivw_res) == 1) c(tmp_ivw_res$Estimate, tmp_ivw_res$StdError, tmp_ivw_res$Pvalue, tmp_ivw_res$SNPs) else c(rep(NA,2),tmp_ivw_res$message,NA)); i <- i + 1
    }
  }
}

mr_ivw_ld_res_outcome_local_chen_supp_met_ratio_list <- list(); i <- 1
for (pheno in names(harmonized_outcome_local_chen_supp_met_ratio_list)){
  tmp_harmonized <- harmonized_outcome_local_chen_supp_met_ratio_list[[pheno]]
  if (length(unique(tmp_harmonized$SNP)) > 1){
    tmp_MRInput <- dat_to_MRInput_yl(tmp_harmonized, 
                                     get_correlations=T, local_correlations=T, local_ld_mat=all_instru_ld_mat, local_ld_mat_dict=all_instru_ld_mat_dict)
    for (out in names(tmp_MRInput)){
      tmp_ivw_res <- tryCatch(MendelianRandomization::mr_ivw(tmp_MRInput[[out]], correl=TRUE),error=function(e) e,warning=function(w) w)
      mr_ivw_ld_res_outcome_local_chen_supp_met_ratio_list[[i]] <- c("local", "ivw_random_ld", "chen_supp_met_ratio", pheno, gsub(paste0(pheno,"."),"",out, fixed=T), if (length(tmp_ivw_res) == 1) c(tmp_ivw_res$Estimate, tmp_ivw_res$StdError, tmp_ivw_res$Pvalue, tmp_ivw_res$SNPs) else c(rep(NA,2),tmp_ivw_res$message,NA)); i <- i + 1
    }
  }
}

##-----------------------------------------
## Combine results
##-----------------------------------------

mr_ivw_ld_res_outcome_server_chen_paper <- as.data.frame(do.call("rbind",mr_ivw_ld_res_outcome_server_chen_paper_list))
mr_ivw_ld_res_outcome_server_chen_supp_met <- as.data.frame(do.call("rbind",mr_ivw_ld_res_outcome_server_chen_supp_met_list))
mr_ivw_ld_res_outcome_server_chen_supp_met_ratio <- as.data.frame(do.call("rbind",mr_ivw_ld_res_outcome_server_chen_supp_met_ratio_list))

mr_ivw_ld_res_outcome_local_chen_paper <- as.data.frame(do.call("rbind",mr_ivw_ld_res_outcome_local_chen_paper_list))
mr_ivw_ld_res_outcome_local_chen_supp_met <- as.data.frame(do.call("rbind",mr_ivw_ld_res_outcome_local_chen_supp_met_list))
mr_ivw_ld_res_outcome_local_chen_supp_met_ratio <- as.data.frame(do.call("rbind",mr_ivw_ld_res_outcome_local_chen_supp_met_ratio_list))

mr_ivw_ld_res_all <- as.data.frame(do.call("rbind",list(mr_ivw_ld_res_outcome_server_chen_paper,
                                                        mr_ivw_ld_res_outcome_server_chen_supp_met,
                                                        mr_ivw_ld_res_outcome_server_chen_supp_met_ratio,
                                                        mr_ivw_ld_res_outcome_local_chen_paper,
                                                        mr_ivw_ld_res_outcome_local_chen_supp_met,
                                                        mr_ivw_ld_res_outcome_local_chen_supp_met_ratio)))

colnames(mr_ivw_ld_res_all) <- c("outcome_source","method","exposure_source","exposure","outcome","beta","se","p","nsnp")

mr_ivw_ld_res <- mr_ivw_ld_res_all %>% 
  filter(!is.na(beta)) %>% 
  mutate(p = as.numeric(p)) %>%
  arrange(p)

##-----------------------------------------
## Save data
##-----------------------------------------

gdata::keep(mr_ivw_ld_res_all, mr_ivw_ld_res,
            sure=T)

save.image("results/MR/mr_results_ivw_ld.RData")
