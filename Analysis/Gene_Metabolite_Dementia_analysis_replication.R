library(dplyr)
library(tidyr)
library(survival)
library(lmtest)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

load("data/merged_data_baseline_hpfs_2023.Rdata")
source("source_covar_list_Gene_Metabolite_Dementia_replication.R")

##-----------------------------------------
## Model
##-----------------------------------------

##-------------------------------
## Interaction for APOE4 
##-------------------------------

model_data <- merged_data_hpfs

main_inter_res_list <- list();  i <- 1
  
for (met in met_selected_hpfs){
  
  print(met)
  tmp_met <- model_data[[met]]
  
  tmp_mod <- tryCatch(coxph(reformulate(c("tmp_met * APOE4_carrier", covar_dem_list[["PGS002280"]]), response=response_subtype), 
                              ties="breslow", data = model_data),error=function(e) e,warning=function(w) w)

  main_inter_res_list[[i]] <- c(met, "dementia", surv_dict$subtype[1], "inter_cox", "m1", "APOE4_carrier", if(length(tmp_mod) > 5) c(summary(tmp_mod)$n, nobs(tmp_mod), summary(tmp_mod)$coefficients["tmp_met",c(1,3,5)], NA, NA, NA, summary(tmp_mod)$coefficients["APOE4_carrier",c(1,3,5)], NA, NA, NA, summary(tmp_mod)$coefficients["tmp_met:APOE4_carrier",c(1,3,5)], NA) else c(tmp_mod$message,rep(NA,17))); i <- i + 1

}

main_inter_res_dat <- as.data.frame(do.call("rbind", main_inter_res_list))
colnames(main_inter_res_dat) <- c("hmdb_id","outcome","subtype","model_type","model","covariates","n","n_events","beta_met","se_met","p_met","beta_gene_1","se_gene_1","p_gene_1","beta_gene_2","se_gene_2","p_gene_2","beta_inter_met_gene_1","se_inter_met_gene_1","p_inter_met_gene_1","beta_inter_met_gene_2","se_inter_met_gene_2","p_inter_met_gene_2","p_LRT")

main_inter_res_dat <- main_inter_res_dat %>%
  group_by(subtype,model) %>%
  mutate(p_fdr_inter_met_gene_2 = p.adjust(p_inter_met_gene_2, "BH"),
         p_fdr_LRT = p.adjust(p_LRT, "BH")) %>%
  pivot_longer(cols = c(p_inter_met_gene_1,p_inter_met_gene_2), names_to = "pvalue_column", values_to = "pvalue") %>%
  mutate(two_inter_fdr = p.adjust(pvalue, method = "fdr")) %>%
  pivot_wider(names_from = "pvalue_column", values_from = c(pvalue, two_inter_fdr))

table(is.na(main_inter_res_dat$beta_inter_met_gene_2))

main_inter_res_dat <- merge(main_inter_res_dat, fdata_filtered_hpfs, by="hmdb_id") 

write.table(main_inter_res_dat, "results/met_dem_cox_baseline_interaction_hpfs_2023.txt", row.names=F, quote=F, sep="\t")
