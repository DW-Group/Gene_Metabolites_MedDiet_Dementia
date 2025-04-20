library(dplyr)
library(tidyr)
library(lmtest)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

load("data/merged_data_baseline_nhs_2023.Rdata")
source("source_covar_list_Gene_Metabolite_Dementia_analysis.R")

##-----------------------------------------
## Model
##-----------------------------------------

##-------------------------------
## Interaction for APOE4 
##-------------------------------

tics_list <- c("av_ntotal","av_nverbl","av_zscre")

model_data <- merged_data %>% filter(!is.na(av_ntotal))

main_inter_res_list <- list();  i <- 1

for (met in met_selected){
  
  print(met)
  tmp_met <- model_data[[met]]
  
  for (tics in tics_list){
    
    tmp_mod <- tryCatch(glm(reformulate(c("tmp_met * APOE4_carrier", covar_tics_list[["PGS002280"]]), response=tics), 
                             data = model_data),error=function(e) e,warning=function(w) w)
  
    main_inter_res_list[[i]] <- c(met, tics, "inter_linear", "m1", "APOE4_carrier", if(length(tmp_mod) > 5) c(sum(!is.na(model_data[[tics]])), nobs(tmp_mod), summary(tmp_mod)$coefficients["tmp_met",c(1,2,4)], NA, NA, NA, summary(tmp_mod)$coefficients["APOE4_carrier",c(1,2,4)], NA, NA, NA, summary(tmp_mod)$coefficients["tmp_met:APOE4_carrier",c(1,2,4)], NA) else c(tmp_mod$message,rep(NA,17))); i <- i + 1

  }
}

main_inter_res_dat <- as.data.frame(do.call("rbind", main_inter_res_list))
colnames(main_inter_res_dat) <- c("hmdb_id","outcome","model_type","model","covariates","n","n_events","beta_met","se_met","p_met","beta_gene_1","se_gene_1","p_gene_1","beta_gene_2","se_gene_2","p_gene_2","beta_inter_met_gene_1","se_inter_met_gene_1","p_inter_met_gene_1","beta_inter_met_gene_2","se_inter_met_gene_2","p_inter_met_gene_2","p_LRT")

main_inter_res_dat <- main_inter_res_dat %>%
  group_by(model, outcome) %>%
  mutate(p_fdr_inter_met_gene_2 = p.adjust(p_inter_met_gene_2, "BH"))

table(is.na(main_inter_res_dat$beta_inter_met_gene_2))

main_inter_res_dat <- merge(main_inter_res_dat, fdata_filtered, by="hmdb_id") 

write.table(main_inter_res_dat, "results/met_tics_interaction.txt", row.names=F, quote=F, sep="\t")