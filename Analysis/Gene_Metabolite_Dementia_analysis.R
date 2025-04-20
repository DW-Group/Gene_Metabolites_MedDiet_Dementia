library(dplyr)
library(tidyr)
library(survival)
library(lmtest)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

load("data/merged_data_baseline_nhs_2023.Rdata")
load("data/genetic/ad_variants_dictionary.RData")
source("source_covar_list_Gene_Metabolite_Dementia_analysis.R")

##-----------------------------------------
## Model
##-----------------------------------------

##-------------------------------
## Main effect
##-------------------------------

model_data <- merged_data

main_res_list <- list();  i <- 1

for (met in met_selected){
  
  print(met)
  
  for (covar in names(covar_dem_list)) {
    
    tmp_mod <- tryCatch(coxph(reformulate(c(eval(met), covar_dem_list[[covar]]), response=response_subtype), ties="breslow", data = model_data),error=function(e) e,warning=function(w) w)
    main_res_list[[i]] <- c(met, "dementia", surv_dict$subtype[1], "main_cox", covar, if(length(tmp_mod) > 5) c(summary(tmp_mod)$n, nobs(tmp_mod), summary(tmp_mod)$coefficients[eval(met),c(1,3,5)]) else c(tmp_mod$message,NA,NA,NA,NA)); i <- i + 1
    
  }
  
}

main_res_dat <- as.data.frame(do.call("rbind", main_res_list))
colnames(main_res_dat) <- c("hmdb_id","outcome","subtype","model_type","covariates","n","n_events","beta","se","p")

main_res_dat <- main_res_dat %>%
  group_by(subtype, model_type) %>%
  mutate(p_fdr = p.adjust(p, "BH"))

table(is.na(main_res_dat$beta))

main_res_dat <- merge(main_res_dat, fdata_filtered, by="hmdb_id") 

write.table(main_res_dat, "results/met_dem_cox_baseline_nhs_2023.txt", row.names=F, quote=F, sep="\t")

##-------------------------------
## Interaction for APOE4 and PRS
##-------------------------------

model_data <- merged_data

main_inter_res_list <- list();  i <- 1

for (met in met_selected){
  
  print(met)
  tmp_met <- model_data[[met]]
  
  tmp_mod1 <- tryCatch(coxph(reformulate(c("tmp_met * APOE4_ncopy_1 + tmp_met * APOE4_ncopy_2", covar_dem_list[["PGS002280"]]), response=response_subtype), 
                            ties="breslow", data = model_data),error=function(e) e,warning=function(w) w)
  tmp_mod2 <- tryCatch(coxph(reformulate(c("tmp_met * APOE4_carrier", covar_dem_list[["PGS002280"]]), response=response_subtype), 
                              ties="breslow", data = model_data),error=function(e) e,warning=function(w) w)                  
  tmp_mod3 <- tryCatch(coxph(reformulate(c("tmp_met * PGS002280", covar_dem_list[["APOE4_2cat"]]), response=response_subtype), 
                              ties="breslow", data = model_data),error=function(e) e,warning=function(w) w)      
  tmp_mod4 <- tryCatch(coxph(reformulate(c("tmp_met * PGS000334", covar_dem_list[["APOE4_2cat"]]), response=response_subtype), 
                              ties="breslow", data = model_data),error=function(e) e,warning=function(w) w)  
  
  tmp_mod_lrt_full <- tryCatch(coxph(reformulate(c("tmp_met * PGS002280 + tmp_met * APOE4_ncopy_1 + tmp_met * APOE4_ncopy_2", covar_dem_list[["PGS002280"]]), response=response_subtype), 
                              ties="breslow", data = model_data),error=function(e) e,warning=function(w) w)   
  tmp_mod_lrt_nested <- tryCatch(coxph(reformulate(c("tmp_met * PGS002280 + tmp_met + APOE4_ncopy_1 + APOE4_ncopy_2", covar_dem_list[["PGS002280"]]), response=response_subtype), 
                                      ties="breslow", data = model_data),error=function(e) e,warning=function(w) w)   
  
  main_inter_res_list[[i]] <- c(met, "dementia", surv_dict$subtype[1], "inter_cox", "m1", "APOE4_2cat", if(length(tmp_mod1) > 5) c(summary(tmp_mod1)$n, nobs(tmp_mod1), summary(tmp_mod1)$coefficients["tmp_met",c(1,3,5)], summary(tmp_mod1)$coefficients["APOE4_ncopy_1",c(1,3,5)], summary(tmp_mod1)$coefficients["APOE4_ncopy_2",c(1,3,5)], summary(tmp_mod1)$coefficients["tmp_met:APOE4_ncopy_1",c(1,3,5)], summary(tmp_mod1)$coefficients["tmp_met:APOE4_ncopy_2",c(1,3,5)], NA) else c(tmp_mod1$message,rep(NA,17))); i <- i + 1
  main_inter_res_list[[i]] <- c(met, "dementia", surv_dict$subtype[1], "inter_cox", "m2", "APOE4_carrier", if(length(tmp_mod2) > 5) c(summary(tmp_mod2)$n, nobs(tmp_mod2), summary(tmp_mod2)$coefficients["tmp_met",c(1,3,5)], NA, NA, NA, summary(tmp_mod2)$coefficients["APOE4_carrier",c(1,3,5)], NA, NA, NA, summary(tmp_mod2)$coefficients["tmp_met:APOE4_carrier",c(1,3,5)], NA) else c(tmp_mod2$message,rep(NA,17))); i <- i + 1
  main_inter_res_list[[i]] <- c(met, "dementia", surv_dict$subtype[1], "inter_cox", "m3", "PGS002280", if(length(tmp_mod3) > 5) c(summary(tmp_mod3)$n, nobs(tmp_mod3), summary(tmp_mod3)$coefficients["tmp_met",c(1,3,5)], NA, NA, NA, summary(tmp_mod3)$coefficients["PGS002280",c(1,3,5)], NA, NA, NA, summary(tmp_mod3)$coefficients["tmp_met:PGS002280",c(1,3,5)], NA) else c(tmp_mod3$message,rep(NA,17))); i <- i + 1
  main_inter_res_list[[i]] <- c(met, "dementia", surv_dict$subtype[1], "inter_cox", "m4", "PGS000334", if(length(tmp_mod4) > 5) c(summary(tmp_mod4)$n, nobs(tmp_mod4), summary(tmp_mod4)$coefficients["tmp_met",c(1,3,5)], NA, NA, NA, summary(tmp_mod4)$coefficients["PGS000334",c(1,3,5)], NA, NA, NA, summary(tmp_mod4)$coefficients["tmp_met:PGS000334",c(1,3,5)], NA) else c(tmp_mod4$message,rep(NA,17))); i <- i + 1
  
  main_inter_res_list[[i]] <- c(met, "dementia", surv_dict$subtype[1], "inter_cox", "lrt_full_apoe4", "APOE4_2cat_PGS002280", if(length(tmp_mod_lrt_full) > 5) c(summary(tmp_mod_lrt_full)$n, nobs(tmp_mod_lrt_full), summary(tmp_mod_lrt_full)$coefficients["tmp_met",c(1,3,5)], summary(tmp_mod_lrt_full)$coefficients["APOE4_ncopy_1",c(1,3,5)], summary(tmp_mod_lrt_full)$coefficients["APOE4_ncopy_2",c(1,3,5)], summary(tmp_mod_lrt_full)$coefficients["tmp_met:APOE4_ncopy_1",c(1,3,5)], summary(tmp_mod_lrt_full)$coefficients["tmp_met:APOE4_ncopy_2",c(1,3,5)], ifelse(length(tmp_mod_lrt_nested) > 5, lrtest(tmp_mod_lrt_full, tmp_mod_lrt_nested)["Pr(>Chisq)"][2,], NA)) else c(tmp_mod_lrt_full$message,rep(NA,17))); i <- i + 1
  main_inter_res_list[[i]] <- c(met, "dementia", surv_dict$subtype[1], "inter_cox", "lrt_full_prs", "APOE4_2cat_PGS002280", if(length(tmp_mod_lrt_full) > 5) c(summary(tmp_mod_lrt_full)$n, nobs(tmp_mod_lrt_full), summary(tmp_mod_lrt_full)$coefficients["tmp_met",c(1,3,5)], NA, NA, NA, summary(tmp_mod_lrt_full)$coefficients["PGS002280",c(1,3,5)], NA, NA, NA, summary(tmp_mod_lrt_full)$coefficients["tmp_met:PGS002280",c(1,3,5)], NA) else c(tmp_mod_lrt_full$message,rep(NA,17))); i <- i + 1
  
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

main_inter_res_dat <- merge(main_inter_res_dat, fdata_filtered, by="hmdb_id") 

write.table(main_inter_res_dat, "results/met_dem_cox_baseline_interaction_nhs_2023.txt", row.names=F, quote=F, sep="\t")

##-------------------------------
## Interaction for ind variants
##-------------------------------

model_data <- merged_data

ind_var_inter_res_list <- list();  i <- 1

for (met in met_selected){
  
  print(met)
  tmp_met <- model_data[[met]]
  
  for (var in ad_dict$name){
    
    tmp_var <- model_data[[var]]
    tmp_eaf <- mean(model_data[[var]], na.rm=T) / 2
    
    tmp_mod <- tryCatch(coxph(reformulate(c("tmp_met * tmp_var + tmp_met * APOE4_ncopy_1 + tmp_met * APOE4_ncopy_2", covar_dem_list[["APOE4_2cat"]]), response=response_subtype), 
                                ties="breslow", data = model_data),error=function(e) e,warning=function(w) w)

    ind_var_inter_res_list[[i]] <- c(met, var, tmp_eaf, "dementia", surv_dict$subtype[1], "inter_cox_ind_var", "m4", "met_var_adjust_APOE4_inter", if(length(tmp_mod) > 5) c(summary(tmp_mod)$n, nobs(tmp_mod), summary(tmp_mod)$coefficients["tmp_met",c(1,3,5)], summary(tmp_mod)$coefficients["tmp_var",c(1,3,5)], summary(tmp_mod)$coefficients["tmp_met:tmp_var",c(1,3,5)]) else c(tmp_mod$message,rep(NA,10))); i <- i + 1
    
  }
}

ind_var_inter_res_dat <- as.data.frame(do.call("rbind", ind_var_inter_res_list))
colnames(ind_var_inter_res_dat) <- c("hmdb_id","variant","EAF","outcome","subtype","model_type","model","covariates","n","n_events","beta_met","se_met","p_met","beta_var","se_var","p_var","beta_inter_met_var","se_inter_met_var","p_inter_met_var")

table(is.na(ind_var_inter_res_dat$beta_met))
table(is.na(ind_var_inter_res_dat$beta_var))

ind_var_inter_res_dat <- merge(ind_var_inter_res_dat, fdata_filtered, by="hmdb_id") 
ind_var_inter_res_dat <- merge(ind_var_inter_res_dat, ad_dict, by.x="variant", by.y="name") 

write.table(ind_var_inter_res_dat, "results/met_dem_cox_baseline_interaction_individual_variants_nhs_2023.txt", row.names=F, quote=F, sep="\t")

##-------------------------------
## Stratified analysis by genetic
##-------------------------------

model_data <- merged_data

strata_res_list <- list();  i <- 1

for (j in 1){
  
  print(surv_dict$subtype[1])
  response_subtype <- surv_dict$response_subtype[1]

  for (met in met_selected){
    
    print(met)
    
    for(subgroup in c("PGS000334_t1","PGS000334_t2","PGS000334_t3",
                      "PGS002280_t1","PGS002280_t2","PGS002280_t3",
                      "APOE4_carrier","APOE4_noncarrier",
                      "APOE4_ncopy_1","APOE4_ncopy_2")){
      
      tmp_data <- model_data %>% filter(!!sym(subgroup)==1); tmp_met <- tmp_data[[met]]
      
      if (subgroup %in% c("PGS000334_t1","PGS000334_t2","PGS000334_t3","PGS002280_t1","PGS002280_t2","PGS002280_t3")){
        tmp_mod <- tryCatch(coxph(reformulate(c("tmp_met",covar_dem_list[["full_covar"]]), response=response_subtype), ties="breslow", data = tmp_data),error=function(e) e,warning=function(w) w)
      } else if (subgroup %in% c("APOE4_carrier","APOE4_noncarrier","APOE4_ncopy_1")){
        tmp_mod <- tryCatch(coxph(reformulate(c("tmp_met",covar_dem_list[["PGS002280"]]), response=response_subtype), ties="breslow", data = tmp_data),error=function(e) e,warning=function(w) w)
      } else if (subgroup %in% c("APOE4_ncopy_2")){
        tmp_mod <- tryCatch(coxph(reformulate(c("tmp_met + blddate + agemo + bmi + actcon + nSES + alcocon + AMED_avg + daykcal + PGS002280 + PC1 + PC2 + PC3 + PC4"), response=response_subtype), ties="breslow", data = tmp_data),error=function(e) e,warning=function(w) w)
      }
   
      strata_res_list[[i]] <- c(met, "dementia", surv_dict$subtype[1], "strata_cox", subgroup, if(length(tmp_mod) > 5) c(summary(tmp_mod)$n, nobs(tmp_mod), summary(tmp_mod)$coefficients["tmp_met",c(1,3,5)]) else c(tmp_mod$message,NA,NA,NA,NA)); i <- i + 1
      
    }
  }
}

strata_res_dat <- as.data.frame(do.call("rbind", strata_res_list))
colnames(strata_res_dat) <- c("hmdb_id", "outcome", "subtype","model_type","model","n","n_events","beta","se","p")

strata_res_dat <- strata_res_dat %>%
  group_by(subtype, model) %>%
  mutate(p_fdr = p.adjust(p, "BH"))

table(is.na(strata_res_dat$beta))

strata_res_dat <- merge(strata_res_dat, fdata_filtered, by="hmdb_id") 

write.table(strata_res_dat, "results/met_dem_cox_baseline_stratified_nhs_2023.txt", row.names=F, quote=F, sep="\t")
