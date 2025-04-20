library(data.table)
library(dplyr)
library(survival)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

load("data/merged_data_baseline_hpfs_2023.Rdata")
load("data/genetic/ad_variants_dictionary.RData")

##-----------------------------------------
## Model
##-----------------------------------------

### Covariates

gene_covar <- c("PC1","PC2","PC3","PC4","platformhuco2","platformillu","platformomni","platformonco")

##-------------------------------
## Gene-dementia association
##-------------------------------

model_data <- merged_data_hpfs

gene_dementia_res_list <- list();  i <- 1

### APOE4

for (j in 1){
  
  print(surv_dict$subtype[j])
  
  response_subtype <- surv_dict$response_subtype[j]
  status_subtype <- surv_dict$status[j]

  tmp_mod1 <- glm(reformulate(c("APOE4_ncopy_1 + APOE4_ncopy_2", gene_covar), response=status_subtype), family=binomial(), data = model_data)
  tmp_mod2 <- glm(reformulate(c("APOE4_ncopy_1 + APOE4_ncopy_2 + agemo", gene_covar), response=status_subtype), family=binomial(), data = model_data)
  tmp_mod3 <- coxph(reformulate(c("APOE4_ncopy_1 + APOE4_ncopy_2", gene_covar), response=response_subtype), ties="breslow", data = model_data)
  tmp_mod4 <- coxph(reformulate(c("APOE4_ncopy_1 + APOE4_ncopy_2 + agemo", gene_covar), response=response_subtype), ties="breslow", data = model_data)
  
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "logistic", "APOE4_ncopy_1", sum(model_data$APOE4_ncopy_1), nobs(tmp_mod1), sum(model_data[[status_subtype]]), summary(tmp_mod1)$coefficients["APOE4_ncopy_1",c(1,2,4)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "logistic_adj_age", "APOE4_ncopy_1", sum(model_data$APOE4_ncopy_1), nobs(tmp_mod2), sum(model_data[[status_subtype]]), summary(tmp_mod2)$coefficients["APOE4_ncopy_1",c(1,2,4)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "cox", "APOE4_ncopy_1", sum(model_data$APOE4_ncopy_1), summary(tmp_mod3)$n, nobs(tmp_mod3), summary(tmp_mod3)$coefficients["APOE4_ncopy_1",c(1,3,5)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "cox_adj_age", "APOE4_ncopy_1", sum(model_data$APOE4_ncopy_1), summary(tmp_mod4)$n, nobs(tmp_mod4), summary(tmp_mod4)$coefficients["APOE4_ncopy_1",c(1,3,5)]); i <- i + 1
  
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "logistic", "APOE4_ncopy_2", sum(model_data$APOE4_ncopy_2), nobs(tmp_mod1), sum(model_data[[status_subtype]]), summary(tmp_mod1)$coefficients["APOE4_ncopy_2",c(1,2,4)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "logistic_adj_age", "APOE4_ncopy_2", sum(model_data$APOE4_ncopy_2), nobs(tmp_mod2), sum(model_data[[status_subtype]]), summary(tmp_mod2)$coefficients["APOE4_ncopy_2",c(1,2,4)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "cox", "APOE4_ncopy_2", sum(model_data$APOE4_ncopy_2), summary(tmp_mod3)$n, nobs(tmp_mod3), summary(tmp_mod3)$coefficients["APOE4_ncopy_2",c(1,3,5)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "cox_adj_age", "APOE4_ncopy_2", sum(model_data$APOE4_ncopy_2), summary(tmp_mod4)$n, nobs(tmp_mod4), summary(tmp_mod4)$coefficients["APOE4_ncopy_2",c(1,3,5)]); i <- i + 1
 
}

### PRS

for (j in 1){
  
  print(surv_dict$subtype[j])
  
  response_subtype <- surv_dict$response_subtype[j]
  status_subtype <- surv_dict$status[j]
  
  tmp_mod1 <- glm(reformulate(c("PGS002280", gene_covar), response=status_subtype), family=binomial(), data = model_data)
  tmp_mod2 <- glm(reformulate(c("PGS002280 + agemo", gene_covar), response=status_subtype), family=binomial(), data = model_data)
  tmp_mod3 <- coxph(reformulate(c("PGS002280", gene_covar), response=response_subtype), ties="breslow", data = model_data)
  tmp_mod4 <- coxph(reformulate(c("PGS002280 + agemo", gene_covar), response=response_subtype), ties="breslow", data = model_data)
  
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "logistic", "PGS002280", NA, nobs(tmp_mod1), sum(model_data[[status_subtype]]), summary(tmp_mod1)$coefficients["PGS002280",c(1,2,4)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "logistic_adj_age", "PGS002280", NA, nobs(tmp_mod2), sum(model_data[[status_subtype]]), summary(tmp_mod2)$coefficients["PGS002280",c(1,2,4)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "cox", "PGS002280", NA, summary(tmp_mod3)$n, nobs(tmp_mod3), summary(tmp_mod3)$coefficients["PGS002280",c(1,3,5)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "cox_adj_age", "PGS002280", NA, summary(tmp_mod4)$n, nobs(tmp_mod4), summary(tmp_mod4)$coefficients["PGS002280",c(1,3,5)]); i <- i + 1

}

### Combine

gene_dementia_res_dat <- as.data.frame(do.call("rbind", gene_dementia_res_list))
colnames(gene_dementia_res_dat) <- c("outcome","subtype","model","gene","EAF_or_n","n","n_events","beta_var","se_var","p_var")

write.table(gene_dementia_res_dat, "results/gene_dementia_hpfs.txt", row.names=F, quote=F, sep="\t")
