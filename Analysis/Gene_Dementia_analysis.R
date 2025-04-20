library(data.table)
library(dplyr)
library(survival)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

load("data/merged_data_baseline_nhs_2023.Rdata")
load("data/genetic/ad_variants_dictionary.RData")

##-----------------------------------------
## Model
##-----------------------------------------

### Covariates

gene_covar <- c("PC1","PC2","PC3","PC4","platformgsa","platformhuco2","platformillu","platformomni","platformonco")
gene_covar_reduced <- c("PC1","PC2","PC3","PC4","platformhuco2","platformillu","platformomni","platformonco")

##-------------------------------
## Gene-dementia association
##-------------------------------

model_data <- merged_data

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

  tmp_mod1 <- glm(reformulate(c("PGS000334", gene_covar), response=status_subtype), family=binomial(), data = model_data)
  tmp_mod2 <- glm(reformulate(c("PGS000334 + agemo", gene_covar), response=status_subtype), family=binomial(), data = model_data)
  tmp_mod3 <- coxph(reformulate(c("PGS000334", gene_covar), response=response_subtype), ties="breslow", data = model_data)
  tmp_mod4 <- coxph(reformulate(c("PGS000334 + agemo", gene_covar), response=response_subtype), ties="breslow", data = model_data)
  
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "logistic", "PGS000334", NA, nobs(tmp_mod1), sum(model_data[[status_subtype]]), summary(tmp_mod1)$coefficients["PGS000334",c(1,2,4)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "logistic_adj_age", "PGS000334", NA, nobs(tmp_mod2), sum(model_data[[status_subtype]]), summary(tmp_mod2)$coefficients["PGS000334",c(1,2,4)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "cox", "PGS000334", NA, summary(tmp_mod3)$n, nobs(tmp_mod3), summary(tmp_mod3)$coefficients["PGS000334",c(1,3,5)]); i <- i + 1
  gene_dementia_res_list[[i]] <- c("dementia", surv_dict$subtype[j], "cox_adj_age", "PGS000334", NA, summary(tmp_mod4)$n, nobs(tmp_mod4), summary(tmp_mod4)$coefficients["PGS000334",c(1,3,5)]); i <- i + 1

}

### Combine

gene_dementia_res_dat <- as.data.frame(do.call("rbind", gene_dementia_res_list))
colnames(gene_dementia_res_dat) <- c("outcome","subtype","model","gene","EAF_or_n","n","n_events","beta_var","se_var","p_var")
gene_dementia_res_dat <- merge(gene_dementia_res_dat, ad_dict, by.x="gene", by.y="name",all.x=T) %>%
  mutate(across(c(n,n_events,beta_var,se_var,p_var), as.numeric))

write.table(gene_dementia_res_dat, "results/gene_dementia_nhs.txt", row.names=F, quote=F, sep="\t")

##-------------------------------
## Gene-TICS association
##-------------------------------

model_data <- merged_data %>% filter(!is.na(av_ntotal))
tics_list <- c("av_ntotal","av_nverbl","av_zscre")

gene_tics_res_list <- list();  i <- 1

### APOE4

for (tics in tics_list){

  tmp_mod1 <- glm(reformulate(c("APOE4_carrier", gene_covar_reduced), response=tics), data = model_data)
  tmp_mod2 <- glm(reformulate(c("APOE4_carrier + agemo", gene_covar_reduced), response=tics), data = model_data)
  
  gene_tics_res_list[[i]] <- c(tics, "linear", "APOE4_carrier", sum(model_data$APOE4_carrier), nobs(tmp_mod1), summary(tmp_mod1)$coefficients["APOE4_carrier",c(1,2,4)]); i <- i + 1
  gene_tics_res_list[[i]] <- c(tics, "linear_adj_age", "APOE4_carrier", sum(model_data$APOE4_carrier), nobs(tmp_mod2), summary(tmp_mod2)$coefficients["APOE4_carrier",c(1,2,4)]); i <- i + 1

}

### PRS

for (tics in tics_list){

  tmp_mod1 <- glm(reformulate(c("PGS002280", gene_covar_reduced), response=tics), data = model_data)
  tmp_mod2 <- glm(reformulate(c("PGS002280 + agemo", gene_covar_reduced), response=tics), data = model_data)
  
  gene_tics_res_list[[i]] <- c(tics, "linear", "PGS002280", NA, nobs(tmp_mod1), summary(tmp_mod1)$coefficients["PGS002280",c(1,2,4)]); i <- i + 1
  gene_tics_res_list[[i]] <- c(tics, "linear_adj_age", "PGS002280", NA, nobs(tmp_mod2), summary(tmp_mod2)$coefficients["PGS002280",c(1,2,4)]); i <- i + 1
  
  tmp_mod1 <- glm(reformulate(c("PGS000334", gene_covar_reduced), response=tics), data = model_data)
  tmp_mod2 <- glm(reformulate(c("PGS000334 + agemo", gene_covar_reduced), response=tics), data = model_data)
  
  gene_tics_res_list[[i]] <- c(tics, "linear", "PGS000334", NA, nobs(tmp_mod1), summary(tmp_mod1)$coefficients["PGS000334",c(1,2,4)]); i <- i + 1
  gene_tics_res_list[[i]] <- c(tics, "linear_adj_age", "PGS000334", NA, nobs(tmp_mod2), summary(tmp_mod2)$coefficients["PGS000334",c(1,2,4)]); i <- i + 1

}

### Combine

gene_tics_res_dat <- as.data.frame(do.call("rbind", gene_tics_res_list))
colnames(gene_tics_res_dat) <- c("outcome","model","gene","EAF_or_n","n","beta_var","se_var","p_var")
gene_tics_res_dat <- merge(gene_tics_res_dat, ad_dict, by.x="gene", by.y="name",all.x=T) %>%
  mutate(across(c(n,beta_var,se_var,p_var), as.numeric))

write.table(gene_tics_res_dat, "results/gene_tics_nhs.txt", row.names=F, quote=F, sep="\t")
