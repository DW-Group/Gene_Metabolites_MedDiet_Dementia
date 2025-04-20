library(data.table)
library(dplyr)
library(survival)
library(lmtest)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

load("data/merged_data_baseline_nhs_2023.Rdata")

##-----------------------------------------
## Model
##-----------------------------------------

##-------------------------------
## Gene-met association
##-------------------------------

model_data <- merged_data

gene_met_res_list <- list();  i <- 1

### APOE4

for (met in met_selected){
  
  print(met)
  
  tmp_met <- model_data[[met]]

  tmp_mod1 <- glm(tmp_met ~ APOE4_ncopy_1 + APOE4_ncopy_2 + PC1 + PC2 + PC3 + PC4 + platformgsa + platformhuco2 + platformillu + platformomni + platformonco, data = model_data)
  tmp_mod2 <- glm(tmp_met ~ APOE4_ncopy_1 + APOE4_ncopy_2 + agemo + PC1 + PC2 + PC3 + PC4 + platformgsa + platformhuco2 + platformillu + platformomni + platformonco, data = model_data)
  
  gene_met_res_list[[i]] <- c("no_age", met, "APOE4_ncopy_1", sum(model_data[!is.na(tmp_met),]$APOE4_ncopy_1), nobs(tmp_mod1), if ("APOE4_ncopy_1" %in% names(summary(tmp_mod1)$coefficients[,1])) summary(tmp_mod1)$coefficients["APOE4_ncopy_1",c(1,2,4)] else c(NA,NA,NA)); i <- i + 1
  gene_met_res_list[[i]] <- c("adj_age", met, "APOE4_ncopy_1", sum(model_data[!is.na(tmp_met),]$APOE4_ncopy_1), nobs(tmp_mod2), if ("APOE4_ncopy_1" %in% names(summary(tmp_mod2)$coefficients[,1])) summary(tmp_mod2)$coefficients["APOE4_ncopy_1",c(1,2,4)] else c(NA,NA,NA)); i <- i + 1
  gene_met_res_list[[i]] <- c("no_age", met, "APOE4_ncopy_2", sum(model_data[!is.na(tmp_met),]$APOE4_ncopy_2), nobs(tmp_mod1), if ("APOE4_ncopy_2" %in% names(summary(tmp_mod1)$coefficients[,1])) summary(tmp_mod1)$coefficients["APOE4_ncopy_2",c(1,2,4)] else c(NA,NA,NA)); i <- i + 1
  gene_met_res_list[[i]] <- c("adj_age", met, "APOE4_ncopy_2", sum(model_data[!is.na(tmp_met),]$APOE4_ncopy_2), nobs(tmp_mod2), if ("APOE4_ncopy_2" %in% names(summary(tmp_mod2)$coefficients[,1])) summary(tmp_mod2)$coefficients["APOE4_ncopy_2",c(1,2,4)] else c(NA,NA,NA)); i <- i + 1

}

### Combine

gene_met_res_dat <- as.data.frame(do.call("rbind", gene_met_res_list))
colnames(gene_met_res_dat) <- c("covar","hmdb_id","gene","EAF_or_n_var","n_model","beta_var","se_var","p_var")
gene_met_res_dat <- merge(gene_met_res_dat, ad_dict, by.x="gene", by.y="name",all.x=T) 
gene_met_res_dat <- merge(gene_met_res_dat, fdata_filtered, by="hmdb_id",all.x=T)  %>%
  mutate(across(c(n_model,beta_var,se_var,p_var), as.numeric))

write.table(gene_met_res_dat, "results/gene_met_nhs_2023.txt", row.names=F, quote=F, sep="\t")
