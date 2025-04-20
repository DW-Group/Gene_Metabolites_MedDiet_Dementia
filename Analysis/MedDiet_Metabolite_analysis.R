library(data.table)
library(dplyr)
library(tidyr)
library(lmtest)
library(writexl)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

load("data/merged_data_baseline_nhs_2023.Rdata")
load("results/met_dem_cox_baseline_nhs_2023_filtered_res.RData")
load("data/genetic/ad_variants_dictionary.RData")

##-----------------------------------------
## Data preparation 
##-----------------------------------------

### Outcome

met_list <- fdata_filtered$hmdb_id

### Exposure

summary(merged_data$AMED_avg)

diet_list <- c("AMED_avg")

### Covariates

source("scripts/source_covar_list_MedDiet_Metabolite_analysis.R")

### Model data

model_data_list <- list("all" = merged_data,
                        "APOE4_carrier" = merged_data %>% filter(APOE4_carrier == 1),
                        "APOE4_noncarrier" = merged_data %>% filter(APOE4_carrier == 0),
                        "APOE4_ncopy_1" = merged_data %>% filter(APOE4_ncopy_1 == 1),
                        "APOE4_ncopy_2" = merged_data %>% filter(APOE4_ncopy_2 == 1),
                        "PGS002280_t1" = merged_data %>% filter(PGS002280_t1 == 1),
                        "PGS002280_t2" = merged_data %>% filter(PGS002280_t2 == 1),
                        "PGS002280_t3" = merged_data %>% filter(PGS002280_t3 == 1),
                        "PGS000334_t1" = merged_data %>% filter(PGS000334_t1 == 1),
                        "PGS000334_t2" = merged_data %>% filter(PGS000334_t2 == 1),
                        "PGS000334_t3" = merged_data %>% filter(PGS000334_t3 == 1))

##-----------------------------------------
## Main and subgroup analysis
##-----------------------------------------

main_subgroup_res_list <- list(); i <- 1

for (diet in diet_list){
  
  print(diet)

  for (met in met_list){
    
    print(met)
    
    for (covar in names(covar_list)){
      
      for (subset in names(model_data_list)){
        
        if (!(grepl("apoe4",covar) & grepl("APOE4",subset)) &
            !(grepl("prs",covar) & grepl("PGS",subset)) &
            !(grepl("prs|apoe4",covar) & grepl("PGS000334",subset))){
          
          tmp_mod <- tryCatch(glm(reformulate(c(diet, covar_list[[covar]]), response=met), data=model_data_list[[subset]]),error=function(e) e,warning=function(w) w)
          
          if (length(tmp_mod) < 10){
            main_subgroup_res_list[[i]] <- c(diet, met, subset, covar, rep(NA,9), tmp_mod$message); i <- i + 1
          } else{
            main_subgroup_res_list[[i]] <- c(diet, met, subset, covar, nobs(tmp_mod),
                                             summary(tmp_mod)$coefficients[diet,],
                                             summary(tmp_mod)$coefficients["agemo",], NA); i <- i + 1
          }
        }
      }
    }
  }
}

main_subgroup_res <- as.data.frame(do.call("rbind",main_subgroup_res_list))
main_subgroup_res[main_subgroup_res=="NaN"] <- NA
colnames(main_subgroup_res) <- c("diet","hmdb_id","subgroup","covar","n","diet_beta","diet_se","diet_t","diet_p","age_beta","age_se","age_t","age_p","message")
main_subgroup_res <- main_subgroup_res %>%
  mutate(across(c(n,diet_beta,diet_se,diet_t,diet_p,age_beta,age_se,age_t,age_p), as.numeric))
main_subgroup_res <- merge(main_subgroup_res,fdata_filtered,by="hmdb_id")

table(main_subgroup_res$diet_p < 0.05)

write_xlsx(main_subgroup_res, "results/diet_gene_met_main_subgroup.xlsx")

##-----------------------------------------
## Interaction analysis
##-----------------------------------------

### APOE and PRS

inter_res_list <- list(); i <- 1

for (diet in diet_list){
  
  print (diet)

  for (met in met_list){
    
    print(met)
      
    tmp_mod1 <- tryCatch(glm(reformulate(c(paste0(diet," * APOE4_carrier"), gene_covar_list), response=met), data=merged_data),error=function(e) e,warning=function(w) w)
    tmp_mod2 <- tryCatch(glm(reformulate(c(paste0(diet," * APOE4_ncopy_1"), paste0(diet," * APOE4_ncopy_2"), gene_covar_list), response=met), data=merged_data),error=function(e) e,warning=function(w) w)
    tmp_mod3 <- tryCatch(glm(reformulate(c(diet,"APOE4_ncopy_1","APOE4_ncopy_2", gene_covar_list), response=met), data=merged_data),error=function(e) e,warning=function(w) w)
    tmp_mod4 <- tryCatch(glm(reformulate(c(paste0(diet," * PGS002280"), covar_list[["apoe4"]]), response=met), data=merged_data),error=function(e) e,warning=function(w) w)
    tmp_mod5 <- tryCatch(glm(reformulate(c(paste0(diet," * PGS000334"), gene_covar_list), response=met), data=merged_data),error=function(e) e,warning=function(w) w)
    
    if (length(tmp_mod1) < 10){
      inter_res_list[[i]] <- c(diet, met, "APOE4_carrier", rep(NA,22), tmp_mod1$message); i <- i + 1
    } else{
      inter_res_list[[i]] <- c(diet, met, "APOE4_carrier", nobs(tmp_mod1),
                               summary(tmp_mod1)$coefficients[diet,],
                               summary(tmp_mod1)$coefficients["APOE4_carrier",],
                               rep(NA,4),
                               summary(tmp_mod1)$coefficients[paste0(diet,":APOE4_carrier"),],
                               rep(NA,4),
                               rep(NA,1),
                               rep(NA,1)); i <- i + 1
    }
    
    if (length(tmp_mod2) < 10){
      inter_res_list[[i]] <- c(diet, met, "APOE4_2cat", rep(NA,22), tmp_mod2$message); i <- i + 1
    } else{
      inter_res_list[[i]] <- c(diet, met, "APOE4_2cat", nobs(tmp_mod2),
                               summary(tmp_mod2)$coefficients[diet,],
                               summary(tmp_mod2)$coefficients["APOE4_ncopy_1",],
                               summary(tmp_mod2)$coefficients["APOE4_ncopy_2",],
                               summary(tmp_mod2)$coefficients[paste0(diet,":APOE4_ncopy_1"),],
                               summary(tmp_mod2)$coefficients[paste0(diet,":APOE4_ncopy_2"),],
                               lrtest(tmp_mod2, tmp_mod3)["Pr(>Chisq)"][2,],
                               rep(NA,1)); i <- i + 1
    }
    if (length(tmp_mod4) < 10){
      inter_res_list[[i]] <- c(diet, met, "PGS002280", rep(NA,22), tmp_mod4$message); i <- i + 1
    } else{
      inter_res_list[[i]] <- c(diet, met, "PGS002280", nobs(tmp_mod4),
                               summary(tmp_mod4)$coefficients[diet,],
                               summary(tmp_mod4)$coefficients["PGS002280",],
                               rep(NA,4),
                               summary(tmp_mod4)$coefficients[paste0(diet,":PGS002280"),],
                               rep(NA,4),
                               rep(NA,1),
                               rep(NA,1)); i <- i + 1
    }
    if (length(tmp_mod5) < 10){
      inter_res_list[[i]] <- c(diet, met, "PGS000334", rep(NA,22), tmp_mod5$message); i <- i + 1
    } else{
      inter_res_list[[i]] <- c(diet, met, "PGS000334", nobs(tmp_mod5),
                               summary(tmp_mod5)$coefficients[diet,],
                               summary(tmp_mod5)$coefficients["PGS000334",],
                               rep(NA,4),
                               summary(tmp_mod5)$coefficients[paste0(diet,":PGS000334"),],
                               rep(NA,4),
                               rep(NA,1),
                               rep(NA,1)); i <- i + 1
    }
  }
}

inter_res <- as.data.frame(do.call("rbind",inter_res_list))
inter_res[inter_res=="NaN"] <- NA
colnames(inter_res) <- c("diet","hmdb_id","var","n",
                         "diet_beta","diet_se","diet_t","diet_p",
                         "var1_beta","var1_se","var1_t","var1_p",
                         "var2_beta","var2_se","var2_t","var2_p",
                         "inter1_beta","inter1_se","inter1_t","inter1_p",
                         "inter2_beta","inter2_se","inter2_t","inter2_p",
                         "lrt_p",
                         "message")
inter_res <- inter_res %>%
  mutate(across(c(n,
                  diet_beta,diet_se,diet_t,diet_p,
                  var1_beta,var1_se,var1_t,var1_p,
                  var2_beta,var2_se,var2_t,var2_p,
                  inter1_beta,inter1_se,inter1_t,inter1_p,
                  inter2_beta,inter2_se,inter2_t,inter2_p,
                  lrt_p), as.numeric))
inter_res <- merge(inter_res, fdata_filtered, by="hmdb_id", all.x=T)

write_xlsx(inter_res, "results/diet_gene_met_inter.xlsx")

### Individual variants

inter_ind_var_res_list <- list(); i <- 1

for (diet in diet_list){
  
  print (diet)
  
  for (met in met_list){
    
    print(met)

    for (var in ad_dict$name){
  
      tmp_mod <- tryCatch(glm(reformulate(c(paste0(diet," * `",var,"`"), covar_list[["apoe4"]]), response=met), data=merged_data),error=function(e) e,warning=function(w) w)

      if (length(tmp_mod) < 10){
        inter_ind_var_res_list[[i]] <- c(diet, met, var, rep(NA,13), tmp_mod$message); i <- i + 1
        } else{
        inter_ind_var_res_list[[i]] <- c(diet, met, var, nobs(tmp_mod),
                                         summary(tmp_mod)$coefficients[diet,],
                                         summary(tmp_mod)$coefficients[paste0("`",var,"`"),],
                                         summary(tmp_mod)$coefficients[paste0(diet,":`",var,"`"),],
                                         NA); i <- i + 1

      }
    }
  }
}

inter_ind_var_res <- as.data.frame(do.call("rbind",inter_ind_var_res_list))
inter_ind_var_res[inter_ind_var_res=="NaN"] <- NA
colnames(inter_ind_var_res) <- c("diet","hmdb_id","var","n",
                                 "diet_beta","diet_se","diet_t","diet_p",
                                 "var_beta","var_se","var_t","var_p",
                                 "inter_beta","inter_se","inter_t","inter_p","message")
inter_ind_var_res <- inter_ind_var_res %>%
  mutate(across(c(n,
                  diet_beta,diet_se,diet_t,diet_p,
                  var_beta,var_se,var_t,var_p,
                  inter_beta,inter_se,inter_t,inter_p), as.numeric))
inter_ind_var_res <- merge(inter_ind_var_res, fdata_filtered, by="hmdb_id")
inter_ind_var_res <- merge(inter_ind_var_res, ad_dict, by.x="var", by.y="name")

write_xlsx(inter_ind_var_res, "results/diet_gene_met_inter_ind_var.xlsx")
