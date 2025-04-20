library(data.table)
library(dplyr)
library(tidyr)
library(lmtest)
library(writexl)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

load("data/merged_data_baseline_hpfs_2023.Rdata")

##-----------------------------------------
## Data preparation 
##-----------------------------------------

### Outcome

met_list <- fdata_filtered_hpfs$hmdb_id

### Exposure

summary(merged_data_hpfs$amedcon)

diet_list <- c("amedcon")

### Covariates

source("scripts/source_covar_list_MedDiet_Metabolite_analysis_replication.R")

### Model data

model_data_list <- list("all" = merged_data_hpfs,
                        "APOE4_carrier" = merged_data_hpfs %>% filter(APOE4_carrier == 1),
                        "APOE4_noncarrier" = merged_data_hpfs %>% filter(APOE4_carrier == 0))

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

main_subgroup_res <- as.data.frame(do.call("rbind",main_subgroup_res_list))
main_subgroup_res[main_subgroup_res=="NaN"] <- NA
colnames(main_subgroup_res) <- c("diet","hmdb_id","subgroup","covar","n","diet_beta","diet_se","diet_t","diet_p","age_beta","age_se","age_t","age_p","message")
main_subgroup_res <- main_subgroup_res %>%
  mutate(across(c(n,diet_beta,diet_se,diet_t,diet_p,age_beta,age_se,age_t,age_p), as.numeric))
main_subgroup_res <- merge(main_subgroup_res,fdata_filtered_hpfs,by="hmdb_id")

table(main_subgroup_res$diet_p < 0.05)

write_xlsx(main_subgroup_res, "results/diet_gene_met_main_subgroup_hpfs.xlsx")