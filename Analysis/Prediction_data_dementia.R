library(data.table)
library(dplyr)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

load("data/merged_data_baseline_nhs_2023.Rdata")
surv_dict <- surv_dict %>% filter(subtype %in% c("all","yr_15"))

##-----------------------------------------
## Predictors
##-----------------------------------------

source("scripts/source_predictiors_list.R")

##-----------------------------------------
## Final data
##-----------------------------------------

tmp <- merged_data %>% 
  subset(select=c("id", "study_ID", surv_dict$time, surv_dict$status, 
                  lifestyle_short_list,
                  apoe4_2cat_list,
                  prs_list))
model_data <- merge(tmp, rf_imputed_data, by="study_ID")

model_data_final <- model_data

##-----------------------------------------
## Split training and test sets
##-----------------------------------------

### Subsets

model_data_final_list <- list("all" = model_data_final,
                              "APOE4_carrier" = model_data_final %>% filter(APOE4_carrier == 1),
                              "APOE4_noncarrier" = model_data_final %>% filter(APOE4_carrier == 0))

### Split

model_data_list <- list()

for (subset in names(model_data_final_list)){
  
  set.seed(123)
  sample_ind <- sample(c(TRUE,FALSE), nrow(model_data_final_list[[subset]]),  
                       replace=TRUE, prob=c(0.6,0.4)) 
  model_data_list[[subset]] <- list("train" = model_data_final_list[[subset]][sample_ind, ],
                                    "test" = model_data_final_list[[subset]][!sample_ind, ])

}

### Standardize covariates (continuous only)

for (subset in names(model_data_final_list)){
  for (var in c("agemo","AMED_avg")){
    tmp_mean <- mean(model_data_list[[subset]][["train"]][[var]])
    tmp_sd <- sd(model_data_list[[subset]][["train"]][[var]])
    model_data_list[[subset]][["train"]][[var]] <- (model_data_list[[subset]][["train"]][[var]] - tmp_mean) / tmp_sd
    model_data_list[[subset]][["test"]][[var]] <- (model_data_list[[subset]][["test"]][[var]] - tmp_mean) / tmp_sd
  }
}

##-----------------------------------------
## Candidate metabolites
##-----------------------------------------

load("results/pred_met_candidate.Rdata")

##-----------------------------------------
## Model data
##-----------------------------------------

### Outcome

surv.y <- list()

for (subset in names(model_data_list)){
  
  surv.y[[subset]] <- list()

  for (j in 1:dim(surv_dict)[1]){
    
    surv.y[[subset]][[surv_dict$subtype[j]]] <- list()
    
    for (dat in c("train","test")){
      
      surv.y[[subset]][[surv_dict$subtype[j]]][[dat]] <- as.matrix(subset(model_data_list[[subset]][[dat]], select=c(surv_dict$time[j],surv_dict$status[j])))
      colnames(surv.y[[subset]][[surv_dict$subtype[j]]][[dat]]) <- c("time","status")

    }
  }
}

### Predictors

source("scripts/source_predictors_final.R")

### Formula and model matrix list

surv_formula_list <- list()
surv_mod.matrix_list <- list()

for (subset in names(model_data_list)){
  
  surv_formula_list[[subset]] <- list()
  surv_mod.matrix_list[[subset]] <- list()
  
  for (j in 1:dim(surv_dict)[1]){
    
    surv_formula_list[[subset]][[surv_dict$subtype[j]]] <- list()
    surv_mod.matrix_list[[subset]][[surv_dict$subtype[j]]] <- list()
    
    if (subset=="all"){
      tmp_covar_list <- covar_all_list[["with_APOE4"]]
      tmp_covar_list[["pred_met_all"]] <- c(lifestyle_APOE4_PRS_list, met_pred_list[["all"]][[surv_dict$subtype[j]]])
    } else{
      tmp_covar_list <- covar_all_list[["no_APOE4"]]
      tmp_covar_list[["pred_met_all"]] <- c(lifestyle_PRS_list,met_pred_list[["all"]][[surv_dict$subtype[j]]])
    }
    
    for (covar in names(tmp_covar_list)){
      
      surv_formula_list[[subset]][[surv_dict$subtype[j]]][[covar]] <- reformulate(tmp_covar_list[[covar]], response=surv_dict$response_subtype[j])
      surv_mod.matrix_list[[subset]][[surv_dict$subtype[j]]][[covar]] <- list()
      
      for (dat in c("train","test")){
        
        surv_mod.matrix_list[[subset]][[surv_dict$subtype[j]]][[covar]][[dat]] <- model.matrix(surv_formula_list[[subset]][[surv_dict$subtype[j]]][[covar]], model_data_list[[subset]][[dat]])[,-1]
        
      }
    }
  }
}

save.image("data/prediction_data.RData")
