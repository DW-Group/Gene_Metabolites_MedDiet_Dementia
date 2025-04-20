library(data.table)
library(dplyr)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

load("data/merged_data_baseline_nhs_2023.Rdata")

##-----------------------------------------
## Data processing
##-----------------------------------------

### Filter on individuals

model_data <- merged_data %>%
  filter(!is.na(av_ntotal))

### Create binary outcomes

model_data <- model_data %>%
  mutate(av_ntotal_t = ntile(av_ntotal, 3),
         av_nverbl_t = ntile(av_nverbl, 3),
         av_zscre_t = ntile(av_zscre, 3)) 

table(model_data$av_ntotal_t);table(model_data$av_ntotal_q)
table(model_data$av_nverbl_t);table(model_data$av_nverbl_q)
table(model_data$av_zscre_t);table(model_data$av_zscre_q)

model_data <- model_data %>%
  mutate(av_ntotal_tb = case_when(av_ntotal_t == 1 ~ 0,
                                      av_ntotal_t == 2 ~ NA,
                                      av_ntotal_t == 3 ~ 1),
         av_nverbl_tb = case_when(av_nverbl_t == 1 ~ 0,
                                      av_nverbl_t == 2 ~ NA,
                                      av_nverbl_t == 3 ~ 1),
         av_zscre_tb = case_when(av_zscre_t == 1 ~ 0,
                                      av_zscre_t == 2 ~ NA,
                                      av_zscre_t == 3 ~ 1))

table(model_data$av_ntotal_tb);table(model_data$av_ntotal_qb)
table(model_data$av_nverbl_tb);table(model_data$av_nverbl_qb)
table(model_data$av_zscre_tb);table(model_data$av_zscre_qb)

##-----------------------------------------
## Predictors
##-----------------------------------------

source("scripts/source_predictiors_list.R")

##-----------------------------------------
## Final data
##-----------------------------------------

tmp <- model_data %>% 
  subset(select=c("id", "study_ID",
                  binary_outcome_list,
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
## Model data
##-----------------------------------------

### Outcome

binary.y <- list()
binary.y.cc <- list()

for (subset in names(model_data_list)){
  
  binary.y[[subset]] <- list()
  binary.y.cc[[subset]] <- list()

  for (outcome in binary_outcome_list){
    
    binary.y[[subset]][[outcome]] <- list()
    binary.y.cc[[subset]][[outcome]] <- list()
    
    for (dat in names(model_data_list[[subset]])){
      
      binary.y[[subset]][[outcome]][[dat]] <- as.matrix(subset(model_data_list[[subset]][[dat]], select=c(outcome)))
      binary.y.cc[[subset]][[outcome]][[dat]] <- as.matrix(na.omit(subset(model_data_list[[subset]][[dat]], select=c(outcome))))
      
    }
  }
}

### Predictors

source("scripts/source_predictors_final.R")
load("results/pred_met_candidate.Rdata")

### Formula and model matrix list

binary_formula_list <- list()
binary_mod.matrix_list <- list()

for (subset in names(model_data_list)){
  
  binary_formula_list[[subset]] <- list()
  binary_mod.matrix_list[[subset]] <- list()
  
  if (subset=="all"){
    tmp_covar_list <- covar_all_list[["with_APOE4"]]
    tmp_covar_list[["pred_met_all"]] <- c(lifestyle_APOE4_PRS_list, met_pred_list[["all"]][["all"]])
  } else{
    tmp_covar_list <- covar_all_list[["no_APOE4"]]
    tmp_covar_list[["pred_met_all"]] <- c(lifestyle_PRS_list, met_pred_list[["all"]][["all"]])
  }
  
  for (outcome in binary_outcome_list){
    
    binary_formula_list[[subset]][[outcome]] <- list()
    binary_mod.matrix_list[[subset]][[outcome]] <- list()
    
    for (covar in names(tmp_covar_list)){
      
      binary_formula_list[[subset]][[outcome]][[covar]] <- reformulate(tmp_covar_list[[covar]], response=outcome)
      binary_mod.matrix_list[[subset]][[outcome]][[covar]] <- list()
      
      for (dat in names(model_data_list[[subset]])){
        
        binary_mod.matrix_list[[subset]][[outcome]][[covar]][[dat]] <- model.matrix(binary_formula_list[[subset]][[outcome]][[covar]], model_data_list[[subset]][[dat]])[,-1]
        
      }
    }
  }
}

save.image("data/prediction_data_tics.RData")
