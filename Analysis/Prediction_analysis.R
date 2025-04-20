library(randomForest)
library(survival)
library(risksetROC)
library(fastshap)
library(dplyr)
library(tibble)
library(writexl)
library(tidyr)

##-----------------------------------------
## Time-dependent AUC and SHAP
##-----------------------------------------

rm(list=ls())

load("data/prediction_data.RData")

subset <- "all"
outcome_selected <- c("all","yr_15")

##------------------------
## Cox model
##------------------------

predictor_selected_list <- list()
predictor_selected_list[["all"]] <- c("pred_met_all",
                                      "lifestyle_APOE4_PRS",
                                      "lifestyle_APOE4")
predictor_selected_list[["APOE4"]] <- c("pred_met_all",
                                       "lifestyle_PRS")

cox_mod_list <- list()
cox_res_list <- list()
for(subset in c("all","APOE4_carrier","APOE4_noncarrier")){
  cox_mod_list[[subset]] <- list()
  cox_res_list[[subset]] <- list()
  for (outcome in outcome_selected){
    cox_mod_list[[subset]][[outcome]] <- list()
    cox_res_list[[subset]][[outcome]] <- list()
    if (grepl("all",subset)){
      predictor_selected <- c(predictor_selected_list[["all"]],"lifestyle_short")
    } else {
      predictor_selected <- c(predictor_selected_list[["APOE4"]],"lifestyle_short")
    }
    for (predictor in predictor_selected){
      cox_mod_list[[subset]][[outcome]][[predictor]] <- coxph(surv_formula_list[[subset]][[outcome]][[predictor]], 
                                                                  ties="breslow", data=model_data_list[[subset]][["train"]])
      cox_res_list[[subset]][[outcome]][[predictor]] <- list(concordance(cox_mod_list[[subset]][[outcome]][[predictor]], newdata=model_data_list[[subset]][["train"]]),
                                                                           concordance(cox_mod_list[[subset]][[outcome]][[predictor]], newdata=model_data_list[[subset]][["test"]]))
    }
  }
}

##------------------------
## Time-dependent AUC
##------------------------

cox_auc_list <- list()
for(subset in c("all","APOE4_carrier","APOE4_noncarrier")){
  cox_auc_list[[subset]] <- list()
  for (outcome in outcome_selected){
    cox_auc_list[[subset]][[outcome]] <- list()
    if (grepl("all",subset)){
      predictor_selected <- c(predictor_selected_list[["all"]],"lifestyle_short")
    } else {
      predictor_selected <- c(predictor_selected_list[["APOE4"]],"lifestyle_short")
    }
    for (predictor in predictor_selected){
      cox_auc_list[[subset]][[outcome]][[predictor]] <- 
        risksetAUC(Stime=model_data_list[[subset]][["test"]][[surv_dict[surv_dict$subtype==outcome,]$time]],
                   status=model_data_list[[subset]][["test"]][[surv_dict[surv_dict$subtype==outcome,]$status]],
                   marker=predict(cox_mod_list[[subset]][[outcome]][[predictor]], newdata=model_data_list[[subset]][["test"]], type="lp"),
                   tmax=max(model_data_list[[subset]][["test"]][[surv_dict[surv_dict$subtype==outcome,]$time]]),
                   plot=F)
    }
  }
}

model_name_dict <- data.frame(predictor = c("lifestyle_short",
                                            "lifestyle_APOE4",
                                            "lifestyle_PRS",
                                            "lifestyle_APOE4_PRS",
                                            "pred_met_all"),
                              anno_predictor = c("Baseline~model",
                                                 "+~italic('APOE4')",
                                                 "+~PRS",
                                                 "+~italic('APOE4')~+~PRS",
                                                 "+~italic('APOE4')~+~PRS~+~Metabolites"))

library(tidyr)
auc_dat_list <- list()
for(subset in c("all","APOE4_carrier","APOE4_noncarrier")){
  auc_dat_list[[subset]] <- list()
  if (grepl("all",subset)){outcome_selected <- c("all","yr_15")} else {outcome_selected <- c("all")}
  for (outcome in outcome_selected){
    if (grepl("all",subset)){
      predictor_selected <- c(predictor_selected_list[["all"]],"lifestyle_short")
    } else {
      predictor_selected <- c(predictor_selected_list[["APOE4"]],"lifestyle_short")
    }
    auc_dat_list[[subset]][[outcome]] <- data.frame(time = cox_auc_list[[subset]][[outcome]][["lifestyle_short"]]$utimes)
    for (predictor in predictor_selected){
      auc_dat_list[[subset]][[outcome]][predictor] <- cox_auc_list[[subset]][[outcome]][[predictor]]$AUC
    }
    auc_dat_list[[subset]][[outcome]] <- auc_dat_list[[subset]][[outcome]] %>% pivot_longer(-time)
    colnames(auc_dat_list[[subset]][[outcome]]) <- c("time","predictor","auc")
    auc_dat_list[[subset]][[outcome]] <- merge(auc_dat_list[[subset]][[outcome]], model_name_dict, by="predictor") %>%
      mutate(anno_predictor = factor(anno_predictor, levels=c("Baseline~model",
                                                              "+~italic('APOE4')",
                                                              "+~PRS",
                                                              "+~italic('APOE4')~+~PRS",
                                                              "+~italic('APOE4')~+~PRS~+~Metabolites")))
  }
}

save(auc_dat_list,file="results/prediction_dementia_auc_by_time.RData")

##------------------------
## SHAP
##------------------------

cox_mod_shap_all <- cox_mod_list[["all"]][["all"]][["pred_met_all"]]
cox_mod_shap_yr_15 <- cox_mod_list[["all"]][["yr_15"]][["pred_met_all"]]

predict_function <- function(object, newdata) {
  predict(object, newdata = newdata, type = "risk")
}

shap_values_all <- fastshap::explain(
  cox_mod_shap_all, 
  X = model_data_list[["all"]][["test"]][, names(cox_mod_shap_all$assign)],
  pred_wrapper = predict_function,
  nsim = 50
)

shap_values_yr_15 <- fastshap::explain(
  cox_mod_shap_yr_15, 
  X = model_data_list[["all"]][["test"]][, names(cox_mod_shap_yr_15$assign)],
  pred_wrapper = predict_function,
  nsim = 50
)

shap_values_all <- as.data.frame(shap_values_all)
shap_values_all$met_shap <- rowSums((dplyr::select(shap_values_all, contains("HMDB"))))
shap_values_all$apoe4 <- rowSums((dplyr::select(shap_values_all, contains("APOE4"))))
shap_values_all$smkk <- rowSums((dplyr::select(shap_values_all, contains("smkk"))))
shap_importance_all <- as.data.frame(colMeans(abs(shap_values_all))) %>%
  rownames_to_column(var = "predictor")
colnames(shap_importance_all) <- c("predictor","overall")

shap_values_yr_15 <- as.data.frame(shap_values_yr_15)
shap_values_yr_15$met_shap <- rowSums((dplyr::select(shap_values_yr_15, contains("HMDB"))))
shap_values_yr_15$apoe4 <- rowSums((dplyr::select(shap_values_yr_15, contains("APOE4"))))
shap_values_yr_15$smkk <- rowSums((dplyr::select(shap_values_yr_15, contains("smkk"))))
shap_importance_yr_15 <- as.data.frame(colMeans(abs(shap_values_yr_15))) %>%
  rownames_to_column(var = "predictor")
colnames(shap_importance_yr_15) <- c("predictor","yr_15")

shap_importance_all <- Reduce(full_join, list(shap_importance_all, shap_importance_yr_15))

shap_importance_filtered <- shap_importance_all %>%
  filter(!grepl("HMDB|APOE4|smkk2|smkk3_4_5", predictor))

ranked_df <- shap_importance_filtered %>%
  dplyr::select(-predictor) %>% 
  mutate(across(everything(), ~ rank(-., ties.method = "first"))) %>% 
  bind_cols(predictor = shap_importance_filtered$predictor, .) %>% 
  arrange(across(everything())) %>% 
  pivot_longer(-predictor, names_to = "interval", values_to = "rank") %>% 
  arrange(interval, rank) %>% 
  pivot_wider(names_from = interval, values_from = predictor) 

write_xlsx(ranked_df,"results/supp/SHAP_rank.xlsx")

gdata::keep(shap_values_all,shap_values_yr_15,
            shap_importance_all,
            model_data_list,
            sure=T)

save.image("results/SHAP_res.Rdata")

##-----------------------------------------
## TICS
##-----------------------------------------

rm(list=ls())

load("data/prediction_data_tics.RData")

subset_selected <- c("all","APOE4_carrier","APOE4_noncarrier")
outcome_tics_selected <- c("av_zscre_tb")

rf_res_list <- list(); i <- 1

for (subset in subset_selected){
  for (outcome in outcome_tics_selected){
      predictor_selected <- c("pred_met_all","lifestyle_APOE4_PRS","lifestyle_APOE4","lifestyle_PRS","lifestyle_short")
      predictor_selected_sub <- predictor_selected[predictor_selected %in% names(binary_formula_list[[subset]][[outcome]])]
      for (predictor in predictor_selected_sub){

        train_data <- model_data_list[[subset]][["train"]]
        test_data <- model_data_list[[subset]][["test"]]
        
        train_data <- train_data[!is.na(train_data[[outcome]]), ]
        test_data <- test_data[!is.na(test_data[[outcome]]), ]
        
        train_data[[outcome]] <- as.factor(train_data[[outcome]])
        test_data[[outcome]] <- as.factor(test_data[[outcome]])
        
        rf_mod <- tryCatch(randomForest(binary_formula_list[[subset]][[outcome]][[predictor]], 
                                        data = train_data, 
                                        importance = TRUE), 
                           warning = function(w) w)
        
        if (length(rf_mod) > 10){
          pred_train <- predict(rf_mod, newdata = model_data_list[[subset]][["train"]], type = "prob")[, 2]
          pred_test <- predict(rf_mod, newdata = model_data_list[[subset]][["test"]], type = "prob")[, 2]
          
          auc_train <- auc(model_data_list[[subset]][["train"]][[outcome]], pred_train)
          auc_test <- auc(model_data_list[[subset]][["test"]][[outcome]], pred_test)
          
          n_train <- sum(!is.na(model_data_list[[subset]][["train"]][[outcome]]))
          n_test <- sum(!is.na(model_data_list[[subset]][["test"]][[outcome]]))
        } else {
          auc_train <- NA; auc_test <- NA
          n_train <- NA; n_test <- NA
        }
        rf_res_list[[i]] <- c("rf", subset, outcome, predictor, auc_train, auc_test, n_train, n_test); i <- i + 1
    }
  }
}

rf_res_dat <- as.data.frame(do.call("rbind", rf_res_list))
colnames(rf_res_dat) <- c("model","subset","subtype","predictor","train_auc","test_auc","train_n","test_n")
rf_res_dat <- rf_res_dat %>%
  mutate(predictor = gsub("lifestyle_PRS", "lifestyle_APOE4_PRS",predictor)) %>%
  mutate(across(c(train_auc,test_auc,train_n,test_n), as.numeric),
         subset_subtype = paste0(subset,"_",subtype)) %>%
  subset(select=c(subset_subtype,predictor,test_auc))

rf_res_out <- reshape(rf_res_dat, idvar = "predictor", timevar = "subset_subtype", direction = "wide")