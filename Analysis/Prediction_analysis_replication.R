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

load("data/prediction_data_hpfs.RData")

subset <- "all"
outcome_selected <- c("all")

##------------------------
## Cox model
##------------------------

predictor_selected <- c("pred_met_all",
                        "lifestyle_APOE4_PRS",
                        "lifestyle_APOE4",
                        "lifestyle_short")

cox_mod_list <- list()
cox_res_list <- list()
for(subset in c("all")){
  cox_mod_list[[subset]] <- list()
  cox_res_list[[subset]] <- list()
  for (outcome in outcome_selected){
    cox_mod_list[[subset]][[outcome]] <- list()
    cox_res_list[[subset]][[outcome]] <- list()
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
for(subset in c("all")){
  cox_auc_list[[subset]] <- list()
  for (outcome in outcome_selected){
    cox_auc_list[[subset]][[outcome]] <- list()
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
                                            "lifestyle_APOE4_PRS",
                                            "pred_met_all"),
                              anno_predictor = c("Baseline~model",
                                                 "+~italic('APOE4')",
                                                 "+~italic('APOE4')~+~PRS",
                                                 "+~italic('APOE4')~+~PRS~+~Metabolites"))

library(tidyr)
auc_dat_list <- list()
for(subset in c("all")){
  auc_dat_list[[subset]] <- list()
  for (outcome in c("all")){
    auc_dat_list[[subset]][[outcome]] <- data.frame(time = cox_auc_list[[subset]][[outcome]][["lifestyle_short"]]$utimes)
    for (predictor in predictor_selected){
      auc_dat_list[[subset]][[outcome]][predictor] <- cox_auc_list[[subset]][[outcome]][[predictor]]$AUC
    }
    auc_dat_list[[subset]][[outcome]] <- auc_dat_list[[subset]][[outcome]] %>% pivot_longer(-time)
    colnames(auc_dat_list[[subset]][[outcome]]) <- c("time","predictor","auc")
    auc_dat_list[[subset]][[outcome]] <- merge(auc_dat_list[[subset]][[outcome]], model_name_dict, by="predictor") %>%
      mutate(anno_predictor = factor(anno_predictor, levels=c("Baseline~model",
                                                              "+~italic('APOE4')",
                                                              "+~italic('APOE4')~+~PRS",
                                                              "+~italic('APOE4')~+~PRS~+~Metabolites")))
  }
}

save(auc_dat_list,file="results/prediction_dementia_auc_by_time_hpfs.RData")

##------------------------
## SHAP
##------------------------

cox_mod_shap_all <- cox_mod_list[["all"]][["all"]][["pred_met_all"]]

predict_function <- function(object, newdata) {
  predict(object, newdata = newdata, type = "risk")
}

shap_values_all <- fastshap::explain(
  cox_mod_shap_all, 
  X = model_data_list[["all"]][["test"]][, names(cox_mod_shap_all$assign)],
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

write_xlsx(ranked_df,"results/supp/SHAP_rank_hpfs.xlsx")

gdata::keep(shap_values_all,
            model_data_list,
            shap_importance_all,
            sure=T)

save.image("results/SHAP_res_hpfs.Rdata")
