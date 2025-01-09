library(dplyr)
library(data.table)
library(caret)
library(pROC)
library(ggplot2)

##-----------------------------------------
## Fig. 3c: RF AMED and metabolites 
##-----------------------------------------

rm(list=ls())

### Load data

load("data/merged_data_baseline_nhs_2023_10072024.Rdata")
load("data/metabolites/rf_imputated_metabolites_nhs_10072024.RData")

merged_sub <- merged_data %>%
  subset(select=c(study_ID,AMED_avg)) %>%
  mutate(AMED_avg_q = ntile(AMED_avg, 4),
         AMED_out = case_when(AMED_avg_q == 1 ~ 0,
                              AMED_avg_q == 4 ~ 1,
                              AMED_avg_q == 2 ~ NA,
                              AMED_avg_q == 3 ~ NA))

rf_reg_dat <- merge(merged_sub,rf_imputed_data,by="study_ID")
rf_reg_sub <- rf_reg_dat %>% 
  filter(!is.na(AMED_out)) %>% 
  mutate(AMED_out = as.factor(AMED_out))

### Split training and test - binary

set.seed(99)
sample_ind <- sample(c(TRUE,FALSE), nrow(rf_reg_sub),  
                     replace=TRUE, prob=c(0.6,0.4)) 
rf_reg_bin_train <- rf_reg_sub[sample_ind,]
rf_reg_bin_test <- rf_reg_sub[!sample_ind,]

### Random forest regression - binary

# trControl <- trainControl(method = "cv",
#                           number = 5,
#                           search = "grid")
# 
# rf_bin_res <- train(rf_reg_bin_train %>% subset(select=grepl("^HMDB",colnames(.))), rf_reg_bin_train$AMED_out,
#                 method = "rf",
#                 metric = "Accuracy",
#                 trControl = trControl)
# save(rf_bin_res, file = "results/diet_met_rf_model_binary_10132024.RData")

load("results/diet_met_rf_model_binary_10132024.RData")

rf_reg_bin_test <- rf_reg_bin_test %>%
  mutate(pred_AMED_out = predict(rf_bin_res, newdata=rf_reg_bin_test, type="prob"))
rf_reg_bin_roc <- roc(rf_reg_bin_test$AMED_out, rf_reg_bin_test$pred_AMED_out[,1])
rf_reg_bin_roc_anno <- paste0("AUC: ",round(rf_reg_bin_roc$auc,2)," (",round(ci(rf_reg_bin_roc)[1],2),", ",round(ci(rf_reg_bin_roc)[3],2),")")

roc_plot_data <- data.frame(TPR = rf_reg_bin_roc$sensitivities,
                            FPR = 1-rf_reg_bin_roc$specificities)

### ROC curve

p_roc <- ggplot() +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               colour = "grey",lty = 2, size=1) +
  geom_line(aes(x=FPR, y=TPR), color="#29608a", data=roc_plot_data, size=2) +
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1)) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1)) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("MedDiet index score \n(top vs. bottom quartile)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=2)) +
  theme(axis.title = element_text(size=20, colour="black"),
        axis.text = element_text(size=20, colour="black"),
        plot.title = element_text(size=20, hjust=0.5, colour="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        legend.background = element_blank(),
        legend.position = c(0.58, 0.12)) +
  guides(colour = guide_legend(keyheight = 1.5)) +
  annotate("text", x = 0.18, y = 0.03, label = rf_reg_bin_roc_anno, hjust=0, vjust = 0, color = "black", size=7) 

pdf("results/figures/fig3c_auc_rf_11202024.pdf", width = 4.5, height = 5, onefile = F)
p_roc
dev.off()
