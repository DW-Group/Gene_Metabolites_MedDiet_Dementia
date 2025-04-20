library(dplyr)
library(data.table)
library(ggplot2)
library(egg)
library(tidyr)
library(caret)
library(pROC)
library(ggpubr)
library(ggrepel)
library(DescTools)
library(grid)
library(gridExtra)
library(readxl)
library(cowplot)

##-----------------------------------------
## Exd Fig. 6: MedDiet results in HPFS
##-----------------------------------------

##-----------------------------------------
## 6a: dementia by AMED
##-----------------------------------------

rm(list=ls())

### Load data

tmp_dat <- fread("scripts/amed_dementia_spline_baseline_hpfs.txt") %>%
  mutate(across(c(AMEDCON,Estimate,Lower,Upper), as.numeric))

### Spline curve

pdf("results/extended/extended_fig6a.pdf", width = 2.5, height = 3, onefile = F)
p_dementia_amed <- ggplot(tmp_dat, aes(x = AMEDCON, y = Estimate)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, linetype=1, fill="#7A8D50") +
  labs(x="MedDiet index score",y="HR (95% CI) of dementia risk") +
  ggtitle("HPFS") +
  geom_line(color="#4F5937", linewidth=1) + 
  geom_hline(yintercept=1, linetype='dashed', alpha=.5) +
  theme_bw() +  
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"), 
    panel.border = element_blank(), 
    axis.ticks = element_line(color = "black") 
  ) +
  annotate("text", x = 0, y = 0.55 ,
           label = "P < 0.001", hjust = 0, vjust = 0,
           color = "black", size = 4)
p_dementia_amed
dev.off()

##-----------------------------------------
## 6b: AMED-dementia by APOE4
##-----------------------------------------

rm(list=ls())

### Load data

subgroup_res <- read.csv("results/amed_APOE4_PRS_subgroup_hpfsfull_baseline.csv") %>%
  filter(variable=="amedcon",
         modelno==1,
         grepl("apoe_3cat",strata)) %>%
  mutate(across(c(beta,see), as.numeric)) %>%
  mutate(hr=exp(beta),
         lower=exp(beta-1.96*see),
         upper=exp(beta+1.96*see)) %>%
  mutate(anno_group = case_when(strata == "apoe_3cat0" ~ "Noncarrier",
                                strata == "apoe_3cat1" ~ "Heterozygote",
                                strata == "apoe_3cat2" ~ "Homozygote"))

inter_res <- read.csv("results/amed_APOE4_interaction_hpfsfull_baseline.csv") %>%
  filter(modelno %in% c(5),
         grepl("int",variable)) 

hr_res <- subgroup_res %>%
  mutate(p_inter = case_when(strata == "apoe_3cat0" ~ NA,
                             strata == "apoe_3cat1" ~ format(round(as.numeric(inter_res[inter_res$variable=="int_apoe3cat1amed",]$ProbChiSq),2),nsmall=2),
                             strata == "apoe_3cat2" ~ format(round(as.numeric(inter_res[inter_res$variable=="int_apoe3cat2amed",]$ProbChiSq),2),nsmall=2))) %>%
  mutate(anno_group = factor(anno_group, levels=c("Noncarrier","Heterozygote","Homozygote")))

### Forest plot

p_forest <- ggplot(data=hr_res, aes(y=anno_group, x=hr, xmin=lower, xmax=upper)) +
  geom_linerange(size=0.5,position=position_dodge(width = 0.5), alpha=0.8) +
  geom_point(aes(fill=anno_group),color="black", shape=21, size=2.5, stroke = 0.5, position=position_dodge(width = 0.5)) +
  labs(title="Dementia risk and MedDiet index\nscore by APOE4 (HPFS)", x='HR (95% CI) of dementia risk', y='') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  scale_y_discrete(limits=rev) +
  scale_fill_manual(values=c("#D4E4B8","#7A8D50","#4F5937"))+
  theme_bw() +
  theme(legend.title= element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(colour = "black",size=9))+
  theme(axis.title.x = element_text(colour = "black", size=8),
        axis.text.y = element_text(colour = "black", size=8),
        axis.text.x = element_text(colour = "black", size=8),
        panel.border = element_blank(), 
        axis.line = element_line(color = "black") 
  ) +
  theme(legend.position="none") 

pdf("results/extended/extended_fig6b.pdf", width = 3, height = 2, onefile = F)
p_forest
dev.off()

##-----------------------------------------
## 6c: RF AMED and metabolites
##-----------------------------------------

rm(list=ls())

### Load data

load("data/merged_data_baseline_hpfs_2023.Rdata")
load("data/metabolites/rf_imputated_metabolites_hpfs.RData")

merged_sub <- merged_data_hpfs %>%
  subset(select=c(study_ID,amedcon)) %>%
  mutate(amedcon_q = ntile(amedcon, 4),
         AMED_out = case_when(amedcon_q == 1 ~ 0,
                              amedcon_q == 4 ~ 1,
                              amedcon_q == 2 ~ NA,
                              amedcon_q == 3 ~ NA))

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
# save(rf_bin_res, file = "results/diet_met_rf_model_binary_hpfs.RData")

load("results/diet_met_rf_model_binary_hpfs.RData")

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
  ggtitle("HPFS: MedDiet index score \n(top vs. bottom quartile)") +
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

pdf("results/extended/extended_fig6c.pdf", width = 4.5, height = 5, onefile = F)
p_roc
dev.off()

##-----------------------------------------
## 6d: subgroup heatmap (NHS + HPFS)
##-----------------------------------------

rm(list=ls())

### Diet results

source("scripts/source_color_superclass.R")
load("data/merged_data_baseline_nhs_2023.Rdata")
load("data/merged_data_baseline_hpfs_2023.Rdata")

# HPFS

diet_res_hpfs <- read_xlsx("results/diet_gene_met_main_subgroup_hpfs.xlsx") %>%
  filter(hmdb_id %in% met_selected) %>%
  mutate(diet_beta = ifelse(is.na(diet_p) | n < 50, NA, diet_beta),
         age_beta = ifelse(is.na(age_p) | n < 50, NA, age_beta)) %>%
  filter(diet=="amedcon",
         (subgroup=="all" & covar=="full") | 
           (subgroup %in% c("APOE4_noncarrier","APOE4_carrier") & covar=="full")) %>%
  mutate(anno_subgroup = case_when(subgroup=="all" ~ "All",
                                   subgroup=="APOE4_noncarrier" ~ "APOE4\nNoncarrier",
                                   subgroup=="APOE4_carrier" ~ "APOE4\nCarrier")) %>%
  mutate(anno_subgroup = factor(anno_subgroup, levels=c("All","APOE4\nNoncarrier","APOE4\nCarrier"))) %>%
  mutate(wins_diet_beta = Winsorize(diet_beta, val=quantile(diet_beta, probs=c(0.01, 0.99), na.rm=T))) %>%
  mutate(diet_fdr_all = p.adjust(diet_p,method="BH"),
         diet_fdr_all_group = case_when(diet_fdr_all <= 0.05 ~ "FDR < 0.05",
                                        diet_fdr_all > 0.05 & diet_p <= 0.05 ~ "p < 0.05",
                                        diet_p > 0.05 | is.na(diet_p) ~ "Others"),
         diet_fdr_all_group = factor(diet_fdr_all_group, levels = c("FDR < 0.05","p < 0.05","Others")),
         sig = case_when(diet_fdr_all <= 0.05 ~ "*",
                         diet_fdr_all > 0.05 ~ "",
                         is.na(diet_fdr_all) ~ ""))
length(unique(diet_res_hpfs$hmdb_id))

# NHS

diet_res_nhs <- read_xlsx("results/diet_gene_met_main_subgroup.xlsx") %>%
  filter(hmdb_id %in% diet_res_hpfs$hmdb_id) %>%
  mutate(diet_beta = ifelse(is.na(diet_p) | n < 50, NA, diet_beta),
         age_beta = ifelse(is.na(age_p) | n < 50, NA, age_beta)) %>%
  filter(diet=="AMED_avg",
         (subgroup=="all" & covar=="full") | 
           (subgroup %in% c("APOE4_noncarrier","APOE4_carrier") & covar=="full")) %>%
  mutate(anno_subgroup = case_when(subgroup=="all" ~ "All",
                                   subgroup=="APOE4_noncarrier" ~ "APOE4\nNoncarrier",
                                   subgroup=="APOE4_carrier" ~ "APOE4\nCarrier")) %>%
  mutate(anno_subgroup = factor(anno_subgroup, levels=c("All","APOE4\nNoncarrier","APOE4\nCarrier"))) %>%
  mutate(wins_diet_beta = Winsorize(diet_beta, val=quantile(diet_beta, probs=c(0.01, 0.99), na.rm=T))) %>%
  mutate(diet_fdr_all = p.adjust(diet_p,method="BH"),
         diet_fdr_all_group = case_when(diet_fdr_all <= 0.05 ~ "FDR < 0.05",
                                        diet_fdr_all > 0.05 & diet_p <= 0.05 ~ "p < 0.05",
                                        diet_p > 0.05 | is.na(diet_p) ~ "Others"),
         diet_fdr_all_group = factor(diet_fdr_all_group, levels = c("FDR < 0.05","p < 0.05","Others")),
         sig = case_when(diet_fdr_all <= 0.05 ~ "*",
                         diet_fdr_all > 0.05 ~ "",
                         is.na(diet_fdr_all) ~ ""))
length(unique(diet_res_nhs$hmdb_id))

tmp_rm_met1 <- diet_res_nhs[is.na(diet_res_nhs$wins_diet_beta),]$hmdb_id
tmp_rm_met2 <- diet_res_hpfs[is.na(diet_res_hpfs$wins_diet_beta),]$hmdb_id

diet_res_nhs <- diet_res_nhs %>% 
  filter(!(hmdb_id %in% c(tmp_rm_met1,tmp_rm_met2))) %>% 
  mutate(diet_beta_nhs = diet_beta) %>% 
  dplyr::select(c(hmdb_id, anno_metabolite_name, super_class_metabolon, anno_subgroup, diet_beta_nhs))
diet_res_hpfs <- diet_res_hpfs %>% 
  filter(hmdb_id %in% diet_res_nhs$hmdb_id,
         !(hmdb_id %in% c(tmp_rm_met1,tmp_rm_met2))) %>% 
  mutate(diet_beta_hpfs = diet_beta) %>% 
  dplyr::select(c(hmdb_id, anno_subgroup, diet_beta_hpfs))
diet_res <- merge(diet_res_nhs, diet_res_hpfs, by=c("hmdb_id","anno_subgroup"))

### Heatmap

diet_res <- diet_res %>% arrange(desc(super_class_metabolon), desc(anno_metabolite_name))
diet_res$super_class_metabolon <- factor(diet_res$super_class_metabolon, levels=names(classcol))
metabolite_name_vec <- unique(diet_res$anno_metabolite_name)
diet_res$anno_metabolite_name <- factor(diet_res$anno_metabolite_name, levels=metabolite_name_vec)

# Annotation

class_met <- ggplot(diet_res, aes(x=anno_metabolite_name, y=1)) +
  geom_tile(aes(fill=super_class_metabolon)) +
  labs(fill="Super class of metabolites") +
  scale_fill_manual(values=classcol) +
  scale_x_discrete(limits=rev) +
  ylab("") +
  theme_bw() +
  theme(panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

# NHS

nhs_heatmap <- ggplot(diet_res, aes(x=anno_metabolite_name, y=anno_subgroup)) +
  geom_tile(aes(fill=diet_beta_nhs)) +
  scale_fill_gradient2(low="#20854EFF",
                       mid="white",
                       high="#7876B1FF",
                       midpoint = 0,
                       na.value="lightgray",
                       guide = "colourbar") +
  scale_x_discrete(limits=rev, position = "top") +
  scale_y_discrete(limits=rev, position = "left") +
  labs(title="NHS",fill=expression(beta), size=12) +
  xlab("") +
  ylab("") +
  theme_bw()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black"),
        panel.border = element_blank(), 
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8, hjust = 1),
        legend.key.size = unit(0.4, 'cm'))

# HPFS

hpfs_heatmap <- ggplot(diet_res, aes(x=anno_metabolite_name, y=anno_subgroup)) +
  geom_tile(aes(fill=diet_beta_hpfs)) +
  scale_fill_gradient2(low="#20854EFF",
                       mid="white",
                       high="#7876B1FF",
                       midpoint = 0,
                       na.value="lightgray",
                       guide = "colourbar") +
  scale_x_discrete(limits=rev, position = "top") +
  scale_y_discrete(limits=rev, position = "left") +
  labs(title="HPFS",fill=expression(beta), size=12) +
  xlab("") +
  ylab("") +
  theme_bw()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8, colour = "black"),
        panel.border = element_blank(),  # Remove all panel borders
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8, hjust = 1),
        legend.key.size = unit(0.4, 'cm'))

pdf("results/extended/extended_fig6d.pdf", width = 9, height = 4, onefile = F)
fig6d_diet <- ggpubr::ggarrange(nhs_heatmap,hpfs_heatmap, class_met, ncol=1, nrow=4,
                                widths = c(30), heights=c(1,1,0.25),
                                align = "v")
fig6d_diet
dev.off()

class_met_legend <- class_met + theme(panel.border = element_rect(size=1.2),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.text.y = element_blank(),
                                      axis.title =element_blank(),
                                      legend.position = "bottom",
                                      legend.text = element_text(colour = "black", size = 14),
                                      legend.title = element_text(colour = "black", size = 14))

pdf("results/extended/extended_fig6d_legend.pdf", width = 5, height = 2, onefile = F)
legend_fig6d <- cowplot::get_plot_component(class_met_legend, 'guide-box', return_all = TRUE)
grid.newpage()
grid.draw(legend_fig6d[[3]])
dev.off()

##-----------------------------------------
## Supp Fig 6: subgroup scatter plot (NHS + HPFS)
##-----------------------------------------

tmp_res <- list()
p_list <- list(); i <- 1
for (subgroup in c("All","APOE4\nNoncarrier","APOE4\nCarrier")){
  tmp_diet_res <- diet_res %>%
    filter(anno_subgroup==subgroup) %>%
    mutate(beta1 = diet_beta_nhs,
           beta2 = diet_beta_hpfs)
  tmp_r2 <- format(round(cor(tmp_diet_res$beta1, tmp_diet_res$beta2, use="complete.obs")[1],2), nsmall = 2)
  tmp_res[[subgroup]] <- tmp_r2
  p_list[[i]] <- ggplot(data=tmp_diet_res, aes(x=beta1, y=beta2)) +
    geom_point(aes(fill=super_class_metabolon), size=2, shape=21, stroke=0.2, alpha=0.8) +
    scale_fill_manual(values=classcol) +
    geom_vline(xintercept=0, linetype='dashed', alpha=.5) +
    geom_hline(yintercept=0, linetype='dashed', alpha=.5) +
    labs(
      x = bquote("NHS" ~ beta["MedDiet"]),
      y = bquote("HPFS" ~ beta["MedDiet"])
    ) +      
    theme_bw() +
    theme(axis.title = element_text(size=10),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),    
          axis.ticks.x = element_line(color="black"),
          axis.ticks.y = element_line(color="black"),    
          panel.border = element_rect(color="black", size=1),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    geom_text_npc(npcx = "left", npcy = "top", label = paste0("r = ",tmp_r2), size=4) +
    geom_text_npc(npcx = "right", npcy = "bottom", label = gsub("\n"," ", subgroup), size=4) +
    theme(legend.position = "none"); i <- i + 1
}

combined_plot <- wrap_plots(p_list, ncol = 3, nrow = 1)
pdf("results/supp/supp_fig_6.pdf", width = 12, height = 4)  
print(combined_plot)
dev.off()
 
