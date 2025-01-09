library(dplyr)
library(data.table)
library(ggplot2)
library(DescTools)
library(ggrepel)
library(readxl)

##-----------------------------------------
## Fig. 3e: Call-out MedDiet-met
##-----------------------------------------

rm(list=ls())

### Diet results

source("scripts/source_color_superclass.R")
load("results/met_dem_cox_baseline_nhs_2023_filtered_res_11112024.RData")
load("data/merged_data_baseline_nhs_2023_11212024.Rdata")

diet_res <- read_xlsx("results/diet_gene_met_main_subgroup_11212024.xlsx") %>%
  mutate(diet_beta = ifelse(is.na(diet_p) | n < 50, NA, diet_beta),
         age_beta = ifelse(is.na(age_p) | n < 50, NA, age_beta)) %>%
  filter(diet=="AMED_avg",
         (subgroup=="all" & covar=="full") | 
           (subgroup %in% c("APOE4_noncarrier","APOE4_carrier","APOE4_ncopy_1") & covar=="prs") | 
           (subgroup %in% c("APOE4_ncopy_2") & covar=="prs_reduced")) %>%
  mutate(anno_subgroup = case_when(subgroup=="all" ~ "All",
                                   subgroup=="APOE4_noncarrier" ~ "Noncarrier",
                                   subgroup=="APOE4_carrier" ~ "Carrier",
                                   subgroup=="APOE4_ncopy_1" ~ "Heterozygotes",
                                   subgroup=="APOE4_ncopy_2" ~ "Homozygotes")) %>%
  mutate(anno_subgroup = factor(anno_subgroup, levels=c("All","Noncarrier","Carrier","Heterozygotes","Homozygotes")))

diet_inter_res <- read_xlsx("results/diet_gene_met_inter_11212024.xlsx") %>% filter(diet=="AMED_avg")
diet_inter_ind_var_res <- read_xlsx("results/diet_gene_met_inter_ind_var_11212024.xlsx")

diet_res_sig_met <- unique(diet_res[diet_res$diet_p < 0.05,]$hmdb_id)
diet_inter_res_sig_met <- unique(diet_inter_res[diet_inter_res$inter1_p < 0.05 | diet_inter_res$inter2_p < 0.05 | diet_inter_res$lrt_p < 0.05,]$hmdb_id)
tmp <- na.omit(intersect(diet_res_sig_met, diet_inter_res_sig_met))
diet_sub <- diet_res %>% filter(hmdb_id %in% tmp)
diet_inter_sub <- diet_inter_res %>% filter(hmdb_id %in% tmp)

### Adjusted residual

for (hmdb_id in c("HMDB0000168", "HMDB0011103", "HMDB0010368", "HMDB0003357")) {
  
  merged_data[[paste0(hmdb_id, "_resid")]] <- NA
  
  for (apoe4_group in unique(merged_data$APOE4)) {
    tmp_data <- subset(merged_data, APOE4 == apoe4_group & !is.na(merged_data[[hmdb_id]]))
    
    tmp_mod <- glm(as.formula(paste(hmdb_id, "~ agemo + blddate + bmi + fast1 + endpoint + caco + hiedu3 + 
                      husbedu2_3_m + phmsstatus3_4 + fhdem + actcon + nSES + marry + 
                      smkk2 + smkk3_4_5 + dep_antidep + htn + sysbp2 + sysbp3 + 
                      sysbpm + hchol + alcocon + daykcal + PGS002280_scaled + PC1 + 
                      PC2 + PC3 + PC4 + platformgsa + platformhuco2 + platformillu + 
                      platformomni + platformonco")), data=tmp_data)
    
    merged_data[merged_data$APOE4 == apoe4_group & !is.na(merged_data[[hmdb_id]]), 
                paste0(hmdb_id, "_resid")] <- residuals(tmp_mod)
    
  }
}

### Scatter plots

merged_data <- merged_data %>% 
  mutate(APOE4_carrier = as.factor(APOE4_carrier),
         APOE4 = as.factor(APOE4))

tmp_met_name <- "HMDB0000168"; tmp_met_name; fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$anno_metabolite_name
merged_data$tmp_met <- merged_data[[paste0(tmp_met_name,"_resid")]]
p1 <- ggplot(data=merged_data, aes(x=AMED_avg, y=tmp_met, group=APOE4)) +
  geom_point(aes(color=APOE4), alpha=0.5, shape=16, size=0.5, show.legend=F)+
  geom_smooth(method = "lm", aes(color=APOE4), size=1)+
  theme_bw()+
  theme(axis.text.x = element_text(size=9, colour = "black"),
        axis.text.y = element_text(size=9, colour = "black"), 
        axis.title.x = element_text(size=11, colour = "black"),
        axis.title.y = element_text(size=11, colour = "black"), 
        legend.position = c(0.35, 0.18), 
        legend.title = element_text(size=8), 
        legend.text = element_text(size=7),
        legend.key.size = unit(0.3, "cm"), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
  scale_color_manual(values=c("tan","tan3", "tan4"), 
                     labels=c("Noncarrier", "Heterozygote", "Homozygote"))+
  guides(color = guide_legend(title = bquote(italic('APOE4')))) +
  labs(x="MedDiet index score", y=fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$anno_metabolite_name) 

tmp_met_name <- "HMDB0003357"; tmp_met_name; fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$anno_metabolite_name
merged_data$tmp_met <- merged_data[[paste0(tmp_met_name,"_resid")]]
p2 <- ggplot(data=merged_data, aes(x=AMED_avg, y=tmp_met, group=APOE4)) +
  geom_point(aes(color=APOE4), alpha=0.5, shape=16, size=0.5, show.legend=F)+
  geom_smooth(method = "lm", aes(color=APOE4), size=1)+
  theme_bw()+
  theme(axis.text.x = element_text(size=9, colour = "black"),
        axis.text.y = element_text(size=9, colour = "black"), 
        axis.title.x = element_text(size=11, colour = "black"),
        axis.title.y = element_text(size=11, colour = "black"), 
        legend.position = c(0.35, 0.18),  
        legend.title = element_text(size=8), 
        legend.text = element_text(size=7),
        legend.key.size = unit(0.3, "cm"), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
  scale_color_manual(values=c("tan","tan3", "tan4"), 
                     labels=c("Noncarrier", "Heterozygote", "Homozygote"))+
  guides(color = guide_legend(title = bquote(italic('APOE4')))) +
  labs(x="MedDiet index score", y=fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$anno_metabolite_name) 

tmp_met_name <- "HMDB0011103"; tmp_met_name; fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$anno_metabolite_name
merged_data$tmp_met <- merged_data[[paste0(tmp_met_name,"_resid")]]
p3 <- ggplot(data=merged_data, aes(x=AMED_avg, y=tmp_met, group=APOE4)) +
  geom_point(aes(color=APOE4), alpha=0.5, shape=16, size=0.5, show.legend=F)+
  geom_smooth(method = "lm", aes(color=APOE4), size=1)+
  theme_bw()+
  theme(axis.text.x = element_text(size=9, colour = "black"),
        axis.text.y = element_text(size=9, colour = "black"), 
        axis.title.x = element_text(size=11, colour = "black"),
        axis.title.y = element_text(size=11, colour = "black"), 
        legend.position = c(0.35, 0.18), 
        legend.title = element_text(size=8), 
        legend.text = element_text(size=7),
        legend.key.size = unit(0.3, "cm"), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
  scale_color_manual(values=c("tan","tan3", "tan4"), 
                     labels=c("Noncarrier", "Heterozygote", "Homozygote"))+
  guides(color = guide_legend(title = bquote(italic('APOE4')))) +
  labs(x="MedDiet index score", y=fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$anno_metabolite_name)

pdf("results/figures/fig3e_call_out_11232024.pdf", width = 6, height = 2.5, onefile = F)
fig3e <- egg::ggarrange(p1, p2, p3, 
                        nrow=1, ncol=3,
                        widths = c(1,1,1), heights=c(1))
fig3e
dev.off()
 