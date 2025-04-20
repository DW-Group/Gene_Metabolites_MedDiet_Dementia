library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(stringr)
library(readxl)

##-----------------------------------------
## Exd Fig. 3: Significant interaction by superclass  
##-----------------------------------------

rm(list=ls())

### Load data

load("results/met_dem_cox_baseline_nhs_2023_filtered_res_r1_03112025.RData")
load("data/genetic/ad_variants_dictionary_r1.RData"); nhs_inter_ind_var_res <- nhs_inter_ind_var_res %>% filter(rsid %in% ad_dict$rsid)

##----------------------------------
## 3a: APOE4 interactions
##----------------------------------

sig_apoe4_inter <- nhs_inter_apoe4_res %>% 
  filter(pvalue_pvalue_p_inter_met_gene_1<0.05|pvalue_pvalue_p_inter_met_gene_2<0.05) %>%
  mutate(fdr_sig = ifelse(two_inter_fdr_pvalue_p_inter_met_gene_1<0.05|two_inter_fdr_pvalue_p_inter_met_gene_2<0.05,1,0))
table(sig_apoe4_inter$fdr_sig)

count_by_super_class <- sig_apoe4_inter %>%
  group_by(super_class_metabolon, fdr_sig) %>%
  count() %>%
  ungroup()

total_counts <- count_by_super_class %>%
  group_by(super_class_metabolon) %>%
  summarise(total_n = sum(n)) %>%
  arrange(desc(total_n))

count_by_super_class <- count_by_super_class %>%
  mutate(super_class_metabolon = factor(count_by_super_class$super_class_metabolon,
                                        levels = total_counts$super_class_metabolon),
         fdr_sig = factor(fdr_sig, levels=c(1,0))
)

pdf("results/extended/extended_fig3a.pdf", width = 7, height = 5, onefile = F)
ggplot(count_by_super_class, aes(y = super_class_metabolon, x = n, fill = factor(fdr_sig))) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = c("1" = "#DAA520", "0" = "#ad4a44"), 
                    labels = c("FDR < 0.05", "P < 0.05 & FDR > 0.05 ")) +
  scale_y_discrete(limits=rev) +
  labs(
    title = " ",
    x = "Number of significant metabolites",
    y = " ",
    fill = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.text = element_text(color = "black"),  # Black legend text
    axis.text.x = element_text(color = "black"),  # Black x-axis text
    axis.text.y = element_text(color = "black"),  # Black y-axis text
    axis.title.x = element_text(color = "black"), # Black x-axis title
    axis.title.y = element_text(color = "black"), # Black y-axis title
    panel.border = element_rect(color = "black", fill = NA),  # Black panel border
    axis.ticks = element_line(color = "black")  # Black axis ticks
  )
dev.off()

##----------------------------------
## 3b: Individual variant interactions
##----------------------------------

sig_var_inter <- nhs_inter_ind_var_res %>% 
  filter(p_inter_met_var<0.05) %>%
  mutate(fdr_sig = ifelse(p_fdr_inter_met_var<0.05,1,0))
table(sig_var_inter$fdr_sig)
sig_var_inter <- sig_var_inter[order(sig_var_inter$fdr_sig, decreasing = T),]
sig_var_inter <- sig_var_inter[!duplicated(sig_var_inter$hmdb_id),]

count_by_super_class <- sig_var_inter %>%
  group_by(super_class_metabolon, fdr_sig) %>%
  count() %>%
  ungroup()

total_counts <- count_by_super_class %>%
  group_by(super_class_metabolon) %>%
  summarise(total_n = sum(n)) %>%
  arrange(desc(total_n))

count_by_super_class <- count_by_super_class %>%
  mutate(super_class_metabolon = factor(count_by_super_class$super_class_metabolon,
                                        levels = total_counts$super_class_metabolon),
         fdr_sig = factor(fdr_sig, levels=c(1,0))
  )

pdf("results/extended/extended_fig3b.pdf", width = 7, height = 5, onefile = F)
ggplot(count_by_super_class, aes(y = super_class_metabolon, x = n, fill = factor(fdr_sig))) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = c("1" = "#DAA520", "0" = "#ad4a44"), 
                    labels = c("FDR < 0.05", "P < 0.05 & FDR > 0.05")) +
  scale_y_discrete(limits=rev) +
  labs(
    title = " ",
    x = "Number of significant metabolites",
    y = " ",
    fill = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.text = element_text(color = "black"),  # Black legend text
    axis.text.x = element_text(color = "black"),  # Black x-axis text
    axis.text.y = element_text(color = "black"),  # Black y-axis text
    axis.title.x = element_text(color = "black"), # Black x-axis title
    axis.title.y = element_text(color = "black"), # Black y-axis title
    panel.border = element_rect(color = "black", fill = NA),  # Black panel border
    axis.ticks = element_line(color = "black")  # Black axis ticks
  )
dev.off()