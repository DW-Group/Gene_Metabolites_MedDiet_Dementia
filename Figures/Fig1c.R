library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)

##-----------------------------------------
## Fig. 1c: Genetic PC plot
##-----------------------------------------

rm(list=ls())

### Data

load("data/merged_data_baseline_nhs_2023.Rdata")
load("results/met_dem_cox_baseline_nhs_2023_filtered_res.RData")

merged_data <- merged_data %>%
  mutate(APOE4_anno_cat = case_when(APOE4 == 0 ~ "Noncarrier",
                                    APOE4 == 1 ~ "Heterozygote",
                                    APOE4 == 2 ~ "Homozygote"),
         APOE4_anno_cat = factor(APOE4_anno_cat, levels=c("Noncarrier","Heterozygote","Homozygote")))
table(merged_data$APOE4_anno_cat)
table(merged_data$APOE4_anno_cat, merged_data$adcase)

### Correlation

pc_corr_list <- c(); i <- 1
for (met_id in met_selected){
  tmp_cor_pc1 <- cor.test(merged_data[[met_id]],merged_data$PC1, use="complete.obs")
  tmp_cor_pc2 <- cor.test(merged_data[[met_id]],merged_data$PC2, use="complete.obs")
  pc_corr_list[[i]] <- c(met_id, sum(!is.na(merged_data[[met_id]])),
                         tmp_cor_pc1$estimate, tmp_cor_pc1$p.value,
                         tmp_cor_pc2$estimate, tmp_cor_pc2$p.value); i <- i + 1
}

pc_corr_dat <- as.data.frame(do.call("rbind",pc_corr_list))
colnames(pc_corr_dat) <- c("hmdb_id","n","corr_pc1","p_pc1","corr_pc2","p_pc2")
pc_corr_dat <- pc_corr_dat %>%
  mutate(across(c(n,corr_pc1,p_pc1,corr_pc2,p_pc2), as.numeric))

pc_corr_dat <- pc_corr_dat %>%
  mutate(corr_pc1_scaled = corr_pc1 * (max(merged_data$PC1) - min(merged_data$PC1)) / 2,
         corr_pc2_scaled = corr_pc2 * (max(merged_data$PC2) - min(merged_data$PC2)) / 2,
         corr_pc1_rank = rank(-abs(corr_pc1)),
         corr_pc2_rank = rank(-abs(corr_pc2)),
         corr_pc1_scaled_rank = rank(-abs(corr_pc1_scaled)),
         corr_pc2_scaled_rank = rank(-abs(corr_pc2_scaled)),
         corr_rank = rank(corr_pc1_rank + corr_pc2_rank),
         corr_scaled_rank = rank(corr_pc1_scaled_rank + corr_pc2_scaled_rank),
         sig_corr_pc1 = ifelse(p_pc1 < 0.05, T, F),
         sig_corr_pc2 = ifelse(p_pc2 < 0.05, T, F),
         sig_corr_both_pc12 = ifelse(p_pc2 < 0.05 & p_pc1 < 0.05, T, F),
         sig_corr_either_pc12 = ifelse(p_pc2 < 0.05 | p_pc1 < 0.05, T, F))

sum(pc_corr_dat$sig_corr_pc1);sum(pc_corr_dat$sig_corr_pc2);sum(pc_corr_dat$sig_corr_either_pc12);sum(pc_corr_dat$sig_corr_both_pc12)

supp_pc_corr_dat <- merge(all_met_dat, pc_corr_dat, by="hmdb_id")

### Plot

pc_corr_dat <- merge(pc_corr_dat, all_met_dat, by="hmdb_id") %>%
  group_by(super_class_metabolon) %>%
  slice_min(corr_rank)

pc_corr <- ggplot() +
  geom_point(data=merged_data, aes(x=PC1, y=PC2, color=APOE4_anno_cat), size=1) +
  scale_color_manual(values = c("tan","tan3","tan4")) +
  labs(color=bquote(italic('APOE4')), x="PC1 (genetic) / Correlation r", y="PC2 (genetic) / Correlation r") +
  new_scale_colour() +
  geom_segment(data = pc_corr_dat, aes(x = 0, y = 0, xend = corr_pc1, yend = corr_pc2, color = super_class_metabolon),
               arrow = arrow(length = unit(0.08, "inches")), size = 0.6, show_guide = FALSE) +
  scale_color_manual(values=classcol) +
  geom_label_repel(
    data=pc_corr_dat,
    aes(x = corr_pc1, y = corr_pc2, label = anno_metabolite_name, color = super_class_metabolon),
    size = 2.5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"), show_guide = FALSE
  ) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(0.5, 'cm'), 
        legend.key.width = unit(0.001, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.background = element_blank(),
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10)
  ) 

pdf("results/figures/fig1c_pc_plot.pdf", width = 6, height = 4, onefile = F)
pc_corr
dev.off()
  