library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)

##-----------------------------------------
## Exd Fig. 4: Metabolites PC and AMED components
##-----------------------------------------

rm(list=ls())

### Load data

load("data/merged_data_baseline_nhs_2023_10072024.Rdata")
load("data/metabolites/rf_imputated_metabolites_nhs_10072024.RData")

rownames(rf_imputed_data) <- rf_imputed_data$study_ID
rf_imputed_data <- rf_imputed_data %>% subset(select=c(-study_ID))

### Calculate PC

pc_met <- prcomp(rf_imputed_data, center=T, scale.=T)
cumsum(summary(pc_met)$importance[2,])[1:10]
pc_met_dat <- as.data.frame(pc_met$x[,1:10])
colnames(pc_met_dat) <- paste0("met_",colnames(pc_met_dat))
pc_met_dat$study_ID <- rownames(pc_met_dat)
pc_met_dat <- merge(pc_met_dat, merged_data, by="study_ID")

table(rownames(rf_imputed_data) == merged_data$study_ID)
pc_met_dist_matrix <- dist(pc_met$x[,1:10])
pc_gene_dist_matrix <- dist(merged_data[,grepl("PC",colnames(merged_data))])

### Correlations between AMED and PCs

amed_list_dat <- data.frame(amed_name = c("AMED_avg","veg_avg","frt_avg","nut_avg","leg_avg","whgrn_avg","fish_avg","mon10_avg","sat10_avg","rmt_avg","etoh10_avg"),
                            anno_amed_name = c("MedDiet index score","Vegetables","Fruits","Nuts","Legumes","Whole grains","Fish","Monounsaturated fat","Saturated fat","Red and processed meat","Alcohol"),
                            component_cat = c("AMED","veg_fruit",
                                              "veg_fruit",
                                              "nut_legume",
                                              "nut_legume",
                                              "wgrains",
                                              "fish",
                                              "fat",
                                              "fat",
                                              "meat",
                                              "alco"))

pc_corr_list <- c(); i <- 1
for (var in amed_list_dat$amed_name){
  tmp_cor_pc1 <- cor.test(pc_met_dat[[var]],pc_met_dat$met_PC1, use="complete.obs")
  tmp_cor_pc2 <- cor.test(pc_met_dat[[var]],pc_met_dat$met_PC2, use="complete.obs")
  pc_corr_list[[i]] <- c(var, sum(!is.na(pc_met_dat[[var]])),
                         tmp_cor_pc1$estimate, tmp_cor_pc1$p.value,
                         tmp_cor_pc2$estimate, tmp_cor_pc2$p.value); i <- i + 1
}

pc_corr_dat <- as.data.frame(do.call("rbind",pc_corr_list))
colnames(pc_corr_dat) <- c("comp","n","corr_pc1","p_pc1","corr_pc2","p_pc2")
pc_corr_dat <- pc_corr_dat %>%
  mutate(across(c(n,corr_pc1,p_pc1,corr_pc2,p_pc2), as.numeric)) %>%
  mutate(corr_pc1_scaled = corr_pc1 * (max(pc_met_dat$met_PC1) - min(pc_met_dat$met_PC1)) / 2,
         corr_pc2_scaled = corr_pc2 * (max(pc_met_dat$met_PC2) - min(pc_met_dat$met_PC2)) / 2)

pc_corr_dat <- merge(pc_corr_dat, amed_list_dat, by.x="comp", by.y="amed_name")

### PC plot

transform_factor <- 150

p_pc <- ggplot(data=pc_met_dat, aes(x=met_PC1, y=met_PC2)) +
  geom_point(size=3, aes(fill=AMED_avg), shape=21, color="white", alpha=0.9) +
  scale_fill_gradient(low = "#ADD8E6", high = "#003366") +
  xlab("PC1 (metabolite)") +
  ylab("PC2 (metabolite)") +
  labs(title="", fill="MedDiet index score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title.align = 0.5,
        legend.position = "right") +
  new_scale_color()  

p_pc <- p_pc +
  geom_segment(data = pc_corr_dat, aes(x = 0, y = 0, xend = corr_pc1 * transform_factor, yend = corr_pc2 * transform_factor),
               arrow = arrow(length = unit(0.08, "inches")), size = 0.8, color = "black") +
  scale_color_manual(values=c("#A8A7A4","black","#D89E5F", "#FAD02E", "#9D5C63", "#2C7BB6", "#7D8C3C", "#D8A8B8")) +
  geom_label_repel(
    data = pc_corr_dat,
    aes(x = corr_pc1 * transform_factor, y = corr_pc2 * transform_factor, label = anno_amed_name, color = component_cat),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"), show_guide = FALSE
  ) +
  theme_cowplot() +
  scale_x_continuous(sec.axis = sec_axis(~ . / transform_factor, name = 'Correlation r with PC1')) +  
  scale_y_continuous(sec.axis = sec_axis(~ . / transform_factor, name = 'Correlation r with PC2'))  

pdf("results/extended/extended_fig4.pdf", width = 6, height = 4.8, onefile = F)
p_pc
dev.off()