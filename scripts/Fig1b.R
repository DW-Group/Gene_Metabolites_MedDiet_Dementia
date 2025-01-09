library(dplyr)
library(data.table)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)

##-----------------------------------------
## Fig. 1b: Metabolites distribution
##-----------------------------------------

rm(list=ls())

### Data

load("data/merged_data_baseline_nhs_2023_10072024.Rdata")
load("results/met_dem_cox_baseline_nhs_2023_filtered_res_10072024.RData")

metdata_no_transform <- read.csv("data/metabolites/metabolomics_data_full_no_transform_nhs.csv")

metdata_no_transform <- metdata_no_transform %>% 
  filter(cohort %in% c("nhs1")) %>%
  mutate(study = case_when(cohort == "nhs1" ~ "NHS"),
         study_ID = paste0(study,"_",id)) %>%
  filter(study_ID %in% merged_data[!is.na(merged_data$APOE4),]$study_ID) %>% 
  subset(select = c("study_ID",met_selected))

### Data for supp table

tmp <- merge(metdata_no_transform, subset(merged_data, select=c(study_ID,APOE4)), by="study_ID")
metdata_no_transform_APOE4_0 <- tmp %>% filter(APOE4==0) %>% subset(select=-c(study_ID,APOE4))
metdata_no_transform_APOE4_1 <- tmp %>% filter(APOE4==1) %>% subset(select=-c(study_ID,APOE4))
metdata_no_transform_APOE4_2 <- tmp %>% filter(APOE4==2) %>% subset(select=-c(study_ID,APOE4))

metdata_no_transform_APOE4_0_norm <- metdata_no_transform_APOE4_0 / rowSums(metdata_no_transform_APOE4_0, na.rm=T)
metdata_no_transform_APOE4_1_norm <- metdata_no_transform_APOE4_1 / rowSums(metdata_no_transform_APOE4_1, na.rm=T)
metdata_no_transform_APOE4_2_norm <- metdata_no_transform_APOE4_2 / rowSums(metdata_no_transform_APOE4_2, na.rm=T)

table(colnames(metdata_no_transform_APOE4_0_norm)==colnames(metdata_no_transform_APOE4_1_norm))
table(colnames(metdata_no_transform_APOE4_0_norm)==colnames(metdata_no_transform_APOE4_2_norm))

abd_met_by_APOE4 <- data.frame(hmdb_id = colnames(metdata_no_transform_APOE4_0_norm),
                               APOE4_noncarrier_mean_abundance = colMeans(metdata_no_transform_APOE4_0_norm, na.rm=T),
                               APOE4_heterozygote_mean_abundance = colMeans(metdata_no_transform_APOE4_1_norm, na.rm=T),
                               APOE4_homozygote_mean_abundance = colMeans(metdata_no_transform_APOE4_2_norm, na.rm=T))

### Data for abundance plot

metdata_no_transform <- subset(metdata_no_transform,select=-c(study_ID))
metdata_no_transform_norm <- metdata_no_transform * 100 / rowSums(metdata_no_transform, na.rm=T)
table(rowSums(metdata_no_transform_norm, na.rm=T))

abd_met <- data.frame(hmdb_id = colnames(metdata_no_transform_norm),
                      mean_abundance = colMeans(metdata_no_transform_norm, na.rm=T)) %>%
  mutate(sd_abundance = apply(metdata_no_transform_norm, 2, sd, na.rm=T),
         cv = sd_abundance / mean_abundance)

summary(abd_met$mean_abundance)
summary(abd_met$cv)

abd_met <- merge(abd_met, all_met_dat, by="hmdb_id") 
abd_met <- abd_met[order(abd_met$mean_abundance, decreasing=T),]

table(abd_met$super_class_metabolon)

### Circle plot - abundance

tmp_circle_all_dat <- abd_met %>%
  arrange(desc(super_class_metabolon), desc(mean_abundance))
tmp_circle_all_dat$super_class_metabolon <- factor(tmp_circle_all_dat$super_class_metabolon, levels=names(classcol))
tmp_circle_all_dat$anno_metabolite_name <- factor(tmp_circle_all_dat$anno_metabolite_name, levels=unique(tmp_circle_all_dat$anno_metabolite_name))

split_super_class <- tmp_circle_all_dat$super_class_metabolon
summary(log10(tmp_circle_all_dat$mean_abundance))
tmp_circle_all_mat <- as.matrix(log10(tmp_circle_all_dat$mean_abundance)) # No Winsorization
summary(tmp_circle_all_mat[,1])
rownames(tmp_circle_all_mat) <- tmp_circle_all_dat$anno_metabolite_name
tmp_circle_all_mat_cv <- as.matrix(log10(tmp_circle_all_dat$cv)) # No Winsorization
rownames(tmp_circle_all_mat_cv) <- tmp_circle_all_dat$anno_metabolite_name

table(rownames(tmp_circle_all_mat) == unique(tmp_circle_all_dat$anno_metabolite_name))
summary(tmp_circle_all_mat)
table(rownames(tmp_circle_all_mat_cv) == unique(tmp_circle_all_dat$anno_metabolite_name))
summary(tmp_circle_all_mat_cv)

col_fun1 <- colorRamp2(c(min(tmp_circle_all_mat),max(tmp_circle_all_mat)), c("white", "#6497b1"))
col_fun2 <- colorRamp2(c(min(tmp_circle_all_mat_cv),max(tmp_circle_all_mat_cv)), c("white", "darkgray"))

pdf("results/figures/fig1b_circle_10132024.pdf", width = 14, height = 14, onefile = F)
circos.heatmap.initialize(tmp_circle_all_mat, split = split_super_class, cluster = FALSE)
circos.heatmap(tmp_circle_all_mat_cv,
               col = col_fun2,
               split = split_super_class,
               cluster = FALSE,
               rownames.side = "outside",
               track.height = 0.1,
               bg.border="black",
               bg.lwd = 1,
               bg.lty = 1)
circos.heatmap(tmp_circle_all_mat,
               col = col_fun1,
               split = split_super_class,
               cluster = FALSE,
               rownames.side = "none",
               track.height = 0.1,
               bg.border="black",
               bg.lwd = 1,
               bg.lty = 1)
circos.heatmap(as.vector(split_super_class),
               col = classcol,
               track.height=0.02)
circos.clear()
dev.off()

pdf("results/figures/fig1b_circle_legend_10132024.pdf", width = 5,height = 5)
lgd_mean_abundance <- Legend(title = expression(bold(paste("log"[10],"(mean abundance, %)",sep=""))),
                             title_gp = gpar(fontsize = 10),
                             direction = "horizontal",
                             col_fun = col_fun1)
lgd_cv <- Legend(title = expression(bold(paste("log"[10],"(coefficient of variation)",sep=""))),
                             title_gp = gpar(fontsize = 10),
                 direction = "horizontal",
                             col_fun = col_fun2)
lgd_super_class <- Legend(title = "Superclass of metabolites",
                          title_gp = gpar(fontsize = 10, fontface = "bold"),
                          at = names(classcol),
                          legend_gp = gpar(fill = classcol))

lgd_list = packLegend(lgd_mean_abundance, lgd_cv, lgd_super_class, max_height = unit(12, "cm"))
draw(lgd_list)
dev.off()
