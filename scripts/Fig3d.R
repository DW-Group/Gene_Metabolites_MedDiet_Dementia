library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(DescTools)
library(grid)
library(gridExtra)
library(readxl)
library(cowplot)

##-----------------------------------------
## Fig. 3d: Subgroup heatmap
##-----------------------------------------

rm(list=ls())

### Diet results

source("scripts/source_color_superclass.R")

diet_res <- read_xlsx("results/diet_gene_met_main_subgroup_11212024.xlsx") %>%
  mutate(diet_beta = ifelse(is.na(diet_p) | n < 50, NA, diet_beta),
         age_beta = ifelse(is.na(age_p) | n < 50, NA, age_beta)) %>%
  filter(diet=="AMED_avg",
         (subgroup=="all" & covar=="full") | 
           (subgroup %in% c("APOE4_noncarrier","APOE4_carrier","APOE4_ncopy_1") & covar=="prs") | 
           (subgroup %in% c("APOE4_ncopy_2") & covar=="prs_reduced")) %>%
  mutate(anno_subgroup = case_when(subgroup=="all" ~ "All",
                                   subgroup=="APOE4_noncarrier" ~ "APOE4\nNoncarrier",
                                   subgroup=="APOE4_carrier" ~ "APOE4\nCarrier",
                                   subgroup=="APOE4_ncopy_1" ~ "APOE4\nHeterozygotes",
                                   subgroup=="APOE4_ncopy_2" ~ "APOE4\nHomozygotes")) %>%
  mutate(anno_subgroup = factor(anno_subgroup, levels=c("All","APOE4\nNoncarrier","APOE4\nCarrier","APOE4\nHeterozygotes","APOE4\nHomozygotes"))) %>%
  mutate(wins_diet_beta = Winsorize(diet_beta, val=quantile(diet_beta, probs=c(0.01, 0.99), na.rm=T))) %>%
  mutate(diet_fdr_all = p.adjust(diet_p,method="BH"),
         diet_fdr_all_group = case_when(diet_fdr_all <= 0.05 ~ "FDR < 0.05",
                                        diet_fdr_all > 0.05 & diet_p <= 0.05 ~ "p < 0.05",
                                        diet_p > 0.05 | is.na(diet_p) ~ "Others"),
         diet_fdr_all_group = factor(diet_fdr_all_group, levels = c("FDR < 0.05","p < 0.05","Others")),
         sig = case_when(diet_fdr_all <= 0.05 ~ "*",
                         diet_fdr_all > 0.05 ~ "",
                         is.na(diet_fdr_all) ~ ""))

tmp_rm_met <- diet_res[is.na(diet_res$wins_diet_beta),]$hmdb_id
diet_res <- diet_res %>%
  filter(!(hmdb_id %in% tmp_rm_met),
         subgroup != "APOE4_carrier")

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

# AMED

amed_strat_heatmap <- ggplot(diet_res %>% filter(diet=="AMED_avg"), aes(x=anno_metabolite_name, y=anno_subgroup)) +
  geom_tile(aes(fill=wins_diet_beta)) +
  scale_fill_gradient2(low="#20854EFF",
                       mid="white",
                       high="#7876B1FF",
                       midpoint = 0,
                       na.value="lightgray",
                       guide = "colourbar") +
  scale_x_discrete(limits=rev, position = "top") +
  scale_y_discrete(limits=rev, position = "left") +
  labs(fill=expression(beta), size=12) +
  xlab("") +
  ylab("") +
  theme_bw()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 2, angle=45, colour = "black", vjust=0, hjust=0),
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

pdf("results/figures/fig3d_diet_heatmap_11252024.pdf", width = 9, height = 2, onefile = F)
fig3d_diet <- ggpubr::ggarrange(amed_strat_heatmap, class_met, ncol=1, nrow=2,
                              widths = c(30), heights=c(1,0.15),
                              align = "v")
fig3d_diet
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
                 legend.title = element_text(colour = "black", size = 14),
                 legend.box.margin = margin(t=200, r=200, b=200, l=200))

pdf("results/figures/fig3d_legend_09202024.pdf", width = 5, height = 2, onefile = F)
legend_fig3d <- cowplot::get_plot_component(class_met_legend, 'guide-box', return_all = TRUE)
grid.newpage()
grid.draw(legend_fig3d[[3]])
dev.off()
