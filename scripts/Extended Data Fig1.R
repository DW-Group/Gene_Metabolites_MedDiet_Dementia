library(dplyr)
library(tidyr)
library(data.table)
library(ggpubr)
library(ggplot2)
library(DescTools)
library(ggrepel)
library(grid)
library(gridExtra)
library(ggh4x)
library(ggtext)
library(stringr)
library(Hmisc)
library(readxl)
library(cowplot)

##-----------------------------------------
## Exd Fig. 1: Manhattan 
##-----------------------------------------

##----------------------------------
## 1a: Genetic PC and metabolites
##----------------------------------

rm(list=ls())

### Data

pc_corr_dat <- read.csv("results/supp/genetic_pc_met_corr.csv")
summary(pc_corr_dat$corr_pc1)
summary(pc_corr_dat$corr_pc2)

pc_corr_dat <- pc_corr_dat %>%
  mutate(sign_pc1 = case_when(corr_pc1 > 0 ~ "Positive",
                              corr_pc1 <= 0 ~ "Inverse"),
         sig_pc1 = case_when(p_pc1 < 0.05 ~ "sig",
                             p_pc1 >= 0.05 ~ "non_sig"),
         abs_pc1 = abs(corr_pc1),
         sign_pc2 = case_when(corr_pc2 > 0 ~ "Positive",
                              corr_pc2 <= 0 ~ "Inverse"),
         sig_pc2 = case_when(p_pc2 < 0.05 ~ "sig",
                             p_pc2 >= 0.05 ~ "non_sig"),
         abs_pc2 = abs(corr_pc2)) %>%
  mutate(sign_pc1 = factor(sign_pc1, levels=c("Positive","Inverse")),
         sign_pc2 = factor(sign_pc2, levels=c("Positive","Inverse")))

source("scripts/source_color_superclass.R")

pc_corr_dat <- pc_corr_dat %>% arrange(desc(super_class_metabolon), desc(metabolite_name))
pc_corr_dat$super_class_metabolon <- factor(pc_corr_dat$super_class_metabolon, levels=names(classcol))
metabolite_name_vec <- unique(pc_corr_dat$metabolite_name)
pc_corr_dat$metabolite_name <- factor(pc_corr_dat$metabolite_name, levels=metabolite_name_vec)

tmp <- as.data.frame(str_split(pc_corr_dat$metabolite_name,"/",simplify=T)) %>% mutate(V1 = capitalize(V1)) 
pc_corr_dat$anno_metabolite_name <- tmp$V1

### Annotation

class_met <- ggplot(pc_corr_dat, aes(x=metabolite_name, y=1)) +
  geom_tile(aes(fill=super_class_metabolon)) +
  labs(fill="Superclass of metabolites") +
  scale_fill_manual(values=classcol) +
  scale_x_discrete(limits=rev) +
  ylab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title =element_blank(),
        legend.position = "right",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14))

pdf("results/extended/extended_fig1a_legend1.pdf", width = 25, height = 3, onefile = F)
legend_exdfig1a1 <- cowplot::get_plot_component(class_met, 'guide-box', return_all = TRUE)
grid.newpage()
grid.draw(legend_exdfig1a1[[1]])
dev.off()

### Manhattan plot: PC1

mhtplot_pc1 <- ggplot(pc_corr_dat, aes(x=metabolite_name, y=-log10(p_pc1), fill=sign_pc1)) +
  geom_point(aes(fill=sign_pc1, size=abs_pc1), shape=21) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  scale_fill_manual(values=c("#E67E22","#34495E")) +
  scale_x_discrete(limits=rev) +
  labs(fill="Direction of correlation",size="Effect size") +
  ylab(expression("-log"[10]*"(P)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +
  annotation_custom(grobTree(textGrob(parse(text=("Correlation~with~genetic~PC1")), x=0.35 ,  y=0.92, hjust=1,
                                      gp=gpar(fontsize=20, fontface="bold")))) +
  theme(legend.position = "right", legend.title = element_text(size=14), legend.text = element_text(size=12), legend.key.size = unit(0.6, 'cm')) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

pdf("results/extended/extended_fig1a_legend2.pdf", width = 25, height = 3, onefile = F)
legend_mhtplot_pc1 <- cowplot::get_plot_component(mhtplot_pc1, 'guide-box', return_all = TRUE)
grid.newpage()
grid.draw(legend_mhtplot_pc1[[1]])
dev.off()

### Manhattan plot: PC2

mhtplot_pc2 <- ggplot(pc_corr_dat, aes(x=metabolite_name, y=-log10(p_pc2), fill=sign_pc2)) +
  geom_point(aes(fill=sign_pc2, size=abs_pc2), shape=21) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  scale_fill_manual(values=c("#E67E22","#34495E")) +
  scale_x_discrete(limits=rev) +
  labs(fill="Correlation",size="Effect size") +
  ylab(expression("-log"[10]*"(P)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +
  annotation_custom(grobTree(textGrob(parse(text=("Correlation~with~genetic~PC2")), x=0.35,  y=0.92, hjust=1,
                                      gp=gpar(fontsize=20, fontface="bold")))) +
  theme(legend.position = "top", legend.title = element_text(size=14), legend.text = element_text(size=12), legend.key.size = unit(0.6, 'cm')) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

pdf("results/extended/extended_fig1a.pdf", width = 11, height = 8, onefile = F)
genetic_pc_met_mhtplots <- ggpubr::ggarrange(mhtplot_pc1, class_met, mhtplot_pc2, ncol=1, nrow=3,
                                             widths = c(30), heights=c(1,0.1,1),
                                             align = "v",
                                             legend = "none")
genetic_pc_met_mhtplots
dev.off()

##----------------------------------
## 1b: APOE4 and metabolites
##----------------------------------

rm(list=ls())

### Data

gene_met_dat <- fread("results/gene_met_nhs_2023_11242024.txt") %>%
  filter(grepl("APOE4_ncopy",gene),
         covar == "adj_age")

summary(gene_met_dat$beta_var)
summary(gene_met_dat$p_var)
unique(gene_met_dat$super_class_metabolon)

gene_met_dat <- gene_met_dat %>%
  mutate(sign = case_when(beta_var > 0 ~ "Positive",
                          beta_var <= 0 ~ "Inverse"),
         sig = case_when(p_var < 0.05 ~ "sig",
                         p_var >= 0.05 ~ "non_sig"),
         abs = abs(beta_var)) %>%
  mutate(sign = factor(sign, levels=c("Positive","Inverse")))

source("scripts/source_color_superclass.R")

gene_met_dat <- gene_met_dat %>% arrange(desc(super_class_metabolon), desc(metabolite_name))
gene_met_dat$super_class_metabolon <- factor(gene_met_dat$super_class_metabolon, levels=names(classcol))
metabolite_name_vec <- unique(gene_met_dat$metabolite_name)
gene_met_dat$metabolite_name <- factor(gene_met_dat$metabolite_name, levels=metabolite_name_vec)

tmp <- as.data.frame(str_split(gene_met_dat$metabolite_name,"/",simplify=T)) %>% mutate(V1 = capitalize(V1)) 
gene_met_dat$anno_metabolite_name <- tmp$V1

gene_met_ncopy_1 <- gene_met_dat %>% filter(gene=="APOE4_ncopy_1")
gene_met_ncopy_2 <- gene_met_dat %>% filter(gene=="APOE4_ncopy_2")

table(gene_met_ncopy_1$anno_metabolite_name == gene_met_ncopy_2$anno_metabolite_name)

### Annotation

class_met <- ggplot(gene_met_ncopy_1, aes(x=metabolite_name, y=1)) +
  geom_tile(aes(fill=super_class_metabolon)) +
  labs(fill="Superclass of metabolites") +
  scale_fill_manual(values=classcol) +
  scale_x_discrete(limits=rev) +
  ylab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title =element_blank(),
        legend.position = "right",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14))

pdf("results/extended/extended_fig1b_legend1.pdf", width = 25, height = 3, onefile = F)
legend_apoe4_met <- cowplot::get_plot_component(class_met, 'guide-box', return_all = TRUE)
grid.newpage()
grid.draw(legend_apoe4_met[[1]])
dev.off()

### Manhattan plot: Heterozygotes

mhtplot_heter <- ggplot(gene_met_ncopy_1, aes(x=metabolite_name, y=-log10(p_var), fill=sign)) +
  geom_point(aes(fill=sign, size=abs), shape=21) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  scale_fill_manual(values=c("#F1C40F","#8E44AD")) +
  scale_x_discrete(limits=rev) +
  labs(fill="Direction of effect",size="Effect size") +
  ylab(expression("-log"[10]*"(P)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +
  annotation_custom(grobTree(textGrob(parse(text=("italic('APOE4')~heterozygosity")), x=0.28,  y=0.92, hjust=1,
                                      gp=gpar(fontsize=20, fontface="bold")))) +
  theme(legend.position = "right", legend.title = element_text(size=14), legend.text = element_text(size=12), legend.key.size = unit(0.6, 'cm')) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

pdf("results/extended/extended_fig1b_legend2.pdf", width = 25, height = 3, onefile = F)
legend_mhtplot_heter <- cowplot::get_plot_component(mhtplot_heter, 'guide-box', return_all = TRUE)
grid.newpage()
grid.draw(legend_mhtplot_heter[[1]])
dev.off()

### Manhattan plot: Homozygotes

mhtplot_homo <- ggplot(gene_met_ncopy_2, aes(x=metabolite_name, y=-log10(p_var), fill=sign)) +
  geom_point(aes(fill=sign, size=abs), shape=21) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  scale_fill_manual(values=c("#F1C40F","#8E44AD")) +
  scale_x_discrete(limits=rev) +
  labs(fill="Direction of effect",size="Effect size") +
  ylab(expression("-log"[10]*"(P)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +
  annotation_custom(grobTree(textGrob(parse(text=("italic('APOE4')~homozygosity")), x=0.28,  y=0.92, hjust=1,
                                      gp=gpar(fontsize=20, fontface="bold")))) +
  theme(legend.position = "right", legend.title = element_text(size=14), legend.text = element_text(size=12), legend.key.size = unit(0.6, 'cm')) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

pdf("results/extended/extended_fig1b.pdf", width = 11, height = 8, onefile = F)
apoe4_met_mhtplots <- ggpubr::ggarrange(mhtplot_heter, class_met, mhtplot_homo, ncol=1, nrow=3,
                                        widths = c(30), heights=c(1,0.1,1),
                                        align = "v",
                                        legend = "none")
apoe4_met_mhtplots
dev.off()