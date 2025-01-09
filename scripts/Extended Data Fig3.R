library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)

##-----------------------------------------
## Exd Fig. 3: Significant interaction by superclass  
##-----------------------------------------

##----------------------------------
## 3a and 3c: TICS and AMED by APOE4 and PRS
##----------------------------------

rm(list=ls())

### Load data

amed_tics <- fread("results/supp/amed_tics_by_gene_11252024.txt") %>%
  filter(!grepl("all",gene)) %>%
  mutate(across(c(beta,se), as.numeric)) %>%
  mutate(lower=beta-1.96*se,
         upper=beta+1.96*se) %>%
  mutate(anno_group = case_when(grepl("PGS",gene) & group==1 ~ "Tertile 1",
                                grepl("PGS",gene) & group==2 ~ "Tertile 2",
                                grepl("PGS",gene) & group==3 ~ "Tertile 3",
                                gene=="APOE4" & group==0 ~ "Noncarrier",
                                gene=="APOE4" & group==1 ~ "Heterozygote",
                                gene=="APOE4" & group==2 ~ "Homozygote"),
         anno_outcome = case_when(tics=="av_ntotal" ~ "Global cognition",
                                  tics=="av_nverbl" ~ "Verbal memory",
                                  tics=="av_zscre" ~ "TICS score")) %>%
  mutate(anno_group = factor(anno_group, levels=c("Noncarrier","Heterozygote","Homozygote","Tertile 1","Tertile 2","Tertile 3")),
         anno_outcome = factor(anno_outcome, levels=c("Global cognition","Verbal memory","TICS score")))

APOE4_res <- amed_tics %>% filter(grepl("APOE4",gene))
PGS002280_res <- amed_tics %>% filter(grepl("PGS002280",gene))
PGS000334_res <- amed_tics %>% filter(grepl("PGS000334",gene))

### Forest plot - APOE4

plot_list <- list()

for (outcome in unique(APOE4_res$anno_outcome)) {
  
  tmp_strat <- APOE4_res %>%
    filter(anno_outcome == outcome)
  
  p <- ggplot(data=tmp_strat, aes(y=anno_group, x=beta, xmin=lower, xmax=upper)) +
    geom_linerange(size=0.7,position=position_dodge(width = 0.5), alpha=0.8) +
    geom_point(aes(fill=anno_group),color="white", shape=21, size=3, stroke = 0.5, position=position_dodge(width = 0.5)) +
    labs(title=paste0(outcome," and MedDiet index\nscore by APOE4"), x='Effect size (95% CI)', y='') +  
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values=c("tan","tan3","tan4"))+
    theme_bw()+
    theme(legend.title= element_blank(),
          plot.title = element_text(size = 10))+
    theme(panel.border = element_rect(color="black", size=0.8),
          axis.title.x = element_text(size=10),
          axis.text.y = element_text(colour = "black"),    
          axis.text.x = element_text(colour = "black", size=8),
          axis.ticks.y = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black")) +
    theme(legend.position="none")
  
  plot_list[[length(plot_list) + 1]] <- p
}

pdf("results/extended/extended_fig3a.pdf", width = 10, height = 2)  # Adjust width and height as needed
grid.arrange(grobs = plot_list, ncol = 3, padding = unit(c(4,4), "cm")) 
dev.off()

### Forest plot - PGS002280

plot_list <- list()

for (outcome in unique(PGS002280_res$anno_outcome)) {
  
  tmp_strat <- PGS002280_res %>%
    filter(anno_outcome == outcome)
  
  p <- ggplot(data=tmp_strat, aes(y=anno_group, x=beta, xmin=lower, xmax=upper)) +
    geom_linerange(size=0.7,position=position_dodge(width = 0.5), alpha=0.8) +
    geom_point(aes(fill=anno_group),color="white", shape=21, size=3, stroke = 0.5, position=position_dodge(width = 0.5)) +
    labs(title=paste0(outcome," and MedDiet index\nscore by PRS (excluding APOE region)"), x='Effect size (95% CI)', y='') +  
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values=c("tan","tan3","tan4"))+
    theme_bw()+
    theme(legend.title= element_blank(),
          plot.title = element_text(size = 10))+
    theme(panel.border = element_rect(color="black", size=0.8),
          axis.title.x = element_text(size=10),
          axis.text.y = element_text(colour = "black"),    
          axis.text.x = element_text(colour = "black", size=8),
          axis.ticks.y = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black")) +
    theme(legend.position="none")
  
  plot_list[[length(plot_list) + 1]] <- p
}

pdf("results/extended/extended_fig3c_1.pdf", width = 10, height = 2)  
grid.arrange(grobs = plot_list, ncol = 3, padding = unit(c(4,4), "cm")) 
dev.off()

### Forest plot - PGS000334

plot_list <- list()

for (outcome in unique(PGS000334_res$anno_outcome)) {
  
  tmp_strat <- PGS000334_res %>%
    filter(anno_outcome == outcome)
  
  p <- ggplot(data=tmp_strat, aes(y=anno_group, x=beta, xmin=lower, xmax=upper)) +
    geom_linerange(size=0.7,position=position_dodge(width = 0.5), alpha=0.8) +
    geom_point(aes(fill=anno_group),color="white", shape=21, size=3, stroke = 0.5, position=position_dodge(width = 0.5)) +
    labs(title=paste0(outcome," and MedDiet index\nscore by PRS (including APOE variants)"), x='Effect size (95% CI)', y='') +  
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values=c("tan","tan3","tan4"))+
    theme_bw()+
    theme(legend.title= element_blank(),
          plot.title = element_text(size = 10))+
    theme(panel.border = element_rect(color="black", size=0.8),
          axis.title.x = element_text(size=10),
          axis.text.y = element_text(colour = "black"),    
          axis.text.x = element_text(colour = "black", size=8),
          axis.ticks.y = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black")) +
    theme(legend.position="none")
  
  plot_list[[length(plot_list) + 1]] <- p
}

pdf("results/extended/extended_fig3c_2.pdf", width = 10, height = 2)  
grid.arrange(grobs = plot_list, ncol = 3, padding = unit(c(4,4), "cm")) 
dev.off()

##----------------------------------
## 3b: dementia and AMED by PRS
##----------------------------------

rm(list=ls())

### Load data

subgroup_res <- read.csv("results/amed_APOE4_PRS_subgroup_nhsfull_baseline.csv") %>%
  filter(variable=="amedcon",
         modelno==2,
         grepl("PGS",strata)) %>%
  mutate(across(c(beta,see), as.numeric)) %>%
  mutate(hr=exp(beta),
         lower=exp(beta-1.96*see),
         upper=exp(beta+1.96*see)) %>%
  mutate(anno_group = case_when(strata %in% c("PGS0003340","PGS0022800") ~ "Tertile 1",
                                strata %in% c("PGS0003341","PGS0022801") ~ "Tertile 2",
                                strata %in% c("PGS0003342","PGS0022802") ~ "Tertile 3")) %>%
  mutate(anno_group = factor(anno_group, levels=c("Tertile 1","Tertile 2","Tertile 3")))

PGS002280_res <- subgroup_res %>% filter(grepl("PGS002280",strata))
PGS000334_res <- subgroup_res %>% filter(grepl("PGS000334",strata))

### Forest plot

p1 <- ggplot(data=PGS002280_res, aes(y=anno_group, x=hr, xmin=lower, xmax=upper)) +
  geom_linerange(size=0.5,position=position_dodge(width = 0.5), alpha=0.8) +
  geom_point(aes(fill=anno_group),color="white", shape=21, size=2.5, stroke = 0.5, position=position_dodge(width = 0.5)) +
  labs(title="Dementia risk and MedDiet index score\nby PRS (excluding APOE region)", x='HR (95% CI) of dementia risk', y='') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  scale_y_discrete(limits=rev) +
  scale_fill_manual(values=c("#D4E4B8","#7A8D50","#4F5937"))+
  theme_bw() +
  theme(legend.title= element_blank(),
        plot.title = element_text(colour = "black",size=9))+
  theme(axis.title.x = element_text(colour = "black", size=8),
        axis.text.y = element_text(colour = "black", size=8),
        axis.text.x = element_text(colour = "black", size=8)
  ) +
  theme(legend.position="none") 

p2 <- ggplot(data=PGS000334_res, aes(y=anno_group, x=hr, xmin=lower, xmax=upper)) +
  geom_linerange(size=0.5,position=position_dodge(width = 0.5), alpha=0.8) +
  geom_point(aes(fill=anno_group),color="white", shape=21, size=2.5, stroke = 0.5, position=position_dodge(width = 0.5)) +
  labs(title="Dementia risk and MedDiet index score\nby PRS (including APOE variants)", x='HR (95% CI) of dementia risk', y='') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  scale_y_discrete(limits=rev) +
  scale_fill_manual(values=c("#D4E4B8","#7A8D50","#4F5937"))+
  theme_bw() +
  theme(legend.title= element_blank(),
        plot.title = element_text(colour = "black",size=9))+
  theme(axis.title.x = element_text(colour = "black", size=8),
        axis.text.y = element_text(colour = "black", size=8),
        axis.text.x = element_text(colour = "black", size=8)
  ) +
  theme(legend.position="none") 

pdf("results/extended/extended_fig3b.pdf", width = 6, height = 2, onefile = F)
tmp_fig <- ggpubr::ggarrange(p1, p2, ncol=2, nrow=1,
                                widths = c(30), heights=c(1,1),
                                align = "v")
tmp_fig
dev.off()