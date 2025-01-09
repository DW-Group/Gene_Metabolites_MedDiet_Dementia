library(dplyr)
library(data.table)
library(ggplot2)

##-----------------------------------------
## Fig. 3b: AMED-dementia by APOE4
##-----------------------------------------

rm(list=ls())

### Load data

subgroup_res <- read.csv("results/amed_APOE4_PRS_subgroup_nhsfull_baseline.csv") %>%
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

inter_res <- read.csv("results/amed_APOE4_interaction_nhsfull_baseline.csv") %>%
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
  labs(title="Dementia risk and MedDiet index\nscore by APOE4", x='HR (95% CI) of dementia risk', y='') +
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

pdf("results/figures/fig3b_dementia_AMED_by_APOE4_11202024.pdf", width = 3, height = 2, onefile = F)
p_forest
dev.off()
