library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(gridExtra)

##-----------------------------------------
## Supp Fig. 2: Forest for APOE4
##-----------------------------------------

rm(list=ls())

### Load data

load("data/merged_data_baseline_nhs_2023.Rdata")
load("results/met_dem_cox_baseline_nhs_2023_filtered_res.RData")
load("data/genetic/ad_variants_dictionary.RData")

### Forest plot

nhs_strat_res <- nhs_strat_res %>%
  mutate(anno_model = case_when(model=="APOE4_noncarrier"~"Noncarrier",
                                model=="APOE4_ncopy_1"~"Heterozygote",
                                model=="APOE4_ncopy_2"~"Homozygote")) %>%
  mutate(anno_model = factor(anno_model, levels=c("Noncarrier","Heterozygote","Homozygote")))

APOE_met_name <- na.omit(nhs_inter_apoe4_res[nhs_inter_apoe4_res$two_inter_fdr_pvalue_p_inter_met_gene_1<0.05 |
                                               nhs_inter_apoe4_res$two_inter_fdr_pvalue_p_inter_met_gene_2<0.05,]$hmdb_id)

plot_list <- list()
dat_list <- list(); i <- 1

for (tmp_met_name in APOE_met_name) {
  
  print(tmp_met_name)
  print(fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$metabolite_name)
  
  tmp_strat <- merge(sig_met_dat, nhs_strat_res, by="hmdb_id", all.x=T, suffixes=c("",".y")) %>%
    filter(hmdb_id == tmp_met_name,
           grepl("APOE4",model),
           model != "APOE4_carrier") %>%
    mutate(hr = exp(beta),
           lower = exp(beta - 1.96*se),
           upper = exp(beta + 1.96*se)) %>%
    mutate(fdr_inter = case_when(model == "APOE4_noncarrier" ~ NA,
                                 model == "APOE4_ncopy_1" ~ round(nhs_inter_apoe4_res[nhs_inter_apoe4_res$hmdb_id == tmp_met_name,]$two_inter_fdr_pvalue_p_inter_met_gene_1,2),
                                 model == "APOE4_ncopy_2" ~ round(nhs_inter_apoe4_res[nhs_inter_apoe4_res$hmdb_id == tmp_met_name,]$two_inter_fdr_pvalue_p_inter_met_gene_2,2)))
  
  dat_list[[i]] <- tmp_strat %>% subset(select=c(anno_metabolite_name_sig,anno_model,hr,lower,upper)); i <- i + 1
  p <- ggplot(data=tmp_strat, aes(y=anno_model, x=hr, xmin=lower, xmax=upper)) +
    geom_linerange(size=0.7,position=position_dodge(width = 0.5), alpha=0.8) +
    geom_point(aes(fill=anno_model),color="white", shape=21, size=3, stroke = 0.5, position=position_dodge(width = 0.5)) +
    labs(title=paste0(unique(tmp_strat$anno_metabolite_name_sig),"\nby APOE4"), x='Hazard ratio (95% CI)', y='') +  
    geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
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

pdf("results/supp/supp_fig_2.pdf", width = 14, height = 20) 
grid.arrange(grobs = plot_list, ncol = 5, padding = unit(c(4,4), "cm")) 
dev.off()
