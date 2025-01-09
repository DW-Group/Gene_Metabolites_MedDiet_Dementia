library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

##-----------------------------------------
## Supp Fig. 4: Forest for ind var
##-----------------------------------------

rm(list=ls())

### Load data

load("data/merged_data_baseline_nhs_2023_10072024.Rdata")
load("results/met_dem_cox_baseline_nhs_2023_filtered_res_11112024.RData")
load("data/genetic/ad_variants_dictionary.RData")
ad_dict <- anno_dat[anno_dat$name %in% nhs_inter_ind_var_res$variant,]

### Forest plot

sig_var <- as.vector(na.omit(unique(nhs_inter_ind_var_res[nhs_inter_ind_var_res$p_fdr_inter_met_var < 0.05,]$variant)))
ind_var_sig <- merge(sig_met_dat, nhs_inter_ind_var_res, by="hmdb_id", all.x=T, suffixes=c("",".y")) %>%
  filter(p_fdr_inter_met_var < 0.05,
         !(hmdb_id %in% c("HMDB0000510","HMDB0011310")))

strata_ind_var_res_dat <- fread("results/met_dem_cox_baseline_stratified_individual_variants_nhs_2023_10102024.txt")

plot_list <- list()
dat_list <- list(); i <- 1

for (i in 1:dim(ind_var_sig)[1]) {
  
  tmp_met_name <- ind_var_sig$hmdb_id[i]
  tmp_var_name <- ind_var_sig$variant[i]
  
  tmp_strat_ind_var <- merge(sig_met_dat, strata_ind_var_res_dat, by="hmdb_id", all.x=T, suffixes=c("",".y")) %>%
    filter(hmdb_id == tmp_met_name,
           variant == tmp_var_name,
           !is.na(beta)) %>%
    mutate(group = case_when(subgroup == "Noncarrier" ~ paste0(alt,alt),
                             subgroup == "Heterozygote" ~ paste0(ref,alt),
                             subgroup == "Homozygote" ~ paste0(ref,ref))) %>%
    mutate(beta = as.numeric(beta),
           se = as.numeric(se),
           hr = exp(beta),
           lower = exp(beta - 1.96*se),
           upper = exp(beta + 1.96*se)) %>%
    mutate(group = factor(group, levels=c(unique(paste0(alt,alt)),unique(paste0(ref,alt)),unique(paste0(ref,ref)))))
  tmp_p_inter <- round(nhs_inter_ind_var_res[nhs_inter_ind_var_res$hmdb_id == tmp_met_name & nhs_inter_ind_var_res$variant == tmp_var_name,]$p_fdr_inter_met_var,2)
  
  dat_list[[i]] <- tmp_strat_ind_var %>% subset(select=c(anno_metabolite_name_sig,rsid,gene,group,hr,lower,upper)); i <- i + 1
  p <- ggplot(data=tmp_strat_ind_var, aes(y=group, x=hr, xmin=lower, xmax=upper)) +
    geom_linerange(size=0.5,position=position_dodge(width = 0.5)) +
    geom_point(aes(fill=group),color="white", shape=21, size=3, stroke = 0.5,position=position_dodge(width = 0.5)) +
    labs(title=paste0(unique(tmp_strat_ind_var$anno_metabolite_name_sig),"\nby ",unique(tmp_strat_ind_var$rsid)," (",unique(tmp_strat_ind_var$gene),")"), x='Hazard ratio (95% CI)', y='') +
    geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values=c("tan","tan3","tan4"), drop=FALSE)+
    theme_bw()+
    theme(legend.title= element_blank(),
          plot.title = element_text(size=10))+
    theme(panel.border = element_rect(color="black", size=0.8),
          axis.title.x = element_text(size=10),
          axis.text.y = element_text(colour = "black", size=8),
          axis.text.x = element_text(colour = "black", size=8),
          axis.ticks.y = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black")) +
    theme(legend.position="none")  
  
  plot_list[[length(plot_list) + 1]] <- p
}

source_dat <- as.data.frame(do.call("rbind",dat_list))

pdf("results/supp/supp_fig_4.pdf", width = 12 , height = 4) 
grid.arrange(grobs = plot_list, ncol = 5, padding = unit(c(4,4), "cm")) 
dev.off()
