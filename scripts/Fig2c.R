library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(grid)
library(DescTools)
library(tidyr)

##-----------------------------------------
## Fig. 2c: Call-out metabolites
##-----------------------------------------

rm(list=ls())

### Data

load("data/merged_data_baseline_nhs_2023_10072024.Rdata")
load("results/met_dem_cox_baseline_nhs_2023_filtered_res_10072024.RData")
load("data/genetic/ad_variants_dictionary.RData")
ad_dict <- anno_dat[anno_dat$name %in% nhs_inter_ind_var_res$variant,]

nhs_strat_res <- nhs_strat_res %>%
  mutate(anno_model = case_when(model=="APOE4_noncarrier"~"Noncarrier",
                                model=="APOE4_ncopy_1"~"Heterozygote",
                                model=="APOE4_ncopy_2"~"Homozygote")) %>%
  mutate(anno_model = factor(anno_model, levels=c("Noncarrier","Heterozygote","Homozygote")))

### Forest plot - APOE4

tmp_met_name <- "HMDB0000043"; tmp_met_name; fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$metabolite_name
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
p1 <- ggplot(data=tmp_strat, aes(y=anno_model, x=hr, xmin=lower, xmax=upper)) +
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
  theme(legend.position="none") +
  geom_text(aes(label=ifelse(is.na(fdr_inter), "", paste0("FDR[~interaction]~'='~",fdr_inter))), 
            position=position_dodge(width = 0.9), 
            parse=T,
            hjust=0,
            vjust=-1.2, size=3) 

tmp_met_name <- "HMDB0010368"; tmp_met_name; fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$metabolite_name
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
p2 <- ggplot(data=tmp_strat, aes(y=anno_model, x=hr, xmin=lower, xmax=upper)) +
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
        axis.text.y = element_blank(),    
        axis.text.x = element_text(colour = "black", size=8),
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black")) +
  theme(legend.position="none") +
  geom_text(aes(label=ifelse(is.na(fdr_inter), "", paste0("FDR[~interaction]~'='~",fdr_inter))), 
            position=position_dodge(width = 0.9), 
            parse=T,
            hjust=0,
            vjust=-1.2, size=3) 

tmp_met_name <- "HMDB0042093"; tmp_met_name; fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$metabolite_name
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
p3 <- ggplot(data=tmp_strat, aes(y=anno_model, x=hr, xmin=lower, xmax=upper)) +
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
        axis.text.y = element_blank(),       
        axis.text.x = element_text(colour = "black", size=8),
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black")) +
  theme(legend.position="none") +
  geom_text(aes(label=ifelse(is.na(fdr_inter), "", paste0("FDR[~interaction]~'='~",fdr_inter))), 
            position=position_dodge(width = 0.9), 
            parse=T,
            vjust=-1.2, size=3) 

### Forest plot - other variants

load("data/merged_data_baseline_nhs_2023_10072024.Rdata")
source("scripts/source_MWAS_covar_list.R")
sig_var <- as.vector(na.omit(unique(nhs_inter_ind_var_res[nhs_inter_ind_var_res$p_fdr_inter_met_var < 0.05,]$variant)))
ind_var_sig <- merge(sig_met_dat, nhs_inter_ind_var_res, by="hmdb_id", all.x=T, suffixes=c("",".y")) %>%
  filter(p_fdr_inter_met_var < 0.05)

# strata_ind_var_res_list <- list();  i <- 1
# 
# for (tmp_var in sig_var){
# 
#   print(tmp_var); print(mean(merged_data[[tmp_var]], na.rm=T)/2)
#   merged_data$cat_geno <- round(merged_data[[tmp_var]])
#   print(mean(merged_data$cat_geno, na.rm=T)/2); print(table(merged_data$cat_geno))
# 
#   for (met in ind_var_sig[ind_var_sig$variant==tmp_var,]$hmdb_id){
# 
#     print(met)
# 
#     tmp_data <- merged_data[merged_data$cat_geno == 0,]; tmp_met <- tmp_data[[met]]
#     tmp_mod1 <- tryCatch(coxph(reformulate(c("tmp_met",covar_dem_list[["full_covar"]]), response="Surv(tad, adcase)"), ties="breslow", data = tmp_data),error=function(e) e,warning=function(w) w)
# 
#     tmp_data <- merged_data[merged_data$cat_geno == 1,]; tmp_met <- tmp_data[[met]]
#     tmp_mod2 <- tryCatch(coxph(reformulate(c("tmp_met",covar_dem_list[["full_covar"]]), response="Surv(tad, adcase)"), ties="breslow", data = tmp_data),error=function(e) e,warning=function(w) w)
# 
#     tmp_data <- merged_data[merged_data$cat_geno == 2,]; tmp_met <- tmp_data[[met]]
#     tmp_mod3 <- tryCatch(coxph(reformulate(c("tmp_met",covar_dem_list[["full_covar"]]), response="Surv(tad, adcase)"), ties="breslow", data = tmp_data),error=function(e) e,warning=function(w) w)
# 
#     strata_ind_var_res_list[[i]] <- c(met, tmp_var, "Noncarrier", if(length(tmp_mod1) > 5) c(summary(tmp_mod1)$n, nobs(tmp_mod1), summary(tmp_mod1)$coefficients["tmp_met",c(1,3,5)]) else c(tmp_mod1$message,NA,NA,NA,NA)); i <- i + 1
#     strata_ind_var_res_list[[i]] <- c(met, tmp_var, "Heterozygote", if(length(tmp_mod2) > 5) c(summary(tmp_mod2)$n, nobs(tmp_mod2), summary(tmp_mod2)$coefficients["tmp_met",c(1,3,5)]) else c(tmp_mod2$message,NA,NA,NA,NA)); i <- i + 1
#     strata_ind_var_res_list[[i]] <- c(met, tmp_var, "Homozygote", if(length(tmp_mod3) > 5) c(summary(tmp_mod3)$n, nobs(tmp_mod3), summary(tmp_mod3)$coefficients["tmp_met",c(1,3,5)]) else c(tmp_mod3$message,NA,NA,NA,NA)); i <- i + 1
# 
#   }
# }
# 
# strata_ind_var_res_dat <- as.data.frame(do.call("rbind", strata_ind_var_res_list))
# colnames(strata_ind_var_res_dat) <- c("hmdb_id", "variant","subgroup","n","n_events","beta","se","p")
# strata_ind_var_res_dat <- merge(strata_ind_var_res_dat, ad_dict, by.x="variant",by.y="name")
# write.table(strata_ind_var_res_dat, "results/met_dem_cox_baseline_stratified_individual_variants_nhs_2023_10102024.txt", row.names=F, quote=F, sep="\t")

strata_ind_var_res_dat <- fread("results/met_dem_cox_baseline_stratified_individual_variants_nhs_2023_10102024.txt")

tmp_met_name <- "HMDB0000898"
tmp_var_name <- "9:107665978_C(/G)"
fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$metabolite_name; ad_dict[ad_dict$name==tmp_var_name,]$gene
table(round(merged_data[[tmp_var_name]]))
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
p4 <- ggplot(data=tmp_strat_ind_var, aes(y=group, x=hr, xmin=lower, xmax=upper)) +
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
  theme(legend.position="none") +
  geom_text_npc(npcx = "left", npcy = "top", label = paste0("FDR[~interaction]~`=`~", tmp_p_inter), parse = T, size=3)

tmp_met_name <- "HMDB0240212"
tmp_var_name <- "21:27473875_C(/T)"
fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$metabolite_name; ad_dict[ad_dict$name==tmp_var_name,]$gene
table(round(merged_data[[tmp_var_name]]))
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
p5 <- ggplot(data=tmp_strat_ind_var, aes(y=group, x=hr, xmin=lower, xmax=upper)) +
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
  theme(legend.position="none") +
  geom_text(label=paste0("FDR[~interaction]~`=`~", tmp_p_inter), x=Inf, y=Inf, hjust=1.1, vjust=2, parse = T, size=3)

tmp_met_name <- "HMDB0007874"
tmp_var_name <- "21:28148191_C(/T)"
fdata_filtered[fdata_filtered$hmdb_id==tmp_met_name,]$metabolite_name; ad_dict[ad_dict$name==tmp_var_name,]$gene
table(round(merged_data[[tmp_var_name]]))
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
p6 <- ggplot(data=tmp_strat_ind_var, aes(y=group, x=hr, xmin=lower, xmax=upper)) +
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
  theme(legend.position="none") +
geom_text_npc(npcx = "left", npcy = "top", label = paste0("FDR[~interaction]~`=`~", tmp_p_inter), parse = T, size=3)

pdf("results/figures/fig2c_call_out_11132024.pdf", width = 8, height = 4.5, onefile = F)
fig2c <- egg::ggarrange(p1, ggplot() + theme_void(), p2, ggplot() + theme_void(), p3,
                        p4, ggplot() + theme_void(), p5, ggplot() + theme_void(), p6,
                                  nrow=2, ncol=5,
                                  widths = c(1,-2,1,-2,1), heights=c(1,1))
fig2c
dev.off()
