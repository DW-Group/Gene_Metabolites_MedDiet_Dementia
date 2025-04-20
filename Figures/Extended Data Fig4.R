library(patchwork)
library(ggplot2)
library(ggpp)
library(DescTools)
  
##-----------------------------------------
## Exd Fig. 4: Comparison of NHS and HPFS interaction results
##-----------------------------------------

rm(list=ls())

source("scripts/source_color_superclass.R")

load("results/met_dem_cox_baseline_nhs_2023_filtered_res.RData")
dementia_inter <- fread("results/met_dem_cox_baseline_interaction_nhs_2023.txt")
hpfs_res_dat <- fread("results/met_dem_cox_baseline_interaction_hpfs_2023.txt")

nhs_res_dat <- merge(subset(nhs_main_res, select=c(metabolite_name,super_class_metabolon)), subset(dementia_inter, select=-super_class_metabolon), by="metabolite_name") %>%
  filter(covariates %in% c("APOE4_carrier"))

met_overlap <- intersect(nhs_res_dat$hmdb_id, hpfs_res_dat$hmdb_id); length(met_overlap)
sig_met_APOE4_inter_dat <- sig_met_dat %>% filter(sig_inter_apoe4==1)
met_overlap <- met_overlap[met_overlap %in% sig_met_APOE4_inter_dat$hmdb_id]; length(met_overlap)

nhs_res_sub <- nhs_res_dat %>% filter(hmdb_id %in% met_overlap)
hpfs_res_sub <- hpfs_res_dat %>% filter(hmdb_id %in% met_overlap, covariates %in% c("APOE4_carrier"))

nhs_res_sub <- nhs_res_sub %>%
  subset(select=c(outcome,hmdb_id,metabolite_name,super_class_metabolon,beta_inter_met_gene_2,se_inter_met_gene_2,p_fdr_inter_met_gene_2))
colnames(nhs_res_sub) <- c("outcome_nhs","hmdb_id","metabolite_name","super_class_metabolon","beta_nhs","se_nhs","p_fdr_nhs")

hpfs_res_sub <- hpfs_res_sub %>%
  subset(select=c(outcome,hmdb_id,beta_inter_met_gene_2,se_inter_met_gene_2,p_fdr_inter_met_gene_2))
colnames(hpfs_res_sub) <- c("outcome_hpfs","hmdb_id","beta_hpfs","se_hpfs","p_fdr_hpfs")

inter_sig_carrier <- merge(nhs_res_sub,hpfs_res_sub,by=c("hmdb_id"),all.x=T,allow.cartesian=TRUE)

inter_sig_carrier <- inter_sig_carrier %>% mutate(flag = ifelse(beta_nhs*beta_hpfs<0,TRUE,FALSE)); table(inter_sig_carrier$flag)
tmp_r2 <- format(round(cor(inter_sig_carrier$beta_nhs, inter_sig_carrier$beta_hpfs, use="complete.obs", method="pearson")[1],2), nsmall = 2); tmp_r2

p_corr <- ggplot(data=inter_sig_carrier, aes(x=beta_nhs, y=beta_hpfs)) +
  geom_point(aes(fill=super_class_metabolon), size=2, shape=21, stroke=0.2, alpha=0.8) +
  scale_fill_manual(values=classcol) +
  geom_vline(xintercept=0, linetype='dashed', alpha=.5) +
  geom_hline(yintercept=0, linetype='dashed', alpha=.5) +
  labs(
    x = bquote(atop("Dementia (NHS)", beta["metabolites ×" * italic("APOE4") * ""])),
    y = bquote(atop("Dementia (HPFS)", beta["metabolites ×" * italic("APOE4") * ""]))
  ) +      
  theme_bw() +
  theme(axis.title = element_text(size=10),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),    
        axis.ticks.x = element_line(color="black"),
        axis.ticks.y = element_line(color="black"),    
        panel.border = element_rect(color="black", size=1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  geom_text_npc(npcx = "right", npcy = "top", label = paste0("r = ",tmp_r2), size=4) +
  theme(legend.position = "none")

pdf("results/extended/extended_fig4.pdf", width = 3, height = 3)
p_corr
dev.off()