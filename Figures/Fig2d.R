library(dplyr)
library(data.table)
library(ggplot2)
library(egg)

##-----------------------------------------
## Fig. 2d: Dementia vs. TICS
##-----------------------------------------

rm(list=ls())

### Data

load("results/met_dem_cox_baseline_nhs_2023_filtered_res.RData")
dementia_inter <- fread("results/met_dem_cox_baseline_interaction_nhs_2023.txt")
tics_inter <- fread("results/met_tics_interaction.txt")

sig_met_APOE4_inter_dat <- sig_met_dat %>%
  filter(sig_inter_apoe4==1)

dementia_inter_sig <- merge(subset(nhs_main_res, select=c(metabolite_name,super_class_metabolon)), subset(dementia_inter, select=-super_class_metabolon), by="metabolite_name") %>%
  filter(hmdb_id %in% sig_met_APOE4_inter_dat$hmdb_id,
         covariates %in% c("APOE4_carrier"))

tics_inter_sig <- merge(subset(nhs_main_res, select=c(metabolite_name,super_class_metabolon)), subset(tics_inter, select=-super_class_metabolon), by="metabolite_name") %>%
  filter(hmdb_id %in% sig_met_APOE4_inter_dat$hmdb_id,
         covariates %in% c("APOE4_carrier"))

### Re-do FDR on TICS

tics_inter_sig <- tics_inter_sig %>%
  group_by(model, outcome) %>%
  mutate(p_fdr_inter_met_gene_1 = p.adjust(p_inter_met_gene_1, "BH"),
         p_fdr_inter_met_gene_2 = p.adjust(p_inter_met_gene_2, "BH"),
         p_fdr_LRT = p.adjust(p_LRT, "BH")) %>%
  pivot_longer(cols = c(p_inter_met_gene_1,p_inter_met_gene_2), names_to = "pvalue_column", values_to = "pvalue") %>%
  mutate(two_inter_fdr = p.adjust(pvalue, "BH")) %>%
  pivot_wider(names_from = "pvalue_column", values_from = c(pvalue, two_inter_fdr))
summary(tics_inter_sig$two_inter_fdr_p_inter_met_gene_1)
summary(tics_inter_sig$two_inter_fdr_p_inter_met_gene_2)

### Correlation plots - Interaction with APOE4 carrier

dementia_inter_sig_carrier <- dementia_inter_sig %>%
  filter(covariates == "APOE4_carrier",
         subtype=="all") %>%
  filter(metabolite_name != "C40:6 PS") %>%
  subset(select=c(outcome,hmdb_id,metabolite_name,super_class_metabolon,beta_inter_met_gene_2,se_inter_met_gene_2,p_fdr_inter_met_gene_2))
colnames(dementia_inter_sig_carrier) <- c("outcome_dementia","hmdb_id","metabolite_name","super_class_metabolon","beta_dementia","se_dementia","p_fdr_dementia")

tics_inter_sig_carrier <- tics_inter_sig %>%
  filter(covariates == "APOE4_carrier") %>%
  filter(metabolite_name != "C40:6 PS") %>%
  subset(select=c(outcome,hmdb_id,metabolite_name,super_class_metabolon,beta_inter_met_gene_2,se_inter_met_gene_2,p_fdr_inter_met_gene_2))
colnames(tics_inter_sig_carrier) <- c("outcome_tics","hmdb_id","metabolite_name","super_class_metabolon","beta_tics","se_tics","p_fdr_tics")

inter_sig_carrier <- merge(dementia_inter_sig_carrier,tics_inter_sig_carrier,by=c("hmdb_id","metabolite_name","super_class_metabolon"),all.x=T)

table(inter_sig_carrier$outcome_tics)

tmp_inter_sig_carrier <- inter_sig_carrier[inter_sig_carrier$outcome_tics=="av_ntotal",] %>% 
  mutate(flag = ifelse(beta_dementia*beta_tics<0,TRUE,FALSE))
tmp_r2 <- format(round(cor(tmp_inter_sig_carrier$beta_dementia, tmp_inter_sig_carrier$beta_tics)[1],2), nsmall = 2)
p1 <- ggplot(data=tmp_inter_sig_carrier, aes(x=beta_tics, y=beta_dementia)) +
  geom_point(aes(fill=super_class_metabolon), size=2, shape=21, stroke=0.2, alpha=0.8) +
  scale_fill_manual(values=classcol) +
  geom_vline(xintercept=0, linetype='dashed', alpha=.5) +
  geom_hline(yintercept=0, linetype='dashed', alpha=.5) +
  labs(x=expression(paste(atop("Global cognition",beta[" metabolites \u00D7" * italic(' APOE4 ') * ""]))),y=expression(paste(atop("Dementia",beta[" metabolites \u00D7" * italic(' APOE4 ') * ""])))) +
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

tmp_inter_sig_carrier <- inter_sig_carrier[inter_sig_carrier$outcome_tics=="av_nverbl",] %>% 
  mutate(flag = ifelse(beta_dementia*beta_tics<0,TRUE,FALSE))
tmp_r2 <- format(round(cor(tmp_inter_sig_carrier$beta_dementia, tmp_inter_sig_carrier$beta_tics)[1],2), nsmall = 2)
p2 <- ggplot(data=tmp_inter_sig_carrier, aes(x=beta_tics, y=beta_dementia)) +
  geom_point(aes(fill=super_class_metabolon), size=2, shape=21, stroke=0.2, alpha=0.8) +
  scale_fill_manual(values=classcol) +
  geom_vline(xintercept=0, linetype='dashed', alpha=.5) +
  geom_hline(yintercept=0, linetype='dashed', alpha=.5) +
  labs(x=expression(paste(atop("Verbal memory",beta[" metabolites \u00D7" * italic(' APOE4 ') * ""]))),y="") +
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

tmp_inter_sig_carrier <- inter_sig_carrier[inter_sig_carrier$outcome_tics=="av_zscre",] %>% 
  mutate(flag = ifelse(beta_dementia*beta_tics<0,TRUE,FALSE))
tmp_r2 <- format(round(cor(tmp_inter_sig_carrier$beta_dementia, tmp_inter_sig_carrier$beta_tics)[1],2), nsmall = 2)
p3 <- ggplot(data=tmp_inter_sig_carrier, aes(x=beta_tics, y=beta_dementia)) +
  geom_point(aes(fill=super_class_metabolon), size=2, shape=21, stroke=0.2, alpha=0.8) +
  scale_fill_manual(values=classcol) +
  geom_vline(xintercept=0, linetype='dashed', alpha=.5) +
  geom_hline(yintercept=0, linetype='dashed', alpha=.5) +
  labs(x=expression(paste(atop("TICS score",beta[" metabolites \u00D7" * italic(' APOE4 ') * ""]))),y="") +
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

pdf("results/figures/fig2d_dementia_tics_carrier.pdf", width = 6.8, height = 2.4, onefile = F)
fig2d <- egg::ggarrange(p1, ggplot() + theme_void(), p2, ggplot() + theme_void(), p3,
                        nrow=1, ncol=5,
                        widths = c(1,-0.1,1,-0.1,1), heights=c(15))
fig2d
dev.off()
