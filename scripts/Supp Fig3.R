library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(ggpubr)

##-----------------------------------------
## Supp Fig. 3: Dementia case vs. death
##-----------------------------------------

rm(list=ls())

### Data

load("results/met_dem_cox_baseline_nhs_2023_filtered_res_10072024.RData")
sens_inter <- fread("results/met_dem_cox_baseline_interaction_nhs_2023_11252024_sensitivity.txt")
dementia_inter <- fread("results/met_dem_cox_baseline_interaction_nhs_2023_10072024.txt")

sig_met_APOE4_inter_dat <- sig_met_dat %>%
  filter(sig_inter_apoe4==1)

sens_inter_sig <- merge(subset(nhs_main_res, select=c(metabolite_name,super_class_metabolon)), subset(sens_inter, select=-super_class_metabolon), by="metabolite_name") %>%
  filter(hmdb_id %in% sig_met_APOE4_inter_dat$hmdb_id,
         covariates %in% c("APOE4_carrier"))

dementia_inter_sig <- merge(subset(nhs_main_res, select=c(metabolite_name,super_class_metabolon)), subset(dementia_inter, select=-super_class_metabolon), by="metabolite_name") %>%
  filter(hmdb_id %in% sig_met_APOE4_inter_dat$hmdb_id,
         covariates %in% c("APOE4_carrier"))

### Correlation plots - Interaction with APOE4 carrier

sens_inter_sig_carrier <- sens_inter_sig %>%
  subset(select=c(outcome,hmdb_id,metabolite_name,super_class_metabolon,beta_inter_met_gene_2))
sens_inter_sig_carrier_case <- sens_inter_sig_carrier %>% filter(outcome=="dementiacase")
sens_inter_sig_carrier_death <- sens_inter_sig_carrier %>% filter(outcome=="dementiadeath")
colnames(sens_inter_sig_carrier_case) <- c("outcome_sens","hmdb_id","metabolite_name","super_class_metabolon","beta_sens_case")
colnames(sens_inter_sig_carrier_death) <- c("outcome_sens","hmdb_id","metabolite_name","super_class_metabolon","beta_sens_death")

dementia_inter_sig_carrier <- dementia_inter_sig %>%
  filter(covariates == "APOE4_carrier",
         subtype=="all") %>%
  subset(select=c(outcome,hmdb_id,metabolite_name,super_class_metabolon,beta_inter_met_gene_2))
colnames(dementia_inter_sig_carrier) <- c("outcome_dementia","hmdb_id","metabolite_name","super_class_metabolon","beta_dementia")

inter_sig_carrier_case <- merge(sens_inter_sig_carrier_case,dementia_inter_sig_carrier,by=c("hmdb_id","metabolite_name","super_class_metabolon"),all.x=T)
inter_sig_carrier_death <- merge(sens_inter_sig_carrier_death,dementia_inter_sig_carrier,by=c("hmdb_id","metabolite_name","super_class_metabolon"),all.x=T)
inter_sig_carrier_both <- merge(sens_inter_sig_carrier_case,sens_inter_sig_carrier_death,by=c("hmdb_id","metabolite_name","super_class_metabolon"),all.x=T)

tmp_r2 <- format(round(cor(inter_sig_carrier_case$beta_dementia, inter_sig_carrier_case$beta_sens_case)[1],2), nsmall = 2)
p1 <- ggplot(data=inter_sig_carrier_case, aes(x=beta_sens_case, y=beta_dementia)) +
  geom_point(aes(fill=super_class_metabolon), size=2, shape=21, stroke=0.2, alpha=0.8) +
  scale_fill_manual(values=classcol) +
  geom_vline(xintercept=0, linetype='dashed', alpha=.5) +
  geom_hline(yintercept=0, linetype='dashed', alpha=.5) +
  labs(x=expression(paste(atop("Self-reported dementia",beta[" metabolites \u00D7" * italic(' APOE4 ') * ""]))),y=expression(paste(atop("Composite dementia",beta[" metabolites \u00D7" * italic(' APOE4 ') * ""])))) +
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

tmp_r2 <- format(round(cor(inter_sig_carrier_death$beta_dementia, inter_sig_carrier_death$beta_sens)[1],2), nsmall = 2)
p2 <- ggplot(data=inter_sig_carrier_death, aes(x=beta_sens_death, y=beta_dementia)) +
  geom_point(aes(fill=super_class_metabolon), size=2, shape=21, stroke=0.2, alpha=0.8) +
  scale_fill_manual(values=classcol) +
  geom_vline(xintercept=0, linetype='dashed', alpha=.5) +
  geom_hline(yintercept=0, linetype='dashed', alpha=.5) +
  labs(x=expression(paste(atop("Dementia death",beta[" metabolites \u00D7" * italic(' APOE4 ') * ""]))),y="") +
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

pdf("results/supp/supp_fig_3.pdf", width = 4.8, height = 2.4, onefile = F)
supp_sens <- egg::ggarrange(p1, ggplot() + theme_void(), p2, 
                            nrow=1, ncol=3,
                            widths = c(1,-0.1,1), heights=c(15))
supp_sens
dev.off()
