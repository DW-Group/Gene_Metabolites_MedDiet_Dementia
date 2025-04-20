library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(DescTools)
library(cowplot)
library(grid)
library(gridExtra)
library(ggh4x)
library(ggtext)

##-----------------------------------------
## Fig. 2a: Significant findings
##-----------------------------------------

rm(list=ls())

### Load data

load("results/met_dem_cox_baseline_nhs_2023_filtered_res.RData")

sig_var <- na.omit(unique(nhs_inter_ind_var_res[nhs_inter_ind_var_res$p_fdr_inter_met_var < 0.05,]$variant))

inter_sig <- merge(sig_met_dat, nhs_inter_apoe4_res, by="hmdb_id", all.x=T, suffixes=c("",".y"))
strat_sig <- merge(sig_met_dat, nhs_strat_res, by="hmdb_id", all.x=T, suffixes=c("",".y")) %>%
  filter(!(model%in% c("PGS002280_t1","PGS002280_t2","PGS002280_t3","APOE4_carrier")))
ind_var_sig <- merge(sig_met_dat, nhs_inter_ind_var_res, by="hmdb_id", all.x=T, suffixes=c("",".y")) %>%
  filter(variant %in% sig_var)

unique(strat_sig$model)
length(unique(inter_sig$hmdb_id)); length(unique(strat_sig$hmdb_id)); length(unique(ind_var_sig$hmdb_id))

### Figure data

strat_sig <- strat_sig %>%
  mutate(sig = case_when(p_fdr < 0.05 ~ "**",
                         p_fdr >= 0.05 &  p_fdr < 0.25 ~ "*",
                         p_fdr >= 0.25 | is.na(p_fdr) ~ ""))
table(strat_sig$sig)
inter_sig <- inter_sig %>%
  mutate(sig_inter_1 = case_when(two_inter_fdr_pvalue_p_inter_met_gene_1 < 0.05 ~ "**",
                                 two_inter_fdr_pvalue_p_inter_met_gene_1 >= 0.05 &  two_inter_fdr_pvalue_p_inter_met_gene_1 < 0.25 ~ "*",
                                 two_inter_fdr_pvalue_p_inter_met_gene_1 >= 0.25 | is.na(two_inter_fdr_pvalue_p_inter_met_gene_1) ~ ""),
         sig_inter_2 = case_when(two_inter_fdr_pvalue_p_inter_met_gene_2 < 0.05 ~ "**",
                                 two_inter_fdr_pvalue_p_inter_met_gene_2 >= 0.05 &  two_inter_fdr_pvalue_p_inter_met_gene_2 < 0.25 ~ "*",
                                 two_inter_fdr_pvalue_p_inter_met_gene_2 >= 0.25 | is.na(two_inter_fdr_pvalue_p_inter_met_gene_2) ~ ""))
table(inter_sig$sig_inter_1)
table(inter_sig$sig_inter_2)
ind_var_sig <- ind_var_sig %>%
  mutate(sig_inter = case_when(p_fdr_inter_met_var < 0.05 ~ "**",
                               p_fdr_inter_met_var >= 0.05 &  p_fdr_inter_met_var < 0.25 ~ "*",
                               p_fdr_inter_met_var | is.na(p_fdr_inter_met_var) >= 0.25 ~ ""))
table(ind_var_sig$sig_inter)
unique(inter_sig$super_class_metabolon)
table(inter_sig$super_class_metabolon)

ordered_super_class_sig <- data.frame(super_class_metabolon=c("Organic acids and derivatives","Organic nitrogen compounds","Organoheterocyclic compounds","Lipids and lipid-like molecules"),
                                  super_class_order=4:1)

strat_sig <- merge(ordered_super_class_sig, strat_sig, by="super_class_metabolon") %>% 
  arrange(super_class_order, desc(anno_metabolite_name_sig))
strat_sig$super_class_metabolon <- factor(strat_sig$super_class_metabolon, levels=names(classcol))
anno_metabolite_name_sig_vec <- unique(strat_sig$anno_metabolite_name_sig)
strat_sig$anno_metabolite_name_sig <- factor(strat_sig$anno_metabolite_name_sig, levels=anno_metabolite_name_sig_vec)

inter_sig <- merge(ordered_super_class_sig, inter_sig, by="super_class_metabolon") %>% 
  arrange(super_class_order, desc(anno_metabolite_name_sig))
inter_sig$super_class_metabolon <- factor(inter_sig$super_class_metabolon, levels=names(classcol))
inter_sig$anno_metabolite_name_sig <- factor(inter_sig$anno_metabolite_name_sig, levels=anno_metabolite_name_sig_vec)

ind_var_sig <- merge(ordered_super_class_sig, ind_var_sig, by="super_class_metabolon") %>% 
  arrange(super_class_order, desc(anno_metabolite_name_sig))
ind_var_sig$super_class_metabolon <- factor(ind_var_sig$super_class_metabolon, levels=names(classcol))
ind_var_sig$anno_metabolite_name_sig <- factor(ind_var_sig$anno_metabolite_name_sig, levels=anno_metabolite_name_sig_vec)

### Annotation

class_met <- ggplot(strat_sig[!duplicated(strat_sig$hmdb_id),], aes(x=1, y=anno_metabolite_name_sig)) +
  geom_tile(aes(fill=super_class_metabolon)) +
  scale_fill_manual(values=classcol) +
  labs(fill="Superclass of metabolites", size=12) +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_rect(color="black", size=0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 14),
        axis.text.x = element_blank(),
        axis.title =element_blank(),
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14)) +
  theme(legend.position="none") 

### Heatmap: stratify by APOE4 and PRS with APOE

table(strat_sig[!duplicated(strat_sig$hmdb_id),]$hmdb_id == ind_var_sig[!duplicated(ind_var_sig$hmdb_id),]$hmdb_id) # 60
table(strat_sig[!duplicated(strat_sig$hmdb_id),]$hmdb_id == strat_sig[!duplicated(strat_sig$hmdb_id),]$hmdb_id) # 60

summary(exp(strat_sig$beta))
summary(Winsorize(exp(strat_sig$beta), val=quantile(exp(strat_sig$beta), probs=c(0.05, 0.95), na.rm=T)))
strat_sig$wins_hr <- Winsorize(exp(strat_sig$beta), val=quantile(exp(strat_sig$beta), probs=c(0.05, 0.95), na.rm=T))
hr_limits_low <- min(strat_sig$wins_hr, na.rm=T)
hr_limits_high <- max(strat_sig$wins_hr, na.rm=T)

table(strat_sig$model)

tmp_strat <- strat_sig %>%
  mutate(gene_group_anno = case_when(grepl("APOE4",model) ~ "italic('APOE4')",
                                     grepl("PGS",model) ~ "PRS")) %>%
  mutate(analysis_group = "Effect~estimates~by~genotype") %>%
  mutate(anno_model = case_when(model == "APOE4_noncarrier" ~ "Noncarrier", 
                                model == "APOE4_ncopy_1" ~ "Heterozygote", 
                                model == "APOE4_ncopy_2" ~ "Homozygote", 
                                model == "PGS000334_t1" ~ "Low", 
                                model == "PGS000334_t2" ~ "Intermediate", 
                                model == "PGS000334_t3" ~ "High")) %>%
  mutate(anno_model = factor(anno_model, levels=c("Noncarrier","Heterozygote","Homozygote",
                                                  "Low","Intermediate","High")))

table(tmp_strat$anno_model)
table(tmp_strat$gene_group_anno)
table(tmp_strat$analysis_group)
  
strat_geno_heatmap <- ggplot(tmp_strat, aes(x=anno_model, y=anno_metabolite_name_sig)) +
  geom_tile(aes(fill=wins_hr)) +
  scale_fill_gradient2(low="#466983",
                       mid="white",
                       high="#8b0000",
                       midpoint = 1,
                       guide = "colourbar",
                       limits=c(hr_limits_low,hr_limits_high),
                       na.value="gray") +
  scale_x_discrete(position = "bottom") +
  geom_text(aes(label=sig), color="black", size=8, nudge_y = -0.3) +
  labs(fill="Hazard ratio", size=12) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14, colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        panel.spacing.x = unit(1,"line"),
        strip.background = element_rect(color="white"),
        ) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text = element_text(size=12), legend.key.size = unit(0.6, 'cm')) +
  facet_nested(~ analysis_group + gene_group_anno, labeller = label_parsed, drop = T, scales="free_x",
               strip = strip_nested(
                 background_x = list(element_rect(color="black",fill = "#466983",size=0.8),
                                     element_rect(color="black",fill = "white",size=0.8),
                                     element_rect(color="black",fill = "white",size=0.8)),
                 text_x = list(element_text(color="white", size=16),
                               element_text(color="black", size=14, face="bold"),
                               element_text(color="black", size=14, face="bold"))))

pdf("results/figures/fig2a_strat_legend.pdf", width = 5, height = 3, onefile = F)
legend_fig2a_strat <- cowplot::get_plot_component(strat_geno_heatmap, 'guide-box', return_all = TRUE)
grid.newpage()
grid.draw(legend_fig2a_strat[[3]])
dev.off()

### Heatmap: interaction for APOE4 and other AD variants

inter_sig <- inter_sig %>%
  mutate(beta_fdr_1 = beta_inter_met_gene_1 * -log10(two_inter_fdr_pvalue_p_inter_met_gene_1),
         beta_fdr_2 = beta_inter_met_gene_2 * -log10(two_inter_fdr_pvalue_p_inter_met_gene_2))
ind_var_sig <- ind_var_sig %>%
  mutate(beta_fdr = beta_inter_met_var * -log10(p_fdr_inter_met_var))

table(inter_sig$covariates)
table(ind_var_sig$covariates)

tmp_inter <- inter_sig %>%
  filter(covariates == "APOE4_2cat_PGS002280") %>%
  subset(select = c(anno_metabolite_name_sig, 
                    beta_inter_met_gene_1, beta_inter_met_gene_2,
                    se_inter_met_gene_1, se_inter_met_gene_2,
                    two_inter_fdr_pvalue_p_inter_met_gene_1, two_inter_fdr_pvalue_p_inter_met_gene_2,
                    beta_fdr_1,beta_fdr_2,
                    sig_inter_1,sig_inter_2))

tmp_inter_ind_var <- ind_var_sig %>%
  mutate(var_gene = paste0(rsid,"-",ref, "~(italic('", gene,"'))")) %>%
  subset(select = c(anno_metabolite_name_sig, 
                    beta_inter_met_var,
                    se_inter_met_var,
                    p_fdr_inter_met_var,
                    beta_fdr,
                    sig_inter,
                    var_gene))

tmp1 <- tmp_inter %>%
  subset(select=c(anno_metabolite_name_sig,beta_inter_met_gene_1,se_inter_met_gene_1,two_inter_fdr_pvalue_p_inter_met_gene_1,beta_fdr_1,sig_inter_1)) %>%
  mutate(inter_term = "Heterozygote",
         gene_group_anno = "italic('APOE4')")
colnames(tmp1) <- c("anno_metabolite_name_sig","Beta","se","FDR","beta_fdr","sig","inter_term","gene_group_anno")

tmp2 <- tmp_inter %>%
  subset(select=c(anno_metabolite_name_sig,beta_inter_met_gene_2,se_inter_met_gene_2,two_inter_fdr_pvalue_p_inter_met_gene_2,beta_fdr_2,sig_inter_2)) %>%
  mutate(inter_term = "Homozygote",
         gene_group_anno = "italic('APOE4')")
colnames(tmp2) <- c("anno_metabolite_name_sig","Beta","se","FDR","beta_fdr","sig","inter_term","gene_group_anno")

tmp3 <- tmp_inter_ind_var %>%
  subset(select=c(anno_metabolite_name_sig,beta_inter_met_var,se_inter_met_var,p_fdr_inter_met_var,beta_fdr,sig_inter,var_gene)) %>%
  mutate(gene_group_anno = "Other~AD/ADRD~variants")
colnames(tmp3) <- c("anno_metabolite_name_sig","Beta","se","FDR","beta_fdr","sig","inter_term","gene_group_anno")

tmp_inter <- as.data.frame(rbind(tmp1, tmp2, tmp3)) %>%
  mutate(wins_beta_fdr = Winsorize(beta_fdr, val=quantile(beta_fdr, probs = c(0.125, 0.875), na.rm=T)),
         analysis_group = "Interaction~of~metabolites~and~genotype")

summary(tmp_inter$beta_fdr)
summary(tmp_inter$wins_beta_fdr)

table(tmp_inter$inter_term)
table(tmp_inter$analysis_group)
table(tmp_inter$gene_group_anno)

inter_geno_heatmap <- ggplot(tmp_inter, aes(x=inter_term, y=anno_metabolite_name_sig)) +
  geom_tile(aes(fill=wins_beta_fdr)) +
  scale_fill_gradient2(low="#84ad2d",
                       mid="white",
                       high="#e88f1a",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks = c(-0.05,0,0.05),
                       na.value="gray") +
  scale_x_discrete(position = "bottom", labels=ggplot2:::parse_safe) +
  geom_text(aes(label=sig), color="black", size=8, nudge_y = -0.3) +
  labs(fill=expression(beta * " \u00D7 -log"[10]*"(FDR)"), size=12) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14, colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        panel.spacing.x = unit(1,"line"),
        strip.background = element_rect(color="white"),
        strip.text.x = element_text(size=12, face ="bold", color="white")) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text = element_text(size=12), legend.key.size = unit(0.6, 'cm')) +
  facet_nested(~ analysis_group + gene_group_anno, labeller = label_parsed, drop = T, scales="free_x", space="free_x",
               strip = strip_nested(
                 background_x = list(element_rect(color="black",fill = "darkgreen",size=0.8),
                                     element_rect(color="black",fill = "white",size=0.8),
                                     element_rect(color="black",fill = "white",size=0.8)),
                 text_x = list(element_text(color="white", size=16, face="bold"),
                               element_text(color="black", size=14, face="bold"),
                               element_text(color="black", size=14, face="bold"))))

pdf("results/figures/fig2a_inter_legend.pdf", width = 5, height = 3, onefile = F)
legend_fig2a_inter <- cowplot::get_plot_component(inter_geno_heatmap, 'guide-box', return_all = TRUE)
grid.newpage()
grid.draw(legend_fig2a_inter[[3]])
dev.off()

blank_fig <- ggplot() + theme_void()

pdf("results/figures/fig2a_sig_findings.pdf", width = 14, height = 21, onefile = F)
fig2a <- egg::ggarrange(class_met, 
                        strat_geno_heatmap + theme(legend.position = "none"), 
                        inter_geno_heatmap + theme(legend.position = "none"),
               blank_fig,
               nrow=1, ncol=4,
               widths = c(0.15,2,2.6,2), heights=c(15))
fig2a
dev.off()
