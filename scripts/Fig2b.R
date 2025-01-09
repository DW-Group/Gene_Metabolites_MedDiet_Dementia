library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)

##-----------------------------------------
## Fig. 2b: Overall findings
##-----------------------------------------

rm(list=ls())

### Data

load("results/met_dem_cox_baseline_nhs_2023_filtered_res_10072024.RData")

tmp1 <- nhs_inter_apoe4_res %>%
  subset(select=c(hmdb_id,metabolite_name,super_class_metabolon,beta_inter_met_gene_1,two_inter_fdr_pvalue_p_inter_met_gene_1)) %>%
  mutate(inter_term = "italic('APOE4')~heterozygote")
colnames(tmp1) <- c("hmdb_id","metabolite_name","super_class_metabolon","Beta","FDR","inter_term")

tmp2 <- nhs_inter_apoe4_res %>%
  subset(select=c(hmdb_id,metabolite_name,super_class_metabolon,beta_inter_met_gene_2,two_inter_fdr_pvalue_p_inter_met_gene_2)) %>%
  mutate(inter_term = "italic('APOE4')~homozygote")
colnames(tmp2) <- c("hmdb_id","metabolite_name","super_class_metabolon","Beta","FDR","inter_term")

inter_all <- as.data.frame(rbind(tmp1, tmp2))  %>%
  mutate(sign = case_when(Beta > 0 ~ "Positive",
                          Beta <= 0 ~ "Inverse"),
         sig = case_when(FDR < 0.05 ~ "sig",
                         FDR >= 0.05 ~ "non_sig"),
         abs_beta = abs(Beta))
inter_all <- inter_all[order(inter_all$FDR),]
inter_all <- inter_all[!duplicated(inter_all$hmdb_id),]
table(inter_all$inter_term)

ind_var_inter_all <- nhs_inter_ind_var_res %>%
  mutate(sign = case_when(beta_inter_met_var > 0 ~ "Positive",
                          beta_inter_met_var <= 0 ~ "Inverse"),
         sig = case_when(p_fdr_inter_met_var < 0.05 ~ "sig",
                         p_fdr_inter_met_var >= 0.05 ~ "non_sig"),
         abs_beta = abs(beta_inter_met_var))
ind_var_inter_all <- ind_var_inter_all[order(ind_var_inter_all$p_fdr_inter_met_var),]
ind_var_inter_all <- ind_var_inter_all[!duplicated(ind_var_inter_all$hmdb_id),]

table(inter_all$hmdb_id %in% ind_var_inter_all$hmdb_id)
inter_all <- merge(subset(ind_var_inter_all,select=c(hmdb_id, metabolite_name, super_class_metabolon)), inter_all, by="hmdb_id", all.x=T, suffixes=c("",".y"))
inter_all$sign <- ifelse(is.na(inter_all$sign),"Positive",inter_all$sign)
inter_all$inter_term <- ifelse(is.na(inter_all$inter_term),"italic('APOE4')~homozygote",inter_all$inter_term)
inter_all$FDR <- ifelse(is.na(inter_all$FDR),1,inter_all$FDR)

table(ind_var_inter_all$hmdb_id %in% ind_var_inter_all$hmdb_id)
ind_var_inter_all <- merge(subset(ind_var_inter_all,select=c(hmdb_id, metabolite_name, super_class_metabolon)), ind_var_inter_all, by="hmdb_id", all.x=T, suffixes=c("",".y"))
ind_var_inter_all$sign <- ifelse(is.na(ind_var_inter_all$sign),"Positive",ind_var_inter_all$sign)
ind_var_inter_all$p_fdr_inter_met_var <- ifelse(is.na(ind_var_inter_all$p_fdr_inter_met_var),1,ind_var_inter_all$p_fdr_inter_met_var)

inter_all$inter_term <- factor(inter_all$inter_term, levels=c("italic('APOE4')~heterozygote","italic('APOE4')~homozygote","Other~AD~variants"))
inter_all$sign <- factor(inter_all$sign, levels=c("Positive","Inverse"))
ind_var_inter_all$sign <- factor(ind_var_inter_all$sign, levels=c("Positive","Inverse"))

table(inter_all$inter_term)

inter_all <- inter_all %>% arrange(desc(super_class_metabolon), desc(metabolite_name))
inter_all$super_class_metabolon <- factor(inter_all$super_class_metabolon, levels=names(classcol))
metabolite_name_vec <- unique(inter_all$metabolite_name)
inter_all$metabolite_name <- factor(inter_all$metabolite_name, levels=metabolite_name_vec)

ind_var_inter_all <- ind_var_inter_all %>% arrange(desc(super_class_metabolon), desc(metabolite_name))
ind_var_inter_all$super_class_metabolon <- factor(ind_var_inter_all$super_class_metabolon, levels=names(classcol))
ind_var_inter_all$metabolite_name <- factor(ind_var_inter_all$metabolite_name, levels=metabolite_name_vec)

table(inter_all$hmdb_id == ind_var_inter_all$hmdb_id)
tmp <- as.data.frame(str_split(inter_all$metabolite_name,"/",simplify=T)) %>% mutate(V1 = capitalize(V1)) 
inter_all$anno_metabolite_name <- tmp$V1

### Annotation

class_met <- ggplot(inter_all, aes(x=metabolite_name, y=1)) +
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
        legend.position = "bottom",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14),
        legend.box.margin = margin(t=200, r=200, b=200, l=200))

pdf("results/figures/fig2_legend_10102024.pdf", width = 25, height = 3, onefile = F)
legend_fig2 <- cowplot::get_plot_component(class_met, 'guide-box', return_all = TRUE)
grid.newpage()
grid.draw(legend_fig2[[3]])
dev.off()

### Manhattan plot: interaction

inter_mhtplot <- ggplot(inter_all, aes(x=metabolite_name, y=-log10(FDR), fill=sign)) +
  geom_point(aes(fill=sign, shape=inter_term), size=2.5) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  scale_fill_manual(values=c("#c42b2b","#7195b0")) +
  scale_shape_manual(values=c(23,22,21),drop = FALSE, labels=ggplot2:::parse_safe) +
  scale_x_discrete(limits=rev) +
  labs(fill="Direction of interaction",shape="Interaction term",size="Effect size") +
  ylab(expression("-log"[10]*"(FDR)")) +
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
  annotation_custom(grobTree(textGrob(parse(text=("Interaction~with~italic('APOE4')")), x=0.98,  y=0.92, hjust=1,
                                      gp=gpar(fontsize=20, fontface="bold")))) +
  theme(legend.position = "top", legend.title = element_text(size=14), legend.text = element_text(size=12), legend.key.size = unit(0.6, 'cm')) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  geom_label_repel(
    data = subset(inter_all, inter_all$anno_metabolite_name %in% c("C18:0 CE","Betaine","C45:0 TAG")),
    aes(label = anno_metabolite_name),
    fill = NA,
    segment.colour = NA,
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    show.legend = F
  )

pdf("results/figures/fig2b_legend_11222024.pdf", width = 25, height = 3, onefile = F)
legend_fig2b <- cowplot::get_plot_component(inter_mhtplot, 'guide-box', return_all = TRUE)
grid.newpage()
grid.draw(legend_fig2b[[4]])
dev.off()

### Manhattan plot: interaction

inter_ind_var_mhtplot <- ggplot(ind_var_inter_all, aes(x=metabolite_name, y=-log10(p_fdr_inter_met_var))) +
  geom_point(aes(fill=sign),shape=21, size=2.5) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  ylim(c(0,3)) +
  scale_fill_manual(values=c("#c42b2b","#7195b0")) +
  scale_x_discrete(limits=rev) +
  labs(fill="Direction of interaction",size="Effect size") +
  ylab(expression("-log"[10]*"(FDR)")) +
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
  theme(legend.position = "none") +
  annotation_custom(grobTree(textGrob("Interaction with other AD variants", x=0.98,  y=0.92, hjust=1,
                                      gp=gpar(fontsize=20)))) +
  theme(plot.title = element_blank())+
  geom_text_repel(
    data = subset(ind_var_inter_all, p_fdr_inter_met_var <= 0.05),
    aes(label = paste0(rsid,"-",ref, "~(italic('", gene,"'))")),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    parse=T
  )

pdf("results/figures/fig2b_overall_findings_11222024.pdf", width = 11, height = 8, onefile = F)
fig2b <- ggpubr::ggarrange(inter_mhtplot, class_met, inter_ind_var_mhtplot, ncol=1, nrow=3,
                           widths = c(30), heights=c(1,0.1,1),
                           align = "v",
                           legend = "none")
fig2b
dev.off()
