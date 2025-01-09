library(dplyr)
library(data.table)
library(readxl)
library(ggplot2)
library(cowplot)
library(ggpubr)

##-----------------------------------------
## Fig. 5c: Forest plot for coloc >= 0.7
##-----------------------------------------

rm(list=ls())

### Data

coloc_hits <- read_excel("results/coloc_res_08072024.xlsx") %>%
  filter(PP4_div_PP3_plus_PP4>=0.7|PP.H4.abf>=0.7)
exposure_dict <- read_excel("results/instruments_selected_sig_final_anno_08062024.xlsx") %>%
  subset(select=c(exposure,exposure_anno)) %>% distinct()
coloc_hits <- merge(coloc_hits, exposure_dict, by="exposure_anno") %>%
  mutate(outcome_anno_exposure = paste0(outcome_anno,"_",exposure))

mr_res <- read_excel("results/mr_results_final_08062024.xlsx") %>%
  mutate(outcome_anno_exposure = paste0(anno_outcome,"_",exposure))

mr_res_selected_anno <- merge(mr_res, coloc_hits, by="outcome_anno_exposure") %>%
  mutate(lower=b-1.96*se,
         upper=b+1.96*se,
         aging_effect=case_when(outcome_anno=="Cognitive performance" & b > 0 ~ "Protective",
                                outcome_anno=="Cognitive performance" & b <= 0 ~ "Adverse",
                                outcome_anno!="Cognitive performance" & b <= 0 ~ "Protective",
                                outcome_anno!="Cognitive performance" & b > 0 ~ "Adverse")) %>%
  arrange(desc(cat_outcome))

fig_dat1 <- mr_res_selected_anno %>%
  filter(outcome_anno == "Dementia") %>%
  mutate(or = exp(b),
         or_lower = exp(lower),
         or_upper = exp(upper))
fig_dat1 <- fig_dat1[!duplicated(fig_dat1$outcome_anno_exposure),]
  
fig_dat2 <- mr_res_selected_anno %>%
  filter(outcome_anno == "Alzheimer's disease") %>%
  mutate(or = exp(b),
         or_lower = exp(lower),
         or_upper = exp(upper))
fig_dat2 <- fig_dat2[!duplicated(fig_dat2$outcome_anno_exposure),]

fig_dat3 <- mr_res_selected_anno %>%
  filter(outcome_anno == "Cognitive performance")
fig_dat3 <- fig_dat3[!duplicated(fig_dat3$outcome_anno_exposure),]

### Plot

forest1 <- ggplot(data=fig_dat1, aes(y=exposure_anno, x=or, xmin=or_lower, xmax=or_upper)) +
  geom_errorbarh(size=0.8,height=0.2, color="#444444")+
  geom_point(aes(color=aging_effect), size=4, stroke = 0.3, position=position_dodge(width = 0.5), shape=15, show.legend=T) +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  scale_color_manual(name="Effect on cognitive health",values=c("#8b0000","#466983")) +
  xlab("OR (95% CI)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(color = "black", size = 14, hjust = 1),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(size = 14, colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        legend.position = "left",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14),
        title=element_text(colour = "black", size = 16)) +
  theme(legend.position = "right", legend.title = element_text(size=14), legend.text = element_text(size=12), legend.key.size = unit(0.6, 'cm')) +
  ggtitle(unique(fig_dat1$outcome_anno)) +
  scale_y_discrete(position = "right",limits=rev)

forest2 <- ggplot(data=fig_dat2, aes(y=exposure_anno, x=or, xmin=or_lower, xmax=or_upper)) +
  geom_errorbarh(size=0.8,height=0.2, color="#444444")+
  geom_point(aes(color=aging_effect), size=4, stroke = 0.3, position=position_dodge(width = 0.5), shape=15, show.legend=F) +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  scale_color_manual(name="Effect on cognitive health",values=c("#8b0000","#466983")) +
  xlab("OR (95% CI)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(color = "black", size = 14, hjust = 1),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(size = 14, colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        legend.position = "left",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14),
        title=element_text(colour = "black", size = 16)) +
  theme(legend.position = "right", legend.title = element_text(size=14), legend.text = element_text(size=12), legend.key.size = unit(0.6, 'cm')) +
  ggtitle(unique(fig_dat2$outcome_anno)) +
  scale_y_discrete(position = "right",limits=rev)

forest3 <- ggplot(data=fig_dat3, aes(y=exposure_anno, x=b, xmin=lower, xmax=upper)) +
  geom_errorbarh(size=0.8,height=0.2, color="#444444")+
  geom_point(aes(color=aging_effect), size=4, stroke = 0.3, position=position_dodge(width = 0.5), shape=15, show.legend=F) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  scale_color_manual(name="Effect on cognitive health",values=c("#8b0000","#466983")) +
  xlab(expression(beta * " coefficient (95% CI)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(color = "black", size = 14, hjust = 1),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(size = 14, colour = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        legend.position = "left",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14),
        title=element_text(colour = "black", size = 16)) +
  theme(legend.position = "right", legend.title = element_text(size=14), legend.text = element_text(size=12), legend.key.size = unit(0.6, 'cm')) +
  ggtitle(unique(fig_dat3$outcome_anno)) +
  scale_y_discrete(position = "right",limits=rev)

pdf("results/figures/fig5c_forest.pdf", width = 10, height = 12)
fig5_forest <- ggpubr::ggarrange(forest1, forest2, forest3, ncol=1, nrow=3,
                                 heights = c(3,6,10),
                                 legend = "none",
                                 align = "v")
fig5_forest
dev.off()

pdf("results/figures/fig5c_forest_legend.pdf", width = 5, height = 2, onefile = F)
legend_fig5c <- cowplot::get_plot_component(forest1, 'guide-box', return_all = TRUE)
grid.newpage()
grid.draw(legend_fig5c[[1]])
dev.off()
