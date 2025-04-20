library(data.table)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

##-----------------------------------------
## Fig. 4a: Time-dependent AUC
##-----------------------------------------

rm(list=ls())

##------------------------
## Load data
##------------------------

load("results/prediction_dementia_auc_by_time.Rdata")

##------------------------
## AUC by time
##------------------------

all <- ggplot(data=auc_dat_list[["all_nocovar"]], aes(x=time, y=auc, color=anno_predictor)) +
  geom_line(alpha=0.8) +
  scale_color_manual(values=brewer.pal(n=8, name="Dark2")[c(1,2,7,3)], labels=ggplot2:::parse_safe) +
  labs(x="Time (months)", y="AUC", title="Overall survival") +
  ylim(c(0.7,0.9)) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        ) 

yr_15 <- ggplot(data=auc_dat_list[["yr_15_nocovar"]], aes(x=time, y=auc, color=anno_predictor)) +
  geom_line(alpha=0.8) +
  scale_color_manual(values=brewer.pal(n=8, name="Dark2")[c(1,2,7,3)], labels=ggplot2:::parse_safe) +
  labs(x="Time (months)", y="AUC", title="15-year survival") +
  ylim(c(0.7,0.9)) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
  )

legend <- all + 
  labs(color="Predictors") +
  theme(legend.title = element_text(color = "black"),
        legend.position = "right",
        legend.text.align = 0) 

legend_fig4a <- cowplot::get_plot_component(legend, 'guide-box', return_all = TRUE)[[1]]

blank_fig <- ggplot() + theme_void()

pdf("results/figures/fig4a_auc_by_time.pdf", width = 8, height = 2.8)
fig4a <- ggpubr::ggarrange(all, yr_15,
                           legend_fig4a,
                           ncol=3, nrow=1,
                           widths = c(1,1,1))
fig4a
dev.off()
