library(ggplot2)
library(RColorBrewer)

##-----------------------------------------
## Exd Fig. 10: Time-dependent AUC (HPFS)
##-----------------------------------------

rm(list=ls())

##------------------------
## Load data
##------------------------

load("results/prediction_dementia_auc_by_time_hpfs.Rdata")

##------------------------
## AUC by time
##------------------------

all <- ggplot(data=auc_dat_list[["all"]][["all"]], aes(x=time, y=auc, color=anno_predictor)) +
  geom_line(alpha=0.8) +
  scale_color_manual(values=brewer.pal(n=8, name="Dark2")[c(1,2,7,3)], labels=ggplot2:::parse_safe) +
  labs(x="Time (months)", y="AUC", title="Overall dementia risk (HPFS)", color = "Predictor") +
  ylim(c(0.6,0.8)) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "right",
  ) 

pdf("results/extended/extended_fig10a.pdf", width = 5, height = 2.8)
all
dev.off()
