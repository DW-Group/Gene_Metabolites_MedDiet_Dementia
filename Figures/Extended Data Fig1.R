library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
  
##-----------------------------------------
## Exd Fig. 1: Overview of HPFS
##-----------------------------------------

##-----------------------------------------
## 1b: dementia risk by genetic
##-----------------------------------------

rm(list=ls())

load("data/merged_data_baseline_hpfs_2023.Rdata")

### APOE4 vs. dementia

merged_data_hpfs <- merged_data_hpfs %>%
  mutate(APOE4 = case_when(APOE4_ncopy_1 == 1 ~ 1,
                           APOE4_ncopy_2 == 1 ~ 2,
                           APOE4_ncopy_1 == 0 & APOE4_ncopy_2 == 0 ~ 0),
         APOE4 = as.factor(APOE4))
table(merged_data_hpfs$APOE4)

fit <- survfit(Surv(tad, adcase) ~ APOE4 , data=merged_data_hpfs)
tmp_genetic <- merged_data_hpfs %>% mutate(APOE4 = relevel(APOE4, ref = "0"))
coxmod_ci1 <- summary(coxph(Surv(tad, adcase) ~ APOE4 , data=tmp_genetic))$conf.int[1,]
coxmod_ci2 <- summary(coxph(Surv(tad, adcase) ~ APOE4 , data=tmp_genetic))$conf.int[2,]
text1 <- paste0(round(coxmod_ci1[1],2)," (",round(coxmod_ci1[3],2),", ",round(coxmod_ci1[4],2),")")
text2 <- paste0(round(coxmod_ci2[1],2)," (",round(coxmod_ci2[3],2),", ",round(coxmod_ci2[4],2),")")
pdf("results/figures/fig1d_cumu_incid_APOE4_unadjusted_all_hpfs.pdf", width = 5, height = 4, onefile = F)
p1_all <- ggsurvplot(fit,  
                     legend.labs = c(0,1,2),
                     conf.int = TRUE,
                     palette = c("tan","tan3","tan4"),
                     fun = "event",
                     title = parse(text=c(italic('APOE4')~"and dementia risk (HPFS)")), 
                     legend.title = parse(text=c(italic('APOE4')~"allele")),
                     xlab = "Time (months)",
                     ylab = "Cumulative incidence", 
                     pval = TRUE)
p1_all$plot <- p1_all$plot +
  ggplot2::annotate(
    "text",
    x = 0, y = 0.45,
    vjust = 0.5, hjust = 0,
    label = paste0("HR (95% CI):\nHeterozygote vs noncarrier: ",text1,"\nHomozygote vs noncarrier: ",text2),
    size = 4
  ) +
  theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 13), 
    legend.key.size = unit(0.5 , "cm")       
  )
p1_all
dev.off()

### PGS002280 tertile vs. dementia

table(merged_data_hpfs$PGS002280_t)
merged_data_hpfs$PGS002280_t  <- factor(merged_data_hpfs$PGS002280_t)

fit <- survfit(Surv(tad, adcase) ~ PGS002280_t, data=merged_data_hpfs)
tmp_genetic <- merged_data_hpfs %>% mutate(PGS002280_t = relevel(PGS002280_t, ref = "1"))
coxmod_ci1 <- summary(coxph(Surv(tad, adcase) ~ PGS002280_t , data=tmp_genetic))$conf.int[1,]
coxmod_ci2 <- summary(coxph(Surv(tad, adcase) ~ PGS002280_t , data=tmp_genetic))$conf.int[2,]
text1 <- paste0(round(coxmod_ci1[1],2)," (",round(coxmod_ci1[3],2),", ",round(coxmod_ci1[4],2),")")
text2 <- paste0(round(coxmod_ci2[1],2)," (",round(coxmod_ci2[3],2),", ",round(coxmod_ci2[4],2),")")
pdf("results/figures/fig1d_cumu_incid_PGS002280_t_unadjusted_all_hpfs.pdf", width = 5, height = 4, onefile = F)
p2_all <- ggsurvplot(fit, 
                 legend.labs = c(1,2,3),
                 conf.int = TRUE,
                 palette = c("tan","tan3","tan4"),
                 fun = "event",
                 title = "ADRD PRS and dementia risk (HPFS)", 
                 legend.title = "PRS tertile", 
                 xlab = "Time (months)",
                 ylab = "Cumulative incidence",
                 pval = TRUE)
p2_all$plot <- p2_all$plot +
  ggplot2::annotate(
    "text",
    x = 0, y = 0.1,
    vjust = 0.5, hjust = 0,
    label = paste0("HR (95% CI):\nT2 vs T1: ",text1,"\nT3 vs T1: ",text2),
    size = 4
  ) +
  theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 13), 
    legend.key.size = unit(0.5 , "cm")       
  )
p2_all
dev.off()

pdf("results/extended/extended_fig1b.pdf", width = 6, height = 9, onefile = F)
arrange_ggsurvplots(list(p1_all,p2_all), print = TRUE,
                    ncol = 1, nrow = 2, surv.plot.height = 0.4)
dev.off()
