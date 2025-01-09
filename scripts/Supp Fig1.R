library(dplyr)
library(data.table)
library(survival)
library(survminer)

##-----------------------------------------
## Supp Fig. 1: Dementia risk by genetics
##-----------------------------------------

rm(list=ls())

load("data/merged_data_baseline_nhs_2023_10072024.Rdata")

### APOE4 vs. dementia

merged_data <- merged_data %>%
  mutate(APOE4 = case_when(APOE4_ncopy_1 == 1 ~ 1,
                           APOE4_ncopy_2 == 1 ~ 2,
                           APOE4_ncopy_1 == 0 & APOE4_ncopy_2 == 0 ~ 0),
         APOE4 = as.factor(APOE4))
table(merged_data$APOE4)

### PGS000334 tertile vs. dementia 

table(merged_data$PGS000334_t)
merged_data$PGS000334_t  <- factor(merged_data$PGS000334_t)

fit <- survfit(Surv(tad, adcase) ~ PGS000334_t, data=merged_data)
tmp_genetic <- merged_data %>% mutate(PGS000334_t = relevel(PGS000334_t, ref = "2"))
coxmod_ci1 <- summary(coxph(Surv(tad, adcase) ~ PGS000334_t , data=tmp_genetic))$conf.int[1,]
coxmod_ci2 <- summary(coxph(Surv(tad, adcase) ~ PGS000334_t , data=tmp_genetic))$conf.int[2,]
text1 <- paste0(round(coxmod_ci1[1],2)," (",round(coxmod_ci1[3],2),", ",round(coxmod_ci1[4],2),")")
text2 <- paste0(round(coxmod_ci2[1],2)," (",round(coxmod_ci2[3],2),", ",round(coxmod_ci2[4],2),")")
pdf("results/supp/supp_fig_1.pdf", width = 5, height = 4, onefile = F)
p <- ggsurvplot(fit, 
                     legend.labs = c(1,2,3),
                     conf.int = TRUE,
                     palette = c("tan","tan3","tan4"),
                     fun = "event",
                     title = "AD PRS and dementia risk", 
                     legend.title = "PRS tertile", 
                     xlab = "Time (months)",
                     ylab = "Cumulative incidence", 
                     pval = TRUE)
p$plot <- p$plot +
  ggplot2::annotate(
    "text",
    x = 0, y = 0.1,
    vjust = 0.5, hjust = 0,
    label = paste0("HR (95% CI):\nT1 vs T2: ",text1,"\nT3 vs T2: ",text2),
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
p
dev.off()