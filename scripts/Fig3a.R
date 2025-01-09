library(dplyr)
library(data.table)
library(ggplot2)
library(egg)
library(tidyr)

##-----------------------------------------
## Fig. 3a: Dementia by AMED
##-----------------------------------------

rm(list=ls())

### Load data

tmp_dat <- fread("scripts/amed_dementia_spline_baseline.txt") %>%
  mutate(across(c(AMEDCON,Estimate,Lower,Upper), as.numeric))

### Spline curve

pdf("results/figures/fig3a_amed_dementia_spline_11202024.pdf", width = 2, height = 3, onefile = F)
p_dementia_amed <- ggplot(tmp_dat, aes(x = AMEDCON, y = Estimate)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, linetype=1, fill="#7A8D50") +
  labs(x="MedDiet index score",y="HR (95% CI) of dementia risk") +
  geom_line(color="#4F5937", linewidth=1) + 
  geom_hline(yintercept=1, linetype='dashed', alpha=.5) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    panel.border = element_blank(),  
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),  
    axis.ticks = element_line(color = "black") 
  ) +
  annotate("text", x = 0, y = 0.55 ,
           label = "P < 0.001", hjust = 0, vjust = 0,
           color = "black", size = 4)
p_dementia_amed
dev.off()

##-----------------------------------------
## Fig. 3a: TICS by AMED
##-----------------------------------------

rm(list=ls())

### Load data

tics_dat <- read.csv("data/dementia/nhs_tics_full.csv")

tics_dat_list <- list("all" = tics_dat,
                      "genetic" = tics_dat %>% filter(!is.na(APOE4)),
                      "APOE4_noncarrier" = tics_dat %>% filter(APOE4==0),
                      "APOE4_carrier" = tics_dat %>% filter(APOE4!=0),
                      "APOE4_heterozygote" = tics_dat %>% filter(APOE4==1),
                      "APOE4_homozygote" = tics_dat %>% filter(APOE4==2))

new_data <- data.frame(
  amedcon = seq(min(tics_dat$amedcon), max(tics_dat$amedcon), length.out = 100),
  agecon = min(tics_dat$agecon),  
  qnSES1 = 0,
  qnSES2 = 0,
  actcc2 = 0,
  actcc3 = 0,
  actcc4 = 0,
  actcc5 = 0,
  smk2 = 0,
  smk3 = 0,
  hbp = 0,
  dprs = 1,
  db = 0,
  anti = 0,
  fhdem = 0,
  marry = 0,
  live_alone = 0,
  nhor2 = 0,
  hiedu2 = 0,
  hiedu3 = 0,
  husbedu2 = 0,
  husbedu3 = 0,
  bmic1 = 0,
  bmic3 = 0,
  bmic4 = 0,
  bmic5 = 0
)

### Fit linear model

linear_model <- lm(av_zscre ~ agecon + amedcon +
                      qnSES1 + qnSES2 +
                      actcc2 + actcc3 + actcc4 + actcc5 +
                      smk2 + smk3 +
                      hbp + dprs + db + anti + fhdem + marry + live_alone + nhor2 +
                      hiedu2 + hiedu3 +
                      husbedu2 + husbedu3 +
                      bmic1 + bmic3 + bmic4 + bmic5, 
                    data = tics_dat)
summary(linear_model)
  
# Predict values and confidence intervals
lm_predictions <- predict(linear_model, newdata = new_data, interval = "confidence")

# Add predictions and CI to new_data
new_data$lm_predicted <- lm_predictions[, "fit"]
new_data$lm_lwr <- lm_predictions[, "lwr"]
new_data$lm_upr <- lm_predictions[, "upr"]

for(subgroup in names(tics_dat_list)){
  linear_model <- lm(av_zscre ~ agecon + amedcon +
                       qnSES1 + qnSES2 +
                       actcc2 + actcc3 + actcc4 + actcc5 +
                       smk2 + smk3 +
                       hbp + dprs + db + anti + fhdem + marry + live_alone + nhor2 +
                       hiedu2 + hiedu3 +
                       husbedu2 + husbedu3 +
                       bmic1 + bmic3 + bmic4 + bmic5
                     # + antihp
                     , 
                     data = tics_dat_list[[subgroup]])

  summary(linear_model)
  
  # Predict values and confidence intervals
  lm_predictions <- predict(linear_model, newdata = new_data, interval = "confidence")
  
  # Add predictions and CI to new_data
  new_data[[paste0(subgroup,"_pred")]] <- lm_predictions[, "fit"]
  new_data[[paste0(subgroup,"_lwr")]] <- lm_predictions[, "lwr"]
  new_data[[paste0(subgroup,"_upr")]] <- lm_predictions[, "upr"]
}

### Plot

pdf("results/figures/fig3a_amed_tics_linear_11262024.pdf", width = 1.5, height = 3, onefile = F)
p_tics_amed <- ggplot(new_data, aes(x = amedcon, y = all_pred)) +
  geom_ribbon(aes(ymin = all_lwr, ymax = all_upr), alpha=0.3, linetype=1, fill="#7A8D50") +
  labs(x="MedDiet index score",y="Predicted TICS score (95% CI)") +
  geom_line(color="#4F5937", linewidth=1) +
  theme_bw() +  
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"), 
    axis.title = element_text(color = "black"), 
    panel.border = element_blank(), 
    plot.margin = margin(5.5, 5.5, 5.5, 5.5), 
    axis.ticks = element_line(color = "black") 
  ) +
  annotate("text", x = 5.5, y = 0,
           label = "P < 0.001", 
           hjust = 0, vjust = 0,
           color = "black", size = 4)
p_tics_amed
dev.off()

pdf("results/figures/fig3a_amed_both_1125024.pdf", width = 5, height = 3, onefile = F)
fig3a <- egg::ggarrange(p_dementia_amed,p_tics_amed,
                        nrow=1, ncol=2,
                        widths = c(1,1), heights=c(1))
fig3a
dev.off()