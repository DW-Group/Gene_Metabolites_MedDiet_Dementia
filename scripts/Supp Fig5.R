library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

##-----------------------------------------
## Supp Fig. 5: TICS by AMED
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

linear_model1 <- lm(av_ntotal ~ agecon + amedcon +
                     qnSES1 + qnSES2 +
                     actcc2 + actcc3 + actcc4 + actcc5 +
                     smk2 + smk3 +
                     hbp + dprs + db + anti + fhdem + marry + live_alone + nhor2 +
                     hiedu2 + hiedu3 +
                     husbedu2 + husbedu3 +
                     bmic1 + bmic3 + bmic4 + bmic5, 
                   data = tics_dat)
linear_model2 <- lm(av_nverbl ~ agecon + amedcon +
                      qnSES1 + qnSES2 +
                      actcc2 + actcc3 + actcc4 + actcc5 +
                      smk2 + smk3 +
                      hbp + dprs + db + anti + fhdem + marry + live_alone + nhor2 +
                      hiedu2 + hiedu3 +
                      husbedu2 + husbedu3 +
                      bmic1 + bmic3 + bmic4 + bmic5, 
                    data = tics_dat)

# Predict values and confidence intervals
lm_predictions1 <- predict(linear_model1, newdata = new_data, interval = "confidence")
lm_predictions2 <- predict(linear_model2, newdata = new_data, interval = "confidence")

# Add predictions and CI to new_data
new_data$lm_predicted1 <- lm_predictions1[, "fit"]
new_data$lm_lwr1 <- lm_predictions1[, "lwr"]
new_data$lm_upr1 <- lm_predictions1[, "upr"]
new_data$lm_predicted2 <- lm_predictions2[, "fit"]
new_data$lm_lwr2 <- lm_predictions2[, "lwr"]
new_data$lm_upr2 <- lm_predictions2[, "upr"]

### Plot

pdf("results/supp/supp_fig_5_1.pdf", width = 3, height = 3, onefile = F)
p_tics_amed <- ggplot(new_data, aes(x = amedcon, y = lm_predicted1)) +
  geom_ribbon(aes(ymin = lm_lwr1, ymax = lm_upr1), alpha=0.3, linetype=1, fill="#7A8D50") +
  labs(x="MedDiet index score",y="Predicted global cognition score (95% CI)") +
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

pdf("results/supp/supp_fig_5_2.pdf", width = 3, height = 3, onefile = F)
p_tics_amed <- ggplot(new_data, aes(x = amedcon, y = lm_predicted2)) +
  geom_ribbon(aes(ymin = lm_lwr2, ymax = lm_upr2), alpha=0.3, linetype=1, fill="#7A8D50") +
  labs(x="MedDiet index score",y="Predicted verbal memory score (95% CI)") +
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
