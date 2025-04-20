library(dplyr)
library(data.table)
library(ggplot2)
library(readxl)

##-----------------------------------------
## Exd Fig. 8: MedDiet and TAG/DAG
##-----------------------------------------

rm(list=ls())

### Diet results

source("scripts/source_color_superclass.R")

diet_res <- read_xlsx("results/diet_gene_met_main_subgroup.xlsx") 
diet_fig <- diet_res %>%
  mutate(diet_beta = ifelse(is.na(diet_p) | n < 50, NA, diet_beta),
         age_beta = ifelse(is.na(age_p) | n < 50, NA, age_beta)) %>%
  filter(diet=="AMED_avg",
         subgroup=="all" & covar=="full",
         class_broad %in% c("Triglycerides","Diglycerides")) %>%
  mutate(n_atoms = as.numeric(gsub("C([0-9]+):[0-9]+ (TAG|DAG)", "\\1", anno_metabolite_name)),
         n_bonds = as.numeric(gsub("C[0-9]+:([0-9]+) (TAG|DAG)", "\\1", anno_metabolite_name))) %>%
  select(c(class_broad,anno_metabolite_name,n_atoms,n_bonds,diet_beta,diet_p)) %>%
  mutate(diet_beta_abs = abs(diet_beta),
         diet_beta_sign_p = -log10(diet_p) * ifelse(diet_beta>0,1,-1))

dag_tag_fig <- ggplot(diet_fig, aes(x = n_atoms, y = n_bonds)) +
  geom_point(aes(size = diet_beta_abs, fill = diet_beta_sign_p), 
             color="black", shape = 21, alpha = 0.8, stroke = 0.2) +
  scale_size_continuous(range = c(2, 10)) +  # Adjust the size of the bubbles
  scale_fill_gradient2(low = "#007F7F", mid = "white", high = "#D85F4D", midpoint = 0) +  # Color scale for diet_p
  labs(x = "Number of carbon atoms",
       y = "Number of double bonds",
       size = "Effect size",
       color = "-log10(P)*sign(beta)") +
  scale_x_continuous(limits = c(26, 60), breaks = seq(26, 58, by = 4)) +  # Increase x-axis density
  scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2)) +   # Increase y-axis density
  facet_wrap(~class_broad, scales = "free") +  # Stratify by class_broad (DAG and TAG)
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size=12),
        axis.text.y = element_text(color = "black", size=12),
        panel.border = element_rect(color = "black", size = 0.6),
        axis.title.x = element_text(size = 14),
        axis.title.y  = element_text(size = 14),
        strip.text = element_text(size = 15),  
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#D3D3D3", color = "black", size = 1))

pdf("results/extended/extended_fig8.pdf", width = 12, height = 6)
dag_tag_fig
dev.off()
