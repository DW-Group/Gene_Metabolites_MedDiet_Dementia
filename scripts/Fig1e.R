library(dplyr)
library(ggplot2)

##-----------------------------------------
## Fig. 1e: diet distribution
##-----------------------------------------

rm(list=ls())

### Data

load("data/merged_data_baseline_nhs_2023_10072024.Rdata")

diet_dist_dat <- data.frame(component = c("Vegetables\n(servings/day)",
                                          "Fruits\n(servings/day)", 
                                          "Nuts\n(servings/day)",
                                          "Legumes\n(servings/day)",
                                          "Whole grains\n(servings/day)", 
                                          "Fish\n(servings/day)",
                                          "Monounsaturated fat\n(10g/day)",
                                          "Saturated fat\n(10g/day)",
                                          "Red and processed meat\n(servings/day)",
                                          "Alcohol\n(10g/day)"),
                            component_cat = c("veg_fruit",
                                              "veg_fruit", 
                                              "nut_legume",
                                              "nut_legume",
                                              "wgrains", 
                                              "fish",
                                              "fat",
                                              "fat",
                                              "meat",
                                              "alco"),
                            comp_mean = c(mean(merged_data$veg_avg, na.rm=T),
                                          mean(merged_data$frt_avg, na.rm=T), 
                                          mean(merged_data$nut_avg, na.rm=T),
                                          mean(merged_data$leg_avg, na.rm=T),
                                          mean(merged_data$whgrn_avg, na.rm=T), 
                                          mean(merged_data$fish_avg, na.rm=T),
                                          mean(merged_data$mon10_avg, na.rm=T),
                                          mean(merged_data$sat10_avg, na.rm=T),
                                          mean(merged_data$rmt_avg, na.rm=T),
                                          mean(merged_data$etoh10_avg, na.rm=T)),
                            comp_median = c(median(merged_data$veg_avg, na.rm=T),
                                          median(merged_data$frt_avg, na.rm=T), 
                                          median(merged_data$nut_avg, na.rm=T),
                                          median(merged_data$leg_avg, na.rm=T),
                                          median(merged_data$whgrn_avg, na.rm=T), 
                                          median(merged_data$fish_avg, na.rm=T),
                                          median(merged_data$mon10_avg, na.rm=T),
                                          median(merged_data$sat10_avg, na.rm=T),
                                          median(merged_data$rmt_avg, na.rm=T),
                                          median(merged_data$etoh10_avg, na.rm=T)))

### Component plot

diet_median_plot <- ggplot(diet_dist_dat) +
  geom_col(aes(x = reorder(component, comp_median), y = comp_median, fill = component_cat), 
           position = "dodge2", show.legend = TRUE, alpha = 0.9) +
  geom_text(aes(x = reorder(component, comp_median), y = comp_median, label = round(comp_median, 1)), 
            hjust = -0.2,  # Adjust vertical position
            color = "black", 
            size = 4.5) +
  labs(y="Median intake") +
  scale_fill_manual(values = c("#A8A7A4", "#D89E5F", "#FAD02E", "#9D5C63", "#2C7BB6", "#7D8C3C", "#D8A8B8")) +
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(color="black",size=15),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "gray12", size = 14, hjust = 1),
    legend.position = "none") +
  theme(
    panel.background = element_rect(color = "white"),
    plot.margin = margin(10, 10, 10, 10)) +
  coord_flip() 

pdf("results/figures/fig1e_diet_median_distribution_10132024.pdf", width = 7, height = 4.5, onefile = F)
diet_median_plot
dev.off()

### Histogram

diet_hist_plot <- ggplot(merged_data, aes(x = AMED_avg)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.5, fill = "#69b3a2", color = "black", alpha = 0.7) +
  geom_density(alpha = 0.2, fill = "#FF6347") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 36, color = "black", margin = margin(t = 5, b = -10)),
    axis.text.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 36, face = "bold") 
  ) +
  labs(title = "MedDiet score", x = "", y = "")

pdf("results/figures/fig1e_diet_histogram_10132024.pdf", width = 7, height = 5, onefile = F)
diet_hist_plot
dev.off()
