library(data.table)
library(dplyr)
library(stringr)
library(readxl)
library(writexl)
library(MendelianRandomization)
library(TwoSampleMR)
library(MRInstruments)
library(MRPRESSO)

rm(list=ls())

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

##-----------------------------------------
## Load results
##-----------------------------------------

load("data/MR/ieu_open_gwas_all_outcomes.RData")

load("results/MR/mr_results.RData")
load("results/MR/mr_results_ivw_ld.RData")
load("results/MR/mr_results_sensitivity.RData")

##-----------------------------------------
## Load instrument info
##-----------------------------------------

### Candidate instruments

chen_supp_met <- read_excel("data/MR/instruments/Chen_NatGenet_2023/Chen_NatGenet_2023_all_identified_associations_metabolites.xlsx") %>%
  subset(select=c("SNP","Metabolite",
                  "BETA","SE","Effect Allele","Non-Effect Allele","FREQ","P*","Effector genes","N",
                  "CHR","Position(hg38)","liftover_output_b37"))
chen_supp_met_ratio <- read_excel("data/MR/instruments/Chen_NatGenet_2023/Chen_NatGenet_2023_all_identified_associations_metabolite_ratio.xlsx") %>%
  subset(select=c("SNP","Metabolite ratios",
                  "BETA","SE","Effect Allele","Non-Effect Allele","FREQ","P*","Effector genes","N",
                  "CHR","POS","liftover_output_b37"))

### All instruments used in the final analysis

load("data/MR/harmonized.RData")

harmonized_outcome_local_chen_paper <- as.data.frame(do.call("rbind",harmonized_outcome_local_chen_paper_list)) %>%
  mutate(exposure_source = "chen_paper", outcome_source = "local")
harmonized_outcome_local_chen_supp_met <- as.data.frame(do.call("rbind",harmonized_outcome_local_chen_supp_met_list)) %>%
  mutate(exposure_source = "chen_supp_met", outcome_source = "local")
harmonized_outcome_local_chen_supp_met_ratio <- as.data.frame(do.call("rbind",harmonized_outcome_local_chen_supp_met_ratio_list)) %>%
  mutate(exposure_source = "chen_supp_met_ratio", outcome_source = "local")

harmonized_outcome_local <- as.data.frame(do.call("rbind.fill", list(harmonized_outcome_local_chen_paper,
                                                                     harmonized_outcome_local_chen_supp_met,
                                                                     harmonized_outcome_local_chen_supp_met_ratio)))

harmonized_outcome_server_chen_paper <- as.data.frame(do.call("rbind",harmonized_outcome_server_chen_paper_list)) %>%
  mutate(exposure_source = "chen_paper", outcome_source = "server")
harmonized_outcome_server_chen_supp_met <- as.data.frame(do.call("rbind",harmonized_outcome_server_chen_supp_met_list)) %>%
  mutate(exposure_source = "chen_supp_met", outcome_source = "server")
harmonized_outcome_server_chen_supp_met_ratio <- as.data.frame(do.call("rbind",harmonized_outcome_server_chen_supp_met_ratio_list)) %>%
  mutate(exposure_source = "chen_supp_met_ratio", outcome_source = "server")

harmonized_outcome_server <- as.data.frame(do.call("rbind.fill", list(harmonized_outcome_server_chen_paper,
                                                                      harmonized_outcome_server_chen_supp_met,
                                                                      harmonized_outcome_server_chen_supp_met_ratio)))

instruments_info <- as.data.frame(rbind.fill(harmonized_outcome_local,harmonized_outcome_server)) %>%
  mutate(id = paste0(outcome, "_", exposure_source, "_", exposure))

##-----------------------------------------
## Filtering annotation
##-----------------------------------------

selected_outcome <- c("wightman_ad",
                      "Any dementia (more controls excluded) || id:finn-b-KRA_PSY_DEMENTIA_EXMORE",
                      "Cognitive performance || id:ebi-a-GCST006572",
                      "Vascular dementia (F5_VASCDEM) || id:finn-b-F5_VASCDEM")

selected_exposure_source <- c("chen_paper","chen_supp_met","chen_supp_met_ratio")

selected_methods <- c("Inverse variance weighted","MR Egger", "Wald ratio")

instruments_info_selected <- instruments_info %>%
  filter(outcome %in% selected_outcome &
           exposure_source %in% selected_exposure_source) %>%
  filter(!(exposure %in% names(harmonized_outcome_local_chen_paper_list) &
             exposure_source %in% c("chen_supp_met","chen_supp_met_ratio"))) %>%
  mutate(exposure_type = case_when(exposure %in% chen_supp_met$Metabolite ~ "met",
                                   exposure %in% chen_supp_met_ratio$`Metabolite ratios` ~ "ratio"))

mr_res_sub <- mr_res_all %>%
  filter(outcome %in% selected_outcome &
         exposure_source %in% selected_exposure_source &
         method %in% selected_methods) %>%
  filter(!(exposure %in% names(harmonized_outcome_local_chen_paper_list) &
             exposure_source %in% c("chen_supp_met","chen_supp_met_ratio"))) %>%
  mutate(exposure_type = case_when(exposure %in% chen_supp_met$Metabolite ~ "met",
                                   exposure %in% chen_supp_met_ratio$`Metabolite ratios` ~ "ratio"))

##-----------------------------------------
## Sensitivity analysis results annotation
##-----------------------------------------

mr_res_sub <- mr_res_sub %>%
  mutate(id = paste0(outcome, "_", exposure_source, "_", exposure),
         unique_id = paste0(outcome, "_", exposure_source, "_", exposure,"_", method))

mr_egger_res <- mr_res_sub %>%
  filter(method == "MR Egger") %>%
  subset(select=c(id, b, se, pval))
colnames(mr_egger_res) <- c("id", "b_mr_egger", "se_mr_egger", "pval_mr_egger")

mr_pleio_res <- mr_pleio_all %>%
  mutate(id = paste0(outcome, "_", exposure_source, "_", exposure)) %>%
  filter(id %in% mr_res_sub$id) %>%
  subset(select=c(id, egger_intercept, se, pval))
colnames(mr_pleio_res) <- c("id", "egger_intercept", "se_pleio", "pval_pleio")

mr_heter_res <- mr_heter_all %>%
  mutate(unique_id = paste0(outcome, "_", exposure_source, "_", exposure,"_", method)) %>%
  filter(unique_id %in% mr_res_sub$unique_id) %>%
  subset(select=c(unique_id, Q, Q_df, Q_pval))

mr_ivw_ld_res <- mr_ivw_ld_res_all %>%
  mutate(id = paste0(outcome, "_", exposure_source, "_", exposure)) %>%
  filter(id %in% mr_res_sub$id) %>%
  subset(select=c(id, beta, se, p, nsnp))
colnames(mr_ivw_ld_res) <- c("id", "b_ivw_ld", "se_ivw_ld", "pval_ivw_ld","nsnp_ivw_ld")

mr_steiger_res <- mr_steiger_all %>%
  mutate(id = paste0(outcome, "_", exposure_source, "_", exposure)) %>%
  filter(id %in% mr_res_sub$id) %>%
  subset(select=c(id, snp_r2.exposure, snp_r2.outcome, correct_causal_direction, steiger_pval))

mr_res_sub_anno <- merge(mr_res_sub, mr_egger_res, by="id", all.x=T)
mr_res_sub_anno <- merge(mr_res_sub_anno, mr_pleio_res, by="id", all.x=T)
mr_res_sub_anno <- merge(mr_res_sub_anno, mr_heter_res, by="unique_id", all.x=T)
mr_res_sub_anno <- merge(mr_res_sub_anno, mr_ivw_ld_res, by="id", all.x=T)
mr_res_sub_anno <- merge(mr_res_sub_anno, mr_steiger_res, by="id", all.x=T)

instruments_selected <- instruments_info_selected %>%
  filter(id %in% mr_res_sub_anno$id) %>%
  mutate(exposure_SNP = paste0(exposure,"_",SNP))

instruments_selected_summary <- instruments_selected %>%
  group_by(id) %>%
  reframe(nsnps_harmonized = n_distinct(SNP),
          snps_harmonized = paste(SNP, collapse=","))

mr_res_sub_anno <- merge(mr_res_sub_anno, instruments_selected_summary, by="id")

##-----------------------------------------
## Filtering
##-----------------------------------------

mr_res_sub_anno_filtered <- mr_res_sub_anno %>% filter(correct_causal_direction) 

mr_res_sub_anno_filtered <- mr_res_sub_anno_filtered %>%
  mutate(keep = case_when(pval_pleio < 0.05 & method == "MR Egger" ~ T,
                          pval_pleio < 0.05 & method == "Inverse variance weighted" ~ F,
                          pval_pleio >= 0.05 & method == "MR Egger" ~ F,
                          pval_pleio >= 0.05 & method == "Inverse variance weighted" ~ T,
                          is.na(pval_pleio) & method == "MR Egger" ~ F,
                          is.na(pval_pleio) & method == "Inverse variance weighted" ~ T,
                          is.na(pval_pleio) & method == "Wald ratio" ~ T))

mr_res_sub_anno_filtered <- mr_res_sub_anno_filtered %>% filter(keep)

##-----------------------------------------
## Multiple testing correction
##-----------------------------------------

### Bonferroni

mr_res_sub_anno_filtered <- mr_res_sub_anno_filtered %>%
  mutate(Bonf_p = 0.05 / n(),
         Bonf_sig = case_when(pval <= Bonf_p ~ TRUE,
                                      pval > Bonf_p ~ FALSE))

### FDR - BH

mr_res_sub_anno_filtered <- mr_res_sub_anno_filtered %>%
  mutate(FDR_BH_p = p.adjust(pval, "BH"),
         FDR_BH_sig = FDR_BH_p <= 0.05)

##-----------------------------------------
## Exposure and outcome annotation 
##-----------------------------------------

### Exposure

chen_met_dict <- read_excel("data/MR/instruments/Chen_NatGenet_2023/Chen_NatGenet_2023_metabolite_dictionary.xlsx")

mr_res_sub_anno_filtered <- merge(mr_res_sub_anno_filtered, chen_met_dict, by.x="exposure", by.y="Metabolites", all.x=T) %>%
  mutate(classification = ifelse(is.na(SUPER_PATHWAY) & grepl(" / ", exposure), "Metabolite ratio", SUPER_PATHWAY))

### Outcome

unique(mr_res_sub_anno_filtered$outcome)
mr_out_dict <- data.frame(outcome = c("Any dementia (more controls excluded) || id:finn-b-KRA_PSY_DEMENTIA_EXMORE",
                                      "wightman_ad",
                                      "Vascular dementia (F5_VASCDEM) || id:finn-b-F5_VASCDEM",
                                      "Cognitive performance || id:ebi-a-GCST006572"),
                          anno_outcome = c("Dementia",
                                           "Alzheimer's disease",
                                           "Vascular dementia",
                                           "Cognitive performance"),
                          cat_outcome = c(rep("Dementia and subtypes",3),
                                          rep("Cognitive function",1)))

detach("package:plyr", unload = TRUE)
outcome_sig_summary <- mr_res_sub_anno_filtered %>% 
  filter(FDR_BH_sig) %>%
  group_by(exposure) %>%
  reframe(n_sig_outcome = n_distinct(outcome),
          sig_outcome = paste(outcome, collapse=","))

mr_res_sub_anno_filtered <- merge(mr_res_sub_anno_filtered, mr_out_dict, by="outcome")
mr_res_sub_anno_filtered <- merge(mr_res_sub_anno_filtered, outcome_sig_summary, by="exposure", all.x=T)
instruments_selected <- merge(instruments_selected, mr_out_dict, by="outcome")