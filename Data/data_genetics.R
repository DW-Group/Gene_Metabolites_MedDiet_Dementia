library(data.table)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggpp)
library(patchwork)

##-----------------------------------------
## PRS 
##-----------------------------------------

rm(list=ls())

apoe4 <- fread("data/genetic/apoe_with_gsa_after_exclusion.txt")

PGS002280_raw_list <- list()
for (pf in c("affy","onco","omni","illu","huco","gsa")){
  PGS002280_raw_list[[pf]] <- fread(paste0("data/genetic/raw_prs/sum_1000G_PGS002280_PRS_",pf,".profile")) %>% 
    mutate(platform=gsub("huco","huco2",pf),
           platform_indID = paste0(platform,"_",IID),
           PGS002280_scaled = scale(SCORESUM))
}
PGS002280_scaled <- as.data.frame(do.call("rbind",PGS002280_raw_list))
PGS002280_scaled <- merge(PGS002280_scaled, apoe4, by="platform_indID", suffixes=c("","_apoe4")) %>% select(c(platform_indID,PGS002280_scaled))

PGS000334_raw_list <- list()
for (pf in c("affy","onco","omni","illu","huco","gsa")){
  PGS000334_raw_list[[pf]] <- fread(paste0("data/genetic/raw_prs/sum_1000G_PGS000334_PRS_",pf,".profile")) %>% 
    mutate(platform=gsub("huco","huco2",pf),
           platform_indID = paste0(platform,"_",IID),
           PGS000334_scaled = scale(SCORESUM))
}
PGS000334_scaled <- as.data.frame(do.call("rbind",PGS000334_raw_list))
PGS000334_scaled <- merge(PGS000334_scaled, apoe4, by="platform_indID", suffixes=c("","_apoe4")) %>% select(c(platform_indID,PGS000334_scaled))

write.table(PGS002280_scaled, "data/genetic/PGS002280_scaled_PRS_1000G_after_exclusion.txt", row.names=F, quote=F)
write.table(PGS000334_scaled, "data/genetic/PGS000334_scaled_PRS_1000G_after_exclusion.txt", row.names=F, quote=F)

##-----------------------------------------
## PC (combined)
##-----------------------------------------

rm(list=ls())

### Projected PCs

pcs_combined_affy <- fread("PCA_projected1KG/Top10PCs_Affy_Studies.txt") %>% 
  mutate(platform = "affy", IID = as.character(IID))
pcs_combined_illu <- fread("PCA_projected1KG/Top10PCs_Illumina_Studies.txt") %>% 
  mutate(platform = "illu", IID = as.character(IID))
pcs_combined_omni <- fread("PCA_projected1KG/Top10PCs_OmniExpress_Studies.txt") %>% 
  mutate(platform = "omni", IID = as.character(IID))
pcs_combined_onco <- fread("PCA_projected1KG/Top10PCs_OncoArray_Studies.txt") %>% 
  mutate(platform = "onco", IID = as.character(IID))
pcs_combined_huco2 <- fread("PCA_projected1KG/Top10PCs_HumanCore2_Studies.txt") %>% 
  mutate(platform = "huco2", IID = as.character(IID))
pcs_combined_gsa <- fread("PCA_projected1KG/Top10PCs_GSA_Studies.txt") %>% 
  mutate(platform = "gsa", IID = as.character(IID))

pcs_combined <- as.data.frame(do.call("rbind",list(pcs_combined_affy,
                                                   pcs_combined_illu,
                                                   pcs_combined_omni,
                                                   pcs_combined_onco,
                                                   pcs_combined_huco2,
                                                   pcs_combined_gsa))) %>%
  mutate(platform_indID = paste0(platform,"_",IID)) %>% select(-c(IID,platform))
colnames(pcs_combined)[grepl("^PC", colnames(pcs_combined))] <- paste0(colnames(pcs_combined)[grepl("^PC", colnames(pcs_combined))], "_comb")

### Save

apoe4 <- fread("data/genetic/apoe_with_gsa_after_exclusion.txt")
pcs_combined <- merge(apoe4, pcs_combined, by="platform_indID") %>% select(c(platform_indID,PC1_comb,PC2_comb,PC3_comb,PC4_comb,PC5_comb,PC6_comb,PC7_comb,PC8_comb,PC9_comb,PC10_comb))
write.table(pcs_combined, "data/genetic/combpcs_after_exclusion.txt", row.names=F, quote=F)

apoe4 <- read.csv("data/genetic/apoe_with_gsa_after_exclusion.csv")
pcs_combined <- fread("data/genetic/combpcs_after_exclusion.txt")
apoe4_comb_pcs <- merge(apoe4, pcs_combined, by="platform_indID") %>% dplyr::select(-c(FID,famID,IID,indID,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10))
write.csv(apoe4_comb_pcs, "data/genetic/apoe_with_gsa_comb_pcs_after_exclusion.csv", row.names=F, quote=F)

##-----------------------------------------
## All genetic data
##-----------------------------------------

apoe <- fread("data/genetic/apoe_with_gsa_after_exclusion.txt")
pcs_comb <- fread("data/genetic/combpcs_after_exclusion.txt")
ad_var <- fread("data/genetic/ad_variants_after_exclusion.txt")
PGS002280 <- fread("data/genetic/PGS002280_scaled_PRS_1000G_after_exclusion.txt"); PGS002280 <- subset(PGS002280, select=-CNT)
PGS000334 <- fread("data/genetic/PGS000334_scaled_PRS_1000G_after_exclusion.txt"); PGS000334 <- subset(PGS000334, select=-CNT)

geno_all <- Reduce(full_join, list(apoe,
                                   PGS002280,PGS000334,
                                   ad_var,pcs_comb))

save(geno_all, file="data/genetic/geno_all.RData")