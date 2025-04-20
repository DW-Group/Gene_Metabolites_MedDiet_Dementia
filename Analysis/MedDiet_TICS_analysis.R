library(splines)
library(dplyr)
library(data.table)

rm(list=ls())

### Load data

tics_dat <- read.csv("data/dementia/nhs_tics_full.csv")
dim(tics_dat); length(unique(tics_dat$id))

load("data/genetic/geno_all.RData")

### Covar list

tics_covar <- c("agecon","hiedu2","hiedu3","qnSES1","qnSES2","actcc2","actcc3","actcc4","actcc5","smk2","smk3",
                "hbp","dprs","anti","db","bmic1","bmic3","bmic4","bmic5",
                "husbedu2","husbedu3","fhdem","nhor2","nhorm","marry","live_alone")
  
### Fit linear model

lm_ntotal <- lm(reformulate(c("amedcon",tics_covar), response="av_ntotal"), data = tics_dat)
lm_nverbl <- lm(reformulate(c("amedcon",tics_covar), response="av_nverbl"), data = tics_dat)
lm_zscre <- lm(reformulate(c("amedcon",tics_covar), response="av_zscre"), data = tics_dat)

summary(lm_ntotal)$coefficients["amedcon",]
summary(lm_nverbl)$coefficients["amedcon",]
summary(lm_zscre)$coefficients["amedcon",]

### Fit linear model by APOE4 and PRS

tics_genetic <- merge(tics_dat,geno_all, by=c("study_ID"), all.x=T, suffixes=c("","_geno")) %>%
  filter(!is.na(APOE4_geno)) %>%
  mutate(PGS002280_t = ntile(residuals(lm(PGS002280 ~ PC1 + PC2 + PC3 + PC4, data = .)), 3),
         PGS000334_t = ntile(residuals(lm(PGS000334 ~ PC1 + PC2 + PC3 + PC4, data = .)), 3))

table(tics_genetic$APOE4_geno)
table(tics_genetic$PGS002280_t)
table(tics_genetic$PGS000334_t)

tics_by_gene_res_list <- list(); i <- 1

for (tics in c("av_ntotal","av_nverbl","av_zscre")){
  
  lm_all <- glm(reformulate(c("amedcon",tics_covar), response=tics), data = tics_dat)
  lm_genetic <- glm(reformulate(c("amedcon",tics_covar), response=tics), data = tics_genetic)
  
  tics_by_gene_res_list[[i]] <- c(tics,"all","all",nobs(lm_all),summary(lm_all)$coefficients["amedcon",][c(1,2,4)]); i <- i + 1
  tics_by_gene_res_list[[i]] <- c(tics,"all","genetic",nobs(lm_genetic),summary(lm_genetic)$coefficients["amedcon",][c(1,2,4)]); i <- i + 1
  
  for (ncopy in c(0,1,2)){
    
    tics_sub_apoe <- tics_genetic %>% filter(APOE4_geno==ncopy)
    lm_apoe4 <- glm(reformulate(c("amedcon",tics_covar), response=tics), data = tics_sub_apoe)
    tics_by_gene_res_list[[i]] <- c(tics,"APOE4",ncopy,nobs(lm_apoe4),summary(lm_apoe4)$coefficients["amedcon",][c(1,2,4)]); i <- i + 1

  }
  
  for (t in 1:3){
    
    tics_sub_prs1 <- tics_genetic %>% filter(PGS002280_t==t)
    tics_sub_prs2 <- tics_genetic %>% filter(PGS000334_t==t)
    
    lm_prs1 <- glm(reformulate(c("amedcon",tics_covar), response=tics), data = tics_sub_prs1)
    lm_prs2 <- glm(reformulate(c("amedcon",tics_covar), response=tics), data = tics_sub_prs2)
    tics_by_gene_res_list[[i]] <- c(tics,"PGS002280",t,nobs(lm_prs1),summary(lm_prs1)$coefficients["amedcon",][c(1,2,4)]); i <- i + 1
    tics_by_gene_res_list[[i]] <- c(tics,"PGS000334",t,nobs(lm_prs2),summary(lm_prs2)$coefficients["amedcon",][c(1,2,4)]); i <- i + 1
    
  }
 
}

tics_by_gene_res_dat <- as.data.frame(do.call("rbind",tics_by_gene_res_list))
colnames(tics_by_gene_res_dat) <- c("tics","gene","group","n","beta","se","p")
  
write.csv(tics_by_gene_res_dat, "results/supp/amed_tics_by_gene.txt", row.names=F, quote=F)
