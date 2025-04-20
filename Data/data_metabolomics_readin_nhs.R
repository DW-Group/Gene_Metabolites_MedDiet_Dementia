library(chanmetab)
library(Biobase)
library(dplyr)
library(tidyr)
library(plyr)
library(stringr)

rm(list=ls())

##---------------------------------
## Load data
##---------------------------------

endpoint_list_all <- c("diabetes","poag","rheumatoid","ibd","colon","prostate","ovarian","breast","racial.diff","exfoliation_glaucoma","stress.method","poag") 
endpoint_list_control <- c("pancreatic","parkinsons","stroke","als")
endpoint_list_substudy <- c("mbs","lvs")

all_data <- merge_metab_data(cohorts = "nhs1",
                             endpoints = endpoint_list_all,
                             methods = c("C8-pos", "HILIC-neg", "HILIC-pos", "C18-neg"),
                             collection_to_use = "first",
                             transformation = "transform_probit_score",
                             impute_cutoff = 0.25,
                             impute_missing_function = "impute_one_half_min",
                             keep_failed_pm_metabolites = F,
                             combine_cohorts = TRUE)

control_data <- merge_metab_data(cohorts = "nhs1",
                                 endpoints = endpoint_list_control,
                                 methods = c("C8-pos", "HILIC-neg", "HILIC-pos", "C18-neg"),
                                 collection_to_use = "first",
                                 controls_only = T,
                                 transformation = "transform_probit_score",
                                 impute_cutoff = 0.25,
                                 impute_missing_function = "impute_one_half_min",
                                 keep_failed_pm_metabolites = F,
                                 combine_cohorts = TRUE)

substudy_data <- merge_metab_data(cohorts = "nhs1",
                                  endpoints = endpoint_list_substudy,
                                  methods = c("C8-pos", "HILIC-neg", "HILIC-pos", "C18-neg"),
                                  collection_to_use = "substudy",
                                  transformation = "transform_probit_score",
                                  impute_cutoff = 0.25,
                                  impute_missing_function = "impute_one_half_min",
                                  keep_failed_pm_metabolites = F,
                                  combine_cohorts = TRUE)

save.image("data/metabolites/metabolomics_data_raw_nhs.RData")

##---------------------------------
## Process data
##---------------------------------

rm(list=ls())

load("data/metabolites/metabolomics_data_raw_nhs.RData")

### Separate data

pdata_all <- pData(all_data$expr_set$all_cohorts) %>%
  mutate(source_data = "all")
pdata_control <- pData(control_data$expr_set$all_cohorts) %>%
  mutate(source_data = "control")
pdata_substudy <- pData(substudy_data$expr_set$all_cohorts) %>%
  mutate(source_data = "substudy")

length(pdata_all$id); length(unique(pdata_all$id)); length(paste0(pdata_all$cohort,"_",pdata_all$id))
length(pdata_control$id); length(unique(pdata_control$id)); length(paste0(pdata_control$cohort,"_",pdata_control$id))
length(pdata_substudy$id); length(unique(pdata_substudy$id)); length(paste0(pdata_substudy$cohort,"_",pdata_substudy$id))

prdata_all <- pData(protocolData(all_data$expr_set$all_cohorts)) %>%
  mutate(source_data = "all")
prdata_control <- pData(protocolData(control_data$expr_set$all_cohorts)) %>%
  mutate(source_data = "control")
prdata_substudy <- pData(protocolData(substudy_data$expr_set$all_cohorts)) %>%
  mutate(source_data = "substudy")

fdata_all <- fData(all_data$expr_set$all_cohorts) %>%
  mutate(source_data = "all")
fdata_control <- fData(control_data$expr_set$all_cohorts) %>%
  mutate(source_data = "control")
fdata_substudy <- fData(substudy_data$expr_set$all_cohorts) %>%
  mutate(source_data = "substudy")

mdata_all <- exprs(all_data$expr_set$all_cohorts)
mdata_control <- exprs(control_data$expr_set$all_cohorts)
mdata_substudy <- exprs(substudy_data$expr_set$all_cohorts)

### Remove duplicated metabolites across methods

fdata <- as.data.frame(rbind.fill(fdata_all, fdata_control, fdata_substudy)) %>%
  mutate(hmdb_id_method = paste0(hmdb_id,"_",method))

length(unique(fdata$hmdb_id)); length(unique(fdata$hmdb_id_method))
tmp1 <- aggregate(method ~ hmdb_id, fdata, n_distinct); table(tmp1$method)
tmp2 <- aggregate(cbind(mean_icc, mean_cv) ~ hmdb_id + method, fdata, mean)
tmp3 <- tmp2[tmp2$hmdb_id %in% tmp1[tmp1$method==2,]$hmdb_id,]
tmp4 <- tmp3 %>% group_by(hmdb_id) %>% slice(which.min(mean_cv)) %>% mutate(hmdb_id_method = paste0(hmdb_id,"_",method))

fdata_sub <- fdata %>%
  filter(!(hmdb_id %in% tmp4$hmdb_id) | (hmdb_id_method %in% tmp4$hmdb_id_method))
table(aggregate(method ~ hmdb_id, fdata_sub, n_distinct)$method)

fdata_sub_agg <- aggregate(cbind(mean_icc, mean_cv) ~ hmdb_id, fdata_sub, mean)

### Remove duplicated samples

pdata <- as.data.frame(rbind.fill(pdata_all, pdata_control, pdata_substudy)) %>%
  mutate(id = substr(id, 1, 6),
         study_id = paste0(cohort,"_",id),
         study_id_endpoint_source_data = paste0(study_id,"_",endpoint,"_",source_data))
length(pdata$id); length(unique(pdata$id)); length(unique(pdata$study_id)); length(unique(pdata$study_id_endpoint_source_data))

tmp1 <- aggregate(source_data ~ study_id, pdata, n_distinct); table(tmp1$source_data)

tmp2_all <- as.data.frame(colSums(!is.na(mdata_all))); colnames(tmp2_all) <- "n"
table(rownames(tmp2_all) == pdata_all$id)
tmp2_all <- tmp2_all %>%
  mutate(id = substr(rownames(tmp2_all), 1, 6),
         cohort = pdata_all$cohort,
         endpoint = pdata_all$endpoint,
         source_data = "all")
tmp2_control <- as.data.frame(colSums(!is.na(mdata_control))); colnames(tmp2_control) <- "n"
table(rownames(tmp2_control) == pdata_control$id)
tmp2_control <- tmp2_control %>%
  mutate(id = substr(rownames(tmp2_control), 1, 6),
         cohort = pdata_control$cohort,
         endpoint = pdata_control$endpoint,
         source_data = "control")
tmp2_substudy <- as.data.frame(colSums(!is.na(mdata_substudy))); colnames(tmp2_substudy) <- "n"
table(rownames(tmp2_substudy) == pdata_substudy$id)
tmp2_substudy <- tmp2_substudy %>%
  mutate(id = substr(rownames(tmp2_substudy), 1, 6),
         cohort = pdata_substudy$cohort,
         endpoint = pdata_substudy$endpoint,
         source_data = "substudy")
tmp2 <- as.data.frame(rbind(tmp2_all, tmp2_control, tmp2_substudy)) %>%
  mutate(study_id = paste0(cohort,"_",id))

tmp3 <- tmp2[tmp2$study_id %in% tmp1[tmp1$source_data > 1,]$study_id,]

tmp4 <- tmp3 %>% group_by(study_id) %>% slice(which.max(n)) %>% mutate(study_id_endpoint_source_data = paste0(study_id,"_",endpoint,"_",source_data))

pdata_sub <- pdata %>%
  filter(!(study_id %in% tmp4$study_id) | (study_id_endpoint_source_data %in% tmp4$study_id_endpoint_source_data))
table(aggregate(source_data ~ study_id, pdata_sub, n_distinct)$source_data)

pdata_all <- pdata_all %>% mutate(id = substr(id, 1, 6),
                                  study_id = paste0(cohort,"_",id),
                                  study_id_endpoint_source_data = paste0(study_id,"_",endpoint,"_",source_data),
                                  keep = !(study_id %in% tmp4$study_id) | (study_id_endpoint_source_data %in% tmp4$study_id_endpoint_source_data))
table(rownames(pdata_all) == colnames(mdata_all))
table(pdata_all$keep)

pdata_control <- pdata_control %>% mutate(id = substr(id, 1, 6),
                                          study_id = paste0(cohort,"_",id),
                                          study_id_endpoint_source_data = paste0(study_id,"_",endpoint,"_",source_data),
                                          keep = !(study_id %in% tmp4$study_id) | (study_id_endpoint_source_data %in% tmp4$study_id_endpoint_source_data))
table(rownames(pdata_control) == colnames(mdata_control))
table(pdata_control$keep)

pdata_substudy <- pdata_substudy %>% mutate(id = substr(id, 1, 6),
                                            study_id = paste0(cohort,"_",id),
                                            study_id_endpoint_source_data = paste0(study_id,"_",endpoint,"_",source_data),
                                            keep = !(study_id %in% tmp4$study_id) | (study_id_endpoint_source_data %in% tmp4$study_id_endpoint_source_data))
table(rownames(pdata_substudy) == colnames(mdata_substudy))
table(pdata_substudy$keep)

### Combine data

mdata_all_t <- data.frame(t(mdata_all)) %>%
  subset(select=gsub("\\*","",fdata_sub[fdata_sub$source_data=="all",]$hmdb_id))
mdata_all_t <- mdata_all_t[pdata_all$keep,]
identical(substr(rownames(mdata_all_t), 1, 6), pdata_sub[pdata_sub$source_data=="all",]$id) # should be TRUE
mdata_all_final <- cbind(mdata_all_t,pdata_sub[pdata_sub$source_data=="all",])

mdata_control_t <- data.frame(t(mdata_control)) %>%
  subset(select=gsub("\\*","",fdata_sub[fdata_sub$source_data=="control",]$hmdb_id))
mdata_control_t <- mdata_control_t[pdata_control$keep,]
identical(substr(rownames(mdata_control_t), 1, 6), pdata_sub[pdata_sub$source_data=="control",]$id) # should be TRUE
mdata_control_final <- cbind(mdata_control_t,pdata_sub[pdata_sub$source_data=="control",])

mdata_substudy_t <- data.frame(t(mdata_substudy)) %>%
  subset(select=gsub("\\*","",fdata_sub[fdata_sub$source_data=="substudy",]$hmdb_id))
mdata_substudy_t <- mdata_substudy_t[pdata_substudy$keep,]
identical(substr(rownames(mdata_substudy_t), 1, 6), pdata_sub[pdata_sub$source_data=="substudy",]$id) # should be TRUE
mdata_substudy_final <- cbind(mdata_substudy_t,pdata_sub[pdata_sub$source_data=="substudy",])

mdata_final <- rbind.fill(mdata_all_final, mdata_control_final, mdata_substudy_final)

table(mdata_final$source_data)
table(mdata_final$cohort)
table(mdata_final$endpoint, mdata_final$caco) 

### Modify ID

mdata_final <- mdata_final %>%
  mutate(id = substr(id, 1, 6),
         study_id <- paste0(cohort,"_",id))

##---------------------------------
## Save data
##---------------------------------

write.csv(mdata_final, "data/metabolites/metabolomics_data_full_nhs.csv", na="", row.names=F, quote=F)

fdata_all <- fdata_sub[fdata_sub$source_data=="all",]; save(fdata_all, file="data/metabolites/fdata_all_nhs.RData")
fdata_control <- fdata_sub[fdata_sub$source_data=="control",]; save(fdata_control, file="data/metabolites/fdata_control_nhs.RData")
fdata_substudy <- fdata_sub[fdata_sub$source_data=="substudy",]; save(fdata_substudy, file="data/metabolites/fdata_substudy_nhs.RData")

save(fdata_sub, file="data/metabolites/fdata_full_nhs.RData")
save(fdata_sub_agg, file="data/metabolites/fdata_full_agg_nhs.RData")

tmp1 <- fdata_sub %>%
  mutate(hmdb_id = gsub("\\*", "", hmdb_id)) %>%
  group_by(hmdb_id) %>%
  summarise(mean_icc_final = mean(mean_icc, na.rm=T),
            mean_cv_final = mean(mean_cv, na.rm=T))
tmp2 <- fdata_sub %>%
  mutate(hmdb_id = gsub("\\*", "", hmdb_id)) %>%
  filter(!duplicated(hmdb_id))
fdata <- tmp2; fdata$mean_cv_final <- tmp1$mean_cv_final; fdata$mean_icc_final <- tmp1$mean_icc_final
save(fdata, file="data/metabolites/fdata_nhs.RData")

save(prdata_all, file="data/metabolites/prdata_nhs.RData")
save(prdata_control, file="data/metabolites/prdata_control_nhs.RData")
save(prdata_substudy, file="data/metabolites/prdata_substudy_nhs.RData")
