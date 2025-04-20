age_covar_list <- c("agemo")
blood_covar_list <- c(age_covar_list,"blddate","bmi","fast1","endpoint","caco")
blood_covar_reduced_list <- c(age_covar_list,"blddate","bmi")
full_covar_list <- c(blood_covar_list,"hiedu3","husbedu2_3_m","phmsstatus3_4","fhdem",
                     "actcon","nSES","marry","smkk2","smkk3_4_5","dep_antidep","htn",
                     "sysbp2","sysbp3","sysbpm","hchol","alcocon","daykcal") 
full_covar_reduced_list <- c(blood_covar_reduced_list,"actcon","nSES","alcocon","daykcal") 
apoe4_covar_list <- c("APOE4_ncopy_1","APOE4_ncopy_2","PC1","PC2","PC3","PC4","platformgsa","platformhuco2","platformillu","platformomni","platformonco")
prs_covar_list <- c("PGS002280","PC1","PC2","PC3","PC4","platformgsa","platformhuco2","platformillu","platformomni","platformonco")
prs_covar_reduced_list <- c("PGS002280","PC1","PC2","PC3","PC4")

covar_list <- list("age" = age_covar_list,
                   "blood" = blood_covar_list,
                   "full" = full_covar_list,
                   "apoe4" = unique(c(full_covar_list,apoe4_2cat_covar_list)),
                   "prs" = unique(c(full_covar_list,prs_covar_list)),
                   "prs_reduced" = unique(c(full_covar_reduced_list,prs_covar_reduced_list)),
                   "apoe4_prs" = unique(c(full_covar_list,prs_covar_list,apoe4_2cat_covar_list)))
gene_covar_list <- unique(c(full_covar_list,"PC1","PC2","PC3","PC4","platformgsa","platformhuco2","platformillu","platformomni","platformonco"))