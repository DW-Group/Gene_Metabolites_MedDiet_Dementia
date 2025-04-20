covar_dem_list <- list()
covar_dem_list[["blood_covar"]] <- c("blddate","agemo","bmi","fast1","strata(endpoint, caco)")
covar_dem_list[["full_covar"]] <- c(covar_dem_list[["blood_covar"]],
                                    "proff2","proff3","fhdem","actcon","nSES","marry","smkk2","smkk3_4_5","dep_antidep","htn","sysbp2","sysbp3","sysbpm","hchol","alcocon","amedcon","daykcal")
covar_dem_list[["APOE4_carrier"]] <- c(covar_dem_list[["full_covar"]],
                                       "APOE4_carrier","PC1","PC2","PC3","PC4","platformhuco2","platformillu","platformomni","platformonco")
covar_dem_list[["PGS002280"]] <- c(covar_dem_list[["full_covar"]],
                                                "PGS002280","PC1","PC2","PC3","PC4","platformhuco2","platformillu","platformomni","platformonco")