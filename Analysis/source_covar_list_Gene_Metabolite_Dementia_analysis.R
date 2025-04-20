covar_dem_list <- list()
covar_dem_list[["blood_covar"]] <- c("blddate","agemo","bmi","fast1","strata(endpoint, caco)")
covar_dem_list[["full_covar"]] <- c(covar_dem_list[["blood_covar"]],
                                    "hiedu3","husbedu2_3_m","phmsstatus3_4","fhdem","actcon","nSES","marry","smkk2","smkk3_4_5","dep_antidep","htn","sysbp2","sysbp3","sysbpm","hchol","alcocon","AMED_avg","daykcal")
covar_dem_list[["APOE4_carrier"]] <- c(covar_dem_list[["full_covar"]],
                                       "APOE4_carrier","PC1","PC2","PC3","PC4","platformgsa","platformhuco2","platformillu","platformomni","platformonco")
covar_dem_list[["APOE4_2cat"]] <- c(covar_dem_list[["full_covar"]],
                                    "APOE4_ncopy_1","APOE4_ncopy_2","PC1","PC2","PC3","PC4","platformgsa","platformhuco2","platformillu","platformomni","platformonco")
covar_dem_list[["PGS002280"]] <- c(covar_dem_list[["full_covar"]],
                                                "PGS002280","PC1","PC2","PC3","PC4","platformgsa","platformhuco2","platformillu","platformomni","platformonco")
covar_dem_list[["APOE4_carrier_PGS002280"]] <- c(covar_dem_list[["full_covar"]],
                                                                "PGS002280","APOE4_carrier","PC1","PC2","PC3","PC4","platformgsa","platformhuco2","platformillu","platformomni","platformonco")
covar_dem_list[["APOE4_2cat_PGS002280"]] <- c(covar_dem_list[["full_covar"]],
                                                            "PGS002280","APOE4_ncopy_1","APOE4_ncopy_2","PC1","PC2","PC3","PC4","platformgsa","platformhuco2","platformillu","platformomni","platformonco")
covar_dem_list[["PGS000334"]] <- c(covar_dem_list[["full_covar"]],
                                                "PGS000334","PC1","PC2","PC3","PC4","platformgsa","platformhuco2","platformillu","platformomni","platformonco")

covar_tics_list <- list()
covar_tics_list[["blood_covar"]] <- c("blddate","agemo","bmi","fast1","endpoint","caco")
covar_tics_list[["full_covar"]] <- c(covar_tics_list[["blood_covar"]],
                                     "hiedu3","husbedu2_3_m","phmsstatus3_4","fhdem","actcon","nSES","marry","smkk2","smkk3_4_5","dep_antidep","htn","sysbp2","sysbp3","sysbpm","hchol","alcocon","AMED_avg","daykcal")
covar_tics_list[["APOE4_carrier"]] <- c(covar_tics_list[["full_covar"]],
                                        "APOE4_carrier","PC1","PC2","PC3","PC4","platformhuco2","platformillu","platformomni","platformonco")