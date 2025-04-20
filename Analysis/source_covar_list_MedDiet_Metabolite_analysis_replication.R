age_covar_list <- c("agemo")
blood_covar_list <- c(age_covar_list,"blddate","bmi","fast1","endpoint","caco")
full_covar_list <- c(blood_covar_list,"proff2","proff3","fhdem",
                     "actcon","nSES","marry","smkk2","smkk3_4_5","dep_antidep","htn",
                     "sysbp2","sysbp3","sysbpm","hchol","alcocon","daykcal") 

covar_list <- list("age" = age_covar_list,
                   "blood" = blood_covar_list,
                   "full" = full_covar_list)
