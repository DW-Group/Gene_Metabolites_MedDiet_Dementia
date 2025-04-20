lifestyle_APOE4_PRS_list <- c(lifestyle_short_list, apoe4_2cat_list, prs_list)

covar_list <- list("lifestyle_short"                = c(lifestyle_short_list),
                   "lifestyle_APOE4"                = c(lifestyle_short_list, apoe4_2cat_list),
                   "lifestyle_APOE4_PRS"            = c(lifestyle_APOE4_PRS_list))

covar_all_list <- list("with_APOE4" = covar_list)
