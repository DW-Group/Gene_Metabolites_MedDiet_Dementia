n_missing <- as.data.frame(rbind(c("finn-b-KRA_PSY_DEMENTIA_EXMORE",5933+166584),
                                 c("finn-b-F5_VASCDEM",881+211508)))
colnames(n_missing) <- c("id.outcome","samplesize.outcome")
n_missing$samplesize.outcome <- as.numeric(n_missing$samplesize.outcome)