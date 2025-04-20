library(dplyr)       
library(data.table)  
library(readxl)       
library(locuscomparer) 

rm(list=ls())

instruments_selected_sig_anno <- read_excel("results/MR/instruments_selected_sig_final_anno.xlsx")

locuscompare_plot_list <- list()

for (accession_id in unique(instruments_selected_sig_anno$accessionId)){
  
  tmp_dat1 <- instruments_selected_sig_anno %>% filter(accessionId==accession_id)
  
  for (out in gsub("6cWYtc","wightman_ad",unique(tmp_dat1$id.outcome))){
    
    tmp_dat2 <- tmp_dat1 %>% filter(id.outcome==gsub("wightman_ad","6cWYtc",out))
    
    for (rsid in unique(tmp_dat2$SNP)){
        
      tmp_out <- paste0("data/MR/colocalization/locuscomparer/locuscomparer_",out,"_",rsid,".tsv")
      tmp_exp <- paste0("data/MR/colocalization/locuscomparer/locuscomparer_chen_",accession_id,"_",rsid,".tsv")
      
      tmp_out_dat <- fread(paste0("data/MR/colocalization/locuscomparer/locuscomparer_",out,"_",rsid,".tsv"))
      tmp_exp_dat <- fread(paste0("data/MR/colocalization/locuscomparer/locuscomparer_chen_",accession_id,"_",rsid,".tsv"))
      tmp_highlight <- ifelse(((rsid %in% tmp_exp_dat$rsid) & (rsid %in% tmp_out_dat$rsid)),list(rsid),list(NULL))
      
      locuscompare_plot_list[[paste0(accession_id,"_",out,"_",rsid)]] <- locuscompare(in_fn1=tmp_out, in_fn2=tmp_exp, snp=tmp_highlight[[1]],
                                                                                      title1=unique(tmp_dat2$anno_outcome), title2=unique(tmp_dat2$exposure_anno),
                                                                                      legend_position="topleft",lz_ylab_linebreak=T)
    }
  } 
}

for (i in names(locuscompare_plot_list)) {
  pdf(paste0("results/MR/figures/mr_locuscompare_",i,".pdf"), width=10,height=5)
  print(locuscompare_plot_list[[i]])
  dev.off()
}