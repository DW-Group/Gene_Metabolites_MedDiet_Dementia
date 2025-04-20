library(locuscomparer)
library(coloc)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(readxl)
library(writexl)
library(dplyr)
library(data.table)
library(stringr)
library(DescTools)

rm(list=ls())

##-----------------------------------------
## Load data
##-----------------------------------------

mr_res_selected_anno_filtered <- read_xlsx("results/MR/mr_results_final.xlsx")
instruments_selected_sig <- read_excel("results/MR/instruments_selected_sig_final.xlsx")

##-----------------------------------------
## Extract outcome summary statistics (server)
##-----------------------------------------

outcome_server_list <- c("finn-b-F5_VASCDEM", 
                         "ebi-a-GCST006572",
                         "finn-b-KRA_PSY_DEMENTIA_EXMORE")

mr_out_dict <- data.frame(outcome = c("Any dementia (more controls excluded) || id:finn-b-KRA_PSY_DEMENTIA_EXMORE",
                                      "wightman_ad",
                                      "Vascular dementia (F5_VASCDEM) || id:finn-b-F5_VASCDEM",
                                      "Cognitive performance || id:ebi-a-GCST006572"),
                          anno_outcome = c("Dementia",
                                           "Alzheimer's disease",
                                           "Vascular dementia",
                                           "Cognitive performance"),
                          cat_outcome = c(rep("Dementia and subtypes",3),
                                          rep("Cognitive function",1)))

instruments_selected_sig_anno <- read_excel("results/MR/instruments_selected_sig_final_anno.xlsx")

for (out in outcome_server_list) {
  
  tmp_dat <- instruments_selected_sig_anno %>% filter(id.outcome==out)
  
  for (rsid in unique(tmp_dat$SNP)){
    
    tmp_chr <- tmp_dat[tmp_dat$SNP==rsid,]$chr.exposure
    tmp_pos1 <- tmp_dat[tmp_dat$SNP==rsid,]$position_b37-5e5
    tmp_pos2 <- tmp_dat[tmp_dat$SNP==rsid,]$position_b37+5e5
    
    tmp_res <- associations(variants=paste0(tmp_chr,":",tmp_pos1,"-",tmp_pos2), proxies=0, id=out, opengwas_jwt=yl_token)
    write.table(tmp_res,paste0("data/MR/colocalization/coloc_",out,"_",rsid,".tsv"), row.names=F, quote=F, sep="\t")
    
  }
  
}

##-----------------------------------------
## Harmonize outcome summary statistics 
##-----------------------------------------

unique(instruments_selected_sig_anno$outcome)

### from server

tmp_file_list1 <- list.files("data/MR/colocalization","coloc_ebi-a-GCST006572",full.names=T)
tmp_file_list2 <- list.files("data/MR/colocalization","coloc_finn-b-KRA_PSY_DEMENTIA_EXMORE",full.names=T)
tmp_file_list3 <- list.files("data/MR/colocalization","coloc_finn-b-F5_VASCDEM",full.names=T)
tmp_file_list <- c(tmp_file_list1,tmp_file_list2,tmp_file_list3)

for (tmp_file in tmp_file_list){
  tmp_dat <- fread(tmp_file) %>%
    mutate(pos_b37=position,
           effect_allele=ea,
           other_allele=nea,
           standard_error=se)
  write.table(tmp_dat,gsub(".tsv","_harmonized.tsv",tmp_file),sep="\t",row.names=F,quote=F)
}

### wightman_ad

tmp_file_list <- list.files("data/MR/colocalization","coloc_wightman_ad",full.names=T)

for (tmp_file in tmp_file_list){
  tmp_dat <- fread(tmp_file) %>%
    mutate(pos_b37=base_pair_location,
           p=p_value)
  write.table(tmp_dat,gsub(".tsv","_harmonized.tsv",tmp_file),sep="\t",row.names=F,quote=F)
}

##-----------------------------------------
## Format for locuscomparer
##-----------------------------------------

rm(list=ls())

instruments_selected_sig_anno <- read_excel("results/MR/instruments_selected_sig_final_anno.xlsx")

for (accession_id in unique(instruments_selected_sig_anno$accessionId)){
  
  tmp_dat1 <- instruments_selected_sig_anno %>% filter(accessionId==accession_id)
  
  for (out in gsub("6cWYtc","wightman_ad",unique(tmp_dat1$id.outcome))){
    
    tmp_dat2 <- tmp_dat1 %>% filter(id.outcome==gsub("wightman_ad","6cWYtc",out))
    
    for (rsid in unique(tmp_dat2$SNP)){
    
      tmp_exp <- fread(paste0("data/MR/colocalization/coloc_chen_",accession_id,"_",rsid,".tsv"), fill=T)
      tmp_out <- fread(paste0("data/MR/colocalization/coloc_",out,"_",rsid,"_harmonized.tsv"), fill=T)
      
      shared_rsid <- tmp_out[tmp_out$rsid %in% tmp_exp$variant_id,]$rsid
      tmp_exp_sub <- tmp_exp %>% filter(variant_id %in% shared_rsid)
      tmp_out_sub <- tmp_out %>% filter(rsid %in% shared_rsid)
        
      low_prob_exp <- ifelse(min(tmp_exp_sub$p_value)==0,0.01,0)
      low_prob_out <- ifelse(min(tmp_out_sub$p)==0,0.01,0)

      tmp_exp_sub <- tmp_exp_sub %>% 
        mutate(rsid = variant_id,
               pval = Winsorize(p_value,val=quantile(p_value,probs=c(low_prob_exp,1)))) %>%
        subset(select=c(rsid,pval))
      tmp_out_sub <- tmp_out_sub %>% 
        mutate(pval = Winsorize(p,val=quantile(p,probs=c(low_prob_out,1)))) %>%
        subset(select=c(rsid,pval))
      
      write.table(tmp_exp_sub, file=paste0("data/MR/colocalization/locuscomparer/locuscomparer_chen_",accession_id,"_",rsid,".tsv"), quote=FALSE, row.names = F, sep='\t')
      write.table(tmp_out_sub, file=paste0("data/MR/colocalization/locuscomparer/locuscomparer_",out,"_",rsid,".tsv"), quote=FALSE, row.names = F, sep='\t')
      
    }
  }
}

##-----------------------------------------
## locuscompare
##-----------------------------------------

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

# gdata::keep(locuscompare_plot_list, sure=T)
# save.image("data/MR/colocalization/locuscomparer_plots.RData")

for (i in names(locuscompare_plot_list)) {
  pdf(paste0("results/MR/figures/mr_locuscompare_",i,".pdf"), width=10,height=5)
  print(locuscompare_plot_list[[i]])
  dev.off()
}

##-----------------------------------------
## Format for coloc
##-----------------------------------------

### Format for coloc input

rm(list=ls())

instruments_selected_sig_anno <- read_excel("results/MR/instruments_selected_sig_final_anno.xlsx")

coloc_res_list <- list(); i <- 1

for (accession_id in unique(instruments_selected_sig_anno$accessionId)){
  
  tmp_dat1 <- instruments_selected_sig_anno %>% filter(accessionId==accession_id)
  
  for (out in gsub("6cWYtc","wightman_ad",unique(tmp_dat1$id.outcome))){
    
    tmp_dat2 <- tmp_dat1 %>% filter(id.outcome==gsub("wightman_ad","6cWYtc",out))
    
    for (rsid in unique(tmp_dat2$SNP)){
      
      tmp_exp <- fread(paste0("data/MR/colocalization/coloc_chen_",accession_id,"_",rsid,".tsv"), fill=T) %>% distinct()
      tmp_out <- fread(paste0("data/MR/colocalization/coloc_",out,"_",rsid,"_harmonized.tsv"), fill=T) %>% distinct()
      
      shared_rsid <- tmp_out[tmp_out$rsid %in% tmp_exp$variant_id,]$rsid
      tmp_exp_sub <- tmp_exp %>% filter(variant_id %in% shared_rsid) %>% arrange(variant_id)
      tmp_out_sub <- tmp_out %>% filter(rsid %in% shared_rsid) %>% arrange(rsid)
      
      tmp_sub <- merge(tmp_exp_sub, tmp_out_sub, by.x="variant_id", by.y="rsid", suffixes=c("_exp","_out"))
      tmp_sub <- tmp_sub %>%
        mutate(beta_out = ifelse(effect_allele_exp == effect_allele_out, beta_out, -1*beta_out)) %>%
        filter((effect_allele_out == effect_allele_exp) | (effect_allele_out == other_allele_exp))
      tmp_sub <- tmp_sub[!duplicated(tmp_sub$variant_id),]

      write.table(tmp_sub, paste0("data/MR/colocalization/formated_chen_",accession_id,"_",out,"_",rsid,".tsv"), quote=F, sep="\t", row.names=F)
     
    }
  }
}

##-----------------------------------------
## Run coloc
##-----------------------------------------

### coloc.abf

rm(list=ls())

instruments_selected_sig_anno <- read_excel("results/MR/instruments_selected_sig_final_anno.xlsx")

coloc_res_list <- list(); i <- 1

for (accession_id in unique(instruments_selected_sig_anno$accessionId)){
  
  tmp_dat1 <- instruments_selected_sig_anno %>% filter(accessionId==accession_id)
  
  for (out in gsub("6cWYtc","wightman_ad",unique(tmp_dat1$id.outcome))){
    
    tmp_dat2 <- tmp_dat1 %>% filter(id.outcome==gsub("wightman_ad","6cWYtc",out))
    
    for (rsid in unique(tmp_dat2$SNP)){
      
      tmp_sub <- fread(paste0("data/MR/colocalization/formated_chen_",accession_id,"_",out,"_",rsid,".tsv"))

      tmp_exp_list <- list()
      tmp_exp_list[["beta"]] <- tmp_sub$beta_exp
      tmp_exp_list[["varbeta"]] <- tmp_sub$standard_error_exp**2
      tmp_exp_list[["snp"]] <- tmp_sub$variant_id
      tmp_exp_list[["position"]] <- tmp_sub$pos_b37
      tmp_exp_list[["type"]] <- "quant"
      tmp_exp_list[["sdY"]] <- 1
      
      tmp_out_list <- list()
      tmp_out_list[["beta"]] <- tmp_sub$beta_out
      tmp_out_list[["varbeta"]] <- tmp_sub$standard_error_out**2
      tmp_out_list[["snp"]] <- tmp_sub$variant_id
      tmp_out_list[["position"]] <- tmp_sub$pos_b37
      if (out == "ebi-a-GCST006572"){
        tmp_out_list[["type"]] <- "quant"
        tmp_out_list[["sdY"]] <- 1
      } else{
        tmp_out_list[["type"]] <- "cc"
      }
      
     tmp_res <- coloc.abf(dataset1=tmp_exp_list,dataset2=tmp_out_list)
     coloc_res_list[[i]] <- c(unique(tmp_dat2$anno_outcome),unique(tmp_dat2$exposure_anno),rsid,
                              as.vector(tmp_res[["summary"]])); i <- i + 1
     
    }
  }
}

coloc_res_dat <- as.data.frame(do.call("rbind",coloc_res_list))
colnames(coloc_res_dat) <- c("outcome_anno","exposure_anno","rsid","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")

coloc_res_dat <- coloc_res_dat %>%
  mutate(PP.H0.abf = as.numeric(PP.H0.abf),
         PP.H1.abf = as.numeric(PP.H1.abf),
         PP.H2.abf = as.numeric(PP.H2.abf),
         PP.H3.abf = as.numeric(PP.H3.abf),
         PP.H4.abf = as.numeric(PP.H4.abf),
         PP3_plus_PP4 = PP.H3.abf + PP.H4.abf,
         PP4_div_PP3_plus_PP4 = PP.H4.abf / (PP3_plus_PP4),
         outcome_exposure_anno = paste0(outcome_anno,"_",exposure_anno))

coloc_res_hits <- coloc_res_dat %>% filter(PP.H4.abf>=0.7|PP4_div_PP3_plus_PP4>=0.7)
length(coloc_res_hits$outcome_exposure_anno); length(unique(coloc_res_hits$outcome_exposure_anno))
range(coloc_res_hits$PP.H4.abf)
range(coloc_res_hits$PP4_div_PP3_plus_PP4)

write_xlsx(coloc_res_dat, "results/MR/colocalization/coloc_res.xlsx")

### coloc.susie

rm(list=ls())

instruments_selected_sig_anno <- read_excel("results/MR/instruments_selected_sig_final_anno.xlsx")
load("~/Documents/coloc_susie_all_unique_rsid_ld_mat.RData")
source("scripts/source_outcome_sample_size.R")

sum(is.na(instruments_selected_sig_anno$samplesize.exposure))
sum(is.na(instruments_selected_sig_anno$samplesize.outcome))
instruments_selected_sig_anno[is.na(instruments_selected_sig_anno$samplesize.outcome),]$outcome
instruments_selected_sig_anno <- merge(instruments_selected_sig_anno, n_missing, by="id.outcome", suffixes=c("","_new"), all.x=T) %>%
  mutate(samplesize.outcome = ifelse(is.na(samplesize.outcome), samplesize.outcome_new, samplesize.outcome))
sum(is.na(instruments_selected_sig_anno$samplesize.outcome))

coloc_susie_dat_list <- list(); i <- 1

for (accession_id in unique(instruments_selected_sig_anno$accessionId)){
  
  tmp_dat1 <- instruments_selected_sig_anno %>% filter(accessionId==accession_id)
  
  for (out in gsub("6cWYtc","wightman_ad",unique(tmp_dat1$id.outcome))){
    
    tmp_dat2 <- tmp_dat1 %>% filter(id.outcome==gsub("wightman_ad","6cWYtc",out))
    
    for (rsid in unique(tmp_dat2$SNP)){
      
      print(i); i <- i + 1
      
      tmp_sub <- fread(paste0("data/MR/colocalization/formated_chen_",accession_id,"_",out,"_",rsid,".tsv"))
      tmp_ld <- ld_mat_list[[paste0(accession_id,"_",out,"_",rsid)]]
      colnames(tmp_ld) <- str_split_fixed(colnames(tmp_ld),"_",3)[,1]
      rownames(tmp_ld) <- str_split_fixed(rownames(tmp_ld),"_",3)[,1]
      tmp_sub <- tmp_sub %>% filter(variant_id %in% colnames(tmp_ld))
      tmp_sub <- tmp_sub[match(colnames(tmp_ld), tmp_sub$variant_id),]
      
      tmp_exp_list <- list()
      tmp_exp_list[["beta"]] <- tmp_sub$beta_exp
      tmp_exp_list[["varbeta"]] <- tmp_sub$standard_error_exp**2
      tmp_exp_list[["snp"]] <- tmp_sub$variant_id
      tmp_exp_list[["position"]] <- tmp_sub$pos_b37
      tmp_exp_list[["type"]] <- "quant"
      tmp_exp_list[["sdY"]] <- 1
      tmp_exp_list[["N"]] <- tmp_dat2[tmp_dat2$SNP==rsid,]$samplesize.exposure
      tmp_exp_list[["LD"]] <- as.matrix(tmp_ld)
      
      tmp_out_list <- list()
      tmp_out_list[["beta"]] <- tmp_sub$beta_out
      tmp_out_list[["varbeta"]] <- tmp_sub$standard_error_out**2
      tmp_out_list[["snp"]] <- tmp_sub$variant_id
      tmp_out_list[["position"]] <- tmp_sub$pos_b37
      if (out == "ebi-a-GCST006572"){
        tmp_out_list[["type"]] <- "quant"
        tmp_out_list[["sdY"]] <- 1
      } else{
        tmp_out_list[["type"]] <- "cc"
      }
      tmp_out_list[["N"]] <- tmp_dat2[tmp_dat2$SNP==rsid,]$samplesize.outcome
      tmp_out_list[["LD"]] <- as.matrix(tmp_ld)
      
      susie_exp <- tryCatch(runsusie(tmp_exp_list, repeat_until_convergence=F),error=function(e) e)
      susie_out <- tryCatch(runsusie(tmp_out_list, repeat_until_convergence=F),error=function(e) e)
      
      coloc_susie_dat_list[[paste0(accession_id,"_",out,"_",rsid)]] <- list("exp"=susie_exp,
                                                                            "out"=susie_out)
      
    }
  }
}

gdata::keep(coloc_susie_dat_list, sure=T)
save.image("~/Documents/coloc_susie_input.RData")

rm(list=ls())
load("~/Documents/coloc_susie_input.RData")

coloc_susie_res_list <- list(); i <- 1

for (unique_id in names(coloc_susie_dat_list)){
  
  tmp_exp_list <- coloc_susie_dat_list[[unique_id]][["exp"]]
  tmp_out_list <- coloc_susie_dat_list[[unique_id]][["out"]]
  
  if (length(tmp_exp_list) > 2 & length(tmp_out_list) > 2){
    
    tmp_res <- coloc.susie(tmp_exp_list,tmp_out_list)
    
    if (!is.data.frame(tmp_res)){
      
      coloc_susie_res_list[[i]] <- tmp_res[["summary"]] %>%
        mutate(unique_id = unique_id); i <- i + 1
      
    }
  }
}

coloc_susie_res_dat <- as.data.frame(do.call("rbind",coloc_susie_res_list))

coloc_susie_res_dat <- coloc_susie_res_dat %>%
  mutate(PP.H0.abf = as.numeric(PP.H0.abf),
         PP.H1.abf = as.numeric(PP.H1.abf),
         PP.H2.abf = as.numeric(PP.H2.abf),
         PP.H3.abf = as.numeric(PP.H3.abf),
         PP.H4.abf = as.numeric(PP.H4.abf),
         PP3_plus_PP4 = PP.H3.abf + PP.H4.abf,
         PP4_div_PP3_plus_PP4 = PP.H4.abf / (PP3_plus_PP4))

instruments_selected_sig_anno <- read_excel("results/MR/instruments_selected_sig_final_anno.xlsx") %>%
  mutate(unique_id = paste0(accessionId,"_",gsub("6cWYtc","wightman_ad",id.outcome),"_",SNP),
         outcome_exposure_anno = paste0(anno_outcome,"_",exposure_anno)) %>%
  subset(select=c(unique_id,outcome_exposure_anno))
coloc_susie_res_dat <- merge(coloc_susie_res_dat,instruments_selected_sig_anno,by="unique_id")

# write_xlsx(coloc_susie_res_dat, "results/MR/colocalization/coloc_susie_res.xlsx")

coloc_res_hits <- coloc_res_dat %>% filter(PP.H4.abf>=0.7|PP4_div_PP3_plus_PP4>=0.7)
coloc_susie_res_hits <- coloc_susie_res_dat %>% filter(PP.H4.abf>=0.7|PP4_div_PP3_plus_PP4>=0.7)

### Call-out locuscompare

rm(list=ls())

instruments_selected_sig_anno <- read_excel("results/MR/instruments_selected_sig_final_anno.xlsx") %>%
  mutate(outcome_exposure_rsid_anno = paste0(anno_outcome,"_",exposure_anno,"_",SNP))
coloc_res_dat <- read_excel("results/MR/colocalization/coloc_res.xlsx")
coloc_res_hits <- coloc_res_dat %>% 
  filter(PP.H4.abf>=0.7|PP4_div_PP3_plus_PP4>=0.7) %>% 
  mutate(outcome_exposure_rsid_anno = paste0(outcome_exposure_anno,"_",rsid))
coloc_res_hits <- merge(coloc_res_hits, instruments_selected_sig_anno, by="outcome_exposure_rsid_anno")

locuscompare_plot_list <- list()

for (i in 1:dim(coloc_res_hits)[1]){
  
  out <- gsub("6cWYtc","wightman_ad",coloc_res_hits$id.outcome[i])
  rsid <- coloc_res_hits$rsid[i]
  accession_id <- coloc_res_hits$accessionId[i]
      
  tmp_out <- paste0("data/MR/colocalization/locuscomparer/locuscomparer_",out,"_",rsid,".tsv")
  tmp_exp <- paste0("data/MR/colocalization/locuscomparer/locuscomparer_chen_",accession_id,"_",rsid,".tsv")
  
  tmp_out_dat <- fread(paste0("data/MR/colocalization/locuscomparer/locuscomparer_",out,"_",rsid,".tsv"))
  tmp_exp_dat <- fread(paste0("data/MR/colocalization/locuscomparer/locuscomparer_chen_",accession_id,"_",rsid,".tsv"))
  tmp_highlight <- ifelse(((rsid %in% tmp_exp_dat$rsid) & (rsid %in% tmp_out_dat$rsid)),list(rsid),list(NULL))
  
  locuscompare_plot_list[[paste0(accession_id,"_",out,"_",rsid)]] <- locuscompare(in_fn1=tmp_out, in_fn2=tmp_exp, snp=tmp_highlight[[1]],
                                                                                  title1=unique(coloc_res_hits$anno_outcome[i]), title2=coloc_res_hits$exposure_anno.x[i],
                                                                                  legend_position="topleft",lz_ylab_linebreak=T)

}

for (i in names(locuscompare_plot_list)) {
  pdf(paste0("results/MR/colocalization/hits_mr_locuscompare_",i,".pdf"), width=10,height=5)
  print(locuscompare_plot_list[[i]])
  dev.off()
}
