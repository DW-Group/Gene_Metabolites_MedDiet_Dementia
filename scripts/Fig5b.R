library(data.table)
library(readxl)
library(ggplot2)
library(circlize)
library(dplyr)
library(ComplexHeatmap)

##-----------------------------------------
## Fig. 5b: Circle plots
##-----------------------------------------

rm(list=ls())

### Data

exposure_dict <- read_excel("results/instruments_selected_sig_final_anno_08062024.xlsx") %>%
  subset(select=c(exposure,exposure_anno)) %>% distinct()
mr_res <- read_excel("results/mr_results_final_08062024.xlsx")

overlaps <- mr_res %>%
  filter(pval<0.05) %>%
  group_by(exposure) %>%
  summarise(n_distinct(outcome))
table(overlaps$`n_distinct(outcome)`)

### Chord diagram

mr_res_sig <- mr_res %>%
  filter(FDR_BH_sig)

chord_dat <- as.data.frame(table(mr_res_sig$exposure,mr_res_sig$anno_outcome)) %>%
  reshape(idvar = "Var1", timevar = "Var2", direction = "wide") %>%
  column_to_rownames("Var1") %>%
  as.matrix()
colnames(chord_dat) <- gsub("Freq.","",colnames(chord_dat))

group <- structure(c(mr_res_sig$classification, "Dementia and subtypes", "Dementia and subtypes", "Dementia and subtypes", "Cognitive function"),
                   names = c(mr_res_sig$exposure, "Dementia", "Alzheimer's disease", "Vascular dementia", "Cognitive performance"))
unique(group)

classcol <- c("#8F7700FF","#20854EFF","#CD534CFF","#EFC000FF","#0072B5FF","#E18727FF","#003C67FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#8F7700FF","#868686FF")
names(classcol) <- c("Dementia and subtypes", "Cognitive function","Metabolite ratio","Amino Acid","Lipid","Cofactors and Vitamins","Nucleotide","Peptide","Xenobiotics","Partially Characterized Molecules","Carbohydrate","Energy","Unknown")

grid.col = structure(c(rep("gray", dim(chord_dat)[1]), c("#C06907B3","#07854780","#9A092F80","#EF8D46B3")),
                     names = c(rownames(chord_dat), colnames(chord_dat)))

pdf("results/figures/fig5b_chord_diagram.pdf", width = 15, height = 15, onefile = F)
circos.par(start.degree = 90,
           track.margin = c(0.01, 0.01),
           track.height = 0.01,
           canvas.xlim = c(-2.2, 1.2), canvas.ylim = c(-1.2, 1))
chordDiagram(t(chord_dat), annotationTrack = "grid", 
             group = group, grid.col=grid.col,  
             preAllocateTracks = 0.5,
             annotationTrackHeight = 0.05,
             big.gap = 5, small.gap = 1)
for(i in unique(group)){
  if(i %in% c("Dementia and subtypes","Cognitive function")){
    tmp_names <- unique(mr_res_sig[mr_res_sig$cat_outcome==i,]$anno_outcome)
  } else{
    tmp_names <- mr_res_sig[mr_res_sig$classification==i,]$exposure
  }
  highlight.sector(tmp_names, 
                   track.index = 1, col = classcol[[i]], niceFacing = TRUE)
}
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + 1.2, sector.name, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.8)
}, bg.border = NA)
circos.clear()
dev.off()

selected_classcol <- classcol[names(classcol) %in% mr_res_sig$classification]
pdf("results/figures/fig5b_chord_diagram.pdf_legend.pdf", width = 5,height = 5)
lgd_super_class <- Legend(title = "Superclass of metabolites",
                          title_gp = gpar(fontsize = 10, fontface = "bold"),
                          at = names(selected_classcol),
                          legend_gp = gpar(fill = selected_classcol))
draw(lgd_super_class)
dev.off()
