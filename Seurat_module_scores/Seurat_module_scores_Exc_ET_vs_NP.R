setwd("/media/hansruedi/3b32f153-0651-4d56-b863-9b070ac7dd55/home/hansruedi/Brain_region_final")
library(ggplot2)
library(Seurat)
Exc_part1 <- readRDS(file = "Ex_raw_set1_2_subset2000.rds")
Exc_part2 <- readRDS(file = "Ex_raw_set2_2_subset2000.rds")
Inh <- readRDS(file = "In_raw_mscores_integrated2_subset2000.rds")
library(Seurat)
Exc_combined <- merge(x = Exc_part1, y = Exc_part2)

Neurons <- merge(x = Exc_combined, y = Inh)

DefaultAssay(Neurons)<-"RNA"
Neurons <- NormalizeData(Neurons, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB")
Marker_genes_full=read.csv("Marker_genes_Exc_NP_ET.csv",header=TRUE)
Marker_genes_full_significant=subset(Marker_genes_full,p_val_adj < 0.05)
Marker_genes=Marker_genes_full_significant[,7:8]
Marker_genes[] <- lapply(Marker_genes, as.character)
Marker_genes_split=split(Marker_genes, Marker_genes$cluster)
Marker_genes_list <- list()              
for (i in 1:length(Marker_genes_split)) {
  name=names(Marker_genes_split)[i]
  Marker_genes_list[[i]]<-Marker_genes_split[[i]]$gene
  names(Marker_genes_list)[i]<-name
}
all_genes=rownames(Neurons)
Marker_genes_list_pruned <- list()              
for (i in 1:length(Marker_genes_list)) {
  name=names(Marker_genes_list)[i]
  Marker_genes_list_pruned[[i]]<-intersect(Marker_genes_list[[i]],all_genes)
  names(Marker_genes_list_pruned)[i]<-name
}
for (i in 1:length(Marker_genes_list_pruned)) {
  name=names(Marker_genes_list_pruned)[i]
  genes.for.scoring <- list(c(Marker_genes_list_pruned[[i]]))
  Neurons <- AddModuleScore(object = Neurons, features = genes.for.scoring, name = name)
}

setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB")
write.csv(Neurons@meta.data,file="Neurons_module_scores_NP_ET.csv")
