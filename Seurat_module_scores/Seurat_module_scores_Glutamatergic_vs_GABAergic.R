setwd("/media/hansruedi/3b32f153-0651-4d56-b863-9b070ac7dd55/home/hansruedi/Brain_region_final")
library(ggplot2)
library(Seurat)
pbmc <- readRDS(file = "In_raw_mscores_integrated2.rds")


addition=read.table("Inh_Metadata_addition_cell_type_highres.txt",header=TRUE,row.names=1,sep="\t")
pbmc=AddMetaData(pbmc,addition)
Idents(pbmc)<-"cell_type_high_resolution"
DefaultAssay(pbmc)<-"RNA"

pbmc_subset <- subset(pbmc, downsample = 2000)
setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB/Data_GABA_Glu_module_scores")
saveRDS(pbmc_subset,file="In_raw_mscores_integrated2_subset2000.rds")

#######################################################################

setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB")
library(Seurat)
pbmc=readRDS("Exc_raw_subset1.rds")

Idents(pbmc)<-"region"
DefaultAssay(pbmc)<-"RNA"

Idents(pbmc)<-"cell_type_high_resolution"


pbmc_subset <- subset(pbmc, downsample = 2000)

setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB/Data_GABA_Glu_module_scores")
saveRDS(pbmc_subset,file="Exc_raw_subset1_subset2000.rds")

#######################################################################

setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB")
library(Seurat)
pbmc=readRDS("Exc_raw_subset2.rds")

Idents(pbmc)<-"region"
DefaultAssay(pbmc)<-"RNA"

Idents(pbmc)<-"cell_type_high_resolution"


pbmc_subset <- subset(pbmc, downsample = 2000)

setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB/Data_GABA_Glu_module_scores")
saveRDS(pbmc_subset,file="Exc_raw_subset2_subset2000.rds")

#######################################################################

setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB/Data_GABA_Glu_module_scores")
library(ggplot2)
library(Seurat)
Exc_part1 <- readRDS(file = "Exc_raw_subset1_subset2000.rds")
Exc_part2 <- readRDS(file = "Exc_raw_subset2_subset2000.rds")
Inh <- readRDS(file = "In_raw_mscores_integrated2_subset2000.rds")
library(Seurat)
Exc_combined <- merge(x = Exc_part1, y = Exc_part2)

Neurons <- merge(x = Exc_combined, y = Inh)

DefaultAssay(Neurons)<-"RNA"
Neurons <- NormalizeData(Neurons, normalization.method = "LogNormalize", scale.factor = 10000)

setwd("/media/hansruedi/3b32f153-0651-4d56-b863-9b070ac7dd55/home/hansruedi/Allen_Institute_06142020")

Marker_genes_full=read.csv("Marker_genes_class_label_only_neurons.csv",header=TRUE)
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

setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB/Data_GABA_Glu_module_scores")

write.csv(Neurons@meta.data,file="Neurons_metadata.csv")
