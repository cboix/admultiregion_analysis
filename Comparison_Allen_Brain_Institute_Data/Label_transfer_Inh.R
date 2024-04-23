library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)


setwd("/media/hansruedi/3b32f153-0651-4d56-b863-9b070ac7dd55/home/hansruedi/Allen_Institute_06142020")
library(Seurat)
pbmc=readRDS("SeuratObject_AllenData.rds")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(pbmc)<-"outlier_call"
pbmc_filtered=subset(pbmc, idents = c("True"), invert = TRUE)
Idents(pbmc_filtered)<-"class_label"
GABA = subset(pbmc_filtered, idents = c("GABAergic"))
Idents(GABA)<-"cell_type_alias_label"

reference=GABA
reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 5000)
reference <- ScaleData(reference)
reference <- RunPCA(reference)


setwd("/media/hansruedi/3b32f153-0651-4d56-b863-9b070ac7dd55/home/hansruedi/Brain_region_final")
pbmc <- readRDS(file = "In_raw_mscores_integrated2.rds")
addition=read.table("Inh_Metadata_addition_cell_type_highres.txt",header=TRUE,row.names=1,sep="\t")
pbmc=AddMetaData(pbmc,addition)
Idents(pbmc)<-"cell_type_high_resolution"
DefaultAssay(pbmc)<-"RNA"

anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

pbmc <- TransferData(
  anchorset = anchors,
  reference = reference,
  query = pbmc,
  refdata = list(
    celltype = "cell_type_alias_label")
)

setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB/Label_transfer")

saveRDS(pbmc,file="In_integrated_batch_3000_module_scores_Allen_predicted.rds")
