setwd("D:/Brain_region_project/Brain_region_final")
library(ggplot2)
library(Seurat)
pbmc=readRDS("Ast_integrated2.rds")

DefaultAssay(pbmc)<-"RNA"

addition=read.table("Ast_Metadata2_addition_cell_type_highres.txt",header=TRUE,row.names=1,sep="\t")
pbmc=AddMetaData(pbmc,addition)
Idents(pbmc)<-"cell_type_high_resolution"

p1=DimPlot(pbmc,cols="Paired",label=TRUE)+NoAxes()
p1

#DropViz data
setwd("D:/DropViz_data")
DropViz=readRDS("FC_PC_HC_TH_merged_ASTROCYTE_human_genes.rds")

#Linnarsson data
setwd("D:/Mousebrain_data")
Linnarsson=readRDS("Astrocytes_selected_regions_human_genes.rds")

#Integration
ifnb.list=list(Mathys=pbmc,Linnarsson=Linnarsson,DropViz=DropViz)

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

rm(pbmc)
rm(DropViz)
rm(Linnarsson)


# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)


# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.1)

setwd("D:/Mouse_data_integration")
saveRDS(immune.combined,file="Astrocytes_integrated_all.rds")