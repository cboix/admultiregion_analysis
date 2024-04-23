setwd("D:/Brain_region_project/Brain_region_final")
library(ggplot2)
library(Seurat)
pbmc <- readRDS(file = "In_raw_mscores_integrated2.rds")


addition=read.table("Inh_Metadata_addition_cell_type_highres.txt",header=TRUE,row.names=1,sep="\t")
pbmc=AddMetaData(pbmc,addition)
Idents(pbmc)<-"cell_type_high_resolution"
DefaultAssay(pbmc)<-"RNA"

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 23
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(nb.cols)

levels(pbmc) <- c("Inh L1-6 LAMP5 CA13","Inh LAMP5 NRG1 (Rosehip)","Inh LAMP5 RELN","Inh PAX6 RELN","Inh VIP ABI3BP","Inh SORCS1 TTN","Inh VIP TSHZ2","Inh RYR3 TSHZ2","Inh VIP THSD7B","Inh SGCD PDE3A","Inh ALCAM TRPM3","Inh PTPRK FAM19A1","Inh VIP CLSTN2","Inh GPC5 RIT2","Inh PVALB SULF1","Inh PVALB HTR4","Inh CUX2 MSR1","Inh PVALB CA8 (Chandelier)","Inh FBN2 EPB41L4A","Inh ENOX2 SPHKAP","Inh L3-5 SST MAFB","Inh L6 SST NPY","Inh MEIS2 FOXP2")
p1=DimPlot(pbmc,cols=mycolors,label=TRUE)+NoAxes()
p1


#human LGN data
setwd("D:/Brain_region_project/Brain_region_final_12TB/AllenBrain_human_LGN/human_LGN_gene_expression_matrices_2018-06-14")
human_LGN=readRDS(file="human_LGN_2.rds")
Idents(human_LGN)="class"
human_LGN=subset(human_LGN,idents="GABAergic")

#mouse cortex
setwd("D:/Brain_region_project/Brain_region_final_12TB/AllenBrain_MouseCortex")
Allen=readRDS("Mouse_cortex_10x_GABAergic_human_genes.rds")
Allen2=subset(Allen,downsample=50000)

#DropViz data
setwd("D:/DropViz_data")
DropViz=readRDS("FC_PC_HC_TH_merged_Neurons_human_genes_class.rds")
Idents(DropViz)="class"
DropViz=subset(DropViz,idents="GABAergic")

#Linnarsson data
setwd("D:/Mousebrain_data")
Linnarsson=readRDS("Neurons_selected_regions_human_genes_class.rds")
Idents(Linnarsson)="class"
Linnarsson=subset(Linnarsson,idents="GABAergic")


#Integration
ifnb.list=list(Mathys=pbmc,Linnarsson=Linnarsson,Allen=Allen2,DropViz=DropViz,human_LGN=human_LGN)

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

rm(pbmc)
rm(Allen)
rm(DropViz)
rm(Linnarsson)
rm(Allen2)
rm(DropViz)


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
saveRDS(immune.combined,file="GABAergic_neurons_integrated_all_class.rds")