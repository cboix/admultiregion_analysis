setwd("D:/Brain_region_project/Brain_region_final")
library(ggplot2)
library(Seurat)
pbmc <- readRDS(file = "Ex_raw_highly_variable_mscores_integrated2_adjusted_umap.rds")

Idents(pbmc)<-"cell_type_high_resolution"
DefaultAssay(pbmc)<-"RNA"

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 32
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(nb.cols)

p1=DimPlot(pbmc,cols=mycolors,label=TRUE)+NoAxes()+NoLegend()
p1

pbmc_small=subset(pbmc,downsample=5000)

#human LGN data
setwd("D:/Brain_region_project/Brain_region_final_12TB/AllenBrain_human_LGN/human_LGN_gene_expression_matrices_2018-06-14")
human_LGN=readRDS(file="human_LGN_2.rds")
Idents(human_LGN)="class"
human_LGN=subset(human_LGN,idents="Glutamatergic")

#mouse cortex
setwd("D:/Brain_region_project/Brain_region_final_12TB/AllenBrain_MouseCortex")
Allen=readRDS("Mouse_cortex_10x_Glutamatergic_200K_human_genes.rds")
Allen2=subset(Allen,downsample=50000)


#DropViz data
setwd("D:/DropViz_data")
DropViz=readRDS("FC_PC_HC_TH_merged_Neurons_human_genes_class.rds")
Idents(DropViz)="class"
DropViz=subset(DropViz,idents="Glutamatergic")
Idents(DropViz)<-"DropViz"
DropViz2=subset(DropViz,downsample=50000)

#Linnarsson data
setwd("D:/Mousebrain_data")
Linnarsson=readRDS("Neurons_selected_regions_human_genes_class.rds")
Idents(Linnarsson)="class"
Linnarsson=subset(Linnarsson,idents="Glutamatergic")

#Integration
ifnb.list=list(Mathys=pbmc_small,Linnarsson=Linnarsson,Allen=Allen2,DropViz=DropViz2,human_LGN=human_LGN)

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

rm(pbmc)
rm(pbmc_small)
rm(Allen)
rm(DropViz)
rm(Linnarsson)
rm(Allen2)
rm(DropViz2)

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
saveRDS(immune.combined,file="Glutamatergic_neurons_integrated_all_class_new.rds")