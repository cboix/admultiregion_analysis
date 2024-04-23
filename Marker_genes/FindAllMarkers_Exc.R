setwd("/media/hansruedi/3b32f153-0651-4d56-b863-9b070ac7dd55/home/hansruedi/Brain_region_final")
library(ggplot2)
library(Seurat)
pbmc <- readRDS(file = "Ex_raw_highly_variable_mscores_integrated2_adjusted_umap.rds")

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 32
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(nb.cols)

p1=DimPlot(pbmc,cols=mycolors,label=TRUE)+NoAxes()
p1

Idents(pbmc)<-"cell_type_high_resolution"
DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,max.cells.per.ident = 1000)

write.csv(pbmc.markers,file="Marker_genes_Exc_cell_type_high_resolution_new.csv")
