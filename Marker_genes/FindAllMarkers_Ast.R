setwd("/media/hansruedi/3b32f153-0651-4d56-b863-9b070ac7dd55/home/hansruedi/Brain_region_final")
library(Seurat)
pbmc=readRDS("Ast_integrated2.rds")

DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


addition=read.table("Ast_Metadata2_addition_cell_type_highres.txt",header=TRUE,row.names=1,sep="\t")
pbmc=AddMetaData(pbmc,addition)
Idents(pbmc)<-"cell_type_high_resolution"

p1=DimPlot(pbmc,cols="Paired",label=TRUE)+NoAxes()
p1

Idents(pbmc)<-"cell_type_high_resolution"
DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers,file="Marker_genes_Ast_all_cells_cell_type_high_resolution.csv")
