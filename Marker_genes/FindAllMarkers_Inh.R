setwd("/media/hansruedi/3b32f153-0651-4d56-b863-9b070ac7dd55/home/hansruedi/Brain_region_final")
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

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,max.cells.per.ident = 1000)
write.csv(pbmc.markers,file="Marker_genes_Inh_cell_type_high_resolution_arranged")

