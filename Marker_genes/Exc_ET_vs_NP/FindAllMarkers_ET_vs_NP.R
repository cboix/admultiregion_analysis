setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB")

library(ggplot2)
library(Seurat)
pbmc <- readRDS(file = "Exc_raw_NP_ET.rds")


#Marker genes
DefaultAssay(pbmc)<-"RNA"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(pbmc)<-"cell_type_high_resolution"

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

setwd("/media/hansruedi/12TB_HD/Brain_region_final_12TB")
write.csv(pbmc.markers,file="Marker_genes_Exc_NP_ET.csv")
