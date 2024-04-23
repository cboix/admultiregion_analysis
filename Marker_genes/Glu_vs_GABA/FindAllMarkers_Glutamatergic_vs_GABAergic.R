setwd("/media/hansruedi/3b32f153-0651-4d56-b863-9b070ac7dd55/home/hansruedi/Allen_Institute_06142020")
data=read.csv("matrix.csv",header=TRUE,row.names=1)
rows=row.names(data)
columns=colnames(data)
metadata=read.csv("metadata.csv")
rownames(metadata)<-metadata$sample_name
library(Seurat)
counts=t(data)
pbmc <- CreateSeuratObject(counts = counts, project = "Allen", min.cells = 0, min.features = 0)
pbmc <- AddMetaData(pbmc, metadata = metadata)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(pbmc)<-"outlier_call"
pbmc_filtered=subset(pbmc, idents = c("True"), invert = TRUE)
Idents(pbmc_filtered)<-"class_label"

pbmc_sub=subset(pbmc_filtered,idents=c("GABAergic","Glutamatergic"))

pbmc.markers <- FindAllMarkers(pbmc_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(pbmc.markers,file="Marker_genes_class_label_only_neurons.csv")
