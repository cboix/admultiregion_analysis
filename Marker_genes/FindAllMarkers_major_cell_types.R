setwd("/media/hansruedi/3b32f153-0651-4d56-b863-9b070ac7dd55/home/hansruedi/Brain_region_final")
library(Seurat)
adata <- ReadH5AD(file = "Brain_region_raw_highly_variable2.h5ad")
pbmc <- CreateSeuratObject(counts = GetAssayData(object = adata, slot = 'counts'), project = "pbmc3k", min.cells = 0, min.features = 0)


addition=read.csv("Brain_region_raw_highly_variable_seurat_metadata_addition_scanpy")
rownames(addition)<-addition$index
pbmc=AddMetaData(pbmc,addition)

Idents(pbmc) <- "region_new"
pbmc=subset(pbmc,idents=c("AG","TH","EC","HC","MT","PFC"))

Idents(pbmc) <- "major_cell_type"
mds=FetchData(pbmc, vars=c("U1_new","U2_new"))

names(mds)[names(mds) == "U1_new"] <- "UMAP_1"
names(mds)[names(mds) == "U2_new"] <- "UMAP_2"
mds=as.matrix(mds)
pbmc[["umap"]] <- CreateDimReducObject(embeddings = mds, key = "UMAP_", assay = DefaultAssay(pbmc))
DimPlot(pbmc,reduction="umap")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,max.cells.per.ident = 5000)

write.csv(pbmc.markers,file="Marker_genes_major_cell_types.csv")