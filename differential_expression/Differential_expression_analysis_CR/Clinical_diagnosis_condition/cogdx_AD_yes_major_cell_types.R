library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)
library(Seurat)

#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")

pbmc=readRDS("Oli_integrated_batch3000_leiden_res02_final.rds")

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"
pbmc[["integrated"]] <- NULL
pbmc[["SCT"]] <- NULL
pbmc[["prediction.score.celltype"]] <- NULL

Idents(pbmc)<-'AD'
pbmc=subset(pbmc,idents=c('yes'))

Idents(pbmc)<-'cogdx'
pbmc=subset(pbmc,idents=c('1','4'))
pbmc.sce <- as.SingleCellExperiment(pbmc)

# remove undetected genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
dim(pbmc.sce)

# remove lowly expressed genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
dim(pbmc.sce)

#sce@colData@listData[["AD"]]=make.names(sce@colData@listData[["AD"]])
#sce@colData@listData[["msex"]]=make.names(sce@colData@listData[["msex"]])

(sce <- prepSCE(pbmc.sce, 
                kid = "cell_type_high_resolution", # subpopulation assignments
                gid = "cogdx",  # group IDs (ctrl/stim)
                sid = "projid",   # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # dr

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids
#t(table(sce$cluster_id, sce$sample_id))

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
assayNames(pb)
#t(head(assay(pb)))

formula = ~group_id+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("group_id","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

rm(pbmc)
rm(sce)
rm(pbmc.sce)
save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cogdx_AD_yes.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogdx_AD_yes")

#Indicate name of the variable analyzed
variable="cogdx"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")

pbmc=readRDS("Opc_integrated_batch_5000_leiden_predicted_final.rds")

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"
pbmc[["integrated"]] <- NULL
pbmc[["SCT"]] <- NULL
pbmc[["prediction.score.celltype"]] <- NULL

Idents(pbmc)<-'AD'
pbmc=subset(pbmc,idents=c('yes'))

Idents(pbmc)<-'cogdx'
pbmc=subset(pbmc,idents=c('1','4'))
pbmc.sce <- as.SingleCellExperiment(pbmc)

# remove undetected genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
dim(pbmc.sce)

# remove lowly expressed genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
dim(pbmc.sce)

#sce@colData@listData[["AD"]]=make.names(sce@colData@listData[["AD"]])
#sce@colData@listData[["msex"]]=make.names(sce@colData@listData[["msex"]])

(sce <- prepSCE(pbmc.sce, 
                kid = "cell_type_high_resolution", # subpopulation assignments
                gid = "cogdx",  # group IDs (ctrl/stim)
                sid = "projid",   # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # dr

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids
#t(table(sce$cluster_id, sce$sample_id))

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
assayNames(pb)
#t(head(assay(pb)))

formula = ~group_id+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("group_id","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

rm(pbmc)
rm(sce)
rm(pbmc.sce)
save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/Opc_cogdx_AD_yes.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogdx_AD_yes")

#Indicate name of the variable analyzed
variable="cogdx"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
#######################################################################

#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")

pbmc=readRDS("Ast_integrated_batch_5000_leiden_predicted_final.rds")

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"
pbmc[["integrated"]] <- NULL
pbmc[["SCT"]] <- NULL
pbmc[["prediction.score.celltype"]] <- NULL

Idents(pbmc)<-'AD'
pbmc=subset(pbmc,idents=c('yes'))

Idents(pbmc)<-'cogdx'
pbmc=subset(pbmc,idents=c('1','4'))

HM=pbmc@meta.data
HM2=HM[,c("projid","cell_type_high_resolution")]

HM3 = HM2 %>% mutate(major_cell_type = case_when(grepl("Exc ",cell_type_high_resolution)~"Exc",
                                                 grepl("CA1 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("CA2, CA3 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("DG granule cells",cell_type_high_resolution)~"Exc",
                                                 grepl("Inh ",cell_type_high_resolution)~"Inh",
                                                 grepl("Ast ",cell_type_high_resolution)~"Ast",
                                                 grepl("Oli ",cell_type_high_resolution)~"Oli",
                                                 grepl("OPC ",cell_type_high_resolution)~"OPC",
                                                 grepl("Mic ",cell_type_high_resolution)~"Mic",
                                                 grepl("CAMs",cell_type_high_resolution)~"CAMs",
                                                 grepl("CPEC",cell_type_high_resolution)~"CPEC",
                                                 grepl("End",cell_type_high_resolution)~"End",
                                                 grepl("Epd",cell_type_high_resolution)~"Epd",
                                                 grepl("Fib",cell_type_high_resolution)~"Fib",
                                                 grepl("Per",cell_type_high_resolution)~"Per",
                                                 grepl("SMC",cell_type_high_resolution)~"SMC",
                                                 grepl("T cells",cell_type_high_resolution)~"T cells"))
HM4=HM3[,-1]
HM5=HM4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,HM5)
Idents(pbmc) <- 'major_cell_type'
pbmc=subset(pbmc, idents="Ast")

pbmc.sce <- as.SingleCellExperiment(pbmc)

# remove undetected genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
dim(pbmc.sce)

# remove lowly expressed genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
dim(pbmc.sce)

#sce@colData@listData[["AD"]]=make.names(sce@colData@listData[["AD"]])
#sce@colData@listData[["msex"]]=make.names(sce@colData@listData[["msex"]])

(sce <- prepSCE(pbmc.sce, 
                kid = "major_cell_type", # subpopulation assignments
                gid = "cogdx",  # group IDs (ctrl/stim)
                sid = "projid",   # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # dr

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids
#t(table(sce$cluster_id, sce$sample_id))

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
assayNames(pb)
#t(head(assay(pb)))

formula = ~group_id+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("group_id","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

rm(pbmc)
rm(sce)
rm(pbmc.sce)
save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_major_cell_type_cogdx_AD_yes.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogdx_AD_yes")

#Indicate name of the variable analyzed
variable="cogdx"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")

pbmc=readRDS("Immune_integrated_batch_5000_predicted_final.rds")

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"
pbmc[["integrated"]] <- NULL
pbmc[["SCT"]] <- NULL
pbmc[["prediction.score.celltype"]] <- NULL

Idents(pbmc)<-'AD'
pbmc=subset(pbmc,idents=c('yes'))

Idents(pbmc)<-'cogdx'
pbmc=subset(pbmc,idents=c('1','4'))

HM=pbmc@meta.data
HM2=HM[,c("projid","cell_type_high_resolution")]

HM3 = HM2 %>% mutate(major_cell_type = case_when(grepl("Exc ",cell_type_high_resolution)~"Exc",
                                                 grepl("CA1 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("CA2, CA3 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("DG granule cells",cell_type_high_resolution)~"Exc",
                                                 grepl("Inh ",cell_type_high_resolution)~"Inh",
                                                 grepl("Ast ",cell_type_high_resolution)~"Ast",
                                                 grepl("Oli ",cell_type_high_resolution)~"Oli",
                                                 grepl("OPC ",cell_type_high_resolution)~"OPC",
                                                 grepl("Mic ",cell_type_high_resolution)~"Mic",
                                                 grepl("CAMs",cell_type_high_resolution)~"CAMs",
                                                 grepl("CPEC",cell_type_high_resolution)~"CPEC",
                                                 grepl("End",cell_type_high_resolution)~"End",
                                                 grepl("Epd",cell_type_high_resolution)~"Epd",
                                                 grepl("Fib",cell_type_high_resolution)~"Fib",
                                                 grepl("Per",cell_type_high_resolution)~"Per",
                                                 grepl("SMC",cell_type_high_resolution)~"SMC",
                                                 grepl("T cells",cell_type_high_resolution)~"T cells"))
HM4=HM3[,-1]
HM5=HM4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,HM5)
Idents(pbmc) <- 'major_cell_type'
pbmc=subset(pbmc, idents="Mic")

pbmc.sce <- as.SingleCellExperiment(pbmc)

# remove undetected genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
dim(pbmc.sce)

# remove lowly expressed genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
dim(pbmc.sce)

#sce@colData@listData[["AD"]]=make.names(sce@colData@listData[["AD"]])
#sce@colData@listData[["msex"]]=make.names(sce@colData@listData[["msex"]])

(sce <- prepSCE(pbmc.sce, 
                kid = "major_cell_type", # subpopulation assignments
                gid = "cogdx",  # group IDs (ctrl/stim)
                sid = "projid",   # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # dr

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids
#t(table(sce$cluster_id, sce$sample_id))

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
assayNames(pb)
#t(head(assay(pb)))

formula = ~group_id+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("group_id","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

rm(pbmc)
rm(sce)
rm(pbmc.sce)
save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_major_cell_type_cogdx_AD_yes.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogdx_AD_yes")

#Indicate name of the variable analyzed
variable="cogdx"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_merge_40_percent_final.rds")

Idents(pbmc)<-"cell_type_high_resolution"
pbmc<-subset(pbmc,idents=c("Exc L2-3 CBLN2 LINC02306", "Exc L3-4 RORB CUX2", "Exc L3-5 RORB PLCH1", "Exc L4-5 RORB GABRG1",    
                           "Exc L4-5 RORB IL1RAPL2", "Exc L5 ET", "Exc L5-6 RORB LINC02196", "Exc L5/6 IT Car3",       
                           "Exc L5/6 NP", "Exc L6 CT", "Exc L6 THEMIS NFIA", "Exc L6b" ))
levels(pbmc)

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"

Idents(pbmc)<-'AD'
pbmc=subset(pbmc,idents=c('yes'))

Idents(pbmc)<-'cogdx'
pbmc=subset(pbmc,idents=c('1','4'))

#add MCT info
DefaultAssay(pbmc)="RNA"

Idents(pbmc)<-"cell_type_high_resolution"

HM=pbmc@meta.data
HM2=HM[,c("projid","cell_type_high_resolution")]

HM3 = HM2 %>% mutate(major_cell_type = case_when(grepl("Exc ",cell_type_high_resolution)~"Exc",
                                                 grepl("CA1 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("CA2, CA3 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("DG granule cells",cell_type_high_resolution)~"Exc",
                                                 grepl("Inh",cell_type_high_resolution)~"Inh",
                                                 grepl("Ast",cell_type_high_resolution)~"Ast",
                                                 grepl("Oli",cell_type_high_resolution)~"Oli",
                                                 grepl("OPC",cell_type_high_resolution)~"OPC",
                                                 grepl("Mic",cell_type_high_resolution)~"Mic",
                                                 grepl("CAMs",cell_type_high_resolution)~"CAMs",
                                                 grepl("CPEC",cell_type_high_resolution)~"CPEC",
                                                 grepl("End",cell_type_high_resolution)~"End",
                                                 grepl("Epd",cell_type_high_resolution)~"Epd",
                                                 grepl("Fib",cell_type_high_resolution)~"Fib",
                                                 grepl("Per",cell_type_high_resolution)~"Per",
                                                 grepl("SMC",cell_type_high_resolution)~"SMC",
                                                 grepl("T cells",cell_type_high_resolution)~"T cells"))
HM4=HM3[,-1]
HM5=HM4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,HM5)
Idents(pbmc) <- 'major_cell_type'
pbmc=subset(pbmc, idents="Exc")
######################################################

pbmc.sce <- as.SingleCellExperiment(pbmc)

# remove undetected genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
dim(pbmc.sce)

# remove lowly expressed genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
dim(pbmc.sce)

#sce@colData@listData[["AD"]]=make.names(sce@colData@listData[["AD"]])
#sce@colData@listData[["msex"]]=make.names(sce@colData@listData[["msex"]])

(sce <- prepSCE(pbmc.sce, 
                kid = "major_cell_type", # subpopulation assignments
                gid = "cogdx",  # group IDs (ctrl/stim)
                sid = "projid",   # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # dr

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids
#t(table(sce$cluster_id, sce$sample_id))

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
assayNames(pb)
#t(head(assay(pb)))

formula = ~group_id+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("group_id","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

rm(pbmc)
rm(sce)
rm(pbmc.sce)

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogdx_AD_yes")

#Indicate name of the variable analyzed
variable="cogdx"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
#######################################################################
#######################################################################

#######################################################################
# MCT_Inh

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("In_integrated_batch_3000_module_scores_predicted_final.rds")
DefaultAssay(pbmc)="RNA"
pbmc[["integrated"]] <- NULL
pbmc[["SCT"]] <- NULL
pbmc[["prediction.score.celltype"]] <- NULL

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"

Idents(pbmc)<-'AD'
pbmc=subset(pbmc,idents=c('yes'))

Idents(pbmc)<-'cogdx'
pbmc=subset(pbmc,idents=c('1','4'))

#add MCT info
DefaultAssay(pbmc)="RNA"

Idents(pbmc)<-"cell_type_high_resolution"

HM=pbmc@meta.data
HM2=HM[,c("projid","cell_type_high_resolution")]

HM3 = HM2 %>% mutate(major_cell_type = case_when(grepl("Exc ",cell_type_high_resolution)~"Exc",
                                                 grepl("CA1 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("CA2, CA3 pyramidal cells",cell_type_high_resolution)~"Exc",
                                                 grepl("DG granule cells",cell_type_high_resolution)~"Exc",
                                                 grepl("Inh",cell_type_high_resolution)~"Inh",
                                                 grepl("Ast",cell_type_high_resolution)~"Ast",
                                                 grepl("Oli",cell_type_high_resolution)~"Oli",
                                                 grepl("OPC",cell_type_high_resolution)~"OPC",
                                                 grepl("Mic",cell_type_high_resolution)~"Mic",
                                                 grepl("CAMs",cell_type_high_resolution)~"CAMs",
                                                 grepl("CPEC",cell_type_high_resolution)~"CPEC",
                                                 grepl("End",cell_type_high_resolution)~"End",
                                                 grepl("Epd",cell_type_high_resolution)~"Epd",
                                                 grepl("Fib",cell_type_high_resolution)~"Fib",
                                                 grepl("Per",cell_type_high_resolution)~"Per",
                                                 grepl("SMC",cell_type_high_resolution)~"SMC",
                                                 grepl("T cells",cell_type_high_resolution)~"T cells"))
HM4=HM3[,-1]
HM5=HM4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,HM5)
Idents(pbmc) <- 'major_cell_type'
pbmc=subset(pbmc, idents="Inh")
######################################################

pbmc.sce <- as.SingleCellExperiment(pbmc)

# remove undetected genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
dim(pbmc.sce)

# remove lowly expressed genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
dim(pbmc.sce)

#sce@colData@listData[["AD"]]=make.names(sce@colData@listData[["AD"]])
#sce@colData@listData[["msex"]]=make.names(sce@colData@listData[["msex"]])

(sce <- prepSCE(pbmc.sce, 
                kid = "major_cell_type", # subpopulation assignments
                gid = "cogdx",  # group IDs (ctrl/stim)
                sid = "projid",   # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # dr

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids
#t(table(sce$cluster_id, sce$sample_id))

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
assayNames(pb)
#t(head(assay(pb)))

formula = ~group_id+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("group_id","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

rm(pbmc)
rm(sce)
rm(pbmc.sce)

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogdx_AD_yes")

#Indicate name of the variable analyzed
variable="cogdx"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
#######################################################################