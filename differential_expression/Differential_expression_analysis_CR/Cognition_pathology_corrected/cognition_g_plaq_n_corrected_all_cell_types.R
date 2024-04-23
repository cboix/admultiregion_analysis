


library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)


##########################################################################################################
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
                gid = "AD",  # group IDs (ctrl/stim)
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

rm(pbmc)
rm(sce)
rm(pbmc.sce)

##############################################
formula = ~cogn_global_lv+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Ast_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Ast_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_raw_set1_final.rds")

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"
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
                gid = "AD",  # group IDs (ctrl/stim)
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

rm(pbmc)
rm(sce)
rm(pbmc.sce)

##############################################
formula = ~cogn_global_lv+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Exc1_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Exc1_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_raw_set2_final.rds")

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"
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
                gid = "AD",  # group IDs (ctrl/stim)
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

rm(pbmc)
rm(sce)
rm(pbmc.sce)


##############################################
formula = ~cogn_global_lv+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Exc2_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Exc2_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_raw_set3_final.rds")

#Add official metadata

setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

DefaultAssay(pbmc)="RNA"
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
                gid = "AD",  # group IDs (ctrl/stim)
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

rm(pbmc)
rm(sce)
rm(pbmc.sce)



##############################################
formula = ~cogn_global_lv+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Exc3_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  name=gsub("/","_",name)
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Exc3_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  name=gsub("/","_",name)
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
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
                gid = "AD",  # group IDs (ctrl/stim)
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

rm(pbmc)
rm(sce)
rm(pbmc.sce)


##############################################
formula = ~cogn_global_lv+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Immune_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Immune_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("In_integrated_batch_3000_module_scores_predicted_final.rds")

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
                gid = "AD",  # group IDs (ctrl/stim)
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

rm(pbmc)
rm(sce)
rm(pbmc.sce)



##############################################
formula = ~cogn_global_lv+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/In_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/In_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
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
#pbmc[["prediction.score.celltype"]] <- NULL
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
                gid = "AD",  # group IDs (ctrl/stim)
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

rm(pbmc)
rm(sce)
rm(pbmc.sce)



##############################################
formula = ~cogn_global_lv+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Oli_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Oli_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
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
                gid = "AD",  # group IDs (ctrl/stim)
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

rm(pbmc)
rm(sce)
rm(pbmc.sce)



##############################################
formula = ~cogn_global_lv+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Opc_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Opc_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
#######################################################################

setwd("G:/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Vasc_integrated_batch_5000_predicted_final.rds")

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
                gid = "AD",  # group IDs (ctrl/stim)
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

rm(pbmc)
rm(sce)
rm(pbmc.sce)



##############################################
formula = ~cogn_global_lv+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Vasc_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Vasc_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
#######################################################################

# MCT

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

Idents(pbmc)<-"cell_type_high_resolution"

HM=pbmc@meta.data
HM2=HM[,c("projid","cell_type_high_resolution")]

HM3 = HM2 %>% mutate(major_cell_type = case_when(grepl("Exc",cell_type_high_resolution)~"Exc",
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
                gid = "AD",  # group IDs (ctrl/stim)
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

rm(pbmc)
rm(sce)
rm(pbmc.sce)



##############################################
formula = ~cogn_global_lv+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Ast_major_cell_type_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Ast_major_cell_type_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

#######################################################################
#######################################################################
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

Idents(pbmc)<-"cell_type_high_resolution"

HM=pbmc@meta.data
HM2=HM[,c("projid","cell_type_high_resolution")]

HM3 = HM2 %>% mutate(major_cell_type = case_when(grepl("Exc",cell_type_high_resolution)~"Exc",
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
                gid = "AD",  # group IDs (ctrl/stim)
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

rm(pbmc)
rm(sce)
rm(pbmc.sce)



##############################################
formula = ~cogn_global_lv+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Mic_major_cell_type_cogn_global_lv.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogn_global_lv")

#Indicate name of the variable analyzed
variable="cogn_global_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~cogng_random_slope+pmi+age_death+plaq_n
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death","plaq_n")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/Mic_major_cell_type_cogng_random_slope.RData")

#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognition_plaq_n_corrected/cogng_random_slope")

#Indicate name of the variable analyzed
variable="cogng_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################

##############################################
