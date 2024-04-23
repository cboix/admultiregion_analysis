


library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)


##########################################################################################################

# MCT

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

#Add CDR score variables
setwd("G:/PFC_429_final_Spring2021/Vanshika")
CRscores=read.csv("PFC427_CDR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

DefaultAssay(pbmc)="RNA"
#pbmc[["integrated"]] <- NULL
#pbmc[["SCT"]] <- NULL
#pbmc[["prediction.score.celltype"]] <- NULL

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
formula = ~gpath_CDR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CDR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))



#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CDR_score")

#Indicate name of the variable analyzed
variable="gpath_CDR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~plaq_n_CDR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CDR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))


#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CDR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CDR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~nft_CDR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CDR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))


#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CDR_score")

#Indicate name of the variable analyzed
variable="nft_CDR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################
##############################################
formula = ~tangles_CDR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles_CDR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))


#Define output directory
setwd("G:/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles_CDR_score")

#Indicate name of the variable analyzed
variable="tangles_CDR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################