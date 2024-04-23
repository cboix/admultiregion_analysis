#Command line code
#conda activate muscat_env
#cd /mnt/f/PFC_429_final_Spring2021/Vanshika/muscat
#R < central_dogma.R --no-save


library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)


#Variables
full_list=list(
  "amyloid",
  "gpath",
  "gpath_3neocort",
  "nft",
  "nft_mf",
  "tangles",
  "plaq_n",
  "plaq_n_mf",
  "plaq_d",
  "plaq_d_mf",
  "bradysc_lv",
  "gaitsc_lv",
  "parksc_lv",
  "arteriol_scler",
  "caa_4gp",
  "cvda_4gp2",
  "ci_num2_gct",
  "ci_num2_gtt",
  "ci_num2_mct",
  "ci_num2_mtt",
  "diabetes_sr_rx_bl",
  "cancer_bl",
  "stroke_bl",
  "headinjrloc_bl",
  "heart_bl",
  "hypertension_bl",
  "tdp_st4",
  "msex",
  "cogn_global_lv",
  "cogn_ep_lv",
  "cogn_po_lv",
  "cogn_ps_lv",
  "cogn_se_lv",
  "cogn_wo_lv",
  "cognwo_random_slope",
  "cognse_random_slope",
  "cognpo_random_slope",
  "cognps_random_slope",
  "cognep_random_slope",
  "cogng_random_slope",
  "plaq_n_CR_score",
  "nft_CR_score",
  "gpath_CR_score"
)

#for(i in (1:length(full_list))){
#  path = paste0("F:/PFC_429_final_Spring2021/Vanshika/muscat/Results/",full_list[i])
#  dir.create(path)
#}

setwd("/mnt/f/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Ast_integrated_batch_5000_leiden_predicted_final.rds")

#Add official metadata

setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("/mnt/f/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

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
formula = ~cogn_global_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cogn_global_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_global_lv")

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
formula = ~cogn_ep_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ep_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cogn_ep_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ep_lv")

#Indicate name of the variable analyzed
variable="cogn_ep_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_po_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_po_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cogn_po_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_po_lv")

#Indicate name of the variable analyzed
variable="cogn_po_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_ps_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ps_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cogn_ps_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ps_lv")

#Indicate name of the variable analyzed
variable="cogn_ps_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_se_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_se_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cogn_se_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_se_lv")

#Indicate name of the variable analyzed
variable="cogn_se_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_wo_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_wo_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cogn_wo_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_wo_lv")

#Indicate name of the variable analyzed
variable="cogn_wo_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognwo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognwo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cognwo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognwo_random_slope")

#Indicate name of the variable analyzed
variable="cognwo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognse_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognse_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cognse_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognse_random_slope")

#Indicate name of the variable analyzed
variable="cognse_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognpo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognpo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cognpo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognpo_random_slope")

#Indicate name of the variable analyzed
variable="cognpo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognps_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognps_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cognps_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognps_random_slope")

#Indicate name of the variable analyzed
variable="cognps_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognep_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognep_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cognep_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognep_random_slope")

#Indicate name of the variable analyzed
variable="cognep_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogng_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cogng_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogng_random_slope")

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
formula = ~plaq_n_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_plaq_n_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_nft_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CR_score")

#Indicate name of the variable analyzed
variable="nft_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_gpath_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CR_score")

#Indicate name of the variable analyzed
variable="gpath_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~amyloid+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("amyloid","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_amyloid.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/amyloid")

#Indicate name of the variable analyzed
variable="amyloid"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_gpath.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

#Indicate name of the variable analyzed
variable="gpath"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_3neocort+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_3neocort","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_gpath_3neocort.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_3neocort")

#Indicate name of the variable analyzed
variable="gpath_3neocort"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_nft.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft")

#Indicate name of the variable analyzed
variable="nft"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_nft_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_mf")

#Indicate name of the variable analyzed
variable="nft_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tangles+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_tangles.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles")

#Indicate name of the variable analyzed
variable="tangles"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_plaq_n.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n")

#Indicate name of the variable analyzed
variable="plaq_n"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_plaq_n_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_mf")

#Indicate name of the variable analyzed
variable="plaq_n_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_plaq_d.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d")

#Indicate name of the variable analyzed
variable="plaq_d"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_plaq_d_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d_mf")

#Indicate name of the variable analyzed
variable="plaq_d_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~bradysc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("bradysc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_bradysc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/bradysc_lv")

#Indicate name of the variable analyzed
variable="bradysc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gaitsc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("gaitsc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_gaitsc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gaitsc_lv")

#Indicate name of the variable analyzed
variable="gaitsc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~parksc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("parksc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_parksc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/parksc_lv")

#Indicate name of the variable analyzed
variable="parksc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~arteriol_scler+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("arteriol_scler","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_arteriol_scler.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/arteriol_scler")

#Indicate name of the variable analyzed
variable="arteriol_scler"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~caa_4gp+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("caa_4gp","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_caa_4gp.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/caa_4gp")

#Indicate name of the variable analyzed
variable="caa_4gp"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cvda_4gp2+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cvda_4gp2","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cvda_4gp2.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cvda_4gp2")

#Indicate name of the variable analyzed
variable="cvda_4gp2"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_ci_num2_gct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gct")

#Indicate name of the variable analyzed
variable="ci_num2_gct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_ci_num2_gtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gtt")

#Indicate name of the variable analyzed
variable="ci_num2_gtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_ci_num2_mct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mct")

#Indicate name of the variable analyzed
variable="ci_num2_mct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_ci_num2_mtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mtt")

#Indicate name of the variable analyzed
variable="ci_num2_mtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~diabetes_sr_rx_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("diabetes_sr_rx_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_diabetes_sr_rx_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/diabetes_sr_rx_bl")

#Indicate name of the variable analyzed
variable="diabetes_sr_rx_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cancer_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cancer_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_cancer_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cancer_bl")

#Indicate name of the variable analyzed
variable="cancer_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~stroke_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("stroke_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_stroke_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/stroke_bl")

#Indicate name of the variable analyzed
variable="stroke_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~headinjrloc_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("headinjrloc_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_headinjrloc_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/headinjrloc_bl")

#Indicate name of the variable analyzed
variable="headinjrloc_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~heart_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("heart_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_heart_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/heart_bl")

#Indicate name of the variable analyzed
variable="heart_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~hypertension_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("hypertension_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_hypertension_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/hypertension_bl")

#Indicate name of the variable analyzed
variable="hypertension_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tdp_st4+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("tdp_st4","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_tdp_st4.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tdp_st4")

#Indicate name of the variable analyzed
variable="tdp_st4"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~msex+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("msex","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_msex.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/msex")

#Indicate name of the variable analyzed
variable="msex"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################



#######################################################################
#######################################################################
#######################################################################

setwd("/mnt/f/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_raw_set1_final.rds")

#Add official metadata

setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("/mnt/f/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

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
formula = ~cogn_global_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cogn_global_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_global_lv")

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
formula = ~cogn_ep_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ep_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cogn_ep_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ep_lv")

#Indicate name of the variable analyzed
variable="cogn_ep_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_po_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_po_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cogn_po_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_po_lv")

#Indicate name of the variable analyzed
variable="cogn_po_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_ps_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ps_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cogn_ps_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ps_lv")

#Indicate name of the variable analyzed
variable="cogn_ps_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_se_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_se_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cogn_se_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_se_lv")

#Indicate name of the variable analyzed
variable="cogn_se_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_wo_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_wo_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cogn_wo_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_wo_lv")

#Indicate name of the variable analyzed
variable="cogn_wo_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognwo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognwo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cognwo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognwo_random_slope")

#Indicate name of the variable analyzed
variable="cognwo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognse_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognse_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cognse_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognse_random_slope")

#Indicate name of the variable analyzed
variable="cognse_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognpo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognpo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cognpo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognpo_random_slope")

#Indicate name of the variable analyzed
variable="cognpo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognps_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognps_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cognps_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognps_random_slope")

#Indicate name of the variable analyzed
variable="cognps_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognep_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognep_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cognep_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognep_random_slope")

#Indicate name of the variable analyzed
variable="cognep_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogng_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cogng_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogng_random_slope")

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
formula = ~plaq_n_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_plaq_n_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_nft_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CR_score")

#Indicate name of the variable analyzed
variable="nft_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_gpath_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CR_score")

#Indicate name of the variable analyzed
variable="gpath_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~amyloid+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("amyloid","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_amyloid.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/amyloid")

#Indicate name of the variable analyzed
variable="amyloid"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_gpath.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

#Indicate name of the variable analyzed
variable="gpath"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_3neocort+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_3neocort","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_gpath_3neocort.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_3neocort")

#Indicate name of the variable analyzed
variable="gpath_3neocort"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_nft.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft")

#Indicate name of the variable analyzed
variable="nft"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_nft_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_mf")

#Indicate name of the variable analyzed
variable="nft_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tangles+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_tangles.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles")

#Indicate name of the variable analyzed
variable="tangles"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_plaq_n.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n")

#Indicate name of the variable analyzed
variable="plaq_n"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_plaq_n_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_mf")

#Indicate name of the variable analyzed
variable="plaq_n_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_plaq_d.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d")

#Indicate name of the variable analyzed
variable="plaq_d"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_plaq_d_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d_mf")

#Indicate name of the variable analyzed
variable="plaq_d_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~bradysc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("bradysc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_bradysc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/bradysc_lv")

#Indicate name of the variable analyzed
variable="bradysc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gaitsc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("gaitsc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_gaitsc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gaitsc_lv")

#Indicate name of the variable analyzed
variable="gaitsc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~parksc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("parksc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_parksc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/parksc_lv")

#Indicate name of the variable analyzed
variable="parksc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~arteriol_scler+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("arteriol_scler","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_arteriol_scler.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/arteriol_scler")

#Indicate name of the variable analyzed
variable="arteriol_scler"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~caa_4gp+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("caa_4gp","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_caa_4gp.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/caa_4gp")

#Indicate name of the variable analyzed
variable="caa_4gp"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cvda_4gp2+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cvda_4gp2","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cvda_4gp2.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cvda_4gp2")

#Indicate name of the variable analyzed
variable="cvda_4gp2"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_ci_num2_gct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gct")

#Indicate name of the variable analyzed
variable="ci_num2_gct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_ci_num2_gtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gtt")

#Indicate name of the variable analyzed
variable="ci_num2_gtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_ci_num2_mct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mct")

#Indicate name of the variable analyzed
variable="ci_num2_mct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_ci_num2_mtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mtt")

#Indicate name of the variable analyzed
variable="ci_num2_mtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~diabetes_sr_rx_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("diabetes_sr_rx_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_diabetes_sr_rx_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/diabetes_sr_rx_bl")

#Indicate name of the variable analyzed
variable="diabetes_sr_rx_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cancer_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cancer_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_cancer_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cancer_bl")

#Indicate name of the variable analyzed
variable="cancer_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~stroke_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("stroke_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_stroke_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/stroke_bl")

#Indicate name of the variable analyzed
variable="stroke_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~headinjrloc_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("headinjrloc_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_headinjrloc_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/headinjrloc_bl")

#Indicate name of the variable analyzed
variable="headinjrloc_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~heart_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("heart_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_heart_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/heart_bl")

#Indicate name of the variable analyzed
variable="heart_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~hypertension_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("hypertension_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_hypertension_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/hypertension_bl")

#Indicate name of the variable analyzed
variable="hypertension_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tdp_st4+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("tdp_st4","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_tdp_st4.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tdp_st4")

#Indicate name of the variable analyzed
variable="tdp_st4"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~msex+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("msex","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc1_msex.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/msex")

#Indicate name of the variable analyzed
variable="msex"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################




#######################################################################
#######################################################################
#######################################################################

setwd("/mnt/f/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_raw_set2_final.rds")

#Add official metadata

setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("/mnt/f/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

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
formula = ~cogn_global_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cogn_global_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_global_lv")

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
formula = ~cogn_ep_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ep_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cogn_ep_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ep_lv")

#Indicate name of the variable analyzed
variable="cogn_ep_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_po_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_po_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cogn_po_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_po_lv")

#Indicate name of the variable analyzed
variable="cogn_po_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_ps_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ps_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cogn_ps_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ps_lv")

#Indicate name of the variable analyzed
variable="cogn_ps_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_se_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_se_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cogn_se_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_se_lv")

#Indicate name of the variable analyzed
variable="cogn_se_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_wo_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_wo_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cogn_wo_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_wo_lv")

#Indicate name of the variable analyzed
variable="cogn_wo_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognwo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognwo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cognwo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognwo_random_slope")

#Indicate name of the variable analyzed
variable="cognwo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognse_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognse_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cognse_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognse_random_slope")

#Indicate name of the variable analyzed
variable="cognse_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognpo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognpo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cognpo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognpo_random_slope")

#Indicate name of the variable analyzed
variable="cognpo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognps_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognps_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cognps_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognps_random_slope")

#Indicate name of the variable analyzed
variable="cognps_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognep_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognep_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cognep_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognep_random_slope")

#Indicate name of the variable analyzed
variable="cognep_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogng_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cogng_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogng_random_slope")

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
formula = ~plaq_n_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_plaq_n_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_nft_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CR_score")

#Indicate name of the variable analyzed
variable="nft_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_gpath_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CR_score")

#Indicate name of the variable analyzed
variable="gpath_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~amyloid+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("amyloid","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_amyloid.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/amyloid")

#Indicate name of the variable analyzed
variable="amyloid"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_gpath.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

#Indicate name of the variable analyzed
variable="gpath"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_3neocort+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_3neocort","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_gpath_3neocort.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_3neocort")

#Indicate name of the variable analyzed
variable="gpath_3neocort"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_nft.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft")

#Indicate name of the variable analyzed
variable="nft"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_nft_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_mf")

#Indicate name of the variable analyzed
variable="nft_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tangles+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_tangles.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles")

#Indicate name of the variable analyzed
variable="tangles"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_plaq_n.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n")

#Indicate name of the variable analyzed
variable="plaq_n"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_plaq_n_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_mf")

#Indicate name of the variable analyzed
variable="plaq_n_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_plaq_d.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d")

#Indicate name of the variable analyzed
variable="plaq_d"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_plaq_d_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d_mf")

#Indicate name of the variable analyzed
variable="plaq_d_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~bradysc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("bradysc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_bradysc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/bradysc_lv")

#Indicate name of the variable analyzed
variable="bradysc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gaitsc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("gaitsc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_gaitsc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gaitsc_lv")

#Indicate name of the variable analyzed
variable="gaitsc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~parksc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("parksc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_parksc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/parksc_lv")

#Indicate name of the variable analyzed
variable="parksc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~arteriol_scler+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("arteriol_scler","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_arteriol_scler.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/arteriol_scler")

#Indicate name of the variable analyzed
variable="arteriol_scler"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~caa_4gp+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("caa_4gp","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_caa_4gp.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/caa_4gp")

#Indicate name of the variable analyzed
variable="caa_4gp"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cvda_4gp2+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cvda_4gp2","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cvda_4gp2.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cvda_4gp2")

#Indicate name of the variable analyzed
variable="cvda_4gp2"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_ci_num2_gct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gct")

#Indicate name of the variable analyzed
variable="ci_num2_gct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_ci_num2_gtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gtt")

#Indicate name of the variable analyzed
variable="ci_num2_gtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_ci_num2_mct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mct")

#Indicate name of the variable analyzed
variable="ci_num2_mct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_ci_num2_mtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mtt")

#Indicate name of the variable analyzed
variable="ci_num2_mtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~diabetes_sr_rx_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("diabetes_sr_rx_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_diabetes_sr_rx_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/diabetes_sr_rx_bl")

#Indicate name of the variable analyzed
variable="diabetes_sr_rx_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cancer_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cancer_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_cancer_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cancer_bl")

#Indicate name of the variable analyzed
variable="cancer_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~stroke_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("stroke_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_stroke_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/stroke_bl")

#Indicate name of the variable analyzed
variable="stroke_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~headinjrloc_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("headinjrloc_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_headinjrloc_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/headinjrloc_bl")

#Indicate name of the variable analyzed
variable="headinjrloc_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~heart_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("heart_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_heart_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/heart_bl")

#Indicate name of the variable analyzed
variable="heart_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~hypertension_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("hypertension_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_hypertension_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/hypertension_bl")

#Indicate name of the variable analyzed
variable="hypertension_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tdp_st4+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("tdp_st4","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_tdp_st4.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tdp_st4")

#Indicate name of the variable analyzed
variable="tdp_st4"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~msex+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("msex","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc2_msex.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/msex")

#Indicate name of the variable analyzed
variable="msex"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################



#######################################################################
#######################################################################
#######################################################################

setwd("/mnt/f/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Exc_raw_set3_final.rds")

#Add official metadata

setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("/mnt/f/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

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
formula = ~cogn_global_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cogn_global_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_global_lv")

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
formula = ~cogn_ep_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ep_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cogn_ep_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ep_lv")

#Indicate name of the variable analyzed
variable="cogn_ep_lv"

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
formula = ~cogn_po_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_po_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cogn_po_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_po_lv")

#Indicate name of the variable analyzed
variable="cogn_po_lv"

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
formula = ~cogn_ps_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ps_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cogn_ps_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ps_lv")

#Indicate name of the variable analyzed
variable="cogn_ps_lv"

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
formula = ~cogn_se_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_se_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cogn_se_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_se_lv")

#Indicate name of the variable analyzed
variable="cogn_se_lv"

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
formula = ~cogn_wo_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_wo_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cogn_wo_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_wo_lv")

#Indicate name of the variable analyzed
variable="cogn_wo_lv"

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
formula = ~cognwo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognwo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cognwo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognwo_random_slope")

#Indicate name of the variable analyzed
variable="cognwo_random_slope"

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
formula = ~cognse_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognse_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cognse_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognse_random_slope")

#Indicate name of the variable analyzed
variable="cognse_random_slope"

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
formula = ~cognpo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognpo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cognpo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognpo_random_slope")

#Indicate name of the variable analyzed
variable="cognpo_random_slope"

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
formula = ~cognps_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognps_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cognps_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognps_random_slope")

#Indicate name of the variable analyzed
variable="cognps_random_slope"

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
formula = ~cognep_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognep_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cognep_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognep_random_slope")

#Indicate name of the variable analyzed
variable="cognep_random_slope"

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
formula = ~cogng_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cogng_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogng_random_slope")

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
formula = ~plaq_n_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_plaq_n_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CR_score"

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
formula = ~nft_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_nft_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CR_score")

#Indicate name of the variable analyzed
variable="nft_CR_score"

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
formula = ~gpath_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_gpath_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CR_score")

#Indicate name of the variable analyzed
variable="gpath_CR_score"

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
formula = ~amyloid+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("amyloid","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_amyloid.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/amyloid")

#Indicate name of the variable analyzed
variable="amyloid"

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
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_gpath.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

#Indicate name of the variable analyzed
variable="gpath"

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
formula = ~gpath_3neocort+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_3neocort","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_gpath_3neocort.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_3neocort")

#Indicate name of the variable analyzed
variable="gpath_3neocort"

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
formula = ~nft+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_nft.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft")

#Indicate name of the variable analyzed
variable="nft"

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
formula = ~nft_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_nft_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_mf")

#Indicate name of the variable analyzed
variable="nft_mf"

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
formula = ~tangles+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_tangles.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles")

#Indicate name of the variable analyzed
variable="tangles"

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
formula = ~plaq_n+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_plaq_n.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n")

#Indicate name of the variable analyzed
variable="plaq_n"

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
formula = ~plaq_n_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_plaq_n_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_mf")

#Indicate name of the variable analyzed
variable="plaq_n_mf"

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
formula = ~plaq_d+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_plaq_d.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d")

#Indicate name of the variable analyzed
variable="plaq_d"

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
formula = ~plaq_d_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_plaq_d_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d_mf")

#Indicate name of the variable analyzed
variable="plaq_d_mf"

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
formula = ~bradysc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("bradysc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_bradysc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/bradysc_lv")

#Indicate name of the variable analyzed
variable="bradysc_lv"

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
formula = ~gaitsc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("gaitsc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_gaitsc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gaitsc_lv")

#Indicate name of the variable analyzed
variable="gaitsc_lv"

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
formula = ~parksc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("parksc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_parksc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/parksc_lv")

#Indicate name of the variable analyzed
variable="parksc_lv"

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
formula = ~arteriol_scler+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("arteriol_scler","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_arteriol_scler.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/arteriol_scler")

#Indicate name of the variable analyzed
variable="arteriol_scler"

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
formula = ~caa_4gp+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("caa_4gp","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_caa_4gp.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/caa_4gp")

#Indicate name of the variable analyzed
variable="caa_4gp"

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
formula = ~cvda_4gp2+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cvda_4gp2","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cvda_4gp2.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cvda_4gp2")

#Indicate name of the variable analyzed
variable="cvda_4gp2"

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
formula = ~ci_num2_gct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_ci_num2_gct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gct")

#Indicate name of the variable analyzed
variable="ci_num2_gct"

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
formula = ~ci_num2_gtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_ci_num2_gtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gtt")

#Indicate name of the variable analyzed
variable="ci_num2_gtt"

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
formula = ~ci_num2_mct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_ci_num2_mct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mct")

#Indicate name of the variable analyzed
variable="ci_num2_mct"

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
formula = ~ci_num2_mtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_ci_num2_mtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mtt")

#Indicate name of the variable analyzed
variable="ci_num2_mtt"

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
formula = ~diabetes_sr_rx_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("diabetes_sr_rx_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_diabetes_sr_rx_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/diabetes_sr_rx_bl")

#Indicate name of the variable analyzed
variable="diabetes_sr_rx_bl"

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
formula = ~cancer_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cancer_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_cancer_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cancer_bl")

#Indicate name of the variable analyzed
variable="cancer_bl"

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
formula = ~stroke_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("stroke_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_stroke_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/stroke_bl")

#Indicate name of the variable analyzed
variable="stroke_bl"

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
formula = ~headinjrloc_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("headinjrloc_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_headinjrloc_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/headinjrloc_bl")

#Indicate name of the variable analyzed
variable="headinjrloc_bl"

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
formula = ~heart_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("heart_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_heart_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/heart_bl")

#Indicate name of the variable analyzed
variable="heart_bl"

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
formula = ~hypertension_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("hypertension_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_hypertension_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/hypertension_bl")

#Indicate name of the variable analyzed
variable="hypertension_bl"

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
formula = ~tdp_st4+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("tdp_st4","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_tdp_st4.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tdp_st4")

#Indicate name of the variable analyzed
variable="tdp_st4"

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
formula = ~msex+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("msex","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Exc3_msex.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/msex")

#Indicate name of the variable analyzed
variable="msex"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  name=gsub("/","_",name)
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################


#######################################################################
#######################################################################
#######################################################################

setwd("/mnt/f/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Immune_integrated_batch_5000_predicted_final.rds")

#Add official metadata

setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("/mnt/f/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

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
formula = ~cogn_global_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cogn_global_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_global_lv")

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
formula = ~cogn_ep_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ep_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cogn_ep_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ep_lv")

#Indicate name of the variable analyzed
variable="cogn_ep_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_po_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_po_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cogn_po_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_po_lv")

#Indicate name of the variable analyzed
variable="cogn_po_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_ps_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ps_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cogn_ps_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ps_lv")

#Indicate name of the variable analyzed
variable="cogn_ps_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_se_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_se_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cogn_se_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_se_lv")

#Indicate name of the variable analyzed
variable="cogn_se_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_wo_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_wo_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cogn_wo_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_wo_lv")

#Indicate name of the variable analyzed
variable="cogn_wo_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognwo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognwo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cognwo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognwo_random_slope")

#Indicate name of the variable analyzed
variable="cognwo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognse_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognse_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cognse_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognse_random_slope")

#Indicate name of the variable analyzed
variable="cognse_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognpo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognpo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cognpo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognpo_random_slope")

#Indicate name of the variable analyzed
variable="cognpo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognps_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognps_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cognps_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognps_random_slope")

#Indicate name of the variable analyzed
variable="cognps_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognep_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognep_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cognep_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognep_random_slope")

#Indicate name of the variable analyzed
variable="cognep_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogng_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cogng_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogng_random_slope")

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
formula = ~plaq_n_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_plaq_n_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_nft_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CR_score")

#Indicate name of the variable analyzed
variable="nft_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_gpath_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CR_score")

#Indicate name of the variable analyzed
variable="gpath_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~amyloid+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("amyloid","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_amyloid.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/amyloid")

#Indicate name of the variable analyzed
variable="amyloid"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_gpath.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

#Indicate name of the variable analyzed
variable="gpath"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_3neocort+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_3neocort","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_gpath_3neocort.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_3neocort")

#Indicate name of the variable analyzed
variable="gpath_3neocort"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_nft.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft")

#Indicate name of the variable analyzed
variable="nft"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_nft_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_mf")

#Indicate name of the variable analyzed
variable="nft_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tangles+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_tangles.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles")

#Indicate name of the variable analyzed
variable="tangles"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_plaq_n.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n")

#Indicate name of the variable analyzed
variable="plaq_n"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_plaq_n_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_mf")

#Indicate name of the variable analyzed
variable="plaq_n_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_plaq_d.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d")

#Indicate name of the variable analyzed
variable="plaq_d"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_plaq_d_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d_mf")

#Indicate name of the variable analyzed
variable="plaq_d_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~bradysc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("bradysc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_bradysc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/bradysc_lv")

#Indicate name of the variable analyzed
variable="bradysc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gaitsc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("gaitsc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_gaitsc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gaitsc_lv")

#Indicate name of the variable analyzed
variable="gaitsc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~parksc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("parksc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_parksc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/parksc_lv")

#Indicate name of the variable analyzed
variable="parksc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~arteriol_scler+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("arteriol_scler","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_arteriol_scler.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/arteriol_scler")

#Indicate name of the variable analyzed
variable="arteriol_scler"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~caa_4gp+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("caa_4gp","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_caa_4gp.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/caa_4gp")

#Indicate name of the variable analyzed
variable="caa_4gp"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cvda_4gp2+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cvda_4gp2","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cvda_4gp2.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cvda_4gp2")

#Indicate name of the variable analyzed
variable="cvda_4gp2"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_ci_num2_gct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gct")

#Indicate name of the variable analyzed
variable="ci_num2_gct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_ci_num2_gtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gtt")

#Indicate name of the variable analyzed
variable="ci_num2_gtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_ci_num2_mct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mct")

#Indicate name of the variable analyzed
variable="ci_num2_mct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_ci_num2_mtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mtt")

#Indicate name of the variable analyzed
variable="ci_num2_mtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~diabetes_sr_rx_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("diabetes_sr_rx_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_diabetes_sr_rx_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/diabetes_sr_rx_bl")

#Indicate name of the variable analyzed
variable="diabetes_sr_rx_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cancer_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cancer_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_cancer_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cancer_bl")

#Indicate name of the variable analyzed
variable="cancer_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~stroke_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("stroke_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_stroke_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/stroke_bl")

#Indicate name of the variable analyzed
variable="stroke_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~headinjrloc_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("headinjrloc_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_headinjrloc_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/headinjrloc_bl")

#Indicate name of the variable analyzed
variable="headinjrloc_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~heart_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("heart_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_heart_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/heart_bl")

#Indicate name of the variable analyzed
variable="heart_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~hypertension_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("hypertension_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_hypertension_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/hypertension_bl")

#Indicate name of the variable analyzed
variable="hypertension_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tdp_st4+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("tdp_st4","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_tdp_st4.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tdp_st4")

#Indicate name of the variable analyzed
variable="tdp_st4"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~msex+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("msex","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Immune_msex.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/msex")

#Indicate name of the variable analyzed
variable="msex"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################



#######################################################################
#######################################################################
#######################################################################

setwd("/mnt/f/PFC_429_final_Spring2021/Data")
pbmc=readRDS("In_integrated_batch_3000_module_scores_predicted_final.rds")

#Add official metadata

setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("/mnt/f/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

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
formula = ~cogn_global_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cogn_global_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_global_lv")

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
formula = ~cogn_ep_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ep_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cogn_ep_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ep_lv")

#Indicate name of the variable analyzed
variable="cogn_ep_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_po_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_po_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cogn_po_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_po_lv")

#Indicate name of the variable analyzed
variable="cogn_po_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_ps_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ps_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cogn_ps_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ps_lv")

#Indicate name of the variable analyzed
variable="cogn_ps_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_se_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_se_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cogn_se_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_se_lv")

#Indicate name of the variable analyzed
variable="cogn_se_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_wo_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_wo_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cogn_wo_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_wo_lv")

#Indicate name of the variable analyzed
variable="cogn_wo_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognwo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognwo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cognwo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognwo_random_slope")

#Indicate name of the variable analyzed
variable="cognwo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognse_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognse_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cognse_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognse_random_slope")

#Indicate name of the variable analyzed
variable="cognse_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognpo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognpo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cognpo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognpo_random_slope")

#Indicate name of the variable analyzed
variable="cognpo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognps_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognps_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cognps_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognps_random_slope")

#Indicate name of the variable analyzed
variable="cognps_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognep_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognep_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cognep_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognep_random_slope")

#Indicate name of the variable analyzed
variable="cognep_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogng_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cogng_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogng_random_slope")

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
formula = ~plaq_n_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_plaq_n_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_nft_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CR_score")

#Indicate name of the variable analyzed
variable="nft_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_gpath_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CR_score")

#Indicate name of the variable analyzed
variable="gpath_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~amyloid+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("amyloid","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_amyloid.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/amyloid")

#Indicate name of the variable analyzed
variable="amyloid"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_gpath.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

#Indicate name of the variable analyzed
variable="gpath"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_3neocort+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_3neocort","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_gpath_3neocort.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_3neocort")

#Indicate name of the variable analyzed
variable="gpath_3neocort"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_nft.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft")

#Indicate name of the variable analyzed
variable="nft"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_nft_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_mf")

#Indicate name of the variable analyzed
variable="nft_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tangles+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_tangles.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles")

#Indicate name of the variable analyzed
variable="tangles"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_plaq_n.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n")

#Indicate name of the variable analyzed
variable="plaq_n"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_plaq_n_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_mf")

#Indicate name of the variable analyzed
variable="plaq_n_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_plaq_d.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d")

#Indicate name of the variable analyzed
variable="plaq_d"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_plaq_d_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d_mf")

#Indicate name of the variable analyzed
variable="plaq_d_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~bradysc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("bradysc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_bradysc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/bradysc_lv")

#Indicate name of the variable analyzed
variable="bradysc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gaitsc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("gaitsc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_gaitsc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gaitsc_lv")

#Indicate name of the variable analyzed
variable="gaitsc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~parksc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("parksc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_parksc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/parksc_lv")

#Indicate name of the variable analyzed
variable="parksc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~arteriol_scler+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("arteriol_scler","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_arteriol_scler.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/arteriol_scler")

#Indicate name of the variable analyzed
variable="arteriol_scler"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~caa_4gp+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("caa_4gp","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_caa_4gp.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/caa_4gp")

#Indicate name of the variable analyzed
variable="caa_4gp"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cvda_4gp2+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cvda_4gp2","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cvda_4gp2.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cvda_4gp2")

#Indicate name of the variable analyzed
variable="cvda_4gp2"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_ci_num2_gct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gct")

#Indicate name of the variable analyzed
variable="ci_num2_gct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_ci_num2_gtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gtt")

#Indicate name of the variable analyzed
variable="ci_num2_gtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_ci_num2_mct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mct")

#Indicate name of the variable analyzed
variable="ci_num2_mct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_ci_num2_mtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mtt")

#Indicate name of the variable analyzed
variable="ci_num2_mtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~diabetes_sr_rx_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("diabetes_sr_rx_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_diabetes_sr_rx_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/diabetes_sr_rx_bl")

#Indicate name of the variable analyzed
variable="diabetes_sr_rx_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cancer_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cancer_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_cancer_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cancer_bl")

#Indicate name of the variable analyzed
variable="cancer_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~stroke_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("stroke_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_stroke_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/stroke_bl")

#Indicate name of the variable analyzed
variable="stroke_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~headinjrloc_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("headinjrloc_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_headinjrloc_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/headinjrloc_bl")

#Indicate name of the variable analyzed
variable="headinjrloc_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~heart_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("heart_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_heart_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/heart_bl")

#Indicate name of the variable analyzed
variable="heart_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~hypertension_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("hypertension_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_hypertension_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/hypertension_bl")

#Indicate name of the variable analyzed
variable="hypertension_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tdp_st4+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("tdp_st4","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_tdp_st4.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tdp_st4")

#Indicate name of the variable analyzed
variable="tdp_st4"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~msex+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("msex","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Inh_msex.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/msex")

#Indicate name of the variable analyzed
variable="msex"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################



#######################################################################
#######################################################################
#######################################################################

setwd("/mnt/f/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Oli_integrated_batch3000_leiden_res02_final.rds")

#Add official metadata

setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("/mnt/f/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

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
formula = ~cogn_global_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cogn_global_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_global_lv")

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
formula = ~cogn_ep_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ep_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cogn_ep_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ep_lv")

#Indicate name of the variable analyzed
variable="cogn_ep_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_po_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_po_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cogn_po_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_po_lv")

#Indicate name of the variable analyzed
variable="cogn_po_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_ps_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ps_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cogn_ps_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ps_lv")

#Indicate name of the variable analyzed
variable="cogn_ps_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_se_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_se_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cogn_se_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_se_lv")

#Indicate name of the variable analyzed
variable="cogn_se_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_wo_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_wo_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cogn_wo_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_wo_lv")

#Indicate name of the variable analyzed
variable="cogn_wo_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognwo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognwo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cognwo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognwo_random_slope")

#Indicate name of the variable analyzed
variable="cognwo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognse_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognse_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cognse_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognse_random_slope")

#Indicate name of the variable analyzed
variable="cognse_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognpo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognpo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cognpo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognpo_random_slope")

#Indicate name of the variable analyzed
variable="cognpo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognps_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognps_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cognps_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognps_random_slope")

#Indicate name of the variable analyzed
variable="cognps_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognep_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognep_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cognep_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognep_random_slope")

#Indicate name of the variable analyzed
variable="cognep_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogng_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cogng_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogng_random_slope")

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
formula = ~plaq_n_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_plaq_n_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_nft_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CR_score")

#Indicate name of the variable analyzed
variable="nft_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_gpath_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CR_score")

#Indicate name of the variable analyzed
variable="gpath_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~amyloid+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("amyloid","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_amyloid.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/amyloid")

#Indicate name of the variable analyzed
variable="amyloid"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_gpath.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

#Indicate name of the variable analyzed
variable="gpath"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_3neocort+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_3neocort","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_gpath_3neocort.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_3neocort")

#Indicate name of the variable analyzed
variable="gpath_3neocort"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_nft.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft")

#Indicate name of the variable analyzed
variable="nft"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_nft_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_mf")

#Indicate name of the variable analyzed
variable="nft_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tangles+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_tangles.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles")

#Indicate name of the variable analyzed
variable="tangles"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_plaq_n.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n")

#Indicate name of the variable analyzed
variable="plaq_n"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_plaq_n_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_mf")

#Indicate name of the variable analyzed
variable="plaq_n_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_plaq_d.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d")

#Indicate name of the variable analyzed
variable="plaq_d"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_plaq_d_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d_mf")

#Indicate name of the variable analyzed
variable="plaq_d_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~bradysc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("bradysc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_bradysc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/bradysc_lv")

#Indicate name of the variable analyzed
variable="bradysc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gaitsc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("gaitsc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_gaitsc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gaitsc_lv")

#Indicate name of the variable analyzed
variable="gaitsc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~parksc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("parksc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_parksc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/parksc_lv")

#Indicate name of the variable analyzed
variable="parksc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~arteriol_scler+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("arteriol_scler","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_arteriol_scler.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/arteriol_scler")

#Indicate name of the variable analyzed
variable="arteriol_scler"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~caa_4gp+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("caa_4gp","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_caa_4gp.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/caa_4gp")

#Indicate name of the variable analyzed
variable="caa_4gp"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cvda_4gp2+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cvda_4gp2","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cvda_4gp2.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cvda_4gp2")

#Indicate name of the variable analyzed
variable="cvda_4gp2"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_ci_num2_gct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gct")

#Indicate name of the variable analyzed
variable="ci_num2_gct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_ci_num2_gtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gtt")

#Indicate name of the variable analyzed
variable="ci_num2_gtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_ci_num2_mct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mct")

#Indicate name of the variable analyzed
variable="ci_num2_mct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_ci_num2_mtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mtt")

#Indicate name of the variable analyzed
variable="ci_num2_mtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~diabetes_sr_rx_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("diabetes_sr_rx_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_diabetes_sr_rx_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/diabetes_sr_rx_bl")

#Indicate name of the variable analyzed
variable="diabetes_sr_rx_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cancer_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cancer_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_cancer_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cancer_bl")

#Indicate name of the variable analyzed
variable="cancer_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~stroke_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("stroke_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_stroke_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/stroke_bl")

#Indicate name of the variable analyzed
variable="stroke_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~headinjrloc_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("headinjrloc_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_headinjrloc_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/headinjrloc_bl")

#Indicate name of the variable analyzed
variable="headinjrloc_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~heart_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("heart_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_heart_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/heart_bl")

#Indicate name of the variable analyzed
variable="heart_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~hypertension_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("hypertension_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_hypertension_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/hypertension_bl")

#Indicate name of the variable analyzed
variable="hypertension_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tdp_st4+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("tdp_st4","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_tdp_st4.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tdp_st4")

#Indicate name of the variable analyzed
variable="tdp_st4"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~msex+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("msex","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Oli_msex.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/msex")

#Indicate name of the variable analyzed
variable="msex"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################



#######################################################################
#######################################################################
#######################################################################

setwd("/mnt/f/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Opc_integrated_batch_5000_leiden_predicted_final.rds")

#Add official metadata

setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("/mnt/f/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

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
formula = ~cogn_global_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cogn_global_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_global_lv")

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
formula = ~cogn_ep_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ep_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cogn_ep_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ep_lv")

#Indicate name of the variable analyzed
variable="cogn_ep_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_po_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_po_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cogn_po_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_po_lv")

#Indicate name of the variable analyzed
variable="cogn_po_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_ps_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ps_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cogn_ps_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ps_lv")

#Indicate name of the variable analyzed
variable="cogn_ps_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_se_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_se_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cogn_se_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_se_lv")

#Indicate name of the variable analyzed
variable="cogn_se_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_wo_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_wo_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cogn_wo_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_wo_lv")

#Indicate name of the variable analyzed
variable="cogn_wo_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognwo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognwo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cognwo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognwo_random_slope")

#Indicate name of the variable analyzed
variable="cognwo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognse_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognse_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cognse_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognse_random_slope")

#Indicate name of the variable analyzed
variable="cognse_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognpo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognpo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cognpo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognpo_random_slope")

#Indicate name of the variable analyzed
variable="cognpo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognps_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognps_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cognps_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognps_random_slope")

#Indicate name of the variable analyzed
variable="cognps_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognep_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognep_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cognep_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognep_random_slope")

#Indicate name of the variable analyzed
variable="cognep_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogng_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cogng_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogng_random_slope")

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
formula = ~plaq_n_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_plaq_n_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_nft_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CR_score")

#Indicate name of the variable analyzed
variable="nft_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_gpath_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CR_score")

#Indicate name of the variable analyzed
variable="gpath_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~amyloid+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("amyloid","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_amyloid.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/amyloid")

#Indicate name of the variable analyzed
variable="amyloid"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_gpath.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

#Indicate name of the variable analyzed
variable="gpath"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_3neocort+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_3neocort","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_gpath_3neocort.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_3neocort")

#Indicate name of the variable analyzed
variable="gpath_3neocort"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_nft.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft")

#Indicate name of the variable analyzed
variable="nft"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_nft_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_mf")

#Indicate name of the variable analyzed
variable="nft_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tangles+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_tangles.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles")

#Indicate name of the variable analyzed
variable="tangles"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_plaq_n.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n")

#Indicate name of the variable analyzed
variable="plaq_n"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_plaq_n_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_mf")

#Indicate name of the variable analyzed
variable="plaq_n_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_plaq_d.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d")

#Indicate name of the variable analyzed
variable="plaq_d"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_plaq_d_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d_mf")

#Indicate name of the variable analyzed
variable="plaq_d_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~bradysc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("bradysc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_bradysc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/bradysc_lv")

#Indicate name of the variable analyzed
variable="bradysc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gaitsc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("gaitsc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_gaitsc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gaitsc_lv")

#Indicate name of the variable analyzed
variable="gaitsc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~parksc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("parksc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_parksc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/parksc_lv")

#Indicate name of the variable analyzed
variable="parksc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~arteriol_scler+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("arteriol_scler","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_arteriol_scler.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/arteriol_scler")

#Indicate name of the variable analyzed
variable="arteriol_scler"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~caa_4gp+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("caa_4gp","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_caa_4gp.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/caa_4gp")

#Indicate name of the variable analyzed
variable="caa_4gp"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cvda_4gp2+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cvda_4gp2","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cvda_4gp2.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cvda_4gp2")

#Indicate name of the variable analyzed
variable="cvda_4gp2"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_ci_num2_gct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gct")

#Indicate name of the variable analyzed
variable="ci_num2_gct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_ci_num2_gtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gtt")

#Indicate name of the variable analyzed
variable="ci_num2_gtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_ci_num2_mct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mct")

#Indicate name of the variable analyzed
variable="ci_num2_mct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_ci_num2_mtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mtt")

#Indicate name of the variable analyzed
variable="ci_num2_mtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~diabetes_sr_rx_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("diabetes_sr_rx_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_diabetes_sr_rx_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/diabetes_sr_rx_bl")

#Indicate name of the variable analyzed
variable="diabetes_sr_rx_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cancer_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cancer_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_cancer_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cancer_bl")

#Indicate name of the variable analyzed
variable="cancer_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~stroke_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("stroke_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_stroke_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/stroke_bl")

#Indicate name of the variable analyzed
variable="stroke_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~headinjrloc_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("headinjrloc_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_headinjrloc_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/headinjrloc_bl")

#Indicate name of the variable analyzed
variable="headinjrloc_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~heart_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("heart_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_heart_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/heart_bl")

#Indicate name of the variable analyzed
variable="heart_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~hypertension_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("hypertension_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_hypertension_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/hypertension_bl")

#Indicate name of the variable analyzed
variable="hypertension_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tdp_st4+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("tdp_st4","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_tdp_st4.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tdp_st4")

#Indicate name of the variable analyzed
variable="tdp_st4"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~msex+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("msex","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/OPC_msex.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/msex")

#Indicate name of the variable analyzed
variable="msex"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################



#######################################################################
#######################################################################
#######################################################################

setwd("/mnt/f/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Vasc_integrated_batch_5000_predicted_final.rds")

#Add official metadata

setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("/mnt/f/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

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
formula = ~cogn_global_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cogn_global_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_global_lv")

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
formula = ~cogn_ep_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ep_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cogn_ep_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ep_lv")

#Indicate name of the variable analyzed
variable="cogn_ep_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_po_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_po_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cogn_po_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_po_lv")

#Indicate name of the variable analyzed
variable="cogn_po_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_ps_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ps_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cogn_ps_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ps_lv")

#Indicate name of the variable analyzed
variable="cogn_ps_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_se_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_se_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cogn_se_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_se_lv")

#Indicate name of the variable analyzed
variable="cogn_se_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_wo_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_wo_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cogn_wo_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_wo_lv")

#Indicate name of the variable analyzed
variable="cogn_wo_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognwo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognwo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cognwo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognwo_random_slope")

#Indicate name of the variable analyzed
variable="cognwo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognse_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognse_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cognse_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognse_random_slope")

#Indicate name of the variable analyzed
variable="cognse_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognpo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognpo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cognpo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognpo_random_slope")

#Indicate name of the variable analyzed
variable="cognpo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognps_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognps_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cognps_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognps_random_slope")

#Indicate name of the variable analyzed
variable="cognps_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognep_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognep_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cognep_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognep_random_slope")

#Indicate name of the variable analyzed
variable="cognep_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogng_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cogng_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogng_random_slope")

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
formula = ~plaq_n_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_plaq_n_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_nft_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CR_score")

#Indicate name of the variable analyzed
variable="nft_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_gpath_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CR_score")

#Indicate name of the variable analyzed
variable="gpath_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~amyloid+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("amyloid","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_amyloid.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/amyloid")

#Indicate name of the variable analyzed
variable="amyloid"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_gpath.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

#Indicate name of the variable analyzed
variable="gpath"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_3neocort+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_3neocort","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_gpath_3neocort.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_3neocort")

#Indicate name of the variable analyzed
variable="gpath_3neocort"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_nft.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft")

#Indicate name of the variable analyzed
variable="nft"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_nft_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_mf")

#Indicate name of the variable analyzed
variable="nft_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tangles+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_tangles.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles")

#Indicate name of the variable analyzed
variable="tangles"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_plaq_n.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n")

#Indicate name of the variable analyzed
variable="plaq_n"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_plaq_n_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_mf")

#Indicate name of the variable analyzed
variable="plaq_n_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_plaq_d.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d")

#Indicate name of the variable analyzed
variable="plaq_d"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_plaq_d_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d_mf")

#Indicate name of the variable analyzed
variable="plaq_d_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~bradysc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("bradysc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_bradysc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/bradysc_lv")

#Indicate name of the variable analyzed
variable="bradysc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gaitsc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("gaitsc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_gaitsc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gaitsc_lv")

#Indicate name of the variable analyzed
variable="gaitsc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~parksc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("parksc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_parksc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/parksc_lv")

#Indicate name of the variable analyzed
variable="parksc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~arteriol_scler+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("arteriol_scler","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_arteriol_scler.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/arteriol_scler")

#Indicate name of the variable analyzed
variable="arteriol_scler"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~caa_4gp+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("caa_4gp","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_caa_4gp.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/caa_4gp")

#Indicate name of the variable analyzed
variable="caa_4gp"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cvda_4gp2+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cvda_4gp2","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cvda_4gp2.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cvda_4gp2")

#Indicate name of the variable analyzed
variable="cvda_4gp2"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_ci_num2_gct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gct")

#Indicate name of the variable analyzed
variable="ci_num2_gct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_ci_num2_gtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gtt")

#Indicate name of the variable analyzed
variable="ci_num2_gtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_ci_num2_mct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mct")

#Indicate name of the variable analyzed
variable="ci_num2_mct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_ci_num2_mtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mtt")

#Indicate name of the variable analyzed
variable="ci_num2_mtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~diabetes_sr_rx_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("diabetes_sr_rx_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_diabetes_sr_rx_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/diabetes_sr_rx_bl")

#Indicate name of the variable analyzed
variable="diabetes_sr_rx_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cancer_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cancer_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_cancer_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cancer_bl")

#Indicate name of the variable analyzed
variable="cancer_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~stroke_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("stroke_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_stroke_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/stroke_bl")

#Indicate name of the variable analyzed
variable="stroke_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~headinjrloc_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("headinjrloc_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_headinjrloc_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/headinjrloc_bl")

#Indicate name of the variable analyzed
variable="headinjrloc_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~heart_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("heart_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_heart_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/heart_bl")

#Indicate name of the variable analyzed
variable="heart_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~hypertension_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("hypertension_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_hypertension_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/hypertension_bl")

#Indicate name of the variable analyzed
variable="hypertension_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tdp_st4+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("tdp_st4","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_tdp_st4.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tdp_st4")

#Indicate name of the variable analyzed
variable="tdp_st4"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~msex+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("msex","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Vasc_msex.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/msex")

#Indicate name of the variable analyzed
variable="msex"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################



#######################################################################
#######################################################################
#######################################################################

# MCT

setwd("/mnt/f/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Ast_integrated_batch_5000_leiden_predicted_final.rds")

#Add official metadata

setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("/mnt/f/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

DefaultAssay(pbmc)="RNA"
pbmc[["integrated"]] <- NULL
pbmc[["SCT"]] <- NULL
pbmc[["prediction.score.celltype"]] <- NULL

Idents(pbmc)<-"cell_type_high_resolution"

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
formula = ~cogn_global_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cogn_global_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_global_lv")

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
formula = ~cogn_ep_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ep_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cogn_ep_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ep_lv")

#Indicate name of the variable analyzed
variable="cogn_ep_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_po_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_po_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cogn_po_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_po_lv")

#Indicate name of the variable analyzed
variable="cogn_po_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_ps_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ps_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cogn_ps_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ps_lv")

#Indicate name of the variable analyzed
variable="cogn_ps_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_se_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_se_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cogn_se_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_se_lv")

#Indicate name of the variable analyzed
variable="cogn_se_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_wo_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_wo_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cogn_wo_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_wo_lv")

#Indicate name of the variable analyzed
variable="cogn_wo_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognwo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognwo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cognwo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognwo_random_slope")

#Indicate name of the variable analyzed
variable="cognwo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognse_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognse_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cognse_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognse_random_slope")

#Indicate name of the variable analyzed
variable="cognse_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognpo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognpo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cognpo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognpo_random_slope")

#Indicate name of the variable analyzed
variable="cognpo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognps_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognps_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cognps_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognps_random_slope")

#Indicate name of the variable analyzed
variable="cognps_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognep_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognep_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cognep_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognep_random_slope")

#Indicate name of the variable analyzed
variable="cognep_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogng_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cogng_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogng_random_slope")

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
formula = ~plaq_n_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_plaq_n_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_nft_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CR_score")

#Indicate name of the variable analyzed
variable="nft_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_gpath_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CR_score")

#Indicate name of the variable analyzed
variable="gpath_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~amyloid+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("amyloid","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_amyloid.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/amyloid")

#Indicate name of the variable analyzed
variable="amyloid"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_gpath.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

#Indicate name of the variable analyzed
variable="gpath"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_3neocort+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_3neocort","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_gpath_3neocort.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_3neocort")

#Indicate name of the variable analyzed
variable="gpath_3neocort"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_nft.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft")

#Indicate name of the variable analyzed
variable="nft"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_nft_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_mf")

#Indicate name of the variable analyzed
variable="nft_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tangles+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_tangles.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles")

#Indicate name of the variable analyzed
variable="tangles"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_plaq_n.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n")

#Indicate name of the variable analyzed
variable="plaq_n"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_plaq_n_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_mf")

#Indicate name of the variable analyzed
variable="plaq_n_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_plaq_d.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d")

#Indicate name of the variable analyzed
variable="plaq_d"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_plaq_d_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d_mf")

#Indicate name of the variable analyzed
variable="plaq_d_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~bradysc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("bradysc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_bradysc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/bradysc_lv")

#Indicate name of the variable analyzed
variable="bradysc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gaitsc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("gaitsc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_gaitsc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gaitsc_lv")

#Indicate name of the variable analyzed
variable="gaitsc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~parksc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("parksc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_parksc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/parksc_lv")

#Indicate name of the variable analyzed
variable="parksc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~arteriol_scler+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("arteriol_scler","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_arteriol_scler.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/arteriol_scler")

#Indicate name of the variable analyzed
variable="arteriol_scler"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~caa_4gp+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("caa_4gp","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_caa_4gp.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/caa_4gp")

#Indicate name of the variable analyzed
variable="caa_4gp"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cvda_4gp2+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cvda_4gp2","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cvda_4gp2.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cvda_4gp2")

#Indicate name of the variable analyzed
variable="cvda_4gp2"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_ci_num2_gct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gct")

#Indicate name of the variable analyzed
variable="ci_num2_gct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_ci_num2_gtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gtt")

#Indicate name of the variable analyzed
variable="ci_num2_gtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_ci_num2_mct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mct")

#Indicate name of the variable analyzed
variable="ci_num2_mct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_ci_num2_mtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mtt")

#Indicate name of the variable analyzed
variable="ci_num2_mtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~diabetes_sr_rx_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("diabetes_sr_rx_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_diabetes_sr_rx_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/diabetes_sr_rx_bl")

#Indicate name of the variable analyzed
variable="diabetes_sr_rx_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cancer_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cancer_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_cancer_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cancer_bl")

#Indicate name of the variable analyzed
variable="cancer_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~stroke_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("stroke_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_stroke_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/stroke_bl")

#Indicate name of the variable analyzed
variable="stroke_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~headinjrloc_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("headinjrloc_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_headinjrloc_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/headinjrloc_bl")

#Indicate name of the variable analyzed
variable="headinjrloc_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~heart_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("heart_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_heart_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/heart_bl")

#Indicate name of the variable analyzed
variable="heart_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~hypertension_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("hypertension_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_hypertension_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/hypertension_bl")

#Indicate name of the variable analyzed
variable="hypertension_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tdp_st4+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("tdp_st4","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_tdp_st4.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tdp_st4")

#Indicate name of the variable analyzed
variable="tdp_st4"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~msex+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("msex","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Ast_MCT_msex.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/msex")

#Indicate name of the variable analyzed
variable="msex"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################



#######################################################################
#######################################################################
#######################################################################

setwd("/mnt/f/PFC_429_final_Spring2021/Data")
pbmc=readRDS("Immune_integrated_batch_5000_predicted_final.rds")

#Add official metadata

setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Add CR score variables
setwd("/mnt/f/PFC_429_final_Spring2021/Jason")
CRscores=read.csv("CR_scores.csv")
df=pbmc@meta.data
df2<-df[,c('projid'),drop=FALSE]
New_metadata=left_join(df2,CRscores,by='projid')
rownames(New_metadata)=rownames(df2)

pbmc=AddMetaData(pbmc,New_metadata)

DefaultAssay(pbmc)="RNA"
pbmc[["integrated"]] <- NULL
pbmc[["SCT"]] <- NULL
pbmc[["prediction.score.celltype"]] <- NULL

Idents(pbmc)<-"cell_type_high_resolution"

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
formula = ~cogn_global_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_global_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cogn_global_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_global_lv")

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
formula = ~cogn_ep_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ep_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cogn_ep_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ep_lv")

#Indicate name of the variable analyzed
variable="cogn_ep_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_po_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_po_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cogn_po_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_po_lv")

#Indicate name of the variable analyzed
variable="cogn_po_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_ps_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_ps_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cogn_ps_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_ps_lv")

#Indicate name of the variable analyzed
variable="cogn_ps_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_se_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_se_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cogn_se_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_se_lv")

#Indicate name of the variable analyzed
variable="cogn_se_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogn_wo_lv+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogn_wo_lv","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cogn_wo_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogn_wo_lv")

#Indicate name of the variable analyzed
variable="cogn_wo_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognwo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognwo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cognwo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognwo_random_slope")

#Indicate name of the variable analyzed
variable="cognwo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognse_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognse_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cognse_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognse_random_slope")

#Indicate name of the variable analyzed
variable="cognse_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognpo_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognpo_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cognpo_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognpo_random_slope")

#Indicate name of the variable analyzed
variable="cognpo_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognps_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognps_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cognps_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognps_random_slope")

#Indicate name of the variable analyzed
variable="cognps_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cognep_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cognep_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cognep_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cognep_random_slope")

#Indicate name of the variable analyzed
variable="cognep_random_slope"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cogng_random_slope+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("cogng_random_slope","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cogng_random_slope.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cogng_random_slope")

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
formula = ~plaq_n_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_plaq_n_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_CR_score")

#Indicate name of the variable analyzed
variable="plaq_n_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_nft_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_CR_score")

#Indicate name of the variable analyzed
variable="nft_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_CR_score+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_CR_score","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_gpath_CR_score.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_CR_score")

#Indicate name of the variable analyzed
variable="gpath_CR_score"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~amyloid+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("amyloid","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_amyloid.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/amyloid")

#Indicate name of the variable analyzed
variable="amyloid"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_gpath.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath")

#Indicate name of the variable analyzed
variable="gpath"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gpath_3neocort+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("gpath_3neocort","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_gpath_3neocort.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gpath_3neocort")

#Indicate name of the variable analyzed
variable="gpath_3neocort"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_nft.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft")

#Indicate name of the variable analyzed
variable="nft"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~nft_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("nft_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_nft_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/nft_mf")

#Indicate name of the variable analyzed
variable="nft_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tangles+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("tangles","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_tangles.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tangles")

#Indicate name of the variable analyzed
variable="tangles"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_plaq_n.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n")

#Indicate name of the variable analyzed
variable="plaq_n"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_n_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_n_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_plaq_n_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_n_mf")

#Indicate name of the variable analyzed
variable="plaq_n_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_plaq_d.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d")

#Indicate name of the variable analyzed
variable="plaq_d"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~plaq_d_mf+pmi+age_death
cd = as.data.frame(colData(pb))
cd2=cd[,c("plaq_d_mf","pmi","age_death")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:4))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_plaq_d_mf.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/plaq_d_mf")

#Indicate name of the variable analyzed
variable="plaq_d_mf"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~bradysc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("bradysc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_bradysc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/bradysc_lv")

#Indicate name of the variable analyzed
variable="bradysc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~gaitsc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("gaitsc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_gaitsc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/gaitsc_lv")

#Indicate name of the variable analyzed
variable="gaitsc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~parksc_lv+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("parksc_lv","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_parksc_lv.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/parksc_lv")

#Indicate name of the variable analyzed
variable="parksc_lv"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~arteriol_scler+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("arteriol_scler","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_arteriol_scler.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/arteriol_scler")

#Indicate name of the variable analyzed
variable="arteriol_scler"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~caa_4gp+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("caa_4gp","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_caa_4gp.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/caa_4gp")

#Indicate name of the variable analyzed
variable="caa_4gp"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cvda_4gp2+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cvda_4gp2","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cvda_4gp2.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cvda_4gp2")

#Indicate name of the variable analyzed
variable="cvda_4gp2"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_ci_num2_gct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gct")

#Indicate name of the variable analyzed
variable="ci_num2_gct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_gtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_gtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_ci_num2_gtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_gtt")

#Indicate name of the variable analyzed
variable="ci_num2_gtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mct+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mct","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_ci_num2_mct.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mct")

#Indicate name of the variable analyzed
variable="ci_num2_mct"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~ci_num2_mtt+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("ci_num2_mtt","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_ci_num2_mtt.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/ci_num2_mtt")

#Indicate name of the variable analyzed
variable="ci_num2_mtt"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~diabetes_sr_rx_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("diabetes_sr_rx_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_diabetes_sr_rx_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/diabetes_sr_rx_bl")

#Indicate name of the variable analyzed
variable="diabetes_sr_rx_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~cancer_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("cancer_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_cancer_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/cancer_bl")

#Indicate name of the variable analyzed
variable="cancer_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~stroke_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("stroke_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_stroke_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/stroke_bl")

#Indicate name of the variable analyzed
variable="stroke_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~headinjrloc_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("headinjrloc_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_headinjrloc_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/headinjrloc_bl")

#Indicate name of the variable analyzed
variable="headinjrloc_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~heart_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("heart_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_heart_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/heart_bl")

#Indicate name of the variable analyzed
variable="heart_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~hypertension_bl+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("hypertension_bl","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_hypertension_bl.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/hypertension_bl")

#Indicate name of the variable analyzed
variable="hypertension_bl"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~tdp_st4+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("tdp_st4","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_tdp_st4.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/tdp_st4")

#Indicate name of the variable analyzed
variable="tdp_st4"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

##############################################
formula = ~msex+pmi+age_death+gpath
cd = as.data.frame(colData(pb))
cd2=cd[,c("msex","pmi","age_death","gpath")]
design=model.matrix(formula,cd2)

pb2=pb[,rownames(design)]

res <- pbDS(pb2,design = design,coef=as.list(2:5))

save.image("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/Mic_MCT_msex.RData")

#Define output directory
setwd("/mnt/f/PFC_429_final_Spring2021/Vanshika/muscat/Results/msex")

#Indicate name of the variable analyzed
variable="msex"

tbl <- res$table[[1]]
names(tbl)

for (i in 1:length(tbl)){
  name=names(tbl)[i]
  k1 <- tbl[[i]]
  write.csv(k1,file=paste0(name,"_",variable,".csv"))
}
##############################################

