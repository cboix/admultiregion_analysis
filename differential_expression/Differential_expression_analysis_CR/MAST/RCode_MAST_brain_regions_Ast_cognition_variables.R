suppressPackageStartupMessages({
  library(ggplot2)
  library(GGally)
  #library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  #library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
})

library(Seurat)
library(dplyr)

#options(mc.cores = detectCores() - 1)
options(mc.cores = 3)
knitr::opts_chunk$set(message = FALSE,error = FALSE,warning = FALSE,cache = FALSE,fig.width=8,fig.height=6)

#freq_expressed <- 0.2
FCTHRESHOLD <- 0.02

BR_list=c("AG","MT","EC","HC","TH")

for(i in 1:length(BR_list)){
  print(BR_list[i])
  
  setwd("G:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
  pbmc = readRDS(paste0(BR_list[i],"_merge_major_cell_type.rds"))
  
  Idents(pbmc) <- "cell_type_high_resolution"
  pbmc=subset(pbmc,idents=c("Ast DPP10","Ast GRM3","Ast LUZP2","Ast DCLK1"))
  
  DefaultAssay(pbmc)="RNA"
  #pbmc[["integrated"]] <- NULL
  #pbmc[["SCT"]] <- NULL
  #pbmc[["prediction.score.celltype"]] <- NULL
  #pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  
  
  #Add official metadata
  
  setwd("G:/PFC_429_final_Spring2021/Vanshika")
  metadata=read.csv("dataset_652_basic_04-23-2020.csv")
  
  ##################################
  #Script modification HM for cognition variables
  PFC427_projids=read.csv("PFC427_projids.csv")
  metadata=left_join(PFC427_projids,metadata,by="projid")
  
  q25=quantile(metadata[,"cogn_global_lv"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cogn_global_lv"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cogn_global_lv"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cogn_global_lv_median =
                                   case_when(cogn_global_lv <= q50 ~ "low", 
                                             cogn_global_lv > q50 ~ "high"))
  
  q25=quantile(metadata[,"cogn_ep_lv"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cogn_ep_lv"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cogn_ep_lv"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cogn_ep_lv_median =
                                   case_when(cogn_ep_lv <= q50 ~ "low", 
                                             cogn_ep_lv > q50 ~ "high"))
  
  q25=quantile(metadata[,"cogn_se_lv"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cogn_se_lv"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cogn_se_lv"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cogn_se_lv_median =
                                   case_when(cogn_se_lv <= q50 ~ "low", 
                                             cogn_se_lv > q50 ~ "high"))
  
  q25=quantile(metadata[,"cogn_po_lv"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cogn_po_lv"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cogn_po_lv"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cogn_po_lv_median =
                                   case_when(cogn_po_lv <= q50 ~ "low", 
                                             cogn_po_lv > q50 ~ "high"))
  
  q25=quantile(metadata[,"cogn_ps_lv"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cogn_ps_lv"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cogn_ps_lv"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cogn_ps_lv_median =
                                   case_when(cogn_ps_lv <= q50 ~ "low", 
                                             cogn_ps_lv > q50 ~ "high"))
  
  q25=quantile(metadata[,"cogn_wo_lv"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cogn_wo_lv"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cogn_wo_lv"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cogn_wo_lv_median =
                                   case_when(cogn_wo_lv <= q50 ~ "low", 
                                             cogn_wo_lv > q50 ~ "high"))
  
  q25=quantile(metadata[,"cogng_random_slope"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cogng_random_slope"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cogng_random_slope"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cogng_random_slope_median =
                                   case_when(cogng_random_slope <= q50 ~ "low", 
                                             cogng_random_slope > q50 ~ "high"))
  
  q25=quantile(metadata[,"cognep_random_slope"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cognep_random_slope"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cognep_random_slope"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cognep_random_slope_median =
                                   case_when(cognep_random_slope <= q50 ~ "low", 
                                             cognep_random_slope > q50 ~ "high"))
  
  q25=quantile(metadata[,"cognse_random_slope"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cognse_random_slope"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cognse_random_slope"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cognse_random_slope_median =
                                   case_when(cognse_random_slope <= q50 ~ "low", 
                                             cognse_random_slope > q50 ~ "high"))
  
  q25=quantile(metadata[,"cognpo_random_slope"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cognpo_random_slope"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cognpo_random_slope"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cognpo_random_slope_median =
                                   case_when(cognpo_random_slope <= q50 ~ "low", 
                                             cognpo_random_slope > q50 ~ "high"))
  
  q25=quantile(metadata[,"cognps_random_slope"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cognps_random_slope"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cognps_random_slope"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cognps_random_slope_median =
                                   case_when(cognps_random_slope <= q50 ~ "low", 
                                             cognps_random_slope > q50 ~ "high"))
  
  q25=quantile(metadata[,"cognwo_random_slope"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(metadata[,"cognwo_random_slope"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(metadata[,"cognwo_random_slope"], prob=c(0.75),na.rm=TRUE)
  metadata = metadata %>% mutate(cognwo_random_slope_median =
                                   case_when(cognwo_random_slope <= q50 ~ "low", 
                                             cognwo_random_slope > q50 ~ "high"))
  
  
  ##################################
  
  data=pbmc@meta.data
  data=data[,"projid",drop=FALSE]
  data2=left_join(data,metadata,by="projid")
  row.names(data2)=row.names(data)
  pbmc=AddMetaData(pbmc,data2)
  
  
  #Add CR score variables
  setwd("G:/PFC_429_final_Spring2021/Vanshika")
  CRscores=read.csv("CR_scores.csv")
  CRscores=na.omit(CRscores)
  
  #Add tangles CR score variables
  setwd("G:/PFC_429_final_Spring2021/Vanshika")
  tanglesCRscores=read.csv("PFC427_tangles_CR_score.csv")
  tanglesCRscores=na.omit(tanglesCRscores)
  
  #Add CDR score variables
  setwd("G:/PFC_429_final_Spring2021/Vanshika")
  CDRscores=read.csv("PFC427_CDR_scores.csv")
  CDRscores=na.omit(CDRscores)
  
  CRscores=left_join(CRscores,tanglesCRscores,by="projid")
  CRscores=left_join(CRscores,CDRscores,by="projid")
  
  q25=quantile(CRscores[,"gpath_CR_score"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(CRscores[,"gpath_CR_score"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(CRscores[,"gpath_CR_score"], prob=c(0.75),na.rm=TRUE)
  CRscores = CRscores %>% mutate(gpath_CR_score_median =
                                     case_when(gpath_CR_score <= q50 ~ "low", 
                                               gpath_CR_score > q50 ~ "high"))
  
  q25=quantile(CRscores[,"plaq_n_CR_score"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(CRscores[,"plaq_n_CR_score"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(CRscores[,"plaq_n_CR_score"], prob=c(0.75),na.rm=TRUE)
  CRscores = CRscores %>% mutate(plaq_n_CR_score_median =
                                     case_when(plaq_n_CR_score <= q50 ~ "low", 
                                               plaq_n_CR_score > q50 ~ "high"))
  
  q25=quantile(CRscores[,"nft_CR_score"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(CRscores[,"nft_CR_score"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(CRscores[,"nft_CR_score"], prob=c(0.75),na.rm=TRUE)
  CRscores = CRscores %>% mutate(nft_CR_score_median =
                                     case_when(nft_CR_score <= q50 ~ "low", 
                                               nft_CR_score > q50 ~ "high"))
  
  q25=quantile(CRscores[,"tangles_CR_score"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(CRscores[,"tangles_CR_score"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(CRscores[,"tangles_CR_score"], prob=c(0.75),na.rm=TRUE)
  CRscores = CRscores %>% mutate(tangles_CR_score_median =
                                     case_when(tangles_CR_score <= q50 ~ "low", 
                                               tangles_CR_score > q50 ~ "high"))
  
  q25=quantile(CRscores[,"gpath_CDR_score"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(CRscores[,"gpath_CDR_score"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(CRscores[,"gpath_CDR_score"], prob=c(0.75),na.rm=TRUE)
  CRscores = CRscores %>% mutate(gpath_CDR_score_median =
                                     case_when(gpath_CDR_score <= q50 ~ "low", 
                                               gpath_CDR_score > q50 ~ "high"))
  
  q25=quantile(CRscores[,"plaq_n_CDR_score"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(CRscores[,"plaq_n_CDR_score"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(CRscores[,"plaq_n_CDR_score"], prob=c(0.75),na.rm=TRUE)
  CRscores = CRscores %>% mutate(plaq_n_CDR_score_median =
                                     case_when(plaq_n_CDR_score <= q50 ~ "low", 
                                               plaq_n_CDR_score > q50 ~ "high"))
  
  q25=quantile(CRscores[,"nft_CDR_score"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(CRscores[,"nft_CDR_score"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(CRscores[,"nft_CDR_score"], prob=c(0.75),na.rm=TRUE)
  CRscores = CRscores %>% mutate(nft_CDR_score_median =
                                     case_when(nft_CDR_score <= q50 ~ "low", 
                                               nft_CDR_score > q50 ~ "high"))
  
  q25=quantile(CRscores[,"tangles_CDR_score"], prob=c(0.25),na.rm=TRUE)
  q50=quantile(CRscores[,"tangles_CDR_score"], prob=c(0.5),na.rm=TRUE)
  q75=quantile(CRscores[,"tangles_CDR_score"], prob=c(0.75),na.rm=TRUE)
  CRscores = CRscores %>% mutate(tangles_CDR_score_median =
                                     case_when(tangles_CDR_score <= q50 ~ "low", 
                                               tangles_CDR_score > q50 ~ "high"))
  
  data=pbmc@meta.data
  data=data[,"projid",drop=FALSE]
  data2=left_join(data,CRscores,by="projid")
  row.names(data2)=row.names(data)
  #data2=data2[,c("projid","gpath_CDR_score","gpath_CDR_score_median")]
  pbmc=AddMetaData(pbmc,data2)
  
  
  #Add major cell type info
  
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
  
  #########################################################
  
  #Downsample 50000 cells
  pbmc=subset(pbmc,downsample=50000)
  
  pbmc.sce <- as.SingleCellExperiment(pbmc)
  rm(pbmc)
  
  # remove undetected genes
  pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
  dim(pbmc.sce)
  
  # remove lowly expressed genes
  pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
  dim(pbmc.sce)
  
  #########################################################
  #########################################################
  
  ## Variables
  
  #This part makes sure that individuals with NA values are excluded
  formula = ~cogn_global_lv_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cogn_global_lv_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cogn_global_lv_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cogn_global_lv.csv")
  
  #########################################################
  #This part makes sure that individuals with NA values are excluded
  formula = ~cogn_ep_lv_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cogn_ep_lv_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cogn_ep_lv_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cogn_ep_lv.csv")
  
  #########################################################
  #This part makes sure that individuals with NA values are excluded
  formula = ~cogn_se_lv_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cogn_se_lv_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cogn_se_lv_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cogn_se_lv.csv")
  
  #########################################################
  #This part makes sure that individuals with NA values are excluded
  formula = ~cogn_po_lv_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cogn_po_lv_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cogn_po_lv_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cogn_po_lv.csv")
  
  #########################################################
  #This part makes sure that individuals with NA values are excluded
  formula = ~cogn_ps_lv_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cogn_ps_lv_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cogn_ps_lv_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cogn_ps_lv.csv")
  
  #########################################################
  #This part makes sure that individuals with NA values are excluded
  formula = ~cogn_wo_lv_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cogn_wo_lv_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cogn_wo_lv_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cogn_wo_lv.csv")
  
  #########################################################
  #This part makes sure that individuals with NA values are excluded
  formula = ~cogng_random_slope_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cogng_random_slope_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cogng_random_slope_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cogng_random_slope.csv")
  
  #########################################################
  #This part makes sure that individuals with NA values are excluded
  formula = ~cognep_random_slope_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cognep_random_slope_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cognep_random_slope_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cognep_random_slope.csv")
  
  #########################################################
  #This part makes sure that individuals with NA values are excluded
  formula = ~cognse_random_slope_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cognse_random_slope_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cognse_random_slope_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cognse_random_slope.csv")
  
  #########################################################
  #This part makes sure that individuals with NA values are excluded
  formula = ~cognpo_random_slope_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cognpo_random_slope_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cognpo_random_slope_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cognpo_random_slope.csv")
  
  #########################################################
  #This part makes sure that individuals with NA values are excluded
  formula = ~cognps_random_slope_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cognps_random_slope_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cognps_random_slope_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cognps_random_slope.csv")
  
  #########################################################
  #This part makes sure that individuals with NA values are excluded
  formula = ~cognwo_random_slope_median + age_death + pmi + msex
  cd = as.data.frame(colData(pbmc.sce))
  cd2=cd[,c("cognwo_random_slope_median","pmi","age_death","msex")]
  design=model.matrix(formula,cd2)
  
  print(colnames(design)[2])
  
  pbmc.sce_sub=pbmc.sce[,rownames(design)]
  
  #SingleCellExperiment to MAST 
  sca = SceToSingleCellAssay(pbmc.sce_sub)
  
  ## Remove invariant genes
  #sca <- sca[sample(which(freq(sca)>0), 6000),]
  
  #Recalculating the cellular detection rate (ngeneson)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  cond<-factor(colData(sca)$cognwo_random_slope_median)
  cond<-relevel(cond,"low")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson + age_death + pmi + msex + (1|projid), sca, exprs_value = 'logcounts')
  
  summaryCond <- summary(zlmCond, doLRT='conditionhigh') 
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionhigh' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionhigh' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
  fcHurdle_All <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle_All, fdr)
  
  setwd(paste0("G:/PFC_429_final_Spring2021/Vanshika/MAST/Results_Brain_Regions/",BR_list[i]))
  write.csv(fcHurdle_All,"Ast_cognwo_random_slope.csv")
  
  #########################################################
  #########################################################
  
}