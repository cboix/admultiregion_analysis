library(dreamlet)
library(muscat)
library(zenith)
library(scater)
library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(purrr)
library(readxl)

setwd("/New12TB_HD1/F/PFC_429_final_Spring2021/Data/Data_Batches/DecontX_data/Grouped_by_cell_class")
pbmc=readRDS("Inh_decontX.rds")

#Add official metadata
setwd("/New12TB_HD1/F/PFC_429_final_Spring2021/Hansruedi/muscat/files")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

metadata = metadata %>% 
  mutate(AD = case_when(niareagansc == 1 ~ 1, 
                        niareagansc == 2 ~ 1,
                        niareagansc == 3 ~ 0,
                        niareagansc == 4 ~ 0))

data=pbmc@meta.data
data=data[,"projid",drop=FALSE]
data2=left_join(data,metadata,by="projid")
row.names(data2)=row.names(data)
pbmc=AddMetaData(pbmc,data2)

#Subsetting individuals without a pathologic diagnosis of AD
Idents(pbmc) <- 'AD'
pbmc=subset(x = pbmc, idents = c("0"))

DefaultAssay(pbmc)="RNA"
pbmc.sce <- as.SingleCellExperiment(pbmc)

###########################################################33
# remove undetected genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
dim(pbmc.sce)

# remove lowly expressed genes
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
dim(pbmc.sce)
##############################
#Alternative
# only keep singlet cells with sufficient reads
#sce <- sce[rowSums(counts(sce) > 0) > 0, ]
#sce <- sce[,colData(sce)$multiplets == 'singlet']

# compute QC metrics
#qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
#ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
#sce <- sce[, !ol]

###############################################
sce=pbmc.sce
# Create pseudobulk data by specifying cluster_id and sample_id
# Count data for each cell type is then stored in the `assay` field
# assay: entry in assayNames(sce) storing raw counts
# cluster_id: variable in colData(sce) indicating cell clusters
# sample_id: variable in colData(sce) indicating sample id for aggregating cells
pb <- aggregateToPseudoBulk(sce,
                            assay = "counts",     
                            cluster_id = "cell_type_high_resolution",  
                            sample_id = "projid",
                            verbose = FALSE)

# one 'assay' per cell type
assayNames(pb)
###############################
###############################
lst = list( test = c("Inh CUX2 MSR1","Inh L3-5 SST MAFB","Inh LAMP5 RELN","Inh L1 PAX6 CA4","Inh PTPRK FAM19A1"), 
            baseline= c("Inh PVALB CA8 (Chandelier)","Inh RYR3 TSHZ2",            
                        "Inh PVALB HTR4","Inh ENOX2 SPHKAP",         
                        "Inh VIP CLSTN2","Inh GPC5 RIT2",            
                        "Inh VIP ABI3BP","Inh PVALB SULF1",          
                        "Inh ALCAM TRPM3","Inh L5-6 SST TH",                    
                        "Inh LAMP5 NRG1 (Rosehip)","Inh L1-6 LAMP5 CA13","Inh L6 SST NPY",            
                        "Inh VIP TSHZ2","Inh SORCS1 TTN","Inh FBN2 EPB41L4A",        
                        "Inh L1-2 PAX6 SCGN","Inh L5-6 PVALB STON2","Inh VIP THSD7B","Inh SGCD PDE3A"))

fit = dreamletCompareClusters( pb, lst, method="fixed")

res = topTable(fit, coef='compare', number=30000)

setwd("/New12TB_HD1/F/PFC_429_final_Spring2021/Hansruedi/Vulnerable_celltype_markers")
write.csv(res,file="Markers_vulnerable_Inh_decontX.csv")

#################################################################################################
#################################################################################################
