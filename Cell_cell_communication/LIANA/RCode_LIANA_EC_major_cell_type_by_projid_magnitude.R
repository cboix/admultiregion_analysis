library(Seurat)
library(dplyr)
library(tidyverse)
library(magrittr)
library(liana)

setwd("F:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc=readRDS("EC_merge_major_cell_type.rds")

Idents(pbmc)<-"cell_type_high_resolution"
pbmc=subset(pbmc,idents=c("Ast DPP10","Ast GRM3","Ast DCLK1","Ast LUZP2","Exc RELN COL5A2","Exc DLC1 SNTG2",
                          "Exc TOX3 TTC6","Exc RELN GPC5","Exc TOX3 INO80D","Exc AGBL1 GPC5","Exc SOX11 NCKAP5",
                          "Exc TOX3 POSTN","Exc COL25A1 SEMA3D","Mic P2RY12","Mic TPT1","CAMs","T cells","OPC GPC5","OPC DOCK5","OPC CEP112"                  
                          ,"Inh GPC5 RIT2","Inh LAMP5 NRG1 (Rosehip)","Inh SGCD PDE3A",              
                          "Inh VIP CLSTN2","Inh RYR3 TSHZ2","Inh SORCS1 TTN",             
                          "Inh L1-6 LAMP5 CA13","Inh CUX2 MSR1","Inh VIP ABI3BP",              
                          "Inh FBN2 EPB41L4A","Inh PVALB HTR4","Inh L3-5 SST MAFB",           
                          "Inh PAX6 RELN","Inh PVALB SULF1","Inh VIP TSHZ2",               
                          "Inh LAMP5 RELN","Inh ALCAM TRPM3","Inh PVALB CA8 (Chandelier)",  
                          "Inh ENOX2 SPHKAP","Inh L6 SST NPY","Inh VIP THSD7B",             
                          "Inh PTPRK FAM19A1","Oli OPALIN","Oli RASGRF1","Per","End","Fib","SMC"                       
                          
))

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

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

seurat_objects=SplitObject(pbmc, split.by = "projid")
rm(pbmc)

for (i in 1:length(seurat_objects)){
  projid=names(seurat_objects)[i]
  pbmc=seurat_objects[[i]]
  
  liana_test <- liana_wrap(pbmc, idents_col = 'major_cell_type')
  
  liana_test <- liana_test %>%
    liana_aggregate(aggregate_how = "magnitude")
  
  output_file=paste0("projid_",projid,".csv")
  setwd("F:/Brain_region_project/LIANA/Results_aggregate_rank_magnitude/EC")
  write.csv(liana_test,file=output_file)
  
}