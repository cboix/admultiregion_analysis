library(Seurat)
library(dplyr)
library(tidyverse)
library(magrittr)
library(liana)

setwd("F:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc=readRDS("PFC_merge_major_cell_type.rds")

Idents(pbmc)<-"cell_type_high_resolution"
pbmc=subset(pbmc,idents=c("Ast DPP10","Ast GRM3","Ast DCLK1","Exc L2-3 CBLN2 LINC02306","Exc L4-5 RORB IL1RAPL2"      
                          ,"Exc L3-4 RORB CUX2","Exc L4-5 RORB GABRG1","Exc L5/6 NP"                 
                          ,"Exc L6 THEMIS NFIA","Exc NRGN","Exc L6 CT"                   
                          ,"Exc L6b","Exc L3-5 RORB PLCH1","Exc L5-6 RORB LINC02196"     
                          ,"Exc L5/6 IT Car3","Exc L5 ET","Mic P2RY12","Mic TPT1",'Mic DUSP1',"CAMs","T cells","OPC GPC5","OPC DOCK5"                  
                          ,"Inh L3-5 SST MAFB","Inh ALCAM TRPM3","Inh CUX2 MSR1","Inh SORCS1 TTN","Inh VIP CLSTN2"              
                          ,"Inh VIP ABI3BP","Inh RYR3 TSHZ2","Inh PTPRK FAM19A1"           
                          ,"Inh VIP TSHZ2","Inh PVALB CA8 (Chandelier)","Inh ENOX2 SPHKAP"            
                          ,"Inh PVALB SULF1","Inh L1-6 LAMP5 CA13","Inh FBN2 EPB41L4A"           
                          ,"Inh PAX6 RELN","Inh VIP THSD7B","Inh PVALB HTR4"              
                          ,"Inh GPC5 RIT2","Inh SGCD PDE3A","Inh LAMP5 NRG1 (Rosehip)"    
                          ,"Inh LAMP5 RELN",'Inh L6 SST NPY',"Oli OPALIN","Oli RASGRF1","Per","End","Fib"                       
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
  setwd("F:/Brain_region_project/LIANA/Results_aggregate_rank_magnitude/PFC")
  write.csv(liana_test,file=output_file)
  
}