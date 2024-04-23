library(Seurat)
library(dplyr)
library(tidyverse)
library(magrittr)
library(liana)

setwd("F:/Brain_region_project/Brain_region_final_12TB/Data_for_Connectome")
pbmc=readRDS("TH_merge_major_cell_type.rds")

Idents(pbmc)<-"cell_type_high_resolution"
pbmc=subset(pbmc,idents=c("Ast DPP10","Ast GRM3","Ast DCLK1","Ast LUZP2","Exc NXPH1 RNF220","Exc VAT1L ERBB4","Exc NRGN",             
                          "Mic P2RY12","Mic TPT1","Mic DUSP1","CAMs","T cells","OPC GPC5","OPC DOCK5","OPC CEP112"                  
                          ,"Inh LAMP5 NRG1 (Rosehip)","Inh MEIS2 FOXP2","Oli OPALIN","Oli RASGRF1","Per","End","Fib","Epd"                       
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
  setwd("F:/Brain_region_project/LIANA/Results_aggregate_rank_magnitude/TH")
  write.csv(liana_test,file=output_file)
  
}
