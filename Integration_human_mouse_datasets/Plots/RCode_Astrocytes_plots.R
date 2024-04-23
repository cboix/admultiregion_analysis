library(ggplot2)
library(Seurat)
setwd("D:/Mouse_data_integration")
pbmc=readRDS(file="Astrocytes_integrated_all.rds")

#write.csv(pbmc@meta.data,file="Astrocytes_integrated_all_metadata.csv")

data=read.csv("Astrocytes_integrated_all_metadata_addition.csv",row.names=1)

pbmc=AddMetaData(pbmc,data)

#saveRDS(pbmc,file="Astrocytes_integrated_all_study.rds")

library(RColorBrewer)
# Define the number of colors you want
#nb.cols <- 23
#mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(nb.cols)

mycolors="Paired"
#levels(pbmc) <- c("Inh L1-6 LAMP5 CA13","Inh LAMP5 NRG1 (Rosehip)","Inh LAMP5 RELN","Inh PAX6 RELN","Inh VIP ABI3BP","Inh SORCS1 TTN","Inh VIP TSHZ2","Inh RYR3 TSHZ2","Inh VIP THSD7B","Inh SGCD PDE3A","Inh ALCAM TRPM3","Inh PTPRK FAM19A1","Inh VIP CLSTN2","Inh GPC5 RIT2","Inh PVALB SULF1","Inh PVALB HTR4","Inh CUX2 MSR1","Inh PVALB CA8 (Chandelier)","Inh FBN2 EPB41L4A","Inh ENOX2 SPHKAP","Inh L3-5 SST MAFB","Inh L6 SST NPY","Inh MEIS2 FOXP2")
p1=DimPlot(pbmc,cols=mycolors,label=TRUE,raster=FALSE)+NoAxes()
p1

p1=DimPlot(pbmc,cols=mycolors,label=FALSE,raster=FALSE)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1


p1=DimPlot(pbmc,cols=mycolors,label=TRUE,raster=FALSE,group.by="species")+NoAxes()
p1

p1=DimPlot(pbmc,cols=mycolors,label=FALSE,raster=FALSE,group.by="species")+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

p1=DimPlot(pbmc,cols=mycolors,label=TRUE,raster=FALSE,group.by="study")+NoAxes()
p1

p1=DimPlot(pbmc,cols=mycolors,label=FALSE,raster=FALSE,group.by="study")+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

library(dplyr)
df=pbmc@meta.data
df2=df[,c("projid","Region")]

df3 = df2 %>% mutate(Region2 = case_when(grepl("Hypothalamus,Thalamus,Midbrain dorsal,Midbrain ventral,Pons,Medulla,Spinal cord",Region)~"Hypothalamus,Thalamus,Midbrain dorsal,Midbrain ventral,Pons,Medulla,Spinal cord",
                                         grepl("Dentate gyrus",Region)~"Dentate gyrus",
                                         grepl("Telencephalon",Region)~"Telencephalon",
                                         is.na(df2$Region)~"not assigned"))
df4=df3[,-1]
df5=df4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,df5)

order=c("Telencephalon","Hypothalamus,Thalamus,Midbrain dorsal,Midbrain ventral,Pons,Medulla,Spinal cord","Dentate gyrus","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#A6CEE3","#1F78B4","#B2DF8A"),label=TRUE,raster=FALSE,group.by="Region2",order=order)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#A6CEE3","#1F78B4","#B2DF8A"),label=FALSE,raster=FALSE,group.by="Region2",order=order)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

df=pbmc@meta.data
df2=df[,c("projid","Tissue")]

df3 = df2 %>% mutate(Tissue2 = case_when(grepl("Amygd",Tissue)~"Amygd",
                                         grepl("CA1",Tissue)~"HC",
                                         grepl("CB",Tissue)~"CB",
                                         grepl("Ctx1",Tissue)~"Cortex",
                                         grepl("Ctx1.5",Tissue)~"Cortex",
                                         grepl("Ctx2",Tissue)~"Cortex",
                                         grepl("Ctx3",Tissue)~"Cortex",
                                         grepl("DentGyr",Tissue)~"HC",
                                         grepl("HC",Tissue)~"HC",
                                         grepl("Hypoth",Tissue)~"Hypoth",
                                         grepl("MBd",Tissue)~"MB",
                                         grepl("MBv",Tissue)~"MB",
                                         grepl("Medulla",Tissue)~"Pons/Medulla",
                                         grepl("OB",Tissue)~"OB",
                                         grepl("Pons",Tissue)~"Pons/Medulla",
                                         grepl("SC",Tissue)~"SC",
                                         grepl("SScortex",Tissue)~"Cortex",
                                         grepl("StriatDor",Tissue)~"Striatum",
                                         grepl("StriatVent",Tissue)~"Striatum",
                                         grepl("Thal",Tissue)~"Thal",
                                         is.na(df2$Tissue)~"not assigned"))
df4=df3[,-1]
df5=df4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,df5)

library(RColorBrewer)
#Define the number of colors you want
#nb.cols <- 23
#mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(nb.cols)
mycolors="Paired"

pbmc$Tissue2=factor(x=pbmc$Tissue2,levels=c("not assigned","Cortex","HC","Amygd","OB","Striatum","CB","Thal","Hypoth","MB","Pons/Medulla","SC"))

p1=DimPlot(pbmc,cols=mycolors,label=TRUE,raster=FALSE,split.by="Tissue2",na.value = "#EAEAEA",ncol=4)+NoAxes()
p1

p1=DimPlot(pbmc,cols=mycolors,label=FALSE,raster=FALSE,split.by="Tissue2",na.value = "#EAEAEA",ncol=4)+NoAxes()+NoLegend()+theme(strip.text=element_text(size=7,color="white"))
p1=AugmentPlot(plot=p1)
p1


df=pbmc@meta.data
df2=df[,c("projid","tissue")]

df3 = df2 %>% mutate(tissue2 = case_when(grepl("FC",tissue)~"Cortex",
                                         grepl("PC",tissue)~"Cortex",
                                         grepl("HC",tissue)~"Cortex",
                                         grepl("TH",tissue)~"TH",
                                         is.na(df2$tissue)~"not assigned"))
df4=df3[,-1]
df5=df4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,df5)

order=c("TH","Cortex","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#1F78B4","#FF7F00"),label=TRUE,raster=FALSE,group.by="tissue2",order=order)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#1F78B4","#FF7F00"),label=FALSE,raster=FALSE,group.by="tissue2",order=order)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

