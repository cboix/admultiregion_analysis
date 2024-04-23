library(ggplot2)
library(Seurat)
setwd("D:/Mouse_data_integration")
pbmc=readRDS(file="GABAergic_neurons_integrated_all_class.rds")

write.csv(pbmc@meta.data,file="GABAergic_neurons_integrated_all_class_metadata.csv")

data=read.csv("GABAergic_neurons_integrated_all_class_metadata_addition.csv",row.names=1)

pbmc=AddMetaData(pbmc,data)

saveRDS(pbmc,file="GABAergic_neurons_integrated_all_class_study.rds")

#Load data
pbmc=readRDS(file="GABAergic_neurons_integrated_all_class_study.rds")
pbmc <- FindClusters(pbmc, resolution = 0.01)

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(nb.cols)

p1=DimPlot(pbmc,cols="Paired",label=TRUE,raster=FALSE)+NoAxes()
p1

p1=DimPlot(pbmc,cols="Paired",label=FALSE,raster=FALSE)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

mycolors="Paired"
p1=DimPlot(pbmc,cols=mycolors,label=TRUE,raster=FALSE,group.by="species")+NoAxes()
p1

p1=DimPlot(pbmc,cols=mycolors,label=FALSE,raster=FALSE,group.by="species")+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

order=c("Bakken","Linnarsson","DropViz","Yao","human_multi_region")
p1=DimPlot(pbmc,cols=mycolors,label=TRUE,raster=FALSE,group.by="study",order=order)+NoAxes()
p1

p1=DimPlot(pbmc,cols=mycolors,label=FALSE,raster=FALSE,group.by="study",order=order)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

p1=DimPlot(pbmc,cols=mycolors,label=TRUE,raster=FALSE,split.by="study",order=order,ncol=2)+NoAxes()
p1

p1=DimPlot(pbmc,cols=mycolors,label=FALSE,raster=FALSE,split.by="study",order=order,ncol=2)+NoAxes()+NoLegend()+theme(strip.text=element_text(size=7,color="white"))
p1=AugmentPlot(plot=p1)
p1


library(dplyr)
df=pbmc@meta.data
df2=df[,c("projid","brain_region")]

df3 = df2 %>% mutate(human_brain_region = case_when(grepl("Angular_gyrus",brain_region)~"Angular gyrus",
                                               grepl("Anterior_thalamus",brain_region)~"Anterior thalamus",
                                               grepl("Entorhinal_cortex",brain_region)~"Entorhinal cortex",
                                               grepl("Hippocampus",brain_region)~"Hippocampus",
                                               grepl("LGN",brain_region)~"LGN",
                                               grepl("Midtemporal_cortex",brain_region)~"Midtemporal cortex",
                                               grepl("Prefrontal_cortex",brain_region)~"Prefrontal cortex",
                                              is.na(df2$brain_region)~"not assigned"))
df4=df3[,-1]
df5=df4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,df5)

order=c("LGN","Angular gyrus","Anterior thalamus","Entorhinal cortex","Hippocampus","Midtemporal cortex","Prefrontal cortex","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#E31A1C","#FB9A99","#33A02C","#B2DF8A","#1F78B4","#A6CEE3","#FDBF6F"),label=TRUE,raster=FALSE,group.by="human_brain_region",order=order)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#E31A1C","#FB9A99","#33A02C","#B2DF8A","#1F78B4","#A6CEE3","#FDBF6F"),label=FALSE,raster=FALSE,group.by="human_brain_region",order=order)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

df=pbmc@meta.data
df2=df[,c("projid","tissue")]

df3 = df2 %>% mutate(tissue2 = case_when(grepl("FC",tissue)~"Cortex",
                                         grepl("PC",tissue)~"Cortex",
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


df=pbmc@meta.data
df2=df[,c("projid","Region")]

df3 = df2 %>% mutate(Region2 = case_when(grepl("Hippocampus,Cortex",Region)~"Cortex",
                                         grepl("Hippocampus",Region)~"Cortex",
                                         grepl("Thalamus",Region)~"Thalamus",
                                         is.na(df2$Region)~"not assigned"))
df4=df3[,-1]
df5=df4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,df5)

order=c("Thalamus","Cortex","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#1F78B4","#FF7F00"),label=TRUE,raster=FALSE,group.by="Region2",order=order)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#1F78B4","#FF7F00"),label=FALSE,raster=FALSE,group.by="Region2",order=order)+NoAxes()+NoLegend()
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
df2=df[,c("projid","region_label")]

df3 = df2 %>% mutate(region_label2 = case_when(grepl("\\<VISp\\>",region_label)~"VISp",
                                         grepl("\\<VISam;VISpm\\>",region_label)~"VISam;VISpm",
                                         grepl("VISl",region_label)~"VISl",
                                         grepl("\\<VIS\\>",region_label)~"VIS",
                                         grepl("TEa;PERI;ECT",region_label)~"TEa;PERI;ECT",
                                         grepl("SSp",region_label)~"SSp",
                                         grepl("RSP",region_label)~"RSP",
                                         grepl("PTLp",region_label)~"PTLp",
                                         grepl("PL;ILA;ORB",region_label)~"PL;ILA;ORB",
                                         grepl("PAR;POST;PRE;SUB;ProS",region_label)~"PAR;POST;PRE;SUB;ProS",
                                         grepl("MOs;FRP",region_label)~"MOs;FRP",
                                         grepl("MOp",region_label)~"MOp",
                                         grepl("HIP",region_label)~"HIP",
                                         grepl("ENT",region_label)~"ENT",
                                         grepl("AUD",region_label)~"AUD",
                                         grepl("\\<SSs;GU;VISC;AIp\\>",region_label)~"SSs;GU;VISC;AIp",
                                         grepl("\\<AId;AIv;AIp\\>",region_label)~"AId;AIv;AIp",
                                         grepl("\\<AId;AIv\\>",region_label)~"AId;AIv",
                                        
                                         grepl("ACA",region_label)~"ACA",
                                         is.na(df2$region_label)~"not assigned"))
df4=df3[,-1]
df5=df4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,df5)

nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(nb.cols)

order=c("VISp","VISam;VISpm","VISl","VIS","TEa;PERI;ECT","SSp","RSP","PTLp","PL;ILA;ORB","PAR;POST;PRE;SUB;ProS","MOs;FRP","MOp","HIP","ENT","AUD","SSs;GU;VISC;AIp","AId;AIv","AId;AIv;AIp","ACA","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#A6CEE3","#5399C6","#3F8EAA","#99CD91","#79C360","#3E9F32","#B89B74","#F47676","#E52829","#F06C45","#FDB762","#FE901E","#ED8F47","#CCAFCA","#9471B4","#825D99","#DDD399","#E0BE6D","#B15928"),label=TRUE,raster=FALSE,group.by="region_label2",order=order)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#A6CEE3","#5399C6","#3F8EAA","#99CD91","#79C360","#3E9F32","#B89B74","#F47676","#E52829","#F06C45","#FDB762","#FE901E","#ED8F47","#CCAFCA","#9471B4","#825D99","#DDD399","#E0BE6D","#B15928"),label=FALSE,raster=FALSE,group.by="region_label2",order=order)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1


df=pbmc@meta.data
df2=df[,c("projid","subclass_label")]

df3 = df2 %>% mutate(subclass_label2 = case_when(grepl("Lamp5",subclass_label)~"Lamp5",
                                               grepl("Vip",subclass_label)~"Vip",
                                               grepl("Pvalb",subclass_label)~"Pvalb",
                                               grepl("\\<Sst Chodl\\>",subclass_label)~"Sst Chodl",
                                               grepl("\\<Sst\\>",subclass_label)~"Sst",
                                               grepl("Sncg",subclass_label)~"Sncg",
                                               is.na(df2$subclass_label)~"not assigned"))
df4=df3[,-1]
df5=df4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,df5)

order=c("Sst Chodl","Sst","Pvalb","Sncg","Vip","Lamp5","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C"),label=TRUE,raster=FALSE,group.by="subclass_label2",order=order)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C"),label=FALSE,raster=FALSE,group.by="subclass_label2",order=order)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1


df=pbmc@meta.data
df2=df[,c("projid","cluster")]

df3 = df2 %>% mutate(cluster2 = case_when(grepl("LGN Inh NTRK1",cluster)~"LGN Inh NTRK1",
                                                 grepl("LGN Inh CTXN3",cluster)~"LGN Inh CTXN3",
                                                 grepl("LGN Inh LAMP5",cluster)~"LGN Inh LAMP5",
                                                 grepl("LGN Inh TRPC4",cluster)~"LGN Inh TRPC4",
                                                 grepl("1",cluster)~"not assigned",
                                                 grepl("2",cluster)~"not assigned",
                                                 grepl("3",cluster)~"not assigned",
                                                 grepl("4",cluster)~"not assigned",
                                                 grepl("5",cluster)~"not assigned",
                                                 is.na(df2$cluster)~"not assigned"))
df4=df3[,-1]
df5=df4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,df5)

order=c("LGN Inh NTRK1","LGN Inh CTXN3","LGN Inh LAMP5","LGN Inh TRPC4","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99"),label=TRUE,raster=FALSE,group.by="cluster2",order=order,pt.size = 2)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99"),label=FALSE,raster=FALSE,group.by="cluster2",order=order,pt.size = 2)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

saveRDS(pbmc,file="GABAergic_neurons_integrated_all_class_study_final.rds")

#Extract metadata 
library(ggplot2)
library(Seurat)
setwd("D:/Mouse_data_integration")
pbmc=readRDS(file="GABAergic_neurons_integrated_all_class_study_final.rds")

write.csv(pbmc@meta.data,file="GABAergic_neurons_integrated_all_class_study_final_metadata.csv")


#Change some of the plots

library(ggplot2)
library(Seurat)
setwd("D:/Mouse_data_integration")

#Load data
pbmc=readRDS(file="GABAergic_neurons_integrated_all_class_study_final.rds")

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(12,"Paired"))(nb.cols)

p1=DimPlot(pbmc,cols="Paired",label=TRUE,raster=FALSE,group.by="seurat_clusters")+NoAxes()
p1

p1=DimPlot(pbmc,cols="Paired",label=FALSE,raster=FALSE,group.by="seurat_clusters")+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

order=c("Bakken","Linnarsson","DropViz","Yao","human_multi_region")
p1=DimPlot(pbmc,cols=mycolors,label=TRUE,raster=FALSE,group.by="study",order=order,pt.size = 1)+NoAxes()
p1

p1=DimPlot(pbmc,cols=mycolors,label=FALSE,raster=FALSE,group.by="study",order=order,pt.size = 2)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

order=c("LGN","Angular gyrus","Anterior thalamus","Entorhinal cortex","Hippocampus","Midtemporal cortex","Prefrontal cortex","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#E31A1C","#FB9A99","#33A02C","#B2DF8A","#1F78B4","#A6CEE3","#FDBF6F"),label=TRUE,raster=FALSE,group.by="human_brain_region",order=order)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#E31A1C","#FB9A99","#33A02C","#B2DF8A","#1F78B4","#A6CEE3","#FDBF6F"),label=FALSE,raster=FALSE,group.by="human_brain_region",order=order,pt.size = 2)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1


order=c("LGN Inh NTRK1","LGN Inh CTXN3","LGN Inh LAMP5","LGN Inh TRPC4","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#A6CEE3","#1F78B4","#33A02C","#B2DF8A","#FB9A99"),label=TRUE,raster=FALSE,group.by="cluster2",order=order,pt.size = 2)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#A6CEE3","#1F78B4","#33A02C","#B2DF8A","#FB9A99"),label=FALSE,raster=FALSE,group.by="cluster2",order=order,pt.size = 2)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

order=c("Bakken","Linnarsson","DropViz","Yao","human_multi_region")
p1=DimPlot(pbmc,cols=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F"),label=TRUE,raster=FALSE,group.by="study",order=order,pt.size = 1)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F"),label=FALSE,raster=FALSE,group.by="study",order=order,pt.size = 2)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1



#Add additional plots

library(ggplot2)
library(Seurat)
setwd("D:/Mouse_data_integration")

#Load data
pbmc=readRDS(file="GABAergic_neurons_integrated_all_class_study_final.rds")

library(dplyr)
df=pbmc@meta.data
df2=df[,c("projid","brain_region")]

df3 = df2 %>% mutate(human_brain_region2 = case_when(grepl("Angular_gyrus",brain_region)~"Angular gyrus",
                                                     grepl("Anterior_thalamus",brain_region)~"Anterior thalamus",
                                                     grepl("Entorhinal_cortex",brain_region)~"Entorhinal cortex",
                                                     grepl("Hippocampus",brain_region)~"Hippocampus",
                                                     grepl("LGN",brain_region)~"not assigned",
                                                     grepl("Midtemporal_cortex",brain_region)~"Midtemporal cortex",
                                                     grepl("Prefrontal_cortex",brain_region)~"Prefrontal cortex",
                                                     is.na(df2$brain_region)~"not assigned"))
df4=df3[,-1]
df5=df4[,-1,drop=FALSE]
pbmc=AddMetaData(pbmc,df5)

order=c("Angular gyrus","Anterior thalamus","Entorhinal cortex","Hippocampus","Midtemporal cortex","Prefrontal cortex","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#E31A1C","#FB9A99","#33A02C","#B2DF8A","#1F78B4","#A6CEE3","#FDBF6F"),label=TRUE,raster=FALSE,group.by="human_brain_region2",order=order)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#E31A1C","#FB9A99","#33A02C","#B2DF8A","#1F78B4","#A6CEE3","#FDBF6F"),label=FALSE,raster=FALSE,group.by="human_brain_region2",order=order)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1


order=c("LGN","Angular gyrus","Anterior thalamus","Entorhinal cortex","Hippocampus","Midtemporal cortex","Prefrontal cortex","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","darkblue"),label=TRUE,raster=FALSE,group.by="human_brain_region",order=order,pt.size=2)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","darkblue"),label=FALSE,raster=FALSE,group.by="human_brain_region",order=order,pt.size=4)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1


order=c("LGN Inh NTRK1","LGN Inh CTXN3","LGN Inh LAMP5","LGN Inh TRPC4","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#6A3D9A","#33A02C","#1F78B4","#FF7F00"),label=TRUE,raster=FALSE,group.by="cluster2",order=order,pt.size = 2)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#6A3D9A","#33A02C","#1F78B4","#FF7F00"),label=FALSE,raster=FALSE,group.by="cluster2",order=order,pt.size = 4)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1



order=c("Thalamus","Cortex","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#1F78B4","#FF7F00"),label=TRUE,raster=FALSE,group.by="Region2",order=order,pt.size=2)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#1F78B4","#FF7F00"),label=FALSE,raster=FALSE,group.by="Region2",order=order,pt.size=2)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1


order=c("TH","Cortex","not assigned")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#1F78B4","#FF7F00"),label=TRUE,raster=FALSE,group.by="tissue2",order=order)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#1F78B4","#FF7F00"),label=FALSE,raster=FALSE,group.by="tissue2",order=order)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1


#Study
#Highlight Bakken
order=c("Bakken","Linnarsson","DropViz","Yao","human_multi_region")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#6A3D9A"),label=TRUE,raster=FALSE,group.by="study",order=order,pt.size = 1)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#6A3D9A"),label=FALSE,raster=FALSE,group.by="study",order=order,pt.size = 4)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1

#Highlight Linnarsson
order=c("Linnarsson","DropViz","Yao","human_multi_region","Bakken")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#33A02C"),label=TRUE,raster=FALSE,group.by="study",order=order,pt.size = 1)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#33A02C"),label=FALSE,raster=FALSE,group.by="study",order=order,pt.size = 1)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1


#Highlight DropViz
order=c("DropViz","Yao","human_multi_region","Bakken","Linnarsson")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#B15928"),label=TRUE,raster=FALSE,group.by="study",order=order,pt.size = 1)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#B15928"),label=FALSE,raster=FALSE,group.by="study",order=order,pt.size = 1)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1


#Highlight Yao
order=c("Yao","human_multi_region","Bakken","Linnarsson","DropViz")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#1F78B4"),label=TRUE,raster=FALSE,group.by="study",order=order,pt.size = 1)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#1F78B4"),label=FALSE,raster=FALSE,group.by="study",order=order)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1


#Highlight human_multi_region
order=c("human_multi_region","Bakken","Linnarsson","DropViz","Yao")
p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#FF7F00"),label=TRUE,raster=FALSE,group.by="study",order=order,pt.size = 1)+NoAxes()
p1

p1=DimPlot(pbmc,cols=c("#EAEAEA","#EAEAEA","#EAEAEA","#EAEAEA","#FF7F00"),label=FALSE,raster=FALSE,group.by="study",order=order)+NoAxes()+NoLegend()
p1=AugmentPlot(plot=p1)
p1