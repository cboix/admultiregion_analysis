# load the reshape2 package
library(reshape2)
library(dplyr)


##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

  data=read.csv("TH.csv")
  data=subset(data,aggregate_rank < 0.05)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  

  table=data2

table_TH=table
##########################################################################################  
  
##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

data=read.csv("HC.csv")
data=subset(data,aggregate_rank < 0.05)
data[,"score"]=-log10(data[,"aggregate_rank"])
data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
data2=data[,c("interaction","score")]


table=data2

table_HC=table
##########################################################################################

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

data=read.csv("EC.csv")
data=subset(data,aggregate_rank < 0.05)
data[,"score"]=-log10(data[,"aggregate_rank"])
data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
data2=data[,c("interaction","score")]


table=data2

table_EC=table
##########################################################################################

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

data=read.csv("AG.csv")
data=subset(data,aggregate_rank < 0.05)
data[,"score"]=-log10(data[,"aggregate_rank"])
data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
data2=data[,c("interaction","score")]


table=data2

table_AG=table
##########################################################################################

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

data=read.csv("MT.csv")
data=subset(data,aggregate_rank < 0.05)
data[,"score"]=-log10(data[,"aggregate_rank"])
data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
data2=data[,c("interaction","score")]


table=data2

table_MT=table
##########################################################################################

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

data=read.csv("PFC.csv")
data=subset(data,aggregate_rank < 0.05)
data[,"score"]=-log10(data[,"aggregate_rank"])
data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
data2=data[,c("interaction","score")]


table=data2

table_PFC=table
##########################################################################################

table_AG=table_AG[,"interaction",drop=FALSE]
table_MT=table_MT[,"interaction",drop=FALSE]
table_PFC=table_PFC[,"interaction",drop=FALSE]
table_EC=table_EC[,"interaction",drop=FALSE]
table_HC=table_HC[,"interaction",drop=FALSE]
table_TH=table_TH[,"interaction",drop=FALSE]

table_AG[,"source"]=sapply(strsplit(table_AG$interaction, "_"), "[", 1)
table_AG[,"target"]=sapply(strsplit(table_AG$interaction, "_"), "[", 2)

table_MT[,"source"]=sapply(strsplit(table_MT$interaction, "_"), "[", 1)
table_MT[,"target"]=sapply(strsplit(table_MT$interaction, "_"), "[", 2)

table_PFC[,"source"]=sapply(strsplit(table_PFC$interaction, "_"), "[", 1)
table_PFC[,"target"]=sapply(strsplit(table_PFC$interaction, "_"), "[", 2)

table_EC[,"source"]=sapply(strsplit(table_EC$interaction, "_"), "[", 1)
table_EC[,"target"]=sapply(strsplit(table_EC$interaction, "_"), "[", 2)

table_HC[,"source"]=sapply(strsplit(table_HC$interaction, "_"), "[", 1)
table_HC[,"target"]=sapply(strsplit(table_HC$interaction, "_"), "[", 2)

table_TH[,"source"]=sapply(strsplit(table_TH$interaction, "_"), "[", 1)
table_TH[,"target"]=sapply(strsplit(table_TH$interaction, "_"), "[", 2)

##########################################################################################
region_list=list(table_AG=table_AG,table_MT=table_MT,table_PFC=table_PFC,table_EC=table_EC,table_HC=table_HC,table_TH=table_TH)

source_list=c("Exc","Inh","Ast","Oli","OPC","Mic","CAMs","T cells","End","Per","Fib","Epd","CPEC")

target_list=c("Exc","Inh","Ast","Oli","OPC","Mic","CAMs","T cells","End","Per","Fib","Epd","CPEC")

for (l in 1:length(region_list)){
  
  for (m in 1:length(region_list)){
  source_region=names(region_list)[l] 
  target_region=names(region_list)[m]
 
for (j in 1:length(source_list)){
  source_ct=source_list[[j]]
  for (k in 1:length(target_list)){
    target_ct=target_list[[k]]
    region_1=region_list[[l]]
    region_2=region_list[[m]]
    region_1=subset(region_1,source==source_ct)
    region_1=subset(region_1,target==target_ct)
    region_2=subset(region_2,source==source_ct)
    region_2=subset(region_2,target==target_ct)
    
    sum=sum(region_1$interaction %in% region_2$interaction)
    n_region_1=nrow(region_1)
    n_inter=sum/n_region_1*100
    n_inter=as.data.frame(n_inter)
    rownames(n_inter)=source_ct
    colnames(n_inter)=target_ct
    
    if (k == 1){
      results_df=n_inter
    }
    if (k > 1){
      results_df=cbind(results_df,n_inter)
    }
    
  }
  if (j == 1){
    results_table=results_df
  }
  if (j > 1){
    results_table=rbind(results_table,results_df)
  }
}

setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Number_overlapping_interactions_between_regions_fraction")
filename_out=paste0(source_region,"_",target_region,".csv")
write.csv(results_table,file=filename_out)
}
}



