# load the reshape2 package
library(reshape2)
library(dplyr)

setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Sign_interaction_table")
table=read.csv("Sign_interaction_table.csv",row.names=1)

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude/TH")

file_list=list.files(pattern="*.csv")

for (i in 1:length(file_list)){
  data=read.csv(file_list[[i]])
  data=subset(data,aggregate_rank < 0.05)
  individual=tools::file_path_sans_ext(file_list[[i]])
  individual=gsub("projid_","",individual)
  individual=paste0("TH_",individual)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- individual
  
    table=left_join(table,data2,by="interaction")
} 
##########################################################################################  
  
##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude/HC")

file_list=list.files(pattern="*.csv")

for (i in 1:length(file_list)){
  data=read.csv(file_list[[i]])
  data=subset(data,aggregate_rank < 0.05)
  individual=tools::file_path_sans_ext(file_list[[i]])
  individual=gsub("projid_","",individual)
  individual=paste0("HC_",individual)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- individual
  
  table=left_join(table,data2,by="interaction")
} 
##########################################################################################

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude/EC")

file_list=list.files(pattern="*.csv")

for (i in 1:length(file_list)){
  data=read.csv(file_list[[i]])
  data=subset(data,aggregate_rank < 0.05)
  individual=tools::file_path_sans_ext(file_list[[i]])
  individual=gsub("projid_","",individual)
  individual=paste0("EC_",individual)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- individual
  
  table=left_join(table,data2,by="interaction")
} 
##########################################################################################

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude/AG")

file_list=list.files(pattern="*.csv")

for (i in 1:length(file_list)){
  data=read.csv(file_list[[i]])
  data=subset(data,aggregate_rank < 0.05)
  individual=tools::file_path_sans_ext(file_list[[i]])
  individual=gsub("projid_","",individual)
  individual=paste0("AG_",individual)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- individual
  
  table=left_join(table,data2,by="interaction")
} 
##########################################################################################

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude/MT")

file_list=list.files(pattern="*.csv")

for (i in 1:length(file_list)){
  data=read.csv(file_list[[i]])
  data=subset(data,aggregate_rank < 0.05)
  individual=tools::file_path_sans_ext(file_list[[i]])
  individual=gsub("projid_","",individual)
  individual=paste0("MT_",individual)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- individual
  
  table=left_join(table,data2,by="interaction")
} 
##########################################################################################

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude/PFC")

file_list=list.files(pattern="*.csv")

for (i in 1:length(file_list)){
  data=read.csv(file_list[[i]])
  data=subset(data,aggregate_rank < 0.05)
  individual=tools::file_path_sans_ext(file_list[[i]])
  individual=gsub("projid_","",individual)
  individual=paste0("PFC_",individual)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- individual
  
  table=left_join(table,data2,by="interaction")
} 
##########################################################################################

setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Sign_interaction_table")
write.csv(table,file="Sign_interactions_score_by_projid.csv")