# load the reshape2 package
library(reshape2)
library(dplyr)

#Table interactions that are significant in at least one region

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

  data=read.csv("AG.csv")
  data=subset(data,aggregate_rank < 0.05)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- "AG"
  table=data2
##########################################################################################  
  
##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")
  
  data=read.csv("MT.csv")
  data=subset(data,aggregate_rank < 0.05)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- "MT"
  
  table=full_join(table,data2,by="interaction")
##########################################################################################

##########################################################################################
  setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")
  
  data=read.csv("PFC.csv")
  data=subset(data,aggregate_rank < 0.05)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- "PFC"
  
  table=full_join(table,data2,by="interaction")
##########################################################################################

##########################################################################################
  setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")
  
  data=read.csv("EC.csv")
  data=subset(data,aggregate_rank < 0.05)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- "EC"
  
  table=full_join(table,data2,by="interaction")
##########################################################################################

##########################################################################################
  setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")
  
  data=read.csv("HC.csv")
  data=subset(data,aggregate_rank < 0.05)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- "HC"
  
  table=full_join(table,data2,by="interaction")
##########################################################################################

##########################################################################################
  setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")
  
  data=read.csv("TH.csv")
  data=subset(data,aggregate_rank < 0.05)
  data[,"score"]=-log10(data[,"aggregate_rank"])
  data[,"interaction"]=paste(data$source,data$target,data$ligand.complex,data$receptor.complex,sep = "_")
  data2=data[,c("interaction","score")]
  names(data2)[names(data2) == "score"] <- "TH"
  
  table=full_join(table,data2,by="interaction")

##########################################################################################

sign_interactions=table[,"interaction",drop=FALSE]

setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Sign_interaction_table")
write.csv(sign_interactions,file="Sign_interaction_table.csv")
write.csv(table,file="Sign_interaction_table_scores.csv")

##########################################################################################
#Subset conserved interactions
table <- subset(table, AG > -log10(0.05) & MT > -log10(0.05) & PFC > -log10(0.05) & EC > -log10(0.05) & HC > -log10(0.05) & TH > -log10(0.05))

write.csv(table,file="Conserved_Sign_interaction_table_scores.csv")