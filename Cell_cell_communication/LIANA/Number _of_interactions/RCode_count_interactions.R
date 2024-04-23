# load the reshape2 package
library(reshape2)
library(dplyr)

source_list=c("Exc","Inh","Ast","Oli","OPC","Mic","CAMs","T cells","End","Per","Fib","Epd","CPEC")

target_list=c("Exc","Inh","Ast","Oli","OPC","Mic","CAMs","T cells","End","Per","Fib","Epd","CPEC")

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

  data=read.csv("AG.csv")
  data=subset(data,aggregate_rank < 0.05)
  
  for (j in 1:length(source_list)){
    source_ct=source_list[[j]]
    data2=subset(data,source == source_ct)
    
    for (k in 1:length(target_list)){
      target_ct=target_list[[k]]
      data3=subset(data2, target == target_ct)
      interactions=nrow(data3)
      interactions=as.data.frame(interactions)
      rownames(interactions)=source_ct
      colnames(interactions)=target_ct
      if (k == 1){
        results=interactions
      }
      if (k > 1){
        results=cbind(results,interactions)
      }
    }
    if (j == 1){
      results_table=results
    }
    if (j > 1){
      results_table=rbind(results_table,results)
    }
  }
  
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Number_of_interactions")
write.csv(results_table,file="Significant_interactions_AG.csv")
##########################################################################################

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

data=read.csv("MT.csv")
data=subset(data,aggregate_rank < 0.05)

for (j in 1:length(source_list)){
  source_ct=source_list[[j]]
  data2=subset(data,source == source_ct)
  
  for (k in 1:length(target_list)){
    target_ct=target_list[[k]]
    data3=subset(data2, target == target_ct)
    interactions=nrow(data3)
    interactions=as.data.frame(interactions)
    rownames(interactions)=source_ct
    colnames(interactions)=target_ct
    if (k == 1){
      results=interactions
    }
    if (k > 1){
      results=cbind(results,interactions)
    }
  }
  if (j == 1){
    results_table=results
  }
  if (j > 1){
    results_table=rbind(results_table,results)
  }
}

setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Number_of_interactions")
write.csv(results_table,file="Significant_interactions_MT.csv")
##########################################################################################


##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

data=read.csv("PFC.csv")
data=subset(data,aggregate_rank < 0.05)

for (j in 1:length(source_list)){
  source_ct=source_list[[j]]
  data2=subset(data,source == source_ct)
  
  for (k in 1:length(target_list)){
    target_ct=target_list[[k]]
    data3=subset(data2, target == target_ct)
    interactions=nrow(data3)
    interactions=as.data.frame(interactions)
    rownames(interactions)=source_ct
    colnames(interactions)=target_ct
    if (k == 1){
      results=interactions
    }
    if (k > 1){
      results=cbind(results,interactions)
    }
  }
  if (j == 1){
    results_table=results
  }
  if (j > 1){
    results_table=rbind(results_table,results)
  }
}

setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Number_of_interactions")
write.csv(results_table,file="Significant_interactions_PFC.csv")
##########################################################################################


##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

data=read.csv("EC.csv")
data=subset(data,aggregate_rank < 0.05)

for (j in 1:length(source_list)){
  source_ct=source_list[[j]]
  data2=subset(data,source == source_ct)
  
  for (k in 1:length(target_list)){
    target_ct=target_list[[k]]
    data3=subset(data2, target == target_ct)
    interactions=nrow(data3)
    interactions=as.data.frame(interactions)
    rownames(interactions)=source_ct
    colnames(interactions)=target_ct
    if (k == 1){
      results=interactions
    }
    if (k > 1){
      results=cbind(results,interactions)
    }
  }
  if (j == 1){
    results_table=results
  }
  if (j > 1){
    results_table=rbind(results_table,results)
  }
}

setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Number_of_interactions")
write.csv(results_table,file="Significant_interactions_EC.csv")
##########################################################################################


##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

data=read.csv("HC.csv")
data=subset(data,aggregate_rank < 0.05)

for (j in 1:length(source_list)){
  source_ct=source_list[[j]]
  data2=subset(data,source == source_ct)
  
  for (k in 1:length(target_list)){
    target_ct=target_list[[k]]
    data3=subset(data2, target == target_ct)
    interactions=nrow(data3)
    interactions=as.data.frame(interactions)
    rownames(interactions)=source_ct
    colnames(interactions)=target_ct
    if (k == 1){
      results=interactions
    }
    if (k > 1){
      results=cbind(results,interactions)
    }
  }
  if (j == 1){
    results_table=results
  }
  if (j > 1){
    results_table=rbind(results_table,results)
  }
}

setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Number_of_interactions")
write.csv(results_table,file="Significant_interactions_HC.csv")
##########################################################################################

##########################################################################################
setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Results_aggregate_rank_magnitude_all_projid")

data=read.csv("TH.csv")
data=subset(data,aggregate_rank < 0.05)

for (j in 1:length(source_list)){
  source_ct=source_list[[j]]
  data2=subset(data,source == source_ct)
  
  for (k in 1:length(target_list)){
    target_ct=target_list[[k]]
    data3=subset(data2, target == target_ct)
    interactions=nrow(data3)
    interactions=as.data.frame(interactions)
    rownames(interactions)=source_ct
    colnames(interactions)=target_ct
    if (k == 1){
      results=interactions
    }
    if (k > 1){
      results=cbind(results,interactions)
    }
  }
  if (j == 1){
    results_table=results
  }
  if (j > 1){
    results_table=rbind(results_table,results)
  }
}

setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Number_of_interactions")
write.csv(results_table,file="Significant_interactions_TH.csv")
##########################################################################################
