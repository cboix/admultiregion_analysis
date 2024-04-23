# load the reshape2 package
library(reshape2)
library(dplyr)
library(lme4)
library(lmerTest)

#metadata
setwd("D:/PFC_427_WS/Vanshika")
meta=read.csv("dataset_652_basic_04-23-2020.csv")
meta=meta[,c("projid","age_death","pmi","msex","gpath")]

setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Sign_interaction_table")
data=read.csv("Sign_interactions_score_by_projid.csv",row.names=1)

datat=t(data)
datat=as.data.frame(datat)
colnames(datat)=datat[1,]
row_names_df_to_remove<-c("interaction")
datat=datat[!(row.names(datat) %in% row_names_df_to_remove),]

colnames(datat)=make.names(colnames(datat))

datat <- datat[,colSums(is.na(datat))<nrow(datat)]

datat <- datat %>% replace(is.na(.), 0)

for (i in 1:length(datat)){
  datat[,i]=as.numeric(datat[,i])
}

region_info=sub("\\_.*","",rownames(datat))
region_info=as.data.frame(region_info)
colnames(region_info)="region"

projid_info=sub(".*_","",rownames(datat))
projid_info=as.data.frame(projid_info)
colnames(projid_info)="projid"

addition=cbind(projid_info,region_info)
addition[,"projid"]=as.numeric(addition[,"projid"])

addition=left_join(addition,meta,by="projid")

data=cbind(addition,datat)

data$region <- relevel(factor(data$region), ref="PFC")

#lme
col_names=colnames(data[7:17256])
for (j in 1:length(col_names)){
  cell_types=col_names[[j]]
  fmla <- as.formula(paste0(cell_types,"~", "region","+","age_death","+","pmi","+","msex","+","(1 | projid)"))
  gpa_mixed = lmer(fmla, data = data)
  results=coef(summary(gpa_mixed))
  results=as.data.frame(results)
  results=results[c("regionAG","regionMT","regionEC","regionHC","regionTH"),]
  results[,"cell_types"]=cell_types
  results[,"region"]=row.names(results)
  
  if (j == 1){
    results_table=results
  }
  if (j > 1){
    results_table=rbind(results_table,results)
  }
}

#Bonferroni correction
results_table[,"p_adj"]=p.adjust(results_table$`Pr(>|t|)`, method = "bonferroni")

setwd("D:/home_wkstn/F/Brain_region_project/LIANA/Analysis_aggregate_rank_magnitude_all_projid/Sign_interaction_table")
write.csv(results_table,file="Results_lmer.csv")