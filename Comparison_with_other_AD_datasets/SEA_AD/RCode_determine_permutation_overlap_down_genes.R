library(dplyr)

celltype_list=c("Exc","Inh","Ast","Oli","OPC","Mic")

region_list=c("AG","EC","HC","MT","PFC","TH")

for (r in 1:length(region_list)){
  region=region_list[[r]]
#####################################
genes_list <- list()

for (i in 1:length(celltype_list)){
  celltype=celltype_list[[i]]

working_directory=paste0("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_SEA_AD/Permutation_test/plaq_n/Random_gene_lists_down_DEGs/",region)

#Select genes considered

setwd(working_directory)
input_file=paste0(celltype,".csv")
data<-read.csv(input_file,header=TRUE,row.names=1)

genes_list[[i]]=data
name=gsub("_","/",celltype)
names(genes_list)[i]=name
}
######################################

#####################################
genes_list_SEAAD <- list()

for (k in 1:length(celltype_list)){
  celltype=celltype_list[[k]]

  working_directory=paste0("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_SEA_AD/Permutation_test/plaq_n/Random_gene_lists_down_DEGs/SEAAD")

  #Select genes considered

  setwd(working_directory)
  input_file=paste0(celltype,".csv")
  data<-read.csv(input_file,header=TRUE,row.names=1)

  genes_list_SEAAD[[k]]=data
  name=gsub("_","/",celltype)
  names(genes_list_SEAAD)[k]=name
}
######################################

###############################

#Determine number of overlapping genes

celltype_list2=c("Exc","Inh","Ast","Oli","OPC","Mic")

for (j in 1:1000){
  iteration=j
for (m in 1:length(celltype_list2)){
  celltype1=celltype_list2[[m]]
  for (n in 1:length(celltype_list2)){
    celltype2=celltype_list2[[n]]
    it_temp1=paste0("iteration_",iteration)
    genes_ct1=genes_list[[celltype1]][[it_temp1]]
    genes_ct2=genes_list_SEAAD[[celltype2]][[it_temp1]]
    overlap=length(intersect(genes_ct1,genes_ct2))
    col_name=paste0(celltype1,"_",celltype2)
    overlap=as.data.frame(overlap)
    colnames(overlap)=col_name
    rownames(overlap)=it_temp1

    if (n==1){
      table=overlap
    }
    if (n>1){
      table=cbind(table,overlap)
    }
  }
  if (m==1){
    results_table=table
  }
  if (m>1){
    results_table=cbind(results_table,table)
  }
}
  if (j==1){
    results_table_final=results_table
  }
  if (j>1){
    results_table_final=rbind(results_table_final,results_table)
  }
}

####################
out_directory=paste0("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_SEA_AD/Permutation_test/plaq_n/Overlap_results_down/",region)
setwd(out_directory)
write.csv(results_table_final,file="Permutation_overlap.csv")
}
