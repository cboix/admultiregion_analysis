library(dplyr)

celltype_list = c("Exc","Inh","Ast","Oli","OPC","Mic")
region_list = c("AG","EC","HC","MT","PFC")

for (i in 1:length(region_list)){
  region=region_list[[i]]

Marker_genes_list <- list()
Marker_genes_list2 <- list()

for (i in 1:length(celltype_list)){
  celltype=celltype_list[[i]]

############################################################
setwd("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/DEGs_Carles")
  
load(file='allmethods.regional_plaq_n.merged.rda')
ct_names=read.csv("Cell_Type_Names.csv")
names(setdflist)
  
# Update the names in the list based on the translation mapping
  for (b in 1:length(setdflist)) {
    if (names(setdflist)[b] %in% ct_names$DEG_CTN) {
      names(setdflist)[b] <- ct_names$LR_CTN[ct_names$DEG_CTN == names(setdflist)[b]]
    }
  }
  
#Select genes considered
DEG_results<-setdflist[[celltype]]
DEG_results=DEG_results[DEG_results$region==region,]
  
#Determine number of DEGs per segment
DEG_significant=subset(DEG_results,col_nm == 1)
DEG_significant=DEG_significant[,"gene",drop=FALSE]
DEG_significant$cluster_id=celltype

Marker_genes=DEG_significant[,1:2]
Marker_genes[] <- lapply(Marker_genes, as.character)
if (nrow(Marker_genes)>0){
  Marker_genes_split=split(Marker_genes, Marker_genes$cluster_id)

  name=names(Marker_genes_split)[1]
  name=gsub("_","/",name)
  Marker_genes_list[[i]]<-Marker_genes_split[[1]]$gene
  names(Marker_genes_list)[i]<-name
}
##############################################################

############################################################
setwd("D:/nb_hans_wkstn2/G/SEA_AD/Analysis/muscat_DecontX/Results/CPS")
input_file2=paste0(celltype,"_CPS.csv")
data2=read.csv(input_file2,row.names=1)

DEG_results<-data2

DEG_significant=subset(DEG_results,p_adj.loc < 0.05)
DEG_significant=subset(DEG_significant,logFC < 0)

Marker_genes=DEG_significant[,1:2]
Marker_genes[] <- lapply(Marker_genes, as.character)
if (nrow(Marker_genes)>0){
  Marker_genes_split=split(Marker_genes, Marker_genes$cluster_id)

  name=names(Marker_genes_split)[1]
  name=gsub("_","/",name)
  Marker_genes_list2[[i]]<-Marker_genes_split[[1]]$gene
  names(Marker_genes_list2)[i]<-name
}
##############################################################

}

#Determine number of overlapping genes

celltype_list2=c("Exc","Inh","Ast","Oli","OPC","Mic")

for (m in 1:length(celltype_list2)){
  celltype1=celltype_list2[[m]]
  for (n in 1:length(celltype_list2)){
    celltype2=celltype_list2[[n]]

    overlap=length(intersect(Marker_genes_list[[celltype1]],Marker_genes_list2[[celltype2]]))
    row_name=paste0(celltype1,"_",celltype2)
    overlap=as.data.frame(overlap)
    rownames(overlap)=row_name

    if (n==1){
      table=overlap
    }
    if (n>1){
      table=rbind(table,overlap)
    }
  }
  if (m==1){
    results_table=table
  }
  if (m>1){
    results_table=rbind(results_table,table)
  }
}

out_directory=paste0("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_SEA_AD/Permutation_test/plaq_n/Overlap_results_down/",region)
setwd(out_directory)
write.csv(results_table,file="Observed_overlap.csv")
}