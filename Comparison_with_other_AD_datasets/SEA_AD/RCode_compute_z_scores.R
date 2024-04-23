library(dplyr)

region_list = c("AG","EC","HC","MT","PFC")

for (i in 1:length(region_list)){
  region=region_list[[i]]

directory=paste0("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_SEA_AD/Permutation_test/plaq_n/Overlap_results_up/",region)
setwd(directory)
observed=read.csv(file="Observed_overlap.csv")
permutation=read.csv(file="Permutation_overlap.csv",check.names=FALSE,row.names=1)
names(observed)[names(observed) == "X"] <- "comparison"

############################
celltype_list = c("Exc","Inh","Ast","Oli","OPC","Mic")

for (i in 1:length(celltype_list)){
  ct1=celltype_list[[i]]
  for (j in 1:length(celltype_list)){
    ct2=celltype_list[[j]]

    comprsn=paste0(ct1,"_",ct2)

    permutation_results <- permutation[,comprsn]
    observed_temp=subset(observed,comparison==comprsn)

    observed_value <- observed_temp$overlap

    z_score <- (observed_value-mean(permutation_results))/sd(permutation_results)

    z_score=as.data.frame(z_score)
    rownames(z_score)=ct2
    colnames(z_score)=ct1

    if (j==1){
    results=z_score
    }
    if (j>1){
    results=rbind(results,z_score)
    }
  }
  if (i==1){
    results_table=results
  }
  if (i>1){
    results_table=cbind(results_table,results)
  }
}

setwd("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_SEA_AD/Permutation_test/plaq_n/zscores_v2")
out_file=paste0("zscores_up_genes_",region,".csv")
write.csv(results_table,file=out_file)

library(ComplexHeatmap)
library(circlize)
library(dplyr)
##########################################################
data=results_table

p1=Heatmap(data,na_col = "#F0F0F0",
           border_gp = gpar(col = "black", lty = 1, lwd = 1),
           rect_gp = gpar(col = "white", lwd = 2),
           width = ncol(data)*unit(8, "mm"),
           height = nrow(data)*unit(8, "mm"),
           cluster_rows = FALSE,cluster_columns = FALSE, col=colorRamp2(c(0, 10, 20),c("#4575b4", "#ffffbf", "#d73027")))
#################################################################
p1
setwd("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_SEA_AD/Permutation_test/plaq_n/zscores_v2")
output_heatmap=paste0("zscores_up_genes_",region,".pdf")
pdf(output_heatmap,width = 20, height = 20)
print(p1)
dev.off()
}