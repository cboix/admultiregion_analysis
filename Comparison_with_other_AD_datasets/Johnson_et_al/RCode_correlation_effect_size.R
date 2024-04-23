library(dplyr)

celltype_list = c("Exc","Inh","Ast","Oli","OPC","Mic")
region_list = c("AG","EC","HC","MT","PFC")

#Select variables
variable_list=list("plaq_n","nft","nrad","cogdxad")

for (k in 1:length(variable_list)) {
variable=variable_list[[k]]

############################################################
setwd("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/DEGs_Carles")
input_file=paste0("allmethods.regional_",variable,".merged.rda")
load(file=input_file)
ct_names=read.csv("Cell_Type_Names.csv")
names(setdflist)

# Update the names in the list based on the translation mapping
for (b in 1:length(setdflist)) {
  if (names(setdflist)[b] %in% ct_names$DEG_CTN) {
    names(setdflist)[b] <- ct_names$LR_CTN[ct_names$DEG_CTN == names(setdflist)[b]]
  }
}

#Loop through regions and cell types
for (r in 1:length(region_list)){
  region=region_list[[r]]

  data_list=list()
  for (c in 1:length(celltype_list)){
  celltype=celltype_list[[c]]
  
  #Select genes considered
  DEG_results<-setdflist[[celltype]]
  DEG_results=DEG_results[DEG_results$region==region,]
  
  #Determine number of DEGs per segment
  DEG_significant <- DEG_results[DEG_results$col_nm %in% c(1, 2), ]
  Sign_genes1=DEG_significant[,c("gene","logFC_nb")]
  data_list[[c]]=Sign_genes1
  names(data_list)[c]=celltype
  }
  
  #load published data
  setwd("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_proteomics/Data")
  Prot=read.csv("Proteomics_AD_vs_control.csv",row.names=1)
  
  DE_results<-Prot
  
  DE_results_significant=subset(DE_results,adj_pval < 0.05)
  sign_genes2=DE_results_significant[,c("gene","effect_size")]
  ########################################
  for (n in 1:length(data_list)) {
    celltype2=names(data_list)[n]
    sign_genes1=data_list[[n]]
    
    combined=left_join(sign_genes1,sign_genes2,by="gene")
    
    OK <- complete.cases(combined$logFC_nb,combined$effect_size)
    x <- combined$logFC_nb[OK]
    y <- combined$effect_size[OK]
    m <- length(x)
    if (m < 3L){
      results_p <- "NA"
      results_p = as.data.frame(results_p)
      rownames(results_p)=region
      colnames(results_p)=celltype2
      
      results_r <- "NA"
      results_r = as.data.frame(results_r)
      rownames(results_r)=region
      colnames(results_r)=celltype2
      
      if (n==1){
        results_df_p=results_p
        results_df_r=results_r
      }
      if (n>1){
        results_df_p=cbind(results_df_p,results_p)
        results_df_r=cbind(results_df_r,results_r)
      }
    } else {
      cor_results=cor.test(combined$logFC_nb,combined$effect_size,alternative = c("two.sided"),method=c("pearson"))
      results_p=as.data.frame(cor_results[["p.value"]])
      rownames(results_p)=region
      colnames(results_p)=celltype2
      
      results_r=as.data.frame(cor_results[["estimate"]])
      rownames(results_r)=region
      colnames(results_r)=celltype2
      
      if (n==1){
        results_df_p=results_p
        results_df_r=results_r
      }
      if (n>1){
        results_df_p=cbind(results_df_p,results_p)
        results_df_r=cbind(results_df_r,results_r)
      }
    }
  }
  if (r==1){
    results_final_p=results_df_p
    results_final_r=results_df_r
  }
  if (r>1){
    results_final_p=rbind(results_final_p,results_df_p)
    results_final_r=rbind(results_final_r,results_df_r)
  }
}

# Reshape the data frame to a vector
p_values <- as.vector(as.matrix(results_final_p))

# Apply BH adjustment to the vector of p-values
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Reshape the adjusted p-values back to a data frame
adjusted_df <- as.data.frame(matrix(adjusted_p_values, nrow = nrow(results_final_p)))

# Assign column names and row names to the adjusted data frame
colnames(adjusted_df) <- colnames(results_final_p)
rownames(adjusted_df) <- rownames(results_final_p)


setwd("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_proteomics/Correlation_EffectSize/Results")
outfile1=paste0("Pval_",variable,".csv")
outfile2=paste0("Correlation_coefficient_",variable,".csv")
outfile3=paste0("Pval_adj_",variable,".csv")
write.csv(results_final_p,file=outfile1)
write.csv(results_final_r,file=outfile2)
write.csv(adjusted_df,file=outfile3)
}
###################



