library(DESeq2)
library(dplyr)

setwd("G:/ROSMAP_RNAseq/ROSMAP_AMP_AD_Jose")
ROSMAP_AMP_AD=readRDS("ROSMAP_AMP-AD.rds")

expression_matrix=ROSMAP_AMP_AD@assays[["data"]]@listData[["counts"]]

#Add official metadata
setwd("G:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

# z-score
metadata_2=na.omit(metadata[,c("projid","age_death","pmi")])
metadata_2$age_death_z=(metadata_2$age_death - mean(metadata_2$age_death))/sd(metadata_2$age_death)
metadata_2$pmi_z=(metadata_2$pmi - mean(metadata_2$pmi))/sd(metadata_2$pmi)
metadata_2=metadata_2[,-c(2,3)]

metadata=left_join(metadata, metadata_2, by="projid")

#Add CDR/CR score variables
setwd("G:/PFC_429_final_Spring2021/Ghada/CR_score_files")
CRscores=read.csv("dataset_652_basic_04-23-2020_CDR_CR_scores.csv")

#Add to coldata
coldata=ROSMAP_AMP_AD@colData
coldata_df=as.data.frame(coldata)
coldata_df=coldata_df[,"projid", drop=FALSE]
coldata_df=left_join(coldata_df,metadata,by="projid")
coldata_df=left_join(coldata_df,CRscores,by="projid")

gene_anno=read.csv("G:/ROSMAP_RNAseq/gene_id_name.csv")
colnames(gene_anno)[1]="Entrez.ID"
gene_anno=gene_anno[!gene_anno$Gene.name=="",]

####################################

variable_list=c(
  "cogn_global_lv", 
  "cogn_ep_lv", 
  "cogn_po_lv",
  "cogn_ps_lv",
  "cogn_se_lv",
  "cogn_wo_lv",
  "cogng_random_slope",
  "cognep_random_slope",
  "cognpo_random_slope",
  "cognps_random_slope",
  "cognse_random_slope",
  "cognwo_random_slope")

for(var in 1:length(variable_list)){
  print(variable_list[var])

#Select projids with complete data
rownames(coldata_df)=coldata_df[,"projid"]
formula = as.formula(paste0("~",variable_list[var],"+","age_death_z","+","pmi_z"))
cd2=coldata_df[,c(variable_list[var],"age_death_z","pmi_z")]
design=model.matrix(formula,cd2)

#Subset expression matrix and coldata
expression_matrix_subset=expression_matrix[,rownames(design)]
coldata_df_subset=coldata_df[rownames(design),]

#DESeq2
dds <- DESeqDataSetFromMatrix(countData = expression_matrix_subset,
                              colData = coldata_df_subset,
                              design = formula)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds,name = variable_list[var])

resOrdered <- res[order(res$pvalue),]
resOrdered <- as.data.frame(resOrdered)

resOrdered=cbind(Entrez.ID=rownames(resOrdered),resOrdered)
resOrdered$Entrez.ID=tools::file_path_sans_ext(resOrdered$Entrez.ID)
resOrdered=left_join(resOrdered,gene_anno,by="Entrez.ID")
resOrdered=na.omit(resOrdered)
#resOrdered=resOrdered[duplicated(resOrdered$Gene.name),]
rownames(resOrdered)=resOrdered[,1]
resOrdered=resOrdered[,-1]

write.csv(resOrdered,file=paste0("G:/PFC_429_final_Spring2021/Ghada/ROSMAP_RNAseq/DESeq2/Results_z/","RNAseq_",variable_list[var],".csv"))
assign(variable_list[var],resOrdered)


#res <- results(dds, name="condition_trt_vs_untrt")
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")


}



