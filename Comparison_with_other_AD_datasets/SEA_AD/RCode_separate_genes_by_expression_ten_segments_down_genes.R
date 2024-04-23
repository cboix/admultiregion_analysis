library(dplyr)

setwd("D:/nb_hans_wkstn2/G/PFC_429_final_Spring2021/Vanshika/muscat_DecontX/Analysis_Results/Evaluate_overlap_DeJager/PermutationTest")
ct_names=read.csv("celltype_order_Mic_Ast_namesR_v3.csv")

setwd("D:/home_wkstn/F/Brain_region_project/Average_expression/Results")

AG=read.csv("AG_average_expression_data_by_celltypes.csv",row.names=1)
EC=read.csv("EC_average_expression_data_by_celltypes.csv",row.names=1)
HC=read.csv("HC_average_expression_data_by_celltypes.csv",row.names=1)
MT=read.csv("MT_average_expression_data_by_celltypes.csv",row.names=1)
PFC=read.csv("PFC_average_expression_data_by_celltypes.csv",row.names=1)
TH=read.csv("TH_average_expression_data_by_celltypes.csv",row.names=1)

setwd("D:/nb_hans_wkstn2/G/PFC_429_final_Spring2021/Sudhagar/Avg_Exp_data")
SEAAD=read.csv("SEAAD_DecontX_average_expression_data_by_celltypes.csv",row.names=1)

rownames(AG)=AG$gene
rownames(EC)=EC$gene
rownames(HC)=HC$gene
rownames(MT)=MT$gene
rownames(PFC)=PFC$gene
rownames(TH)=TH$gene
rownames(SEAAD)=SEAAD$gene

AG=AG[-1]
EC=EC[-1]
HC=HC[-1]
MT=MT[-1]
PFC=PFC[-1]
TH=TH[-1]
SEAAD=SEAAD[-1]

AG_cols=colnames(AG)
EC_cols=colnames(EC)
HC_cols=colnames(HC)
MT_cols=colnames(MT)
PFC_cols=colnames(PFC)
TH_cols=colnames(TH)
SEAAD_cols=colnames(SEAAD)

AG_cols=as.data.frame(AG_cols)
EC_cols=as.data.frame(EC_cols)
HC_cols=as.data.frame(HC_cols)
MT_cols=as.data.frame(MT_cols)
PFC_cols=as.data.frame(PFC_cols)
TH_cols=as.data.frame(TH_cols)
SEAAD_cols=as.data.frame(SEAAD_cols)

colnames(AG_cols)="celltypeR"
colnames(EC_cols)="celltypeR"
colnames(HC_cols)="celltypeR"
colnames(MT_cols)="celltypeR"
colnames(PFC_cols)="celltypeR"
colnames(TH_cols)="celltypeR"
colnames(SEAAD_cols)="celltypeR"

AG_cols=left_join(AG_cols,ct_names,by="celltypeR")
EC_cols=left_join(EC_cols,ct_names,by="celltypeR")
HC_cols=left_join(HC_cols,ct_names,by="celltypeR")
MT_cols=left_join(MT_cols,ct_names,by="celltypeR")
PFC_cols=left_join(PFC_cols,ct_names,by="celltypeR")
TH_cols=left_join(TH_cols,ct_names,by="celltypeR")
SEAAD_cols=left_join(SEAAD_cols,ct_names,by="celltypeR")

# Replace the column names
colnames(AG) <- AG_cols$celltype_file
colnames(EC) <- EC_cols$celltype_file
colnames(HC) <- HC_cols$celltype_file
colnames(MT) <- MT_cols$celltype_file
colnames(PFC) <- PFC_cols$celltype_file
colnames(TH) <- TH_cols$celltype_file
colnames(SEAAD) <- SEAAD_cols$celltype_file

average_expression_list=list(AG=AG,EC=EC,HC=HC,MT=MT,PFC=PFC,TH=TH)

#############################################

celltype_list = c("Exc","Inh","Ast","Oli","OPC","Mic")
region_list = c("AG","EC","HC","MT","PFC","TH")

for (r in 1:length(region_list)){
  region=region_list[[r]]

for (c in 1:length(celltype_list)){
  celltype=celltype_list[[c]]

##########################################################################################################################
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

All_genes=DEG_results[,"gene",drop=FALSE]

average_expression=average_expression_list[[region]]
average_expression=average_expression[,celltype,drop=FALSE]
colnames(average_expression) = "average_expression"
average_expression$gene=rownames(average_expression)

All_genes_expression=left_join(All_genes,average_expression,by="gene")

q10=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.1))
q20=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.2))
q30=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.3))
q40=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.4))
q50=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.5))
q60=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.6))
q70=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.7))
q80=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.8))
q90=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.9))

All_genes_expression = All_genes_expression %>% mutate(expression_quartile =
                                   case_when(average_expression <= q10 ~ "segment1",
                                             average_expression > q10 & average_expression <= q20 ~ "segment2",
                                             average_expression > q20 & average_expression <= q30 ~ "segment3",
                                             average_expression > q30 & average_expression <= q40 ~ "segment4",
                                             average_expression > q40 & average_expression <= q50 ~ "segment5",
                                             average_expression > q50 & average_expression <= q60 ~ "segment6",
                                             average_expression > q60 & average_expression <= q70 ~ "segment7",
                                             average_expression > q70 & average_expression <= q80 ~ "segment8",
                                             average_expression > q80 & average_expression <= q90 ~ "segment9",
                                             average_expression > q90 ~ "segment10")
)

segment1=subset(All_genes_expression, expression_quartile=="segment1")
segment2=subset(All_genes_expression, expression_quartile=="segment2")
segment3=subset(All_genes_expression, expression_quartile=="segment3")
segment4=subset(All_genes_expression, expression_quartile=="segment4")
segment5=subset(All_genes_expression, expression_quartile=="segment5")
segment6=subset(All_genes_expression, expression_quartile=="segment6")
segment7=subset(All_genes_expression, expression_quartile=="segment7")
segment8=subset(All_genes_expression, expression_quartile=="segment8")
segment9=subset(All_genes_expression, expression_quartile=="segment9")
segment10=subset(All_genes_expression, expression_quartile=="segment10")


gene_lists=list(segment1=segment1$gene,
                segment2=segment2$gene,
                segment3=segment3$gene,
                segment4=segment4$gene,
                segment5=segment5$gene,
                segment6=segment6$gene,
                segment7=segment7$gene,
                segment8=segment8$gene,
                segment9=segment9$gene,
                segment10=segment10$gene)

#Determine number of DEGs per segment
DEG_significant=subset(DEG_results,col_nm == 1)

nDEGs_segment1=nrow(DEG_significant[DEG_significant$gene %in% gene_lists[["segment1"]], ])
nDEGs_segment2=nrow(DEG_significant[DEG_significant$gene %in% gene_lists[["segment2"]], ])
nDEGs_segment3=nrow(DEG_significant[DEG_significant$gene %in% gene_lists[["segment3"]], ])
nDEGs_segment4=nrow(DEG_significant[DEG_significant$gene %in% gene_lists[["segment4"]], ])
nDEGs_segment5=nrow(DEG_significant[DEG_significant$gene %in% gene_lists[["segment5"]], ])
nDEGs_segment6=nrow(DEG_significant[DEG_significant$gene %in% gene_lists[["segment6"]], ])
nDEGs_segment7=nrow(DEG_significant[DEG_significant$gene %in% gene_lists[["segment7"]], ])
nDEGs_segment8=nrow(DEG_significant[DEG_significant$gene %in% gene_lists[["segment8"]], ])
nDEGs_segment9=nrow(DEG_significant[DEG_significant$gene %in% gene_lists[["segment9"]], ])
nDEGs_segment10=nrow(DEG_significant[DEG_significant$gene %in% gene_lists[["segment10"]], ])
##########################################################################################################

####################################################################################
#Determine number of iterations
for (m in 1:1000){

  seg1_genes=sample(gene_lists[["segment1"]], nDEGs_segment1, replace = FALSE)
  seg2_genes=sample(gene_lists[["segment2"]], nDEGs_segment2, replace = FALSE)
  seg3_genes=sample(gene_lists[["segment3"]], nDEGs_segment3, replace = FALSE)
  seg4_genes=sample(gene_lists[["segment4"]], nDEGs_segment4, replace = FALSE)
  seg5_genes=sample(gene_lists[["segment5"]], nDEGs_segment5, replace = FALSE)
  seg6_genes=sample(gene_lists[["segment6"]], nDEGs_segment6, replace = FALSE)
  seg7_genes=sample(gene_lists[["segment7"]], nDEGs_segment7, replace = FALSE)
  seg8_genes=sample(gene_lists[["segment8"]], nDEGs_segment8, replace = FALSE)
  seg9_genes=sample(gene_lists[["segment9"]], nDEGs_segment9, replace = FALSE)
  seg10_genes=sample(gene_lists[["segment10"]], nDEGs_segment10, replace = FALSE)

  all_segments=c(seg1_genes,seg2_genes,seg3_genes,seg4_genes,seg5_genes,seg6_genes,seg7_genes,seg8_genes,seg9_genes,seg10_genes)

  all_segments=as.data.frame(all_segments)
  col_name=paste0("iteration_",m)
  colnames(all_segments)=col_name

  if (m==1){
    table=all_segments
  }
  if (m>1){
    table=cbind(table,all_segments)
  }
}

out_directory=paste0("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_SEA_AD/Permutation_test/plaq_n/Random_gene_lists_down_DEGs/",region)
setwd(out_directory)
output_file=paste0(celltype,".csv")
write.csv(table,file=output_file)
####################################################################################

}
}
#################################################
##############################################
##########################################


for (c in 1:length(celltype_list)){
    celltype=celltype_list[[c]]

##########################################################################################################################
#Select genes considered

setwd("D:/nb_hans_wkstn2/G/SEA_AD/Analysis/muscat_DecontX/Results/CPS")
input_file2=paste0(celltype,"_CPS.csv")
data2=read.csv(input_file2,row.names=1)

#data2=subset(data2, p_adj.loc < 0.05)
#data2_up=subset(data2, logFC > 0)
#data2_down=subset(data2, logFC < 0)

DEG_results<-data2

All_genes=DEG_results[,"gene",drop=FALSE]

average_expression=SEAAD
average_expression=average_expression[,celltype,drop=FALSE]
colnames(average_expression) = "average_expression"
average_expression$gene=rownames(average_expression)

All_genes_expression=left_join(All_genes,average_expression,by="gene")

q10=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.1))
q20=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.2))
q30=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.3))
q40=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.4))
q50=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.5))
q60=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.6))
q70=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.7))
q80=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.8))
q90=quantile(All_genes_expression$average_expression,na.rm = T,prob = c(0.9))

All_genes_expression = All_genes_expression %>% mutate(expression_quartile =
                                                         case_when(average_expression <= q10 ~ "segment1",
                                                                   average_expression > q10 & average_expression <= q20 ~ "segment2",
                                                                   average_expression > q20 & average_expression <= q30 ~ "segment3",
                                                                   average_expression > q30 & average_expression <= q40 ~ "segment4",
                                                                   average_expression > q40 & average_expression <= q50 ~ "segment5",
                                                                   average_expression > q50 & average_expression <= q60 ~ "segment6",
                                                                   average_expression > q60 & average_expression <= q70 ~ "segment7",
                                                                   average_expression > q70 & average_expression <= q80 ~ "segment8",
                                                                   average_expression > q80 & average_expression <= q90 ~ "segment9",
                                                                   average_expression > q90 ~ "segment10")
)

segment1_SEAAD=subset(All_genes_expression, expression_quartile=="segment1")
segment2_SEAAD=subset(All_genes_expression, expression_quartile=="segment2")
segment3_SEAAD=subset(All_genes_expression, expression_quartile=="segment3")
segment4_SEAAD=subset(All_genes_expression, expression_quartile=="segment4")
segment5_SEAAD=subset(All_genes_expression, expression_quartile=="segment5")
segment6_SEAAD=subset(All_genes_expression, expression_quartile=="segment6")
segment7_SEAAD=subset(All_genes_expression, expression_quartile=="segment7")
segment8_SEAAD=subset(All_genes_expression, expression_quartile=="segment8")
segment9_SEAAD=subset(All_genes_expression, expression_quartile=="segment9")
segment10_SEAAD=subset(All_genes_expression, expression_quartile=="segment10")


gene_lists_SEAAD=list(segment1=segment1_SEAAD$gene,
                segment2=segment2_SEAAD$gene,
                segment3=segment3_SEAAD$gene,
                segment4=segment4_SEAAD$gene,
                segment5=segment5_SEAAD$gene,
                segment6=segment6_SEAAD$gene,
                segment7=segment7_SEAAD$gene,
                segment8=segment8_SEAAD$gene,
                segment9=segment9_SEAAD$gene,
                segment10=segment10_SEAAD$gene)

#Determine number of DEGs per segment
DEG_significant=subset(DEG_results,p_adj.loc < 0.05)
DEG_significant=subset(DEG_significant,logFC < 0)

nDEGs_segment1_SEAAD=nrow(DEG_significant[DEG_significant$gene %in% gene_lists_SEAAD[["segment1"]], ])
nDEGs_segment2_SEAAD=nrow(DEG_significant[DEG_significant$gene %in% gene_lists_SEAAD[["segment2"]], ])
nDEGs_segment3_SEAAD=nrow(DEG_significant[DEG_significant$gene %in% gene_lists_SEAAD[["segment3"]], ])
nDEGs_segment4_SEAAD=nrow(DEG_significant[DEG_significant$gene %in% gene_lists_SEAAD[["segment4"]], ])
nDEGs_segment5_SEAAD=nrow(DEG_significant[DEG_significant$gene %in% gene_lists_SEAAD[["segment5"]], ])
nDEGs_segment6_SEAAD=nrow(DEG_significant[DEG_significant$gene %in% gene_lists_SEAAD[["segment6"]], ])
nDEGs_segment7_SEAAD=nrow(DEG_significant[DEG_significant$gene %in% gene_lists_SEAAD[["segment7"]], ])
nDEGs_segment8_SEAAD=nrow(DEG_significant[DEG_significant$gene %in% gene_lists_SEAAD[["segment8"]], ])
nDEGs_segment9_SEAAD=nrow(DEG_significant[DEG_significant$gene %in% gene_lists_SEAAD[["segment9"]], ])
nDEGs_segment10_SEAAD=nrow(DEG_significant[DEG_significant$gene %in% gene_lists_SEAAD[["segment10"]], ])
##########################################################################################################

####################################################################################
#Determine number of iterations
for (m in 1:1000){

  seg1_genes2=sample(gene_lists_SEAAD[["segment1"]], nDEGs_segment1_SEAAD, replace = FALSE)
  seg2_genes2=sample(gene_lists_SEAAD[["segment2"]], nDEGs_segment2_SEAAD, replace = FALSE)
  seg3_genes2=sample(gene_lists_SEAAD[["segment3"]], nDEGs_segment3_SEAAD, replace = FALSE)
  seg4_genes2=sample(gene_lists_SEAAD[["segment4"]], nDEGs_segment4_SEAAD, replace = FALSE)
  seg5_genes2=sample(gene_lists_SEAAD[["segment5"]], nDEGs_segment5_SEAAD, replace = FALSE)
  seg6_genes2=sample(gene_lists_SEAAD[["segment6"]], nDEGs_segment6_SEAAD, replace = FALSE)
  seg7_genes2=sample(gene_lists_SEAAD[["segment7"]], nDEGs_segment7_SEAAD, replace = FALSE)
  seg8_genes2=sample(gene_lists_SEAAD[["segment8"]], nDEGs_segment8_SEAAD, replace = FALSE)
  seg9_genes2=sample(gene_lists_SEAAD[["segment9"]], nDEGs_segment9_SEAAD, replace = FALSE)
  seg10_genes2=sample(gene_lists_SEAAD[["segment10"]], nDEGs_segment10_SEAAD, replace = FALSE)

  all_segments=c(seg1_genes2,seg2_genes2,seg3_genes2,seg4_genes2,seg5_genes2,seg6_genes2,seg7_genes2,seg8_genes2,seg9_genes2,seg10_genes2)

  all_segments=as.data.frame(all_segments)
  col_name=paste0("iteration_",m)
  colnames(all_segments)=col_name

  if (m==1){
    table=all_segments
  }
  if (m>1){
    table=cbind(table,all_segments)
  }
}

setwd("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_SEA_AD/Permutation_test/plaq_n/Random_gene_lists_down_DEGs/SEAAD")
output_file=paste0(celltype,".csv")
write.csv(table,file=output_file)
####################################################################################

}
