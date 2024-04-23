library(ComplexHeatmap)
library(circlize)
library(dplyr)

setwd("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_proteomics/Correlation_EffectSize/Results")

pval=read.csv("Pval_adj_plaq_n.csv",row.names=1)
correlation_coefficient=read.csv("Correlation_coefficient_plaq_n.csv",row.names=1)

desired_order <- c("PFC", "MT", "AG", "EC", "HC")
celltypes <- c("Exc","Inh","Ast", "Oli", "OPC", "Mic")
# Reorder the rows of the data frame based on the desired order
pval2 <- pval[desired_order,celltypes]
correlation_coefficient2 <- correlation_coefficient[desired_order,celltypes]
#########################################################################################
#########################################################################################
#########################################################################################


gb = textGrob("*")
gb_w = convertWidth(grobWidth(gb), "mm")
gb_h = convertHeight(grobHeight(gb), "mm")

##########################################################
data=correlation_coefficient2
data2=pval2

h1=Heatmap(data,na_col = "#F0F0F0",
           border_gp = gpar(col = "black", lty = 1, lwd = 1),
           rect_gp = gpar(col = "white", lwd = 2),
           width = ncol(data)*unit(8, "mm"), 
           height = nrow(data)*unit(4, "mm"),
           cluster_rows = FALSE,cluster_columns = FALSE, col=colorRamp2(c(-1, 0, 1),c("#2166AC", "white", "#B2182B")), cell_fun = function(j, i, x, y, w, h, fill) {
             if(data2[i, j] < 0.05 & data2[i, j] >= 0.01) {
               grid.text("*", x, y - gb_h*0.5 + gb_w*0.5)
             } 
             if(data2[i, j] < 0.01 & data2[i, j] >= 0.001) {
               grid.text("**", x, y - gb_h*0.5 + gb_w*0.5)
             }
             if(data2[i, j] < 0.001 & data2[i, j] >= 0.0001) {
               grid.text("***", x, y - gb_h*0.5 + gb_w*0.5)
             }
             if(data2[i, j] < 0.0001) {
               grid.text("****", x, y - gb_h*0.5 + gb_w*0.5)
             }
             })
h1
#################################################################


#################################################################
setwd("D:/Projects/Large_scale_snRNA_seq/Brain_region_project/Brain_region_final/Evaluate_overlap_proteomics/Correlation_EffectSize/Results")

pdf("Correlation_effect_size_plaq_n_neurons_glia.pdf",height=4,width=8)
h1
dev.off()

#######################################################