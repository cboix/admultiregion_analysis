library(corrplot)
library(dplyr)
library(RColorBrewer)

setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Results")

data=read.csv("Ast_metadata.csv")

colnames(data)

library(ggplot2)
######################################################################
#################################################
# Assuming your data frame is named 'data'
summary_data <- data %>%
  group_by(region, projid) %>%   # Group by 'region' and 'projid'
  summarize(gpath_up = mean(gpath_up_genes1))

df_list <- split(summary_data, summary_data$region)

for (i in 1:length(df_list)){
  data1=df_list[[i]]
  region=names(df_list)[i]
  data2=data1[,c("projid","gpath_up")]
  names(data2)[names(data2) == "gpath_up"] <- region
  
  if (i==1){
    results_df=data2
  }
  if (i>1){
    results_df=left_join(results_df,data2,by="projid")
  }
}

cor_data=results_df[-1]
M<-cor(cor_data, use="pairwise.complete.obs")
testRes = cor.mtest(cor_data, conf.level = 0.95)

pvals=testRes$p

# Flatten the matrix into a vector
p_values_vector <- as.vector(pvals)

# Adjust p-values for the entire matrix
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "BH")

# Reshape the result back into a matrix
adjusted_p_values_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(pvals))

rownames(adjusted_p_values_matrix) <- rownames(pvals)
colnames(adjusted_p_values_matrix) <- colnames(pvals)

mat=-log10(adjusted_p_values_matrix)
mat[is.infinite(mat)] <- 0
mat <- round(mat, digits = 1)

################################################

setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Analysis/Analysis_averaged_by_individual/Correlation_between_regions/Plots_significance")

pdf(file = "Ast_gpath_up.pdf",   # The directory you want to save the file in
    width = 3, # The width of the plot in inches
    height = 3) # The height of the plot in inches

corrplot(mat, is.corr=FALSE, method="color",addgrid.col = 'white',type="lower",
         col=brewer.pal(n=8, name="Oranges"),col.lim = c(0, 24),
         diag = FALSE, addCoef.col = 'white',tl.srt = 45, tl.col = 'black',cl.ratio=0.3)

dev.off()
######################################################################

#Down genes

#################################################
# Assuming your data frame is named 'data'
summary_data <- data %>%
  group_by(region, projid) %>%   # Group by 'region' and 'projid'
  summarize(gpath_down = mean(gpath_down_genes1))

df_list <- split(summary_data, summary_data$region)

for (i in 1:length(df_list)){
  data1=df_list[[i]]
  region=names(df_list)[i]
  data2=data1[,c("projid","gpath_down")]
  names(data2)[names(data2) == "gpath_down"] <- region
  
  if (i==1){
    results_df=data2
  }
  if (i>1){
    results_df=left_join(results_df,data2,by="projid")
  }
}

cor_data=results_df[-1]
M<-cor(cor_data, use="pairwise.complete.obs")
testRes = cor.mtest(cor_data, conf.level = 0.95)

pvals=testRes$p

# Flatten the matrix into a vector
p_values_vector <- as.vector(pvals)

# Adjust p-values for the entire matrix
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "BH")

# Reshape the result back into a matrix
adjusted_p_values_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(pvals))

rownames(adjusted_p_values_matrix) <- rownames(pvals)
colnames(adjusted_p_values_matrix) <- colnames(pvals)

mat=-log10(adjusted_p_values_matrix)
mat[is.infinite(mat)] <- 0
mat <- round(mat, digits = 1)

################################################

setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Analysis/Analysis_averaged_by_individual/Correlation_between_regions/Plots_significance")

pdf(file = "Ast_gpath_down.pdf",   # The directory you want to save the file in
    width = 3, # The width of the plot in inches
    height = 3) # The height of the plot in inches

corrplot(mat, is.corr=FALSE, method="color",addgrid.col = 'white',type="lower",
         col=brewer.pal(n=8, name="Oranges"),col.lim = c(0, 24),
         diag = FALSE, addCoef.col = 'white',tl.srt = 45, tl.col = 'black',cl.ratio=0.3)

dev.off()
######################################################################
######################################################################
######################################################################
######################################################################

#Mic
setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Results")

data=read.csv("Mic_metadata.csv")

colnames(data)

library(ggplot2)
######################################################################
#################################################
# Assuming your data frame is named 'data'
summary_data <- data %>%
  group_by(region, projid) %>%   # Group by 'region' and 'projid'
  summarize(gpath_up = mean(gpath_up_genes1))

df_list <- split(summary_data, summary_data$region)

for (i in 1:length(df_list)){
  data1=df_list[[i]]
  region=names(df_list)[i]
  data2=data1[,c("projid","gpath_up")]
  names(data2)[names(data2) == "gpath_up"] <- region
  
  if (i==1){
    results_df=data2
  }
  if (i>1){
    results_df=left_join(results_df,data2,by="projid")
  }
}

cor_data=results_df[-1]
M<-cor(cor_data, use="pairwise.complete.obs")
testRes = cor.mtest(cor_data, conf.level = 0.95)

pvals=testRes$p

# Flatten the matrix into a vector
p_values_vector <- as.vector(pvals)

# Adjust p-values for the entire matrix
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "BH")

# Reshape the result back into a matrix
adjusted_p_values_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(pvals))

rownames(adjusted_p_values_matrix) <- rownames(pvals)
colnames(adjusted_p_values_matrix) <- colnames(pvals)

mat=-log10(adjusted_p_values_matrix)
mat[is.infinite(mat)] <- 0
mat <- round(mat, digits = 1)

################################################

setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Analysis/Analysis_averaged_by_individual/Correlation_between_regions/Plots_significance")

pdf(file = "Mic_gpath_up.pdf",   # The directory you want to save the file in
    width = 3, # The width of the plot in inches
    height = 3) # The height of the plot in inches

corrplot(mat, is.corr=FALSE, method="color",addgrid.col = 'white',type="lower",
         col=brewer.pal(n=8, name="Oranges"),col.lim = c(0, 24),
         diag = FALSE, addCoef.col = 'white',tl.srt = 45, tl.col = 'black',cl.ratio=0.3)

dev.off()
######################################################################

#Down genes

#################################################
# Assuming your data frame is named 'data'
summary_data <- data %>%
  group_by(region, projid) %>%   # Group by 'region' and 'projid'
  summarize(gpath_down = mean(gpath_down_genes1))

df_list <- split(summary_data, summary_data$region)

for (i in 1:length(df_list)){
  data1=df_list[[i]]
  region=names(df_list)[i]
  data2=data1[,c("projid","gpath_down")]
  names(data2)[names(data2) == "gpath_down"] <- region
  
  if (i==1){
    results_df=data2
  }
  if (i>1){
    results_df=left_join(results_df,data2,by="projid")
  }
}

cor_data=results_df[-1]
M<-cor(cor_data, use="pairwise.complete.obs")
testRes = cor.mtest(cor_data, conf.level = 0.95)

pvals=testRes$p

# Flatten the matrix into a vector
p_values_vector <- as.vector(pvals)

# Adjust p-values for the entire matrix
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "BH")

# Reshape the result back into a matrix
adjusted_p_values_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(pvals))

rownames(adjusted_p_values_matrix) <- rownames(pvals)
colnames(adjusted_p_values_matrix) <- colnames(pvals)

mat=-log10(adjusted_p_values_matrix)
mat[is.infinite(mat)] <- 0
mat <- round(mat, digits = 1)

################################################

setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Analysis/Analysis_averaged_by_individual/Correlation_between_regions/Plots_significance")

pdf(file = "Mic_gpath_down.pdf",   # The directory you want to save the file in
    width = 3, # The width of the plot in inches
    height = 3) # The height of the plot in inches

corrplot(mat, is.corr=FALSE, method="color",addgrid.col = 'white',type="lower",
         col=brewer.pal(n=8, name="Oranges"),col.lim = c(0, 24),
         diag = FALSE, addCoef.col = 'white',tl.srt = 45, tl.col = 'black',cl.ratio=0.3)

dev.off()
######################################################################
######################################################################
######################################################################
######################################################################


#Oli
setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Results")

data=read.csv("Oli_metadata.csv")

colnames(data)

library(ggplot2)
######################################################################
#################################################
# Assuming your data frame is named 'data'
summary_data <- data %>%
  group_by(region, projid) %>%   # Group by 'region' and 'projid'
  summarize(gpath_up = mean(gpath_up_genes1))

df_list <- split(summary_data, summary_data$region)

for (i in 1:length(df_list)){
  data1=df_list[[i]]
  region=names(df_list)[i]
  data2=data1[,c("projid","gpath_up")]
  names(data2)[names(data2) == "gpath_up"] <- region
  
  if (i==1){
    results_df=data2
  }
  if (i>1){
    results_df=left_join(results_df,data2,by="projid")
  }
}

cor_data=results_df[-1]
M<-cor(cor_data, use="pairwise.complete.obs")
testRes = cor.mtest(cor_data, conf.level = 0.95)

pvals=testRes$p

# Flatten the matrix into a vector
p_values_vector <- as.vector(pvals)

# Adjust p-values for the entire matrix
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "BH")

# Reshape the result back into a matrix
adjusted_p_values_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(pvals))

rownames(adjusted_p_values_matrix) <- rownames(pvals)
colnames(adjusted_p_values_matrix) <- colnames(pvals)

mat=-log10(adjusted_p_values_matrix)
mat[is.infinite(mat)] <- 0
mat <- round(mat, digits = 1)

################################################

setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Analysis/Analysis_averaged_by_individual/Correlation_between_regions/Plots_significance")

pdf(file = "Oli_gpath_up.pdf",   # The directory you want to save the file in
    width = 3, # The width of the plot in inches
    height = 3) # The height of the plot in inches

corrplot(mat, is.corr=FALSE, method="color",addgrid.col = 'white',type="lower",
         col=brewer.pal(n=8, name="Oranges"),col.lim = c(0, 24),
         diag = FALSE, addCoef.col = 'white',tl.srt = 45, tl.col = 'black',cl.ratio=0.3)

dev.off()
######################################################################

#Down genes

#################################################
# Assuming your data frame is named 'data'
summary_data <- data %>%
  group_by(region, projid) %>%   # Group by 'region' and 'projid'
  summarize(gpath_down = mean(gpath_down_genes1))

df_list <- split(summary_data, summary_data$region)

for (i in 1:length(df_list)){
  data1=df_list[[i]]
  region=names(df_list)[i]
  data2=data1[,c("projid","gpath_down")]
  names(data2)[names(data2) == "gpath_down"] <- region
  
  if (i==1){
    results_df=data2
  }
  if (i>1){
    results_df=left_join(results_df,data2,by="projid")
  }
}

cor_data=results_df[-1]
M<-cor(cor_data, use="pairwise.complete.obs")
testRes = cor.mtest(cor_data, conf.level = 0.95)

pvals=testRes$p

# Flatten the matrix into a vector
p_values_vector <- as.vector(pvals)

# Adjust p-values for the entire matrix
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "BH")

# Reshape the result back into a matrix
adjusted_p_values_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(pvals))

rownames(adjusted_p_values_matrix) <- rownames(pvals)
colnames(adjusted_p_values_matrix) <- colnames(pvals)

mat=-log10(adjusted_p_values_matrix)
mat[is.infinite(mat)] <- 0
mat <- round(mat, digits = 1)

################################################

setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Analysis/Analysis_averaged_by_individual/Correlation_between_regions/Plots_significance")

pdf(file = "Oli_gpath_down.pdf",   # The directory you want to save the file in
    width = 3, # The width of the plot in inches
    height = 3) # The height of the plot in inches

corrplot(mat, is.corr=FALSE, method="color",addgrid.col = 'white',type="lower",
         col=brewer.pal(n=8, name="Oranges"),col.lim = c(0, 24),
         diag = FALSE, addCoef.col = 'white',tl.srt = 45, tl.col = 'black',cl.ratio=0.3)

dev.off()
######################################################################
######################################################################
######################################################################
######################################################################


#OPC
setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Results")

data=read.csv("OPC_metadata.csv")

colnames(data)

library(ggplot2)
######################################################################
#################################################
# Assuming your data frame is named 'data'
summary_data <- data %>%
  group_by(region, projid) %>%   # Group by 'region' and 'projid'
  summarize(gpath_up = mean(gpath_up_genes1))

df_list <- split(summary_data, summary_data$region)

for (i in 1:length(df_list)){
  data1=df_list[[i]]
  region=names(df_list)[i]
  data2=data1[,c("projid","gpath_up")]
  names(data2)[names(data2) == "gpath_up"] <- region
  
  if (i==1){
    results_df=data2
  }
  if (i>1){
    results_df=left_join(results_df,data2,by="projid")
  }
}

cor_data=results_df[-1]
M<-cor(cor_data, use="pairwise.complete.obs")
testRes = cor.mtest(cor_data, conf.level = 0.95)

pvals=testRes$p

# Flatten the matrix into a vector
p_values_vector <- as.vector(pvals)

# Adjust p-values for the entire matrix
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "BH")

# Reshape the result back into a matrix
adjusted_p_values_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(pvals))

rownames(adjusted_p_values_matrix) <- rownames(pvals)
colnames(adjusted_p_values_matrix) <- colnames(pvals)

mat=-log10(adjusted_p_values_matrix)
mat[is.infinite(mat)] <- 0
mat <- round(mat, digits = 1)

################################################

setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Analysis/Analysis_averaged_by_individual/Correlation_between_regions/Plots_significance")

pdf(file = "OPC_gpath_up.pdf",   # The directory you want to save the file in
    width = 3, # The width of the plot in inches
    height = 3) # The height of the plot in inches

corrplot(mat, is.corr=FALSE, method="color",addgrid.col = 'white',type="lower",
         col=brewer.pal(n=8, name="Oranges"),col.lim = c(0, 24),
         diag = FALSE, addCoef.col = 'white',tl.srt = 45, tl.col = 'black',cl.ratio=0.3)

dev.off()
######################################################################

#Down genes

#################################################
# Assuming your data frame is named 'data'
summary_data <- data %>%
  group_by(region, projid) %>%   # Group by 'region' and 'projid'
  summarize(gpath_down = mean(gpath_down_genes1))

df_list <- split(summary_data, summary_data$region)

for (i in 1:length(df_list)){
  data1=df_list[[i]]
  region=names(df_list)[i]
  data2=data1[,c("projid","gpath_down")]
  names(data2)[names(data2) == "gpath_down"] <- region
  
  if (i==1){
    results_df=data2
  }
  if (i>1){
    results_df=left_join(results_df,data2,by="projid")
  }
}

cor_data=results_df[-1]
M<-cor(cor_data, use="pairwise.complete.obs")
testRes = cor.mtest(cor_data, conf.level = 0.95)

pvals=testRes$p

# Flatten the matrix into a vector
p_values_vector <- as.vector(pvals)

# Adjust p-values for the entire matrix
adjusted_p_values_vector <- p.adjust(p_values_vector, method = "BH")

# Reshape the result back into a matrix
adjusted_p_values_matrix <- matrix(adjusted_p_values_vector, nrow = nrow(pvals))

rownames(adjusted_p_values_matrix) <- rownames(pvals)
colnames(adjusted_p_values_matrix) <- colnames(pvals)

mat=-log10(adjusted_p_values_matrix)
mat[is.infinite(mat)] <- 0
mat <- round(mat, digits = 1)

################################################

setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Analysis/Analysis_averaged_by_individual/Correlation_between_regions/Plots_significance")

pdf(file = "OPC_gpath_down.pdf",   # The directory you want to save the file in
    width = 3, # The width of the plot in inches
    height = 3) # The height of the plot in inches

corrplot(mat, is.corr=FALSE, method="color",addgrid.col = 'white',type="lower",
         col=brewer.pal(n=8, name="Oranges"),col.lim = c(0, 24),
         diag = FALSE, addCoef.col = 'white',tl.srt = 45, tl.col = 'black',cl.ratio=0.3)

dev.off()
######################################################################
######################################################################
######################################################################
######################################################################
