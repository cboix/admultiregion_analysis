library(dplyr)
library(gridExtra)
setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Results")

data=read.csv("Mic_metadata.csv")

colnames(data)

###############################################################################
#Load metadata
setwd("F:/PFC_429_final_Spring2021/Vanshika")
metadata=read.csv("dataset_652_basic_04-23-2020.csv")

#Determine mean module scores
summary_data <- data %>%
  group_by(region, projid) %>%   # Group by 'region' and 'projid'
  summarize(average_gpath_up = mean(gpath_up_genes1))

summary_data=left_join(summary_data,metadata,by="projid")

summary_data_down <- data %>%
  group_by(region, projid) %>%   # Group by 'region' and 'projid'
  summarize(average_gpath_down = mean(gpath_down_genes1))

summary_data_down=left_join(summary_data_down,metadata,by="projid")
########################################################################

library(ggplot2)

#One single plot

library(RColorBrewer)

custom_colors <- data.frame(
  region = unique(data$region),
  color = c("#8DD3C7", "#FFED6F", "#BEBADA", "#80B1D3", "#FDB462", "#B3DE69") # Add colors of your choice
)

#Up module
gpath=ggplot(data = summary_data, aes(x = gpath, y = average_gpath_up, color = region)) +
  geom_smooth(method = "loess", se = TRUE, level=0.95, linetype = 0, show.legend = FALSE, fill = "#D3D3D3")+
  geom_smooth(method = "loess", se = FALSE, level=0.95, linewidth =1.5, show.legend = TRUE)+
  labs(x = "Global AD pathology", y = "Module score") +
  scale_color_manual(values = custom_colors$color) +
  theme(aspect.ratio=1)+
  theme_bw()

#Down module
gpath_down=ggplot(data = summary_data_down, aes(x = gpath, y = average_gpath_down, color = region)) +
  geom_smooth(method = "loess", se = TRUE, level=0.95, linetype = 0, show.legend = FALSE, fill = "#D3D3D3")+
  geom_smooth(method = "loess", se = FALSE, level=0.95, linewidth =1.5, show.legend = TRUE)+
  labs(x = "Global AD pathology", y = "Module score") +
  scale_color_manual(values = custom_colors$color) +
  theme(aspect.ratio=1)+
  theme_bw()
######################################

###################

setwd("F:/Brain_region_project/Module_scores_DEGs/Module_scores_PFC427DEGs/Analysis/Analysis_averaged_by_individual/Plots")

pdf("Mic_including_TH_Legend.pdf", width = 4, height = 8)  # Adjust width and height as needed
grid.arrange(gpath, 
             gpath_down,
             ncol = 1)  # Adjust ncol to specify the number of columns
dev.off()

##################################

