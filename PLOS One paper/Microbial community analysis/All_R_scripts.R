#2D NMDS ANALYSIS

library(vegan)
pc = read.csv("OTU_NMDS.txt", header= TRUE, sep="\t")

#make community matrix - extract columns with abundance information
com = pc[,8:ncol(pc)]

#turn abundance data frame into a matrix
m_com = as.matrix(com)

set.seed(123)
nmds = metaMDS(m_com, distance = "bray")

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))

#add columns to data frame 
data.scores$Sample_ID = pc$Sample_ID
data.scores$Site = pc$Site
data.scores$Dist = pc$Dist
data.scores$Day = pc$Day
data.scores$Spot = pc$Spot
data.scores$Depth = pc$Depth
data.scores$Sub_flow = pc$Sub_flow
 
head(data.scores)

library(ggplot2)
library(viridis)

plot_2D = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(size = 4, aes(colour = Sub_flow))+ 
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
    axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
    legend.text = element_text(size = 12, face ="bold", colour ="black"), 
    legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
    axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
    legend.title = element_text(size = 14, colour = "black", face = "bold"), 
    panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
    legend.key=element_blank()) + 
    labs(x = "NMDS1", colour = "Flow", y = "NMDS2") 
 
plot_2D

#ANOSIM ANALYSIS

#Find unique groups in data
groups = unique(pc[c("Sub_flow")])
group_list = as.list(groups$Sub_flow)

#Calculate pairwise ANOSIM statistics
for (i in group_list) {
    for (j in group_list) {
    if (i==j) next
    sub_df1 = subset(pc, subset=(Sub_flow==i))
    sub_df2 = subset(pc, subset=(Sub_flow==j))
    sub_df = rbind(sub_df1, sub_df2)
    nums = sub_df[,8:ncol(sub_df)]
    nums_matrix = as.matrix(nums)
    set.seed(123)
    ano = anosim(nums_matrix, sub_df$Sub_flow, distance = "bray", permutations = 9999)
    Rvalue = ano$statistic
    pvalue = ano$signif
    cat(i," vs ",j ,": R-value: ", Rvalue, " P-value: ", pvalue, '\n')

    }
}

#3D NMDS ANALYSIS

set.seed(123)
nmds = metaMDS(m_com, k=3, distance = "bray")

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))

#add columns to data frame 
data.scores$Sample_ID = pc$Sample_ID
data.scores$Site = pc$Site
data.scores$Dist = pc$Dist
data.scores$Day = pc$Day
data.scores$Spot = pc$Spot
data.scores$Depth = pc$Depth
data.scores$Sub_flow = pc$Sub_flow

#Load up libraries
library(plotly)

#Generate plot
plot_3D <- plot_ly(data.scores, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, color = ~Sub_flow, marker = list(size = 5), text = ~paste('<br>Site:', Site, '<br>Distance from inflow:', Dist, '<br>Exposure:', Sub_flow)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'NMDS1'),
                     yaxis = list(title = 'NMDS2'),
                     zaxis = list(title = 'NMDS3')))
#Show plot
plot_3D
