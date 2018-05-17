#Set working directory
setwd("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/Reletive_Abundance/")

#load R data and librarys
load("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/VESTA_Final_base.RData")
#install.packages("reshape")
library("reshape")
library("ggplot2")
library("plyr")
library("RColorBrewer")

#calculate distance matrix of all samples
d = distance(rar_pro, method='bray')

#melt matrix to data frame
melt_d <- as.data.frame(melt(as.matrix(d)))
melt_d$X1 <- factor(melt_d$X1, levels = c("BHF", "SOBEC", "VESTA1", "VESTA2", "VESTA3", "VESTA4", "VESTA5", "VESTA6", "VESTA7", "VESTA8", "VESTA9", "VESTA10", "VESTA12", "VESTA13", "VESTA14"))
melt_d$X2 <- factor(melt_d$X2, levels = c("BHF", "SOBEC", "VESTA1", "VESTA2", "VESTA3", "VESTA4", "VESTA5", "VESTA6", "VESTA7", "VESTA8", "VESTA9", "VESTA10", "VESTA12", "VESTA13", "VESTA14"))
factor(melt_d$X2)

#draw heatmap
heatmap <- ggplot(data = melt_d, aes(x=X1, y=X2, fill=value)) + 
  geom_tile() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_fill_gradient(low = "steelblue", high = "white") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(size = 12)) +
  labs(fill = "Distance")
heatmap
ggsave(filename = "VESTA_Product_heatmap.pdf",plot=heatmap, width = 7, height = 6, dpi = 120)

