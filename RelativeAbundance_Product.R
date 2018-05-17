#Set working directory
setwd("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/Reletive_Abundance/")

#load R data and librarys
load("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/VESTA_Final_base.RData")
#install.packages("reshape")
library("reshape")
library("ggplot2")
library("plyr")
library("RColorBrewer")

#Agglomerate Taxa 
rar_sam_trans  = transform_sample_counts(rar_pro, function(x) 100 * x / sum(x) )
glom_Phylum <- tax_glom(rar_sam_trans, taxrank = "Phylum")
glom_Class <- tax_glom(rar_sam_trans, taxrank = "Class")
glom_Order <- tax_glom(rar_sam_trans, taxrank = "Order")
glom_Family <- tax_glom(rar_sam_trans, taxrank = "Family")
glom_Genus <- tax_glom(rar_sam_trans, taxrank = "Genus")
glom_Phylum #38 Phylums
glom_Class #84 Classes
glom_Order #164 Orders
glom_Family #301 Families
glom_Genus #482 Genuses

#Shotgun
melt_Class1 <- read.csv(file = "Product_16s_shotgun.txt", sep = "\t", header = TRUE)
#sort sample
melt_Class$Sample <- factor(melt_Class$Sample, levels = c("BHF", "SOBEC", "VESTA1", "VESTA2", "VESTA3", "VESTA4", "VESTA5", "VESTA6", "VESTA7", "VESTA8", "VESTA9", "VESTA10", "VESTA12", "VESTA13", "VESTA14"))
#plot relative abundance within all sample types
abun_plot_Class <- ggplot(data=melt_Class1, aes(x=Type, y=Abundance, fill=factor(Class))) +
  geom_bar(stat="identity",position = "fill") +
  scale_fill_manual(values = c("#d9d9d9", "#eff3ff", "#9ecae1", "#6baed6", "#2171b5", "#edf8e9", "#c7e9c0", "#74c476", "#41ab5d", "#ffffb2", "#FEEAB2", "#fed976", "#feb24c", "#f16913")) +
  theme(axis.text.x=element_text(size=12,color="black",angle=45,hjust=1), axis.text.y=element_text(size=12,color="black"), axis.title=element_text(size=16,face="bold"),text=element_text(size=16)) +
  guides(fill=guide_legend(title="Class"))
abun_plot_Class
ggsave(filename = "Relative_abundance_shot.pdf",plot=abun_plot_Class, width = 5, height = 6, dpi = 120)
write.csv(melt_Class, file = "melt_class.csv")


## By Class
others_Class <- names(sort(taxa_sums(glom_Class), decreasing = TRUE)[13:84])
glom_Class <- merge_taxa(glom_Class, others_Class,1)
melt_Class <- psmelt(glom_Class)
#convert Class to a character vector from a factor 
str(melt_Class$Class)
melt_Class$Class <- as.character(melt_Class$Class)
#replace all NAs with Other
melt_Class$Class[is.na(melt_Class$Class)] <- "Other"
# convert Class back to a factor
melt_Class$Class <- as.factor(melt_Class$Class)
levels(factor(melt_Class$Class))
#sort taxa
melt_Class$Class <- factor(melt_Class$Class, levels = c("Other", "GemmatimonadetesPH", "Deltaproteobacteria", "Flavobacteria", "Bacilli", "AcidobacteriaPH", "ActinobacteriaPH", "TM7PH", "SR1PH", "Betaproteobacteria", "Gammaproteobacteria", "Sphingobacteria", "Alphaproteobacteria"))
summary(melt_Class)
#sort sample
melt_Class$Sample <- factor(melt_Class$Sample, levels = c("BHF", "SOBEC", "VESTA1", "VESTA2", "VESTA3", "VESTA4", "VESTA5", "VESTA6", "VESTA7", "VESTA8", "VESTA9", "VESTA10", "VESTA12", "VESTA13", "VESTA14"))
#plot relative abundance within all sample types
abun_plot_Class <- ggplot(data=melt_Class, aes(x=Sample, y=Abundance, fill=factor(Class))) +
  geom_bar(stat="identity",position = "fill") +
  scale_fill_manual(values = c("#d9d9d9", "#eff3ff", "#9ecae1", "#6baed6", "#2171b5", "#edf8e9", "#c7e9c0", "#74c476", "#41ab5d", "#ffffb2", "#fed976", "#feb24c", "#f16913")) +
  theme(axis.text.x=element_text(size=12,color="black",angle=45,hjust=1), axis.text.y=element_text(size=12,color="black"), axis.title=element_text(size=16,face="bold"),text=element_text(size=16)) +
  guides(fill=guide_legend(title="Class"))
abun_plot_Class
ggsave(filename = "Relative_abundance_all_Product.pdf",plot=abun_plot_Class, width = 8, height = 6, dpi = 120)


rar_pro_actinos <- subset_taxa(rar_pro, Class=="Alphaproteobacteria")
glom_Family <- tax_glom(rar_pro_actinos, taxrank = "Family")
