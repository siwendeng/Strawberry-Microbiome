#Set working directory
setwd("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/Reletive_Abundance/")

#load R data and librarys
load("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/VESTA_Final_base.RData")
#install.packages("reshape")
library("reshape")
library("ggplot2")
library("plyr")
library("RColorBrewer")
library("colorspace")

#Agglomerate Taxa 
rar_sam_trans  = transform_sample_counts(rar_sam, function(x) 100 * x / sum(x) )
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

#####
#Class level sort based on root taxa abundance
## By Class
others_Class <- names(sort(taxa_sums(subset_samples(glom_Class, SampleType=="root")), decreasing = TRUE)[13:84])
glom_Class <- merge_taxa(glom_Class, others_Class,1)
glom_Class
tax_table(glom_Class)
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
melt_Class$Class <- factor(melt_Class$Class, levels = c("Other", "unclassified", "Anaerolineae", "Verrucomicrobiae", "Opitutae", "TM7PH", "Deltaproteobacteria", "Flavobacteria", "Betaproteobacteria", "Sphingobacteria", "Alphaproteobacteria", "ActinobacteriaPH", "Gammaproteobacteria"))
summary(melt_Class)
#sort sample type
melt_Class$SampleType <- factor(melt_Class$SampleType, levels = c("soil", "rhizosphere","root"))
#plot relative abundance within all sample types
abun_plot_Class <- ggplot(data=melt_Class, aes(x=Treatment, y=Abundance, fill=factor(Class))) +
  geom_bar(stat="identity",position = "fill") +
  facet_grid(SampleType~Timepoint) + 
  scale_fill_manual(values = c("#d9d9d9", "#eff3ff", "#9ecae1", "#6baed6", "#2171b5", "#edf8e9", "#c7e9c0", "#74c476", "#41ab5d", "#ffffb2", "#fed976", "#feb24c", "#f16913")) +
  theme(axis.text.x=element_text(size=12,color="black",angle=90), axis.text.y=element_text(size=12,color="black"), axis.title=element_text(size=16,face="bold"),text=element_text(size=16)) +
  guides(fill=guide_legend(title="Class"))
abun_plot_Class
ggsave(filename = "Relative_abundance_all_sample_types_Class_sortbyroot.pdf",plot=abun_plot_Class, width = 8, height = 6)



## By Phylum
others_Phylum <- names(sort(taxa_sums(glom_Phylum), decreasing = TRUE)[11:38])
glom_Phylum <- merge_taxa(glom_Phylum, others_Phylum,1)
melt_Phylum <- psmelt(glom_Phylum)
#convert Phylum to a character vector from a factor 
str(melt_Phylum$Phylum)
melt_Phylum$Phylum <- as.character(melt_Phylum$Phylum)
#replace all NAs with Other
melt_Phylum$Phylum[is.na(melt_Phylum$Phylum)] <- "Other"
# convert Phylum back to a factor
melt_Phylum$Phylum <- as.factor(melt_Phylum$Phylum)
#sort taxa
melt_Phylum$Phylum <- factor(melt_Phylum$Phylum, levels = c("Other", "TM7", "Verrucomicrobia", "Gemmatimonadetes", "Acidobacteria", "Chloroflexi", "Bacteroidetes", "Firmicutes", "Actinobacteria", "Proteobacteria"))
summary(melt_Phylum)
#sort sample type
melt_Phylum$SampleType <- factor(melt_Phylum$SampleType, levels = c("soil", "rhizosphere","root"))
#plot relative abundance within all sample types
abun_plot_Phylum <- ggplot(data=melt_Phylum, aes(x=Treatment, y=Abundance, fill=factor(Phylum))) +
  geom_bar(stat="identity",position = "fill") +
  facet_wrap(~SampleType) +
  scale_fill_manual(values = c("grey", (brewer.pal(9,"Spectral")))) +
  theme(axis.text.x=element_text(size=12,color="black",angle=90), axis.text.y=element_text(size=12,color="black"), axis.title=element_text(size=16,face="bold"),text=element_text(size=16)) +
  guides(fill=guide_legend(title="Phylum"))
abun_plot_Phylum
ggsave(filename = "Relative_abundance_all_sample_types_Phylum.jpg",plot=abun_plot_Phylum, width = 8, height = 6, dpi = 120)

## By Class
others_Class <- names(sort(taxa_sums(glom_Class), decreasing = TRUE)[16:84])
# others_Class <- names(sort(taxa_sums(subset_samples(glom_Class, SampleType=="root")), decreasing = TRUE)[16:84])
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
melt_Class$Class <- factor(melt_Class$Class, levels = c("Other", "AcidobacteriaPH", "Anaerolineae", "TM7PH", "Verrucomicrobiae", "GemmatimonadetesPH", "Deltaproteobacteria", "Opitutae", "Flavobacteria", "Chloracidobacteria", "Sphingobacteria", "Betaproteobacteria", "ActinobacteriaPH", "Alphaproteobacteria", "Bacilli", "Gammaproteobacteria"))
summary(melt_Class)
#sort sample type
melt_Class$SampleType <- factor(melt_Class$SampleType, levels = c("soil", "rhizosphere","root"))
#plot relative abundance within all sample types
abun_plot_Class <- ggplot(data=melt_Class, aes(x=Treatment, y=Abundance, fill=factor(Class))) +
  geom_bar(stat="identity",position = "fill") +
  facet_wrap(~SampleType) + 
  scale_fill_manual(values = c("grey", (brewer.pal(5,"Blues")), (brewer.pal(5,"Greens")), (brewer.pal(5,"YlOrRd")))) +
  theme(axis.text.x=element_text(size=12,color="black",angle=90), axis.text.y=element_text(size=12,color="black"), axis.title=element_text(size=16,face="bold"),text=element_text(size=16)) +
  guides(fill=guide_legend(title="Class"))
abun_plot_Class
ggsave(filename = "Relative_abundance_all_sample_types_Class.jpg",plot=abun_plot_Class, width = 8, height = 6, dpi = 120)

## By Order
others_Order <- names(sort(taxa_sums(glom_Order), decreasing = TRUE)[16:164])
glom_Order <- merge_taxa(glom_Order, others_Order,1)
melt_Order <- psmelt(glom_Order)
#convert Order to a character vector from a factor 
str(melt_Order$Order)
melt_Order$Order <- as.character(melt_Order$Order)
#replace all NAs with Other
melt_Order$Order[is.na(melt_Order$Order)] <- "Other"
# convert Order back to a factor
melt_Order$Order <- as.factor(melt_Order$Order)
levels(factor(melt_Order$Order))
#sort taxa
melt_Order$Order <- factor(melt_Order$Order, levels = c("Other", "unclassified", "Rhodospirillales", "Sphingomonadales", "Myxococcales", "Flavobacteriales", "Pseudomonadales", "Rhizobiales", "Methylophilales", "Alteromonadales", "Burkholderiales", "Oceanospirillales", "Bacillales", "Sphingobacteriales", "Actinomycetales", "Xanthomonadales"))
summary(melt_Order)
#sort sample type
melt_Order$SampleType <- factor(melt_Order$SampleType, levels = c("soil", "rhizosphere","root"))
#plot relative abundance within all sample types
abun_plot_Order <- ggplot(data=melt_Order, aes(x=Treatment, y=Abundance, fill=factor(Order))) +
  geom_bar(stat="identity",position = "fill") +
  facet_wrap(~SampleType) + 
  scale_fill_manual(values = c("grey", (brewer.pal(5,"Blues")), (brewer.pal(5,"Greens")), (brewer.pal(5,"YlOrRd")))) +
  theme(axis.text.x=element_text(size=12,color="black",angle=90), axis.text.y=element_text(size=12,color="black"), axis.title=element_text(size=16,face="bold"),text=element_text(size=16)) +
  guides(fill=guide_legend(title="Order"))
abun_plot_Order
ggsave(filename = "Relative_abundance_all_sample_types_Order.jpg",plot=abun_plot_Order, width = 8, height = 6, dpi = 120)

## By Family
others_Family <- names(sort(taxa_sums(glom_Family), decreasing = TRUE)[16:301])
# #sort based on root taxa sum
# others_Family <- names(sort(taxa_sums(subset_samples(glom_Family, SampleType=="root")), decreasing = TRUE)[16:301])
glom_Family <- merge_taxa(glom_Family, others_Family,1)
melt_Family <- psmelt(glom_Family)
#convert Family to a character vector from a factor 
str(melt_Family$Family)
melt_Family$Family <- as.character(melt_Family$Family)
#replace all NAs with Other
melt_Family$Family[is.na(melt_Family$Family)] <- "Other"
# convert Family back to a factor
melt_Family$Family <- as.factor(melt_Family$Family)
levels(factor(melt_Family$Family))
#sort taxa
melt_Family$Family <- factor(melt_Family$Family, levels=c("Other", "unclassified", "Bacillaceae", "Paenibacillaceae", "Micrococcaceae", "Sphingomonadaceae", "Sinobacteraceae", "Pseudomonadaceae", "SphingobacterialesOR", "Flexibacteraceae", "Methylophilaceae", "Streptomycetaceae", "Alteromonadaceae", "BurkholderialesOR", "Comamonadaceae", "FCPT525"))
summary(melt_Family)
#sort sample type
melt_Family$SampleType <- factor(melt_Family$SampleType, levels = c("soil", "rhizosphere","root"))
#plot relative abundance within all sample types
abun_plot_Family <- ggplot(data=melt_Family, aes(x=Treatment, y=Abundance, fill=factor(Family))) +
  geom_bar(stat="identity",position = "fill") +
  facet_wrap(~SampleType) + 
  scale_fill_manual(values = c("grey", (brewer.pal(5,"Blues")), (brewer.pal(5,"Greens")), (brewer.pal(5,"YlOrRd")))) +
  theme(axis.text.x=element_text(size=12,color="black",angle=90), axis.text.y=element_text(size=12,color="black"), axis.title=element_text(size=16,face="bold"),text=element_text(size=16)) +
  guides(fill=guide_legend(title="Family"))
abun_plot_Family
ggsave(filename = "Relative_abundance_all_sample_types_Family.jpg",plot=abun_plot_Family, width = 8, height = 6, dpi = 120)

## By Genus
others_Genus <- names(sort(taxa_sums(glom_Genus), decreasing = TRUE)[16:482])
# #sort based on root taxa sum
# others_Genus <- names(sort(taxa_sums(subset_samples(glom_Genus, SampleType=="root")), decreasing = TRUE)[16:482])
glom_Genus <- merge_taxa(glom_Genus, others_Genus,1)
melt_Genus <- psmelt(glom_Genus)
#convert Genus to a character vector from a factor 
str(melt_Genus$Genus)
melt_Genus$Genus <- as.character(melt_Genus$Genus)
#replace all NAs with Other
melt_Genus$Genus[is.na(melt_Genus$Genus)] <- "Other"
# convert Genus back to a factor
melt_Genus$Genus <- as.factor(melt_Genus$Genus)
levels(factor(melt_Genus$Genus))
#sort taxa
melt_Genus$Genus <- factor(melt_Genus$Genus, levels=c("Other", "unclassified", "Bacillus", "Paenibacillus", "Arthrobacter", "SphingobacterialesOR", "BrucellaceaeFA", "Steroidobacter", "OceanospirillalesOR", "Cytophaga", "Pseudomonas", "Flavobacterium", "Methylotenera", "Streptomyces", "Cellvibrio", "FCPT525FA"))
summary(melt_Genus)
#sort sample type
melt_Genus$SampleType <- factor(melt_Genus$SampleType, levels = c("soil", "rhizosphere","root"))
#plot relative abundance within all sample types
abun_plot_Genus <- ggplot(data=melt_Genus, aes(x=Treatment, y=Abundance, fill=factor(Genus))) +
  geom_bar(stat="identity",position = "fill") +
  facet_wrap(~SampleType) + 
  scale_fill_manual(values = c("grey", (brewer.pal(5,"Blues")), (brewer.pal(5,"Greens")), (brewer.pal(5,"YlOrRd")))) +
  theme(axis.text.x=element_text(size=12,color="black",angle=90), axis.text.y=element_text(size=12,color="black"), axis.title=element_text(size=16,face="bold"),text=element_text(size=16)) +
  guides(fill=guide_legend(title="Genus"))
abun_plot_Genus
ggsave(filename = "Relative_abundance_all_sample_types_Genus.jpg",plot=abun_plot_Genus, width = 8, height = 6, dpi = 120)

save.image("Relative_abundance.RData")
