#Set working directory
setwd("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/Beta_Diversity/")

#load R data and source itag_diversity scripts
load("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/VESTA_Final_base.RData")
#install.packages("reshape")
library("reshape")
library("ggplot2")
library("RColorBrewer")

# TreatmentByTimepoint <- as.factor(with(sample_data(rar_sam),interaction(Treatment,Timepoint)))
# sample_data(rar_sam)$TreatmentByTimepoint <- TreatmentByTimepoint
# sample_data(rar_sam) <- as.data.frame(sample_data(rar_sam))

SampleByTreatment <- as.factor(with(sample_data(rar_sam),interaction(SampleType,Treatment)))
sample_data(rar_sam)$SampleByTreatment <- SampleByTreatment
sample_data(rar_sam) <- as.data.frame(sample_data(rar_sam))

SampleByTimepoint <- as.factor(with(sample_data(rar_sam),interaction(SampleType,Timepoint)))
sample_data(rar_sam)$SampleByTimepoint <- SampleByTimepoint
sample_data(rar_sam) <- as.data.frame(sample_data(rar_sam))


## Beta Diversity
orduniw = ordinate(rar_sam, "PCoA", "unifrac", weighted=TRUE)
# ordununiw = ordinate(rar_sam, "PCoA", "unifrac", weighted=FALSE)
ordbray = ordinate(rar_sam, "PCoA", "bray")

#sort sample type
sample_data(rar_sam)$SampleType <- factor(sample_data(rar_sam)$SampleType, levels = c("soil", "rhizosphere","root"))
sample_data(rar_sam)$SampleByTimepoint <- factor(sample_data(rar_sam)$SampleByTimepoint, levels = c("soil.1", "soil.2","soil.3","soil.4","rhizosphere.1","rhizosphere.2", "rhizosphere.3", "rhizosphere.4","root.1", "root.2","root.3","root.4"))
# sample_data(rar_sam)$TreatmentByTimepoint <- factor(sample_data(rar_sam)$TreatmentByTimepoint, levels = c("Control.1", "Control.2","Control.3","Control.4","VESTA.1","VESTA.2", "VESTA.3", "VESTA.4"))
sample_data(rar_sam)$SampleByTreatment <- factor(sample_data(rar_sam)$SampleByTreatment, levels = c("soil.Control", "soil.VESTA", "rhizosphere.Control", "rhizosphere.VESTA", "root.Control", "root.VESTA"))

## Plot
#unifrac weighted
unifrac_wei_plot<- plot_ordination(rar_sam, orduniw, color="SampleByTimepoint", shape="Treatment") + 
  geom_point(colour="black",size=3.5) +
  geom_point(size=3) +
  scale_color_manual(values=c("#E6C3B2", "#D49576", "#C1693D", "#87492B", "#EEE690", "#E4D84E", "#C6B81E", "#847B14", "#BCDCC1","#89C191", "#58A663", "#3E7445")) +
  scale_fill_manual(values=c("soil.1", "soil.2","soil.3","soil.4","rhizosphere.1","rhizosphere.2", "rhizosphere.3", "rhizosphere.4","root.1", "root.2","root.3","root.4")) +
  ggtitle("PCoA on weighted-UniFrac distance") + 
  theme(axis.text.x=element_text(size=12,color="black",angle=90),
        axis.text.y=element_text(size=12,color="black"),
        axis.title=element_text(size=12,face="bold"),
        text=element_text(size=12))

unifrac_wei_plot
ggsave(filename = "PCoA_AllSamples_unifrac_weighted.pdf",plot = unifrac_wei_plot,width=6,height=5)

# #unifrac unweighted
# unifrac_unwei_plot<- plot_ordination(rar_sam, ordununiw, color="Treatment", shape="SampleType") + 
#   geom_point(size=3) +
#   scale_color_manual(values=c("#FFBB00","#20948B")) +
#   scale_fill_manual(values=c("Control", "VESTA")) +
#   ggtitle("PCoA on unweighted-UniFrac distance") + 
#   theme(axis.text.x=element_text(size=12,color="black",angle=90), axis.text.y=element_text(size=12,color="black"), axis.title=element_text(size=12,face="bold"),text=element_text(size=12))
# unifrac_unwei_plot
# ggsave(filename = "PCoA_AllSamples_unifrac_unweighted.jpg",plot = unifrac_unwei_plot,width=8,height=6)

#bay-curtis
bray_plot<- plot_ordination(rar_sam, ordbray, color="SampleByTreatment") + 
  geom_point(colour="black",size=3.5) +
  geom_point(size=3) +
  scale_color_manual(values=c("#CE8764","#87492B","#E1D337","#847B14","#9ACAA1","#478550")) +
  scale_fill_manual(values=c("Control", "VESTA")) +
  ggtitle("PCoA on bay-curtis distance") + 
  theme(axis.text.x=element_text(size=12,color="black",angle=90),
        axis.text.y=element_text(size=12,color="black"),
        axis.title=element_text(size=12,face="bold"),
        text=element_text(size=12))
bray_plot
ggsave(filename = "PCoA_AllSamples_bray.pdf",plot = bray_plot,width=6,height=5)

save.image("PCoA.RData")
