#To plot the CAPs results for UniFrac and Bray Curtis distances

# Sets the main working directory containing your OTU table and mapping file.
setwd("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/Beta_Diversity/")

#load R data and source itag_diversity scripts
load("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/VESTA_Final_base.RData")
library(reshape)
library(vegan)
library(ggplot2)

#To perform the Canonical Analysis of Principal Coordinates:
#NOTE: The "inertia" in a data set is analogous to the variance.

#### 
CAP_table<-read.table("CAPS_table_Timepoint.txt",header = T,sep="\t")
CAP_table$Timepoint<-factor(CAP_table$Timepoint,levels=c("1","2","3","4"))
CAP_table$Factor<-factor(CAP_table$Factor,levels=c("SampleType","Treatment","Replicate"))

#Final method split by distance type 
CAP_table_bray<-subset(CAP_table,Distance=="Bray")
CAP_table_unifrac<-subset(CAP_table,Distance=="UniFrac")

CAPS_plot_bray<-ggplot(data=CAP_table_bray,aes(x=Factor,y=PercVariance,fill=log10(pvalue))) + 
  geom_bar(stat="identity")+
  facet_grid(~Timepoint) +
  theme(axis.text.x=element_text(size=16,color="black",angle=90), axis.text.y=element_text(size=16,color="black"), axis.title=element_text(size=16,face="bold"),text=element_text(size=16)) +
  scale_fill_gradient(low="black", high="grey")
CAPS_plot_bray
ggsave(filename = "CAPs_plot_bray_byTP.pdf",plot=CAPS_plot_bray,width=4.5,height=5.5)

CAPS_plot_unifrac<-ggplot(data=CAP_table_unifrac,aes(x=Factor,y=PercVariance,fill=log10(pvalue))) + 
  geom_bar(stat="identity")+
  facet_grid(~Timepoint) +
  theme(axis.text.x=element_text(size=16,color="black",angle=90), axis.text.y=element_text(size=16,color="black"), axis.title=element_text(size=16,face="bold"),text=element_text(size=16)) +
  scale_fill_gradient(low="black", high="grey")
CAPS_plot_unifrac
ggsave(filename = "CAPs_plot_wUF_byTP.pdf",plot=CAPS_plot_unifrac,width=4.5,height=5.5)
