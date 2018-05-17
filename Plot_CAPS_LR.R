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
CAP_table<-read.table("CAPS_LR_1.txt",header = T,sep="\t")
CAP_table$SampleType <- factor(CAP_table$SampleType, levels = c("Soil", "Rhizo", "Root"))


CAPS_plot_LR <- ggplot(CAP_table, aes(x=Days, y=Variance, colour=SampleType)) +
  geom_point() + 
  #scale_color_manual(values=c("#87492B","#847B14","#478550")) + # Darker
  scale_color_manual(values=c("#CE8764","#E1D337","#9ACAA1")) + # Lighter
  geom_smooth(method=lm,  se=FALSE) 
CAPS_plot_LR
ggsave(filename = "CAPS_LR.pdf",plot=CAPS_plot_LR,width=5,height=3.5)

#linear regression
Soil_sub <- subset(CAP_table, SampleType == "Soil")
Rhizo_sub <- subset(CAP_table, SampleType == "Rhizo")
Root_sub <- subset(CAP_table, SampleType == "Root")

fit_Soil <- lm(Variance ~ Days, data = Soil_sub)
summary(fit_Soil)
fit_Rhizo <- lm(Variance ~ Days, data = Rhizo_sub)
summary(fit_Rhizo)
fit_Root <- lm(Variance ~ Days, data = Root_sub)
summary(fit_Root)
