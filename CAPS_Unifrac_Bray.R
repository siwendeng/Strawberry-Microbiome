#This script runs the CAPS analysis of the Bray-Curtis and Unifrac distances.
#Set working directory
setwd("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/Beta_Diversity/")

#load R data and source itag_diversity scripts
load("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/VESTA_Final_base.RData")
library(reshape)
library(vegan)

#get Bray-Curtis and Unifrac distance from Phyloseq

BrayDist =phyloseq::distance(rar_sam, "bray")
UnifracDist = phyloseq::distance(rar_sam,"wunifrac")
env_current <- data.frame(sample_data(rar_sam))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)

## CAPs analysis of Bray and UNIFRAC distance for all samples
CAPS_Unifrac<-capscale(UnifracDist~SampleType+Treatment+Timepoint+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
#UniFrac
# Model: capscale(formula = UnifracDist ~ SampleType + Treatment + Timepoint + Replicate, data = env_current)
#              Df SumOfSqs        F Pr(>F)    
# SampleType   2  2.88441 108.4102  0.001 ***
#   Treatment    1  0.54363  40.8648  0.001 ***
#   Timepoint    3  0.29536   7.4008  0.001 ***
#   Replicate    5  0.06521   0.9804  0.438    
# Residual   122  1.62299         

CAPS_Bray<-capscale(BrayDist~SampleType+Treatment+Timepoint+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")
#Bray
# Model: capscale(formula = BrayDist ~ SampleType + Treatment + Timepoint + Replicate, data = env_current)
#              Df SumOfSqs       F Pr(>F)    
# SampleType   2  13.5589 68.5415  0.001 ***
#   Treatment    1   4.0932 41.3828  0.001 ***
#   Timepoint    3   1.6578  5.5868  0.001 ***
#   Replicate    5   0.4650  0.9402  0.516    
# Residual   122  12.0670     

##----------------------------
#For soil
BrayDist =phyloseq::distance(rar_soil, "bray")
UnifracDist = phyloseq::distance(rar_soil,"wunifrac")
env_current <- data.frame(sample_data(rar_soil))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Timepoint+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Timepoint+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")

#For rhizosphere
BrayDist =phyloseq::distance(rar_rhizo, "bray")
UnifracDist = phyloseq::distance(rar_rhizo,"wunifrac")
env_current <- data.frame(sample_data(rar_rhizo))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Timepoint+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Timepoint+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")

#For root
BrayDist =phyloseq::distance(rar_root, "bray")
UnifracDist = phyloseq::distance(rar_root,"wunifrac")
env_current <- data.frame(sample_data(rar_root))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Timepoint+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Timepoint+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")

##------------------------------
#For TP1
BrayDist =phyloseq::distance(rar_1, "bray")
UnifracDist = phyloseq::distance(rar_1,"wunifrac")
env_current <- data.frame(sample_data(rar_1))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
env_current$SampleType = as.factor(env_current$SampleType)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~SampleType+Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~SampleType+Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")

#For TP2
BrayDist =phyloseq::distance(rar_2, "bray")
UnifracDist = phyloseq::distance(rar_2,"wunifrac")
env_current <- data.frame(sample_data(rar_2))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~SampleType+Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~SampleType+Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")

#For TP3
BrayDist =phyloseq::distance(rar_3, "bray")
UnifracDist = phyloseq::distance(rar_3,"wunifrac")
env_current <- data.frame(sample_data(rar_3))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~SampleType+Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~SampleType+Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")

#For TP4
BrayDist =phyloseq::distance(rar_4, "bray")
UnifracDist = phyloseq::distance(rar_4,"wunifrac")
env_current <- data.frame(sample_data(rar_4))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~SampleType+Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~SampleType+Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")

##--------------------
##subset soil
rar_soil_1 = subset_samples(rar_soil, Timepoint == "1")
rar_soil_2 = subset_samples(rar_soil, Timepoint == "2")
rar_soil_3 = subset_samples(rar_soil, Timepoint == "3")
rar_soil_4 = subset_samples(rar_soil, Timepoint == "4")

##subset rhizo
rar_rhizo_1 = subset_samples(rar_rhizo, Timepoint == "1")
rar_rhizo_2 = subset_samples(rar_rhizo, Timepoint == "2")
rar_rhizo_3 = subset_samples(rar_rhizo, Timepoint == "3")
rar_rhizo_4 = subset_samples(rar_rhizo, Timepoint == "4")

##subset root
rar_root_1 = subset_samples(rar_root, Timepoint == "1")
rar_root_2 = subset_samples(rar_root, Timepoint == "2")
rar_root_3 = subset_samples(rar_root, Timepoint == "3")
rar_root_4 = subset_samples(rar_root, Timepoint == "4")

#For soil 1
BrayDist =phyloseq::distance(rar_soil_1, "bray")
UnifracDist = phyloseq::distance(rar_soil_1,"wunifrac")
env_current <- data.frame(sample_data(rar_soil_1))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")
#For soil 2
BrayDist =phyloseq::distance(rar_soil_2, "bray")
UnifracDist = phyloseq::distance(rar_soil_2,"wunifrac")
env_current <- data.frame(sample_data(rar_soil_2))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")
#For soil 3
BrayDist =phyloseq::distance(rar_soil_3, "bray")
UnifracDist = phyloseq::distance(rar_soil_3,"wunifrac")
env_current <- data.frame(sample_data(rar_soil_3))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")
#For soil 4
BrayDist =phyloseq::distance(rar_soil_4, "bray")
UnifracDist = phyloseq::distance(rar_soil_4,"wunifrac")
env_current <- data.frame(sample_data(rar_soil_4))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")

#For rhizo 1
BrayDist =phyloseq::distance(rar_rhizo_1, "bray")
UnifracDist = phyloseq::distance(rar_rhizo_1,"wunifrac")
env_current <- data.frame(sample_data(rar_rhizo_1))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")
#For rhizo 2
BrayDist =phyloseq::distance(rar_rhizo_2, "bray")
UnifracDist = phyloseq::distance(rar_rhizo_2,"wunifrac")
env_current <- data.frame(sample_data(rar_rhizo_2))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")
#For rhizo 3
BrayDist =phyloseq::distance(rar_rhizo_3, "bray")
UnifracDist = phyloseq::distance(rar_rhizo_3,"wunifrac")
env_current <- data.frame(sample_data(rar_rhizo_3))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")
#For rhizo 4
BrayDist =phyloseq::distance(rar_rhizo_4, "bray")
UnifracDist = phyloseq::distance(rar_rhizo_4,"wunifrac")
env_current <- data.frame(sample_data(rar_rhizo_4))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")

#For root 1
BrayDist =phyloseq::distance(rar_root_1, "bray")
UnifracDist = phyloseq::distance(rar_root_1,"wunifrac")
env_current <- data.frame(sample_data(rar_root_1))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")
#For root 2
BrayDist =phyloseq::distance(rar_root_2, "bray")
UnifracDist = phyloseq::distance(rar_root_2,"wunifrac")
env_current <- data.frame(sample_data(rar_root_2))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")
#For root 3
BrayDist =phyloseq::distance(rar_root_3, "bray")
UnifracDist = phyloseq::distance(rar_root_3,"wunifrac")
env_current <- data.frame(sample_data(rar_root_3))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")
#For root 4
BrayDist =phyloseq::distance(rar_root_4, "bray")
UnifracDist = phyloseq::distance(rar_root_4,"wunifrac")
env_current <- data.frame(sample_data(rar_root_4))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#CAPS
CAPS_Unifrac<-capscale(UnifracDist~Treatment+Replicate,data=env_current)
anova(CAPS_Unifrac,by="terms")
CAPS_Bray<-capscale(BrayDist~Treatment+Replicate,data=env_current)
anova(CAPS_Bray,by="terms")
