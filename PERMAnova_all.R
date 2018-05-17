#Set working directory
setwd("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/Beta_Diversity/")

#load R data and source itag_diversity scripts
load("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/VESTA_Final_base.RData")
library(reshape)
library(vegan)

#Root
Root_OTU <- as.data.frame(t(as.data.frame(otu_table(rar_root))))
env_current <- data.frame(sample_data(rar_root))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#Root_dist<-vegdist(Root_OTU,"bray")
adonis(Root_OTU ~ Treatment*Timepoint,data=env_current,permutations = 9999)

#Soil
Soil_OTU <- as.data.frame(t(as.data.frame(otu_table(rar_soil))))
env_current <- data.frame(sample_data(rar_soil))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#Soil_dist<-vegdist(Soil_OTU,"bray")
adonis(Soil_OTU ~ Treatment*Timepoint,data=env_current,permutations = 9999)

#Rhizosphere
Rhizo_OTU <- as.data.frame(t(as.data.frame(otu_table(rar_rhizo))))
env_current <- data.frame(sample_data(rar_rhizo))
env_current$Replicate = as.factor(env_current$Replicate)
env_current$Timepoint = as.factor(env_current$Timepoint)
#Rhizo_dist<-vegdist(Rhizo_OTU,"bray")
adonis(Rhizo_OTU ~ Treatment*Timepoint,data=env_current,permutations = 9999)