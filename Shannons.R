#Set working directory
setwd("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/Alpha_Diversity/")

#load R data and source itag_diversity scripts
load("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/VESTA_Final_base.RData")
#install.packages("reshape")
library("reshape")
library("ggplot2")

## Alpha Diversity (Shannon's index accounts for both abundance and evenness of the species present.)
Diversity_table <- estimate_richness(rar_sam,  measures=c("Shannon"))
Diversity_table <- merge(Diversity_table, sample_data(rar), by="row.names")
Diversity_table$SampleType <- factor(Diversity_table$SampleType, levels = c("soil", "rhizosphere", "root"))

# Plot
Shannons_plot_density<-ggplot(data=Diversity_table,aes(y=Shannon,x=Treatment,colour=Treatment)) + 
  geom_boxplot() +
  scale_color_manual(values=c("#FFBB00","#20948B")) +
  scale_fill_manual(values=c("Control", "VESTA")) +
  facet_grid(SampleType ~ Timepoint,scales = "free_y",space = "free_y") +
  labs(y = "Alpha Diversity Measures (Shannon)") +
  theme(axis.text.x=element_text(size=12,color="black",angle=90),
        axis.text.y=element_text(size=12,color="black"),
        axis.title=element_text(size=12,face="bold"),
        text=element_text(size=12))
Shannons_plot_density

ggsave(filename = "Shannons_plot_ALL.pdf",plot=Shannons_plot_density, width = 4, height = 6, dpi = 120)

library(agricolae)
#For generating the Tukeys comparison of sample type and treatment
X<-with(Diversity_table,interaction(SampleType,Treatment))
Diversity_table$Timepoint <- as.factor(Diversity_table$Timepoint)
Diversity_table$Replicate <- as.factor(Diversity_table$Replicate)

anova_SampleType_Treatment<-aov(data=Diversity_table,formula = Shannon~X)
summary(anova_SampleType_Treatment)
TukeyHSD(anova_SampleType_Treatment)
tukeyresult<-HSD.test(anova_SampleType_Treatment,"X",group=T)
tukeyresult$groups
#                   trt    means  M
# 1 soil.Control        5.798428  a
# 2 rhizosphere.Control 5.578588 ab
# 3 soil.VESTA          5.241065 bc
# 4 rhizosphere.VESTA   4.907169  c
# 5 root.Control        3.814820  d
# 6 root.VESTA          3.677880  d

#for generating the Anova stats for sample type treatement and timepoint
anova_SampleType_Treatement_Timepoint<-aov(data=Diversity_table,formula = Shannon~SampleType+Treatment+Timepoint+Replicate)
summary(anova_SampleType_Treatement_Timepoint)

#             Df Sum Sq Mean Sq F value   Pr(>F)    
# SampleType    2  82.59   41.30 269.743  < 2e-16 ***
#   Treatment     1   6.68    6.68  43.651 9.54e-10 ***
#   Timepoint     1   3.30    3.30  21.555 8.40e-06 ***
#   Replicate     1   0.16    0.16   1.043    0.309    
# Residuals   128  19.60    0.15                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

save.image("Shannon.RData")
