#Set working directory
setwd("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/Alpha_Diversity/")

#load R data and source itag_diversity scripts
load("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/VESTA_Final_base.RData")
#install.packages("reshape")
library("reshape")
library("ggplot2")

## Alpha Diversity (Observed's index only accounts for abundance)
Diversity_table <- estimate_richness(rar_sam,  measures=c("Observed"))
Diversity_table <- merge(Diversity_table, sample_data(rar), by="row.names")
Diversity_table$SampleType <- factor(Diversity_table$SampleType, levels = c("soil", "rhizosphere", "root"))

# Plot
Observed_plot_density<-ggplot(data=Diversity_table,aes(y=Observed,x=Treatment,colour=Treatment)) + 
  geom_boxplot() +
  scale_color_manual(values=c("#FFBB00","#20948B")) +
  scale_fill_manual(values=c("Control", "VESTA")) +
  facet_grid(SampleType ~ Timepoint, scales = "free_y", space = "free_y") +
  labs(y = "Alpha Diversity Measures (Observed)") +
  theme(axis.text.x=element_text(size=12,color="black",angle=90), 
        axis.text.y=element_text(size=12,color="black"), 
        axis.title=element_text(size=12,face="bold"),
        text=element_text(size=12))
Observed_plot_density

ggsave(filename = "Observed_plot_ALL.pdf",plot=Observed_plot_density, width = 4, height = 6, dpi = 120)

library(agricolae)
#For generating the Tukeys comparison of sample type and treatment
X<-with(Diversity_table,interaction(SampleType,Treatment))
anova_SampleType_Treatment<-aov(data=Diversity_table,formula = Observed~X)
summary(anova_SampleType_Treatment)
TukeyHSD(anova_SampleType_Treatment)
tukeyresult<-HSD.test(anova_SampleType_Treatment,"X",group=T)
tukeyresult$groups
#                   trt    means  M
# 1 soil.Control        1635.7368 a
# 2 rhizosphere.Control 1401.9167 b
# 3 soil.VESTA          1218.2500 c
# 4 rhizosphere.VESTA   1048.0000 d
# 5 root.Control         535.0833 e
# 6 root.VESTA           419.7083 f

#for generating the Anova stats for sample type treatement and timepoint
anova_SampleType_Treatement_Timepoint<-aov(data=Diversity_table,formula = Observed~SampleType+Treatment+Timepoint+Replicate)
summary(anova_SampleType_Treatement_Timepoint)

#               Df Sum Sq Mean Sq F value   Pr(>F)    
#   SampleType    2 22468574 11234287 824.926  < 2e-16 ***
#   Treatment     1  2757147  2757147 202.455  < 2e-16 ***
#   Timepoint     1   343425   343425  25.217 1.68e-06 ***
#   Replicate     1       96       96   0.007    0.933    
# Residuals   128  1743173    13619                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

save.image("Observed.RData")
