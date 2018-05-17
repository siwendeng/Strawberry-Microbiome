#Set working directory
setwd("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/Indicator/")

#load R data and source itag_diversity scripts
load("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/VESTA_Final_base.RData")
#install.packages("reshape")
library("reshape")
library("ggplot2")
library("VennDiagram")
library("gplots")
library("indicspecies")
library("data.table")


#remove otus with no read count
rar_root <- prune_taxa(taxa_sums(rar_root) > 0, rar_root)
rar_vesta_sub <- subset_samples(rar_pro, Timepoint == "2")
rar_vesta_sub <- prune_taxa(taxa_sums(rar_vesta_sub) > 0, rar_vesta_sub)
rar_vesta_all <- prune_taxa(taxa_sums(rar_pro) > 0, rar_pro)
#transpose otu table of root subset
root_otu = as.data.frame(t(otu_table(rar_root)))
#get metadata for root subset
root_meta = as.data.frame(sample_data(rar_root))
#indicator species analysis
n = ncol(root_meta)
N = n+1
b = cbind(root_meta,root_otu)
wetpt = multipatt(b[,N:ncol(b)], b$Treatment, duleg = TRUE, control = how(nperm=999)) #a list of control values describing properties of the permutation design, as returned by a call to how
write.csv(wetpt$sign, "Indicator.csv")
mergeindi = cbind(wetpt$sign,tax_table(rar_root),otu_table(rar_root))
mergeindi = subset(mergeindi,p.value<=0.01)
# write.csv(mergeindi, "merge_Indicator.csv")

#Indicator for control and VESTA enriched
Control = subset(mergeindi,s.Control == 1)
#write.csv(Control, "Control_merge_Indicator.csv")
Control_enrich = row.names(Control)
VESTA = subset(mergeindi,s.VESTA == 1)
#write.csv(VESTA, "VESTA_merge_Indicator.csv")
VESTA_enrich = row.names(VESTA)

#VESTA
rar_vesta_sub <- row.names(otu_table(rar_vesta_sub))
rar_vesta_all <- row.names(otu_table(rar_vesta_all))
#Venn all
venn_all <- venn(list(Control_enrich = Control_enrich, VESTA_enrich = VESTA_enrich, VESTA_Pro = rar_vesta_sub))

#ggsave(filename = "venn_all.pdf", venn_all, width = 8, height = 6)

str(venn_all)
VESTA_enrich_VESTA_Pro <- data.frame(attr(venn_all, "intersections")[3])
VESTA_enrich_VESTA_Pro$Type <- "VESTA_enrich_VESTA_Pro"
VESTA_enrich_VESTA_Pro$alpha <- "1"
colnames(VESTA_enrich_VESTA_Pro)[1] <- "OTU_ID"
Control_enrich_VESTA_Pro <- data.frame(attr(venn_all, "intersections")[4])
Control_enrich_VESTA_Pro$Type <- "Control_enrich_VESTA_Pro"
Control_enrich_VESTA_Pro$alpha <- "1"
colnames(Control_enrich_VESTA_Pro)[1] <- "OTU_ID"
VESTA_Pro <- data.frame(attr(venn_all, "intersections")[5])
VESTA_Pro$Type <- "VESTA_Pro"
VESTA_Pro$alpha <- "0.05"
colnames(VESTA_Pro)[1] <- "OTU_ID"
str(VESTA_enrich_VESTA_Pro)
str(Control_enrich_VESTA_Pro)
str(VESTA_Pro)
all <- rbind(VESTA_enrich_VESTA_Pro, Control_enrich_VESTA_Pro, VESTA_Pro)
str(all)

#VESTA OTU
rar_vesta_sub <- subset_samples(rar_pro, Timepoint == "2")
rar_vesta_sub <- prune_taxa(taxa_sums(rar_vesta_sub) > 0, rar_vesta_sub)
VESTA_OTU <- data.frame(otu_table(rar_vesta_sub))
str(VESTA_OTU)
colnames(VESTA_OTU) <- "Abundance"
setDT(VESTA_OTU, keep.rownames = TRUE)[]
colnames(VESTA_OTU)[1] <- "OTU_ID"
rownames(VESTA_OTU) <- NULL

#Merge
nrow(all) == nrow(VESTA_OTU)
Final <- merge(all, VESTA_OTU, by="OTU_ID")
Final$Type <- factor(Final$Type, levels = c("VESTA_Pro", "Control_enrich_VESTA_Pro", "VESTA_enrich_VESTA_Pro"))
levels(Final$Type)
write.csv(Final, "Final_Rand_Abundance_new_new.csv")

#Plot
p <- ggplot(Final, aes(x=reorder(OTU_ID, -Abundance), y= Abundance)) + 
  geom_point(aes(colour = Type, alpha = alpha)) +
  scale_color_manual(values = c("#8EB4E3", "#FFBB00","#20948B")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_line(),
        axis.line.y=element_line()) +
  geom_vline(xintercept = 382, colour = "#FFBB00",linetype="dotted") +
  geom_vline(xintercept = 262, colour = "#20948B", linetype="dotted") +
  geom_vline(xintercept = 305, colour = "#8EB4E3", linetype="dotted")
p
ggsave(filename = "rank_abundance.pdf",plot=p,width=5.5,height=2.5)
