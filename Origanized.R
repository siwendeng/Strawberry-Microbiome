library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
#Define a default theme for ggplot graphics.
theme_set(theme_bw())
library("scales")
library("grid")
library("DESeq2")
library("ape")
vignette("phyloseq_basics")
vignette("phyloseq_analysis")

setwd("~/Documents/Colmen-DerrLab/Projects/VESTA/VESTA_Strawberry/Redo_Data_analysis/itagger/R/Phyloseq/")

##import biom and map file
biom_file = "st_otu_table.biom"
map_file = "VESTA_Strawberry_metadata_Three_VESTAs.txt"
tree= read.tree("st_otu_itags.phy")
otutax = import_biom(biom_file,parseFunction = parse_taxonomy_greengenes)
sam = import_qiime_sample_data(map_file)

#use merge_phyloseq to build up the vesta_st object
vesta_st = merge_phyloseq(otutax, sam, tree)
vesta_st

#Take a look of the data
ntaxa(vesta_st)
taxa_names(vesta_st)[1:10]
nsamples(vesta_st)
sample_names(vesta_st)
rank_names(vesta_st)
sample_variables(vesta_st)
tax_table(vesta_st)[1:5, 1:4]
phy_tree(vesta_st)

##randomly root your tree -- or you could put an outgroup in and then trim it after tree is rooted
tree_Q = root(tree_Q, 1, resolve.root = T)

###Merge
Bushman = merge_phyloseq(biomot, tree_Q,bmsd)

Bu_rare=rarefy_even_depth(Bushman, sample.size = min(sample_sums(Bushman)),
                  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

###alpha_diversity
plot_richness(Bu_rare,x="SampleType") + geom_boxplot()

##Rank level
rank_names(Bushman)
get_taxa_unique(Bushman, "Phylum")
colnames(tax_table(Bushman)) <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")


##Normalize with Treatment
diagdds = phyloseq_to_deseq2(Bushman, ~ Treatment)
diagdds = estimateSizeFactors(diagdds)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)
dim(diagvst)

Bushman0 = Bushman
otu_table(Bushman0) <- otu_table(diagvst, taxa_are_rows = TRUE)
Bushman0

##Normalize with SampleType
diagdds2 = phyloseq_to_deseq2(Bushman, ~ SampleType)
diagdds2 = estimateSizeFactors(diagdds2)
diagdds2 = estimateDispersions(diagdds2)
diagvst2 = getVarianceStabilizedData(diagdds2)
dim(diagvst2)

Bushman1 = Bushman
otu_table(Bushman1) <- otu_table(diagvst2, taxa_are_rows = TRUE)
Bushman1

##Top 10 abundance phylum and Genus,full 10,SampleType
Bushman_genus<-tax_glom(Bushman1, taxrank="Phylum")
Phylum10 = names(sort(taxa_sums(Bushman_genus), TRUE)[1:10])
ent11 = prune_taxa(Phylum10, Bushman1)
plot_bar(ent11, x="SampleType",y="Abundance", facet_grid = ~Phylum,fill = "Genus")

##Top 10 abundance phylum and Genus,full 10,Treatment
Bushman_genus<-tax_glom(Bushman0, taxrank="Phylum")
Phylum10 = names(sort(taxa_sums(Bushman_genus), TRUE)[1:10])
ent10 = prune_taxa(Phylum10, Bushman0)
plot_bar(ent10, x="Treatment",y="Abundance", facet_grid = ~Phylum,fill = "Genus")

###heatmap for treatment
plot_heatmap(ent10, taxa.order = "Phylum", taxa.label = "Genus", sample.label = "Treatment", sample.order = "Treatment")
plot_heatmap(ent10, taxa.order = "Kingdom", taxa.label = "Phylum", sample.label = "Treatment", sample.order = "Treatment")

###heatmap for sampletype
plot_heatmap(ent11, taxa.order = "Phylum", taxa.label = "Genus", sample.label = "SampleType", sample.order = "SampleType")
plot_heatmap(ent11, taxa.order = "Kingdom", taxa.label = "Phylum", sample.label = "SampleType", sample.order = "SampleType")


##PCoA
set.seed(711L)
###GP.chl = subset_taxa(Bushman, Phylum == "Actinobacteria")
###GP.chl = subset_taxa(Bushman,Kingdom=="Bacteria")
GP.chl = subset_taxa(Bushman)
GP.chl = prune_samples(names(which(sample_sums(GP.chl) >= 20)), GP.chl)
# remove the samples that have less than 20 total reads
plot_ordination(Bushman1, ordinate(GP.chl, "MDS"), color = "SampleType") + geom_point(size = 5)


##Network
set.seed(711L)
Bushman_rare= subset_samples(Bu_rare, !is.na(Bushman))
###plot_net(Bushman_rare, maxdist = 0.4, point_label = "SampleType")
plot_net(Bushman, maxdist = 0.3, color = "SampleType", shape="Treatment")

##################subset_samples(kostic1,Timepoint%in%c("14"))
###DESeq2 Fold_change Treatment
kostic = subset_samples(ent100,Treatment != "Product")  ##top 100 abundance OTU
##kostic = subset_samples(Bushman,Treatment != "Product") ##all data
diagdds = phyloseq_to_deseq2(kostic, ~ Treatment)
dds <- estimateSizeFactors(diagdds)
norm=counts(dds, normalized=TRUE)

diagdds = estimateSizeFactors(diagdds)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)
dim(diagvst)

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
head(res,10)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#Tree,input tree first
par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
GP <- prune_taxa(taxa_sums(Bushman) > 0, Bushman)
GP.chl <- subset_taxa(GP, Genus == "Bacillus")
plot_tree(GP.chl, color = "SampleType", shape = "Family", label.tips = "Genus", 
          size = "abundance")
#plot_tree(GP, color = "SampleType", shape = "Genotype", label.tips = "Timepoint", size = "abundance")

###Get unique genus for Actinobacteria
GP <- prune_taxa(taxa_sums(Bushman) > 0, Bushman)
GP.chl <- subset_taxa(GP, Phylum == "Actinobacteria")
###GP.chl <- subset_taxa(GP, Genus == "unclassified")
get_taxa_unique(GP.chl, "Genus")


