# load packages
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("scales")
library("grid")
library("DESeq2")
library("ape")
library("reshape2")
library("vegan")

# set working directory
setwd("~/Downloads/VESTA_Strawberry/Phyloseq/")

# import biom and map file
biom_file = "final_otu_table.biom"
map_file = "VESTA_Final_metadata.txt"
tree= read.tree("VESTA.phylo.tre")
otutax = import_biom(biom_file, parseFunction=parse_taxonomy_greengenes)
sam = import_qiime_sample_data(map_file)
otutax
# randomly root tree
is.rooted(tree)
tree = root(tree, 1, resolve.root = T)
is.rooted(tree)

# use merge_phyloseq to build up the phyloseq object
all = merge_phyloseq(otutax, sam, tree)
all

# # total reads per sample
# readsumsdf = data.frame(nreads = sort(taxa_sums(all), TRUE), sorted = 1:ntaxa(all), type = "OTUs")
# readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(all), TRUE), sorted = 1:nsamples(all), type = "Samples"))
# title = "Total number of reads"
# p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
# p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#filtering
filter = prune_samples(sample_sums(all)>=1000, all)
filter = prune_taxa(taxa_sums(filter)>=5, filter)
filter = filter_taxa((filter), function(x) sum(x > 5) > 5, TRUE)
filter

##define a random number seed so that others can reproduce your subsetted data
##set.seed(28132)
#rarefy sample

rar = rarefy_even_depth(filter)
rar # strawberry and product

# subset sample
## product only
rar_pro = subset_samples(rar, Treatment == "Product")
## strawberry only
rar_sam = subset_samples(rar, Treatment != "Product")
## by treatment
rar_control = subset_samples(rar, Treatment == "Control")
rar_vesta = subset_samples(rar, Treatment == "VESTA")
## by sample type
rar_root = subset_samples(rar_sam, SampleType == "root")
rar_rhizo = subset_samples(rar_sam, SampleType == "rhizosphere")
rar_soil = subset_samples(rar_sam, SampleType == "soil")
## by timepoint
rar_1 = subset_samples(rar_sam, Timepoint == "1")
rar_2 = subset_samples(rar_sam, Timepoint == "2")
rar_3 = subset_samples(rar_sam, Timepoint == "3")
rar_4 = subset_samples(rar_sam, Timepoint == "4")

save.image("VESTA_Final_base.RData")
