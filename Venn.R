#Set working directory
setwd("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/Venn/")

#load R data and source itag_diversity scripts
load("/Users/siwendeng/Downloads/VESTA_Strawberry/Phyloseq/VESTA_Final_base.RData")
#install.packages("reshape")
library("reshape")
library("ggplot2")
library("VennDiagram")
library("gplots")

# rar_control <- prune_taxa(taxa_sums(rar_control) > 0, rar_control)
# rar_vesta <- prune_taxa(taxa_sums(rar_vesta) > 0, rar_vesta)
# rar_control_OTU <- row.names(otu_table(rar_control))
# rar_vesta_OTU <- row.names(otu_table(rar_vesta))
# venn(list(rar_control_OTU,rar_vesta_OTU))

#control
rar_control <- prune_taxa(taxa_sums(rar_control) > 0, rar_control)
rar_control_otu <- row.names(otu_table(rar_control))
#treated
rar_vesta <- prune_taxa(taxa_sums(rar_vesta) > 0, rar_vesta)
rar_vesta_otu <- row.names(otu_table(rar_vesta))

venn_1 <- venn(list(control = rar_control_otu, treated = rar_vesta_otu, VESTA = rar_vesta_sub))


#VESTA
rar_vesta_sub <- subset_samples(rar_pro, Timepoint == "2")
rar_vesta_sub <- prune_taxa(taxa_sums(rar_vesta_sub) > 0, rar_vesta_sub)
rar_vesta_sub <- row.names(otu_table(rar_vesta_sub))

#Root
rar_root <- subset_samples(rar_root, Treatment == "VESTA")
rar_root <- prune_taxa(taxa_sums(rar_root) > 0, rar_root)
rar_root_otu <- row.names(otu_table(rar_root))

#Rhizo
rar_rhizo <- subset_samples(rar_rhizo, Treatment == "VESTA")
rar_rhizo <- prune_taxa(taxa_sums(rar_rhizo) > 0, rar_rhizo)
rar_rhizo_otu <- row.names(otu_table(rar_rhizo))

#soil
rar_soil <- subset_samples(rar_soil, Treatment == "VESTA")
rar_soil <- prune_taxa(taxa_sums(rar_soil) > 0, rar_soil)
rar_soil_otu <- row.names(otu_table(rar_soil))

#Venn all
venn_all <- venn(list(Root = rar_root_otu, Rhizo = rar_rhizo_otu, Soil = rar_soil_otu, VESTA = rar_vesta_sub))
ggsave(filename = "venn_all.pdf", venn_all, width = 8, height = 6)


venn_all <- venn_all[1:16, 1:5]
venn_all
dev.off()
venn.plot <- draw.quad.venn(
  area1 = 200,
  area2 = 10,
  area3 = 35,
  area4 = 4,
  n12 = 19,
  n13 = 754,
  n14 = 115,
  n23 = 0,
  n24 = 2,
  n34 = 3,
  n123 = 3,
  n124 = 45,
  n134 = 8,
  n234 = 1688,
  n1234 = 430,
  category = c("First", "Second", "Third", "Fourth"),
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
);

#root, control vs treated
rar_root_nt <- subset_samples(rar_root, Treatment == "Control")
rar_root_nt <- prune_taxa(taxa_sums(rar_root_nt) > 0, rar_root_nt)
rar_root_vt <- subset_samples(rar_root, Treatment == "VESTA")
rar_root_vt <- prune_taxa(taxa_sums(rar_root_vt) > 0, rar_root_vt)
rar_root_nt <- row.names(otu_table(rar_root_nt))
rar_root_vt <- row.names(otu_table(rar_root_vt))
venn_root <- venn(list(root_control = rar_root_nt, root_treated = rar_root_vt))
venn_root <- venn_root[1:4, 1:3]
venn_root
#     num root_control root_treated
# 00    0            0            0
# 01  396            0            1
# 10  594            1            0
# 11 1189            1            1

dev.off()
draw.pairwise.venn(
  area1 = 594+1189,
  area2 = 396+1189,
  cross.area = 1189,
  category = c("Root_Control", "Root_Treated"),
  fill = c("#FFBB00","#20948B"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = c(-30, 150)
)


#rhizosphere, control vs treated
rar_rhizo_nt <- subset_samples(rar_rhizo, Treatment == "Control")
rar_rhizo_nt <- prune_taxa(taxa_sums(rar_rhizo_nt) > 0, rar_rhizo_nt)
rar_rhizo_vt <- subset_samples(rar_rhizo, Treatment == "VESTA")
rar_rhizo_vt <- prune_taxa(taxa_sums(rar_rhizo_vt) > 0, rar_rhizo_vt)
rar_rhizo_nt <- row.names(otu_table(rar_rhizo_nt))
rar_rhizo_vt <- row.names(otu_table(rar_rhizo_vt))
venn_rhizo <- venn(list(rhizo_control = rar_rhizo_nt,rhizo_treated = rar_rhizo_vt))
venn_rhizo <- venn_rhizo[1:4, 1:3]
venn_rhizo

#     num rhizo_control rhizo_treated
# 00    0             0             0
# 01  161             0             1
# 10  524             1             0
# 11 2378             1             1

dev.off()
draw.pairwise.venn(
  area1 = 524+2378,
  area2 = 161+2378,
  cross.area = 2378,
  category = c("Rhizo_Control", "Rhizo_Treated"),
  fill = c("#FFBB00","#20948B"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = c(-30, 150)
)

#soil, control vs treated
rar_soil_nt <- subset_samples(rar_soil, Treatment == "Control")
rar_soil_nt <- prune_taxa(taxa_sums(rar_soil_nt) > 0, rar_soil_nt)
rar_soil_vt <- subset_samples(rar_soil, Treatment == "VESTA")
rar_soil_vt <- prune_taxa(taxa_sums(rar_soil_vt) > 0, rar_soil_vt)
rar_soil_nt <- row.names(otu_table(rar_soil_nt))
rar_soil_vt <- row.names(otu_table(rar_soil_vt))
venn_soil <- venn(list(soil_control = rar_soil_nt, soil_treated = rar_soil_vt))
venn_soil <- venn_soil[1:4, 1:3]
venn_soil
#     num soil_control soil_treated
# 00    0            0            0
# 01  219            0            1
# 10  420            1            0
# 11 2399            1            1

dev.off()
draw.pairwise.venn(
  area1 = 420+2399,
  area2 = 219+2399,
  cross.area = 2399,
  category = c("Soil_Control", "Soil_Treated"),
  fill = c("#FFBB00","#20948B"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = c(-30, 150)
)

