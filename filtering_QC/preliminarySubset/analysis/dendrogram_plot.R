# Plot unrooted Neighbor-Joining Tree

library(ape)
library(hierfstat)
library(adegenet)

# "Whitelist" loci (1,000 random loci)
# Pull in genepop object (with file suffix updated--Stacks writes as ".genepop", but needs to be ".gen")
prelimWL_genind <- read.genepop("populations.snps.gen")
# Number of loci
length(prelimWL_genind@loc.n.all)
# Create distance matrix, based off of SNP matrix (individuals in rows, alleles in columns)
D_WL <- dist(tab(prelimWL_genind))
# Create Neighbor-Joining Tree, based off of SNP matrix
tree_WL <- nj(D_WL)
plot(tree_WL, type="unrooted", edge.w=2, main="1,000 loci")
plot(tree_WL, type="phylogram", edge.w=2, main="1,000 loci")

# All loci
# Pull in genepop object (with file suffix updated--Stacks writes as ".genepop", but needs to be ".gen")
prelimAll_genind <- read.genepop("populations.snps.gen")
# Number of loci
length(prelimAll_genind@loc.n.all)
# Create distance matrix, based off of SNP matrix (individuals in rows, alleles in columns)
D <- dist(tab(prelimAll_genind))
# Create Neighbor-Joining Tree, based off of SNP matrix
tree <- nj(D)
plot(tree, type="unrooted", edge.w=2, main="47,266 loci")
plot(tree, type="phylogram", edge.w=2, main="47,266 loci")
