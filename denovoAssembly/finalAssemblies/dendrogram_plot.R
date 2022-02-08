# Plot unrooted Neighbor-Joining Tree

# Load libraries
library(ape)
library(hierfstat)
library(adegenet)
# Working directory set to /RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/

# QUAC de novo final assembly (DNFA)----
# Pull in genepop object (with file suffix updated--Stacks writes as ".genepop", but needs to be ".gen")
QUAC_DNFA_genind <- read.genepop("QUAC/output/populations_Sum/populations.snps.gen")
# Number of loci
length(QUAC_DNFA_genind@loc.n.all)
# Create distance matrix, based off of SNP matrix (individuals in rows, alleles in columns)
D <- dist(tab(QUAC_DNFA_genind))
# Create Neighbor-Joining Tree, based off of SNP matrix
QUAC.DNFA.tree <- nj(D)

popmap <- read.table("QUAC/QUAC_popmap", skip = 1)
popmap[,2]
QUAC.DNFA.tree$tip.label <- popmap[,2]

plot(QUAC.DNFA.tree, type="unrooted", edge.w=2, main="QUAC de novo: 5,208 loci", show.tip.label=FALSE)
plot(QUAC.DNFA.tree, type="phylogram", edge.w=2, main="QUAC de novo: 5,208 loci")
plot(QUAC.DNFA.tree, type="phylogram", edge.w=2, main="QUAC de novo: 5,208 loci", show.tip.label=FALSE)

# QUBO de novo final assembly (DNFA)----
# Pull in genepop object (with file suffix updated--Stacks writes as ".genepop", but needs to be ".gen")
QUBO_DNFA_genind <- read.genepop("QUBO/output/populations_Sum/populations.snps.gen")
# Number of loci
length(QUBO_DNFA_genind@loc.n.all)
# Create distance matrix, based off of SNP matrix (individuals in rows, alleles in columns)
D <- dist(tab(QUBO_DNFA_genind))
# Create Neighbor-Joining Tree, based off of SNP matrix
QUBO.DNFA.tree <- nj(D)

popmap <- read.table("QUAC/QUBO_popmap", skip = 1)
popmap[,2]
QUBO.DNFA.tree$tip.label <- popmap[,2]

plot(QUBO.DNFA.tree, type="unrooted", edge.w=2, main="QUBO de novo: 4,758 loci", show.tip.label=FALSE)
plot(QUBO.DNFA.tree, type="phylogram", edge.w=2, main="QUBO de novo: 4,758 loci")
plot(QUBO.DNFA.tree, type="phylogram", edge.w=2, main="QUBO de novo: 4,758 loci", show.tip.label=FALSE)