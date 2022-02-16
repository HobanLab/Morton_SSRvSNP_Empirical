# Garden and Wild loci comparison

library(adegenet)
library(hierfstat)
library(parallel)
library(doParallel)

# Using the QUAC Reference dataset, as the generated the largest number of polymorphic loci (compared to de novo and Hybrid)
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_sum/")

# Pull in genepop object (with file suffix updated--Stacks writes as ".genepop", but needs to be ".gen")
QUAC.genind <- read.genepop("populations.snps.gen")
# This genind object has incorrect population assignments (populations are somehow sample names instead)
# Read in the Stacks popmap values, and use these to replace @pop values (using the pops accessor)
pop(QUAC.genind) <- factor(read.table("../../../QUAC_popmap", header=FALSE)[,2])

# Genind tab slot contains a matrix of allele counts
# Total number of loci in assembly (after R80 filtering) is 14,033. But, this includes monomorphic loci
# Genpop file only includes polymorphic loci, of which there are 6,361
nLoc(QUAC.genind)
# Each locus contains two alleles, which leads to 6,361 * 2 = 12,722 columns of the sample x allele matrix
ncol(QUAC.genind@tab)

# Create a genpop object from genind
QUAC.genpop <- genind2genpop(QUAC.genind)
QUAC.genpop@tab[1:6,1:6]
# Separate garden and wild populations
QUAC.genpop.garden <- QUAC.genpop[1,]
QUAC.genpop.wild <- QUAC.genpop[2:6,]

# Before we determine the genetic capture of wild populations,
# we have to see whether there are alleles that are present only in the garden samples
length(which(colSums(QUAC.genpop.wild@tab) == 0))
# There are--563 alleles

# Remove the missing alleles from the wild matrix
names(which(colSums(QUAC.genpop.wild@tab) != 0))
test <- QUAC.genpop.wild[, loc=(names(which(colSums(QUAC.genpop.wild@tab) != 0)))]
test <- QUAC.genpop.wild[, loc="2385_35.01"]

test <- QUAC.genpop.wild@tab[,names(which(colSums(QUAC.genpop.wild@tab) != 0))]

locNames(QUAC.genpop.wild)

unique(alleles(QUAC.genpop.wild))
QUAC.genpop.wild@tab


# 
nAll(QUAC.genpop[2:6,, drop=TRUE])
length(nAll(QUAC.genpop[1:3,, drop=TRUE]))

length(nAll(QUAC.genpop))

which(QUAC.genpop.wild@tab == 0)

length(which(QUAC.genpop.wild@tab == 0))
length(which(QUAC.genpop[2,]@tab == 0))
length(which(QUAC.genpop[3,]@tab == 0))
length(which(QUAC.genpop[4,]@tab == 0))
length(which(QUAC.genpop[5,]@tab == 0))
length(which(QUAC.genpop[6,]@tab == 0))







# Alleles that are missing from garden populations
which(QUAC.genpop.garden@tab == 0)
# Alleles that are missing from wild populations
which(QUAC.genpop.wild@tab == 0)

length(match(which(QUAC.genpop.wild@tab == 0), which(QUAC.genpop.garden@tab == 0)))

length(which(QUAC.genpop.wild@tab == 0))

# Gardens are missing 565 alleles
length(which(QUAC.genpop.garden@tab == 0))
# So, gardens capture ~95.6% of alleles
(length(which(QUAC.genpop.garden@tab != 0))/(nLoc(QUAC.genpop.garden)*2))*100
