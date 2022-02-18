# %%% EX SITU SAMPLING %%%

library(adegenet)
library(hierfstat)
library(parallel)
library(doParallel)

# %%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUERCUS ACERIFOLIA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%

# Using the QUAC Reference dataset (aligned using GSNAP) 
# This generated the largest number of polymorphic loci (compared to de novo and Hybrid analyses)
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_sum/")

# Pull in genepop object (with the file suffix updated--Stacks writes this as ".genepop", but it needs to be ".gen")
QUAC.genind <- read.genepop("populations.snps.gen")
# Read in the Stacks popmap values, and use these to replace @pop values (using the pops accessor; original pop names are incorrect)
pop(QUAC.genind) <- factor(read.table("../../../QUAC_popmap2", header=FALSE)[,2])

# Genind/genpop tab slot contains a matrix of allele counts
# Total number of loci in assembly (after Stacks populations -R80 filtering) is 14,033. But, this includes monomorphic loci
# Genpop file only includes polymorphic loci, of which there are 6,361
nLoc(QUAC.genind)
# Each locus contains two alleles, which leads to 6,361 * 2 = 12,722 columns of the sample x allele matrix
ncol(QUAC.genind@tab)

# GENETIC CAPTURE OF WILD POPULATIONS----
# Create a genpop object from genind, collapsing samples based on their populations
QUAC.genpop <- genind2genpop(QUAC.genind)
# Separate garden and wild populations, dropping alleles that are absent from each dataset
QUAC.genpop.garden <- QUAC.genpop[1,drop=TRUE]
QUAC.genpop.wild <- QUAC.genpop[2:6,drop=TRUE]
# 12,157 alleles across garden samples; 12,159 across wild populations
ncol(QUAC.genpop.garden@tab) ; ncol(QUAC.genpop.wild@tab)
# There are 563 alleles that are unique to the garden population (i.e. present in gardens, but not in the wild) 
ncol(QUAC.genpop@tab) - ncol(QUAC.genpop.wild@tab)
# Using which to match allele names between garden and wild matrices, gardens capture 95.3% of wild alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% colnames(QUAC.genpop.wild@tab)))/length(colnames(QUAC.genpop.wild@tab))) * 100

# I'm not sure why colnames has to be used, or what the line below is comparing between the allele matrices
# But, it doesn't seem to be correct
# length(which(QUAC.genpop.garden@tab %in% QUAC.genpop.wild@tab))/ncol(QUAC.genpop.wild@tab) * 100

# ALLELE CATEGORIES----
# Categorize wild alleles (based on Sean's 4 categories), then determine how many of each allele category gardens are capturing

QUAC.genpop.wild@tab[,1:6]
colSums(QUAC.genpop.wild@tab[,1:6])

QUAC.wildAlleleSums <- colSums(QUAC.genpop.wild@tab)
sum(colSums(QUAC.genpop.wild@tab))

# For each allele (matrix  columnn), calculate its frequency

# Sean and Kaylee's code use a denominator of population size * 2 (for diploids; i.e. the number of chromosomes in the population)
# These fractions demonstrate the frequency of finding an allele when drawing out a single individual
QUAC.wildAlleleFreqs <- (colSums(QUAC.genpop.wild@tab)/(nInd(QUAC.genind)*2))*100
sum(QUAC.wildAlleleFreqs)
# Max frequency is 49.5%; min is 0.25%
max(QUAC.wildAlleleFreqs)
min(QUAC.wildAlleleFreqs)

# This calculation, on the other hand, gives you the frequency of drwaing an allele when pulling from the entire allele pool
# QUAC.wildAlleleFreqs <- (colSums(QUAC.genpop.wild@tab)/sum(colSums(QUAC.genpop.wild@tab))*100)
