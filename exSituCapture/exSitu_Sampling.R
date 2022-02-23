# %%% EX SITU SAMPLING %%%

library(adegenet)
library(hierfstat)
library(parallel)
library(doParallel)

# %%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUERCUS ACERIFOLIA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%

# READ IN AND PROCESS GENEPOP FILE----
# Using the QUAC Reference dataset (aligned using GSNAP) 
# This generated the largest number of polymorphic loci (compared to de novo and Hybrid analyses)
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_sum/")

# Pull in genepop object (with the file suffix updated--Stacks writes this as ".genepop", but it needs to be ".gen")
QUAC.genind <- read.genepop("populations.snps.gen")
# Read in the Stacks popmap values, and use these to replace @pop values 
# (using the pops accessor; original pop names are incorrect)
pop(QUAC.genind) <- factor(read.table("../../../QUAC_popmap2", header=FALSE)[,2])

# Genind/genpop tab slot contains a matrix of allele counts
# Total loci in assembly (after Stacks populations -R80 filtering) is 14,033. But, this includes monomorphic loci
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
(length(which(colnames(QUAC.genpop.garden@tab) %in% colnames(QUAC.genpop.wild@tab)))/length(colnames(QUAC.genpop.wild@tab)))*100

# ALLELE CATEGORIES----
# Categorize wild alleles, then determine how many of each allele category gardens are capturing
# These fractions demonstrate the frequency of finding an allele when drawing out a single individual
QUAC.wildAlleleFreqs <- (colSums(QUAC.genpop.wild@tab)/(nInd(QUAC.genind)*2))*100
# Max frequency is 49.5%; min is 0.25%
cat(paste("Max wild allele frequency: ",max(QUAC.wildAlleleFreqs),"\n","Min wild allele frequency: ",min(QUAC.wildAlleleFreqs)))
# Gardens capture:
# 100% of very common alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs > 10)])))/length(which(QUAC.wildAlleleFreqs > 10)))*100
# 100% of very common alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs > 5)])))/length(which(QUAC.wildAlleleFreqs > 5)))*100
# 77.4% of low frequency alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs < 10 & QUAC.wildAlleleFreqs > 1)])))/length(which(QUAC.wildAlleleFreqs < 10 & QUAC.wildAlleleFreqs > 1)))*100
# 100% of rare alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs < 1)])))/length(which(QUAC.wildAlleleFreqs < 1)))*100

# %%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUERCUS BOYNTONII %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%

# READ IN AND PROCESS GENEPOP FILE----
# Using the QUBO Reference dataset (aligned using GSNAP) 
# This generated the largest number of polymorphic loci (compared to de novo and Hybrid analyses)
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP/output/populations_sum/")

# Pull in genepop object (with the file suffix updated--Stacks writes this as ".genepop", but it needs to be ".gen")
QUBO.genind <- read.genepop("populations.snps.gen")
# Read in the Stacks popmap values, and use these to replace @pop values 
# (using the pops accessor; original pop names are incorrect)
pop(QUBO.genind) <- factor(read.table("../../../QUBO_popmap2", header=FALSE)[,2])

# Genind/genpop tab slot contains a matrix of allele counts
# Total loci in assembly (after Stacks populations -R80 filtering) is 29,964. But, this includes monomorphic loci
# Genpop file only includes polymorphic loci, of which there are 13,821
nLoc(QUBO.genind)
# Each locus contains two alleles, which leads to 13,821 * 2 = 27,642 columns of the sample x allele matrix
ncol(QUBO.genind@tab)

# GENETIC CAPTURE OF WILD POPULATIONS----
# Create a genpop object from genind, collapsing samples based on their populations
QUBO.genpop <- genind2genpop(QUBO.genind)
# Separate garden and wild populations, dropping alleles that are absent from each dataset
QUBO.genpop.garden <- QUBO.genpop[1,drop=TRUE]
QUBO.genpop.wild <- QUBO.genpop[2:12,drop=TRUE]
# 25,817 alleles across garden samples; 26,711 across wild populations
ncol(QUBO.genpop.garden@tab) ; ncol(QUBO.genpop.wild@tab)
# There are 931 alleles that are unique to the garden population (i.e. present in gardens, but not in the wild) 
ncol(QUBO.genpop@tab) - ncol(QUBO.genpop.wild@tab)
# Using which to match allele names between garden and wild matrices, gardens capture 93.2% of wild alleles
(length(which(colnames(QUBO.genpop.garden@tab) %in% colnames(QUBO.genpop.wild@tab)))/length(colnames(QUBO.genpop.wild@tab)))*100

# ALLELE CATEGORIES----
# Categorize wild alleles, then determine how many of each allele category gardens are capturing
# These fractions demonstrate the frequency of finding an allele when drawing out a single individual
QUBO.wildAlleleFreqs <- (colSums(QUBO.genpop.wild@tab)/(nInd(QUBO.genind)*2))*100
# Max frequency is 53.0%; min is 0.28%
cat(paste("Max wild allele frequency: ",max(QUBO.wildAlleleFreqs),"\n","Min wild allele frequency: ",min(QUBO.wildAlleleFreqs)))
# Gardens capture:
# 99.99281% of very common alleles
(length(which(colnames(QUBO.genpop.garden@tab) %in% names(QUBO.wildAlleleFreqs[which(QUBO.wildAlleleFreqs > 10)])))/length(which(QUBO.wildAlleleFreqs > 10)))*100
# 99.98578% of very common alleles
(length(which(colnames(QUBO.genpop.garden@tab) %in% names(QUBO.wildAlleleFreqs[which(QUBO.wildAlleleFreqs > 5)])))/length(which(QUBO.wildAlleleFreqs > 5)))*100
# 71.9% of very common alleles
(length(which(colnames(QUBO.genpop.garden@tab) %in% names(QUBO.wildAlleleFreqs[which(QUBO.wildAlleleFreqs < 10 & QUBO.wildAlleleFreqs > 1)])))/length(which(QUBO.wildAlleleFreqs < 10 & QUBO.wildAlleleFreqs > 1)))*100
# 94.9% of rare alleles
(length(which(colnames(QUBO.genpop.garden@tab) %in% names(QUBO.wildAlleleFreqs[which(QUBO.wildAlleleFreqs < 1)])))/length(which(QUBO.wildAlleleFreqs < 1)))*100
