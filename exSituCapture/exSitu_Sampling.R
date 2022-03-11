# %%% EX SITU SAMPLING %%%

library(adegenet)
library(hierfstat)
library(parallel)
library(doParallel)
library(RColorBrewer)

# %%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUERCUS ACERIFOLIA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%

# READ IN AND PROCESS GENEPOP FILE----
# Set the filepath variable below to specify which populations directory so source the genpop file from
# (Note that Stacks will write this files with the suffix ".genepop", but it needs to be ".gen")
# Using the QUAC Reference dataset (aligned using GSNAP), and the "summary" populations 
# This dataset uses an -R 80 parameter, meaning loci had to be present in 80% of all samples to be retained
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_sum/"
setwd(genpop.filePath)
# Pull in genepop object 
QUAC.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Read in the Stacks popmap values, and use these to replace @pop values 
# (using the pops accessor; this is necessary because original pop names are incorrect)
pop(QUAC.genind) <- factor(read.table("../../../QUAC_popmap2", header=FALSE)[,2])
# Genind/genpop tab slot contains a matrix of allele counts
# Total loci in assembly (after Stacks populations -R80 filtering) is 14,033. But, this includes monomorphic loci
# Genpop file only includes polymorphic loci, of which there are 6,361
nLoc(QUAC.genind)
# Each locus contains two alleles, which leads to 6,361 * 2 = 12,722 columns of the sample x allele matrix
ncol(QUAC.genind@tab)

# GENETIC CAPTURE OF WILD POPULATIONS----
# %%% All alleles %%%
# How many garden individuals
QUAC.garden.inds <- length(which(pop(QUAC.genind)=="garden"))
# Row in QUAC.genind@tab matrix where wild individuals start
QUAC.wild.index <- QUAC.garden.inds + 1
# Garden allele columns that are non-empty (present)
garden.alleles.present <- which(colSums(QUAC.genind@tab[1:QUAC.garden.inds,], na.rm = TRUE) != 0)
# Wild allele columns present
wild.alleles.present <- which(colSums(QUAC.genind@tab[QUAC.wild.index:nInd(QUAC.genind),], na.rm = TRUE) != 0)
# Number of captured wild alleles
length(which(garden.alleles.present %in% wild.alleles.present))
# Gardens capture 95.35% of wild alleles
(length(which(garden.alleles.present %in% wild.alleles.present))/length(wild.alleles.present))*100

# %%% Singletons/doubletons removed %%%
# Create allele matrix with WILD singletons/doubletons removed
QUAC.genind.SDR <- QUAC.genind@tab[,which(colSums(QUAC.genind@tab[wild.index:nrow(QUAC.genind.SDR),], na.rm=TRUE) >= 3)]
# Garden allele columns that are non-empty (present)
garden.alleles.present.SDR <- which(colSums(QUAC.genind.SDR[1:QUAC.garden.inds,], na.rm = TRUE) != 0)
# Wild allele columns present
wild.alleles.present.SDR <- which(colSums(QUAC.genind.SDR[QUAC.wild.index:nrow(QUAC.genind.SDR),], na.rm = TRUE) != 0)
# Number of captured wild alleles
length(which(garden.alleles.present.SDR %in% wild.alleles.present.SDR))
# Gardens capture 94.38% of wild alleles, with wild singletons/doubletons removed
(length(which(garden.alleles.present.SDR %in% wild.alleles.present.SDR))/length(wild.alleles.present.SDR))*100

# Old approach, using names
# Create a genpop object from genind, collapsing samples based on their populations
QUAC.genpop <- genind2genpop(QUAC.genind)
# Separate garden and wild populations, dropping alleles that are absent from each dataset
QUAC.genpop.garden <- QUAC.genpop[1,drop=TRUE]
QUAC.genpop.wild <- QUAC.genpop[2:6,drop=TRUE]
# 12,157 alleles across garden samples; 12,159 across wild populations
ncol(QUAC.genpop.garden@tab) ; ncol(QUAC.genpop.wild@tab)
# Using which to match allele names between garden and wild matrices, gardens capture 95.3% of wild alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% colnames(QUAC.genpop.wild@tab)))/length(colnames(QUAC.genpop.wild@tab)))*100

# ALLELE CATEGORIES----
# %%% All alleles %%%
# Categorize wild alleles, then determine how many of each allele category gardens are capturing
# These fractions demonstrate the frequency of finding an allele when drawing out a single individual
QUAC.wildAlleleFreqs <- (colSums(QUAC.genpop.wild@tab)/(nInd(QUAC.genind)*2))*100
# Max frequency is 49.5%; min is 0.25%
cat(paste("Max wild allele frequency: ",max(QUAC.wildAlleleFreqs),"\n","Min wild allele frequency: ",min(QUAC.wildAlleleFreqs)))
# Gardens capture:
# 100% of very common alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs > 10)])))/length(which(QUAC.wildAlleleFreqs > 10)))*100
# 100% of common alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs > 5)])))/length(which(QUAC.wildAlleleFreqs > 5)))*100
# 77.4% of low frequency alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs < 10 & QUAC.wildAlleleFreqs > 1)])))/length(which(QUAC.wildAlleleFreqs < 10 & QUAC.wildAlleleFreqs > 1)))*100
# 100% of rare alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs < 1)])))/length(which(QUAC.wildAlleleFreqs < 1)))*100

# %%% Removing singletons/doubletons %%%
# Find alleles which are only present once or twice total in the entire wild sample set
# Singletons
QUAC.genpop.wild.SDRmat <- QUAC.genpop.wild@tab[,-which(colSums(QUAC.genpop.wild@tab, na.rm=TRUE) == 1)]
# Doubletons
QUAC.genpop.wild.SDRmat <- QUAC.genpop.wild.SDRmat[,-which(colSums(QUAC.genpop.wild.SDRmat, na.rm=TRUE) == 2)]
# Categorize wild alleles, then determine how many of each allele category gardens are capturing
QUAC.wildAlleleFreqs.SDR <- (colSums(QUAC.genpop.wild.SDRmat)/(nrow(QUAC.genind.wild@tab)*2))*100
# Max frequency is 49.5%; min is 0.25%
cat(paste("Max wild allele frequency: ",max(QUAC.wildAlleleFreqs.SDR),"\n","Min wild allele frequency: ",min(QUAC.wildAlleleFreqs.SDR)))
# Gardens capture:
# 100% of very common alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs.SDR[which(QUAC.wildAlleleFreqs.SDR > 10)])))/length(which(QUAC.wildAlleleFreqs.SDR > 10)))*100
# 100% of common alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs.SDR[which(QUAC.wildAlleleFreqs.SDR > 5)])))/length(which(QUAC.wildAlleleFreqs.SDR > 5)))*100
# 77.4% of low frequency alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs.SDR[which(QUAC.wildAlleleFreqs.SDR < 10 & QUAC.wildAlleleFreqs.SDR > 1)])))/length(which(QUAC.wildAlleleFreqs.SDR < 10 & QUAC.wildAlleleFreqs.SDR > 1)))*100
# 100% of rare alleles
(length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs.SDR[which(QUAC.wildAlleleFreqs.SDR < 1)])))/length(which(QUAC.wildAlleleFreqs.SDR < 1)))*100

# DRAFT: ALLELE CATEGORIES NEW----
# 1. Get allele categories to work without matching strings
# 2. Calculate number of individuals needed to capture 95% of allelic diversity (w/ and w/o singletons/doubletons)
# 3. Optimize sampling curve calculations; then, test with increased reps
# 4. Calculate diversity metrics: allelic richness, heterozygosity, and possibly estimated effective population size

# Subsetting wild samples using seppop/repool
QUAC.genind.wild <- seppop(QUAC.genind)[2:6]
# (Using repool will subset the resulting genind object to only include loci present in the combined populations)
QUAC.genind.wild <- repool(QUAC.genind.wild$porterMt,QUAC.genind.wild$magazineMt,QUAC.genind.wild$pryorMt,
                           QUAC.genind.wild$sugarloaf_midlandPeak,QUAC.genind.wild$kessler_shaleBarrenRidge)

# Subsetting wild samples using repool_new function from Sean
repool_new <- function(genind_obj,vect_pops){
  genind_obj_sep <- seppop(genind_obj)
  genind_obj_merge <- genind_obj_sep[[vect_pops[1]]]
  for (i in 1:(length(vect_pops)-1)) genind_obj_merge<-repool(genind_obj_merge,genind_obj_sep[[vect_pops[i+1]]])
  genind_obj_merge
}
QUAC.genind.wild.NEW <- repool_new(QUAC.genind,2:6)

# NOTE: there are differences in allele matrics depending on how genind/genpop objects are constructed
ncol(QUAC.genind@tab)
ncol(QUAC.genpop@tab)
# Without na.rm argument, colSums don't match across genind and genpop objects
identical(colSums(QUAC.genind@tab),colSums(QUAC.genpop@tab))
identical(colSums(QUAC.genind@tab, na.rm = TRUE),colSums(QUAC.genpop@tab, na.rm = TRUE))
# Subset objects have NAs removed
ncol(QUAC.genind.wild@tab)
ncol(QUAC.genind.wild.NEW@tab)
ncol(QUAC.genpop.wild@tab)


# Calculate wild allele frequencies
QUAC.wildAlleleFreqs <- (colSums(QUAC.genind@tab[wild.index:nInd(QUAC.genind),], na.rm = TRUE)/(nInd(QUAC.genind)*2))*100
# Calculate allelic capture for different categories, based on garden.alleles.present vector from before
# Wild allele frequences that are very commmon (greater than 10%)
which(QUAC.wildAlleleFreqs > 10)
# Vector of TRUE/FALSE for garden alleles being included within the very common wild alleles
garden.alleles.present %in% which(QUAC.wildAlleleFreqs > 10)
# Garden alleles that ARE included within very common wild alleles (i.e. get TRUE values)
which(garden.alleles.present %in% which(QUAC.wildAlleleFreqs > 10))
# NUMBER of garden alleles included within very common wild alleles
length(which(garden.alleles.present %in% which(QUAC.wildAlleleFreqs > 10)))
# Number of garden alleles included within very common wild alleles, divided by number of very common alleles
length(which(garden.alleles.present %in% which(QUAC.wildAlleleFreqs > 10)))/length(which(QUAC.wildAlleleFreqs > 10))*100
# Common alleles
length(which(garden.alleles.present %in% which(QUAC.wildAlleleFreqs > 5)))/length(which(QUAC.wildAlleleFreqs > 5))*100
# Low frequency alleles
length(which(garden.alleles.present %in% which(QUAC.wildAlleleFreqs < 10 & QUAC.wildAlleleFreqs > 1)))/length(which(QUAC.wildAlleleFreqs < 10 & QUAC.wildAlleleFreqs > 1))*100
# Rare alleles
length(which(garden.alleles.present %in% which(QUAC.wildAlleleFreqs < 1)))/length(which(QUAC.wildAlleleFreqs < 1))*100

# RESAMPLING----

# Create a matrix of wild individuals versus alleles
QUAC.genind.wild <- seppop(QUAC.genind)[2:6]
QUAC.genind.wild <- repool(QUAC.genind.wild$porterMt,QUAC.genind.wild$magazineMt,QUAC.genind.wild$pryorMt,
                           QUAC.genind.wild$sugarloaf_midlandPeak,QUAC.genind.wild$kessler_shaleBarrenRidge)
# Processing duplicates
# Only one duplicate, in QUAC Wild data set
which(rownames(QUAC.genind.wild@tab) == "QUAC_W_DUP_SH_Q2121")
QUAC.genind.wild@tab[65:66,1:4] # QUAC_W_SH_Q2121 and QUAC_W_DUP_SH_Q2121
length(which(QUAC.genind.wild@tab[65,] != QUAC.genind.wild@tab[66,]))
# 874/12,159 alleles (7.19%) are different between these duplicate samples
# For now, just removing duplicate sample
QUAC.genind.wild@tab <- QUAC.genind.wild@tab[-66,]

# Function for measuring wild allelic capture of samples
# The two arguments are the vector of frequencies of wild alleles, and the sample of the wild genind object
get.allele.cat <- function(freq.vector, sample.mat){
  # Calculate percentages
  # Very common alleles (greater than 10%)
  v_com <- (length(which(names(which(colSums(sample.mat, na.rm = TRUE)!=0)) %in% names(freq.vector[which(freq.vector > 10)])))/length(which(freq.vector > 10)))*100
  # Common alleles (greater than 5%)
  com <- (length(which(names(which(colSums(sample.mat, na.rm = TRUE)!=0)) %in% names(freq.vector[which(freq.vector > 5)])))/length(which(freq.vector > 5)))*100
  # Low frequency alleles (between 1% and 10%)
  low_freq <- (length(which(names(which(colSums(sample.mat, na.rm = TRUE)!=0)) %in% names(freq.vector[which(freq.vector < 10 & freq.vector > 1)])))/length(which(freq.vector < 10 & freq.vector > 1)))*100
  # Rare alleles (less than 1%)
  rare <- (length(which(names(which(colSums(sample.mat, na.rm = TRUE)!=0)) %in% names(freq.vector[which(freq.vector < 1)])))/length(which(freq.vector < 1)))*100
  # Concatentate values to a vector, and return
  return(c(v_com,com,low_freq,rare))
}

# Build a 3D array, with rows being sample numbers and columns being allele frequency categories
# 3rd dimension is replicates. nrow is number of individuals-1 because sample doesn't work for vectors of length 1
num_reps <- 10 # In Hoban et al. 2020 (where QUBO is considered), num_reps = 75,000
list_allele_cat <- c("v_com","com","low_freq","rare")
samplingResults <- array(dim=c(nrow(QUAC.genind.wild@tab)-1,length(list_allele_cat),num_reps))
colnames(samplingResults) <- list_allele_cat

# For each replicate (which is the third dimension, in the samplingResults array)...
for(i in 1:num_reps){
  # loop through sampling from 2 to the maximum number of wild individuals
  for(j in 2:nrow(QUAC.genind.wild@tab)){
    # Create a sample of the wild allele matrix, of "j" size
    samp <- QUAC.genind.wild@tab[sample(nrow(QUAC.genind.wild@tab), size=j, replace = FALSE),]
    # Calculate how many alleles of each category that sample captures,
    # and place those percentages into the row of the samplingResults array
    samplingResults[j-1,,i] <- get.allele.cat(QUAC.wildAlleleFreqs,samp)
  }
}

# Plotting
v.com_means <- apply(samplingResults[,1,], 1, mean)
v.com_sd <- apply(samplingResults[,1,], 1, sd)

com_means <- apply(samplingResults[,2,], 1, mean)
com_sd <- apply(samplingResults[,2,], 1, sd)

lowfr_means <- apply(samplingResults[,3,], 1, mean)
lowfr_sd <- apply(samplingResults[,3,], 1, sd)

rare_means <- apply(samplingResults[,4,], 1, mean)
rare_sd <- apply(samplingResults[,4,], 1, sd)

plotColors <- brewer.pal(n=4, name="Dark2")
plot(v.com_means, ylim=c(0,110), col=plotColors[1], pch=16, main="QUAC: Resampling Curve (10 replicates)",
     xlab="Number of samples", ylab="Mean percent captured")
points(com_means, col=plotColors[2], pch=16)
points(lowfr_means, col=plotColors[3], pch=16)
points(rare_means, col=plotColors[4], pch=16)
legend("bottomright", inset = 0.05, legend = c("Very common: >10%","Common: >5%","Low frequency: 1--10%", "Rare: <1%"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2)

# Zoomed in 
plot(v.com_means, ylim=c(0,110), xlim=c(0,5), col=plotColors[1], pch=16, main="QUAC: Resampling Curve (10 replicates)",
     xlab="Number of samples", ylab="Mean percent captured")
points(com_means, col=plotColors[2], pch=16)
points(lowfr_means, col=plotColors[3], pch=16)
points(rare_means, col=plotColors[4], pch=16)
legend("bottomright", inset = 0.05, legend = c("Very common: >10%","Common: >5%","Low frequency: 1--10%", "Rare: <1%"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2)

# %%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUERCUS BOYNTONII %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%

# READ IN AND PROCESS GENEPOP FILE----
# Using the QUBO Reference dataset (aligned using GSNAP) 
# This generated the largest number of polymorphic loci (compared to de novo and Hybrid analyses)
# This dataset uses an -R 80 parameter, meaning loci had to be present in 80% of all samples to be retained
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
