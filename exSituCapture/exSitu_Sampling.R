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
# Make variables for number of garden and wild individuals
QUAC.garden.inds <- length(which(pop(QUAC.genind)=="garden")); QUAC.wild.inds <- nInd(QUAC.genind)-QUAC.garden.inds
# Row in QUAC.genind@tab matrix where wild individuals start
QUAC.wild.index <- QUAC.garden.inds + 1

# GENETIC CAPTURE OF WILD POPULATIONS----
# %%% All alleles %%%
# Garden allele columns that are non-empty (present)
QUAC.garden.alleles.present <- which(colSums(QUAC.genind@tab[1:QUAC.garden.inds,], na.rm = TRUE) != 0)
# Wild allele columns present
QUAC.wild.alleles.present <- which(colSums(QUAC.genind@tab[QUAC.wild.index:nInd(QUAC.genind),], na.rm = TRUE) != 0)
# Number of captured wild alleles
length(which(QUAC.garden.alleles.present %in% QUAC.wild.alleles.present))
# Gardens capture 95.35% of wild alleles
(length(which(QUAC.garden.alleles.present %in% QUAC.wild.alleles.present))/length(QUAC.wild.alleles.present))*100

# %%% Singletons/doubletons removed %%%
# Create allele matrix with wild singletons/doubletons removed
# Take all rows of the genind tab object, but only the columns with sums greater than or equal to 3
QUAC.genind.SDR <- QUAC.genind@tab[,which(colSums(QUAC.genind@tab, na.rm=TRUE) >= 3)]
# Garden allele columns that are non-empty (present)
QUAC.garden.alleles.present.SDR <- which(colSums(QUAC.genind.SDR[1:QUAC.garden.inds,], na.rm = TRUE) != 0)
# Wild allele columns present
QUAC.wild.alleles.present.SDR <- which(colSums(QUAC.genind.SDR[QUAC.wild.index:nrow(QUAC.genind.SDR),], na.rm = TRUE) != 0)
# Number of captured wild alleles
length(which(QUAC.garden.alleles.present.SDR %in% QUAC.wild.alleles.present.SDR))
# Gardens capture 95.35% of wild alleles, with wild singletons/doubletons removed
(length(which(QUAC.garden.alleles.present.SDR %in% QUAC.wild.alleles.present.SDR))/length(QUAC.wild.alleles.present.SDR))*100
# So, there are no singleton/doubleton alleles (which(colSums(QUAC.genind@tab, na.rm=TRUE) <= 2) == 0)

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
# These fractions demonstrate the frequency of finding an allele when drawing out a single individual by
# taking the colSums of the matrix for wild individuals and dividing it by the number of individuals * 2 (diploids)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate wild allele frequencies
# Subset the rows to just look at wild samples, and the columns to only look at alleles present in wild samples
QUAC.wildAlleleFreqs <- (colSums(QUAC.genind@tab[QUAC.wild.index:nInd(QUAC.genind),QUAC.wild.alleles.present], na.rm = TRUE)/(QUAC.wild.inds*2))*100
# Calculate allelic capture for different categories, based on garden.alleles.present vector from before
# The lines below step through the buildup to the final calculation of the garden genetic capture (for very common alleles)
# Wild allele frequencies that are very common (greater than 10%)
which(QUAC.wildAlleleFreqs > 10)
# Vector of TRUE/FALSE for garden alleles being included within the very common wild alleles
QUAC.garden.alleles.present %in% which(QUAC.wildAlleleFreqs > 10)
# Garden alleles that ARE included within very common wild alleles (i.e. get TRUE values)
which(QUAC.garden.alleles.present %in% which(QUAC.wildAlleleFreqs > 10))
# NUMBER of garden alleles included within very common wild alleles
length(which(QUAC.garden.alleles.present %in% which(QUAC.wildAlleleFreqs > 10)))
# Number of garden alleles included within very common wild alleles, divided by number of very common alleles
length(which(QUAC.garden.alleles.present %in% which(QUAC.wildAlleleFreqs > 10)))/length(which(QUAC.wildAlleleFreqs > 10))*100
# Common alleles
length(which(QUAC.garden.alleles.present %in% which(QUAC.wildAlleleFreqs > 5)))/length(which(QUAC.wildAlleleFreqs > 5))*100
# Low frequency alleles
length(which(QUAC.garden.alleles.present %in% which(QUAC.wildAlleleFreqs < 10 & QUAC.wildAlleleFreqs > 1)))/length(which(QUAC.wildAlleleFreqs < 10 & QUAC.wildAlleleFreqs > 1))*100
# Rare alleles
length(which(QUAC.garden.alleles.present %in% which(QUAC.wildAlleleFreqs < 1)))/length(which(QUAC.wildAlleleFreqs < 1))*100

# %%% Singletons/doubletons removed %%%
# Calculate wild allele frequencies
# Take the colSums of the matrix for wild individuals and divide it by the number of individuals * 2 (diploids)
QUAC.wildAlleleFreqs.SDR <- (colSums(QUAC.genind.SDR[QUAC.wild.index:nrow(QUAC.genind.SDR),QUAC.wild.alleles.present.SDR], na.rm = TRUE)/(QUAC.wild.inds*2))*100
# Calculate allelic capture for different categories, based on garden.alleles.present vector from before
# The steps below show the buildup of calculating the capture in gardens (for very common alleles)
# Wild allele frequencies that are very common (greater than 10%)
which(QUAC.wildAlleleFreqs.SDR > 10)
# Vector of TRUE/FALSE for garden alleles being included within the very common wild alleles
QUAC.garden.alleles.present.SDR %in% which(QUAC.wildAlleleFreqs.SDR > 10)
# Garden alleles that ARE included within very common wild alleles (i.e. get TRUE values)
which(QUAC.garden.alleles.present.SDR %in% which(QUAC.wildAlleleFreqs.SDR > 10))
# NUMBER of garden alleles included within very common wild alleles
length(which(QUAC.garden.alleles.present.SDR %in% which(QUAC.wildAlleleFreqs.SDR > 10)))
# Number of garden alleles included within very common wild alleles, divided by number of very common alleles
length(which(QUAC.garden.alleles.present.SDR %in% which(QUAC.wildAlleleFreqs.SDR > 10)))/length(which(QUAC.wildAlleleFreqs.SDR > 10))*100
# Common alleles
length(which(QUAC.garden.alleles.present.SDR %in% which(QUAC.wildAlleleFreqs.SDR > 5)))/length(which(QUAC.wildAlleleFreqs.SDR > 5))*100
# Low frequency alleles
length(which(QUAC.garden.alleles.present.SDR %in% which(QUAC.wildAlleleFreqs.SDR < 10 & QUAC.wildAlleleFreqs.SDR > 1)))/length(which(QUAC.wildAlleleFreqs.SDR < 10 & QUAC.wildAlleleFreqs.SDR > 1))*100
# Rare alleles
length(which(QUAC.garden.alleles.present.SDR %in% which(QUAC.wildAlleleFreqs.SDR < 1)))/length(which(QUAC.wildAlleleFreqs.SDR < 1))*100

# # Old approach, using QUAC.genpop object
# QUAC.wildAlleleFreqs <- (colSums(QUAC.genpop.wild@tab)/(QUAC.wild.inds*2))*100
# # Max frequency is 49.5%; min is 0.25%
# cat(paste("Max wild allele frequency: ",max(QUAC.wildAlleleFreqs),"\n","Min wild allele frequency: ",min(QUAC.wildAlleleFreqs)))
# # Gardens capture:
# # 100% of very common alleles
# (length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs > 10)])))/length(which(QUAC.wildAlleleFreqs > 10)))*100
# # 100% of common alleles
# (length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs > 5)])))/length(which(QUAC.wildAlleleFreqs > 5)))*100
# # 77.4% of low frequency alleles
# (length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs < 10 & QUAC.wildAlleleFreqs > 1)])))/length(which(QUAC.wildAlleleFreqs < 10 & QUAC.wildAlleleFreqs > 1)))*100
# # 100% of rare alleles
# (length(which(colnames(QUAC.genpop.garden@tab) %in% names(QUAC.wildAlleleFreqs[which(QUAC.wildAlleleFreqs < 1)])))/length(which(QUAC.wildAlleleFreqs < 1)))*100

# # Demonstration of different alleles being retained based on different separation techniques
# # Subsetting wild samples using seppop/repool
# QUAC.genind.wild <- seppop(QUAC.genind)[2:6]
# # (Using repool will subset the resulting genind object to only include loci present in the combined populations)
# QUAC.genind.wild <- repool(QUAC.genind.wild$porterMt,QUAC.genind.wild$magazineMt,QUAC.genind.wild$pryorMt,
#                            QUAC.genind.wild$sugarloaf_midlandPeak,QUAC.genind.wild$kessler_shaleBarrenRidge)
# 
# # Subsetting wild samples using repool_new function from Sean
# repool_new <- function(genind_obj,vect_pops){
#   genind_obj_sep <- seppop(genind_obj)
#   genind_obj_merge <- genind_obj_sep[[vect_pops[1]]]
#   for (i in 1:(length(vect_pops)-1)) genind_obj_merge<-repool(genind_obj_merge,genind_obj_sep[[vect_pops[i+1]]])
#   genind_obj_merge
# }
# QUAC.genind.wild.NEW <- repool_new(QUAC.genind,2:6)
# 
# # NOTE: there are differences in allele matrices depending on how genind/genpop objects are constructed
# ncol(QUAC.genind@tab)
# ncol(QUAC.genpop@tab)
# # Without na.rm argument, colSums don't match across genind and genpop objects
# identical(colSums(QUAC.genind@tab),colSums(QUAC.genpop@tab))
# identical(colSums(QUAC.genind@tab, na.rm = TRUE),colSums(QUAC.genpop@tab, na.rm = TRUE))
# # Subset objects have NAs removed
# ncol(QUAC.genind.wild@tab)
# ncol(QUAC.genind.wild.NEW@tab)
# ncol(QUAC.genpop.wild@tab)

# RESAMPLING----
# Create a matrix of only wild individuals with present alleles
QUAC.wild.mat <- QUAC.genind@tab[QUAC.wild.index:nInd(QUAC.genind),QUAC.wild.alleles.present]
# Processing duplicates
# Only one duplicate, in QUAC Wild data set
which(rownames(QUAC.wild.mat) == "QUAC_W_DUP_SH_Q2121")
QUAC.wild.mat[65:66,1:4] # QUAC_W_SH_Q2121 and QUAC_W_DUP_SH_Q2121
length(which(QUAC.wild.mat[65,] != QUAC.wild.mat[66,]))
# 874/12,159 alleles (7.19%) are different between these duplicate samples
# For now, just removing duplicate sample
QUAC.wild.mat <- QUAC.wild.mat[-66,]

# Function for measuring wild allelic capture of samples
# The two arguments are the vector of frequencies of wild alleles, and the sample of the wild genind object
get.allele.cat <- function(freq.vector, sample.mat){
  # Total alleles
  # Determine how many alleles in the sample (i.e. greater than 0) are found in the frequency vector 
  total <- length(which(colSums(sample.mat, na.rm = TRUE) > 0) %in% freq.vector)/length(freq.vector)*100
  # Very common alleles (greater than 10%)
  v_com <- length(which(which(colSums(sample.mat, na.rm = TRUE) > 0) %in% which(freq.vector > 10)))/length(which(freq.vector > 10))*100
  # Common alleles (greater than 5%)
  com <- length(which(which(colSums(sample.mat, na.rm = TRUE) > 0) %in% which(freq.vector > 5)))/length(which(freq.vector > 5))*100
  # Low frequency alleles (between 1% and 10%)
  low_freq <- length(which(which(colSums(sample.mat, na.rm = TRUE) > 0) %in% which(freq.vector < 10 & freq.vector > 1)))/length(which(freq.vector < 10 & freq.vector > 1))*100
  # Rare alleles (less than 1%)
  rare <- length(which(which(colSums(sample.mat, na.rm = TRUE) > 0) %in% which(freq.vector < 1)))/length(which(freq.vector < 1))*100
  # Concatentate values to a vector, and return
  return(c(total,v_com,com,low_freq,rare))
}

# # Old function for measuring wild allelic capture of samples, matching strings
# # The two arguments are the vector of frequencies of wild alleles, and the sample of the wild genind object
# get.allele.cat.OLD <- function(freq.vector, sample.mat){
#   # Calculate percentages
#   # Total alleles
#   total <- (length(which(names(which(colSums(sample.mat, na.rm = TRUE)!=0)) %in% names(freq.vector)))/length(freq.vector))*100
#   # Very common alleles (greater than 10%)
#   v_com <- (length(which(names(which(colSums(sample.mat, na.rm = TRUE)!=0)) %in% names(freq.vector[which(freq.vector > 10)])))/length(which(freq.vector > 10)))*100
#   # Common alleles (greater than 5%)
#   com <- (length(which(names(which(colSums(sample.mat, na.rm = TRUE)!=0)) %in% names(freq.vector[which(freq.vector > 5)])))/length(which(freq.vector > 5)))*100
#   # Low frequency alleles (between 1% and 10%)
#   low_freq <- (length(which(names(which(colSums(sample.mat, na.rm = TRUE)!=0)) %in% names(freq.vector[which(freq.vector < 10 & freq.vector > 1)])))/length(which(freq.vector < 10 & freq.vector > 1)))*100
#   # Rare alleles (less than 1%)
#   rare <- (length(which(names(which(colSums(sample.mat, na.rm = TRUE)!=0)) %in% names(freq.vector[which(freq.vector < 1)])))/length(which(freq.vector < 1)))*100
#   # Concatentate values to a vector, and return
#   return(c(total,v_com,com,low_freq,rare))
# }

# Build a 3D array, with rows being sample numbers and columns being allele frequency categories
# 3rd dimension is replicates. nrow is number of individuals-1 because sample doesn't work for vectors of length 1
num_reps <- 5 # In Hoban et al. 2020 (where QUBO is considered), num_reps = 75,000
list_allele_cat <- c("tot","v_com","com","low_freq","rare")
samplingResults <- array(dim=c(nrow(QUAC.wild.mat)-1,length(list_allele_cat),num_reps))
colnames(samplingResults) <- list_allele_cat

# For each replicate (which is the third dimension, in the samplingResults array)...
for(i in 1:num_reps){
  # loop through sampling from 2 to the maximum number of wild individuals
  for(j in 2:nrow(QUAC.wild.mat)){
    # Create a sample of the wild allele matrix, of "j" size
    samp <- QUAC.wild.mat[sample(nrow(QUAC.wild.mat), size=j, replace = FALSE),]
    # Calculate how many alleles of each category that sample captures,
    # and place those percentages into the row of the samplingResults array
    samplingResults[j-1,,i] <- get.allele.cat(QUAC.wildAlleleFreqs,samp)
  }
}

# RECREATING ABOVE RESAMPLING CODE USING APPLY FUNCTIONS----

# First, trying to recreate the results of the inner loop (j), which samples the wild allele matrix
# and then calculates the captured allele frequencies of that sample
samp <- QUAC.wild.mat[sample(nrow(QUAC.wild.mat), size=2, replace = FALSE),]

# 2 approaches: use a lambda function within some kind of apply function,
# or create a new function and then try to apply that over the wild allele matrix

# Lambda function
# Works for 2, but not multiple values
NEW.samp <- apply(QUAC.wild.mat,2, function(x) sample(x, size=2, replace = FALSE))

NEW.samp <- apply(QUAC.wild.mat,2, function(x) sample(x, size=seq(from=2, to=nrow(x)), replace = FALSE))
# The 2 calls below still only use 2 as the sample size
NEW.samp <- apply(QUAC.wild.mat,2, function(x) sample(x, size=2:96, replace = FALSE))

NEW.samp <- apply(QUAC.wild.mat,2, function(x) sample(x, size=seq(from=2, to=nrow(QUAC.wild.mat)), replace = FALSE))

# New function
resample.test <- function(wild.matrix, sample.number, allele.frequencies){
  # Create a sample of the wild allele matrix, of "j" size
  # wild.sample <- wild.matrix[sample(nrow(wild.matrix), size=sample.number, replace = FALSE)]
  wild.sample <- wild.matrix[sample(wild.matrix, size=sample.number, replace = FALSE),]
  # Calculate how many alleles of each category that sample captures,
  # and output those percentages
  output <- get.allele.cat(allele.frequencies, wild.sample)
}

# Single example
test <- resample.test(QUAC.wild.mat, sample.number = 2, allele.frequencies = QUAC.wildAlleleFreqs)
# Again, this doesn't work, because I'm trying to pass a vector to the size argument of 'sample'...
test.resamplingMat <- apply(QUAC.wild.mat, 2, resample.test, sample.number=seq(from=2, to=nrow(QUAC.wild.mat)), allele.frequencies=QUAC.wildAlleleFreqs)
# Attempts using replicate function
replicate(96, sample(QUAC.wild.mat, size = 2:nrow(QUAC.wild.mat), replace = FALSE), simplify=FALSE)
replicate(96, sample(QUAC.wild.mat, replace = FALSE), simplify=FALSE)
yo <- replicate(96, QUAC.wild.mat[sample(nrow(QUAC.wild.mat), size=2:nrow(QUAC.wild.mat), replace = FALSE),], simplify=FALSE)
# Variable samples per matrix, but not ordered...size argument is itself called from sample?
yo <- replicate(96, QUAC.wild.mat[sample(nrow(QUAC.wild.mat), sample(2:nrow(QUAC.wild.mat), 1L), replace = FALSE),], simplify=FALSE)
# 2 samples per matrix
yo <- replicate(96, QUAC.wild.mat[sample(nrow(QUAC.wild.mat), 2:nrow(QUAC.wild.mat), replace = FALSE),], simplify=FALSE)
# 96 samples per matrix
yo <- replicate(96, QUAC.wild.mat[sample(nrow(QUAC.wild.mat), 96, replace = FALSE),], simplify=FALSE)

# (tapply/lambda function example)
input <- 1:10
grouping <- rep(letters[1:2], 5)
tapply(input, grouping, function(x) sd(x)/sqrt(length(x)))

# Plotting 
total_means <- apply(samplingResults[,1,], 1, mean)
total_sd <- apply(samplingResults[,1,], 1, sd)

v.com_means <- apply(samplingResults[,2,], 1, mean)
v.com_sd <- apply(samplingResults[,2,], 1, sd)

com_means <- apply(samplingResults[,3,], 1, mean)
com_sd <- apply(samplingResults[,3,], 1, sd)

lowfr_means <- apply(samplingResults[,4,], 1, mean)
lowfr_sd <- apply(samplingResults[,4,], 1, sd)

rare_means <- apply(samplingResults[,5,], 1, mean)
rare_sd <- apply(samplingResults[,5,], 1, sd)

# plotColors <- brewer.pal(n=5, name="Dark2")
plotColors <- c("red","firebrick","darkorange3","coral","deeppink4")
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend("bottomright", inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.5)

# Zoomed in 
plot(v.com_means, ylim=c(0,110), xlim=c(0,5), col=plotColors[1], pch=16, main="QUAC: Resampling Curve (10 replicates)",
     xlab="Number of samples", ylab="Mean percent captured")
points(com_means, col=plotColors[2], pch=16)
points(lowfr_means, col=plotColors[3], pch=16)
points(rare_means, col=plotColors[4], pch=16)
legend("bottomright", inset = 0.05, legend = c("Very common: >10%","Common: >5%","Low frequency: 1--10%", "Rare: <1%"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2)

# How many samples reach 95% genetic capture threshold?
# Extract just the results for all alleles
samplingResults[,1,]
# Below line averages across reps (supposedly? Seems like MARGIN should be 2, not 1...)
apply(samplingResults[,1,],1,mean)
# 60 samples captures 95% of total genetic diversity

# Heterozygosity in wild and garden samples----
# Pull in genpop object 
QUAC.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Assign populations
pop(QUAC.genind) <- factor(read.table("../../../QUAC_popmap3", header=FALSE)[,2])

Hs(QUAC.genind)
mean(summary(QUAC.genind)[7]$Hexp)
sd(summary(QUAC.genind)[7]$Hexp)

QUAC_hexp_mean_df <- as.data.frame(rbind(mean(QUAC_hexp[[1]][,1]), mean(QUAC_hexp[[2]][,1])))
QUAC_hexp_mean_df$pop_type <- NA
QUAC_hexp_mean_df$pop_type <- c("Garden", "Wild")


# Barplot for expected heterozygosity, microsattelite markers (sourced from 3_QUAC_garden_wild_comparison.R)
barplot(QUAC_hexp_mean_df[,1], beside = TRUE, 
        ylim = c(0,1), col = c("darkgreen", "darkseagreen1"),
        names = c("Garden", "Wild"), 
        main = "Heterozygosity: SNPs", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0)

# Barplot for expected heterozygosity, SNP markers
barplot(Hs(QUAC.genind), beside = TRUE, 
        ylim = c(0,1), col = c("darkgreen", "darkseagreen1"),
        names = c("Garden", "Wild"), 
        main = "Heterozygosity: Microsatellites", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0)


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
