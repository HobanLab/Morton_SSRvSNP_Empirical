# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% RESAMPLING WITH MICROSATELLITES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script generates Resampling arrays, and plots their results, 
# for Quercus acerifolia (QUAC; optimized Stacks de novo assembly, m 7, M/n 4, gt-alpha 0.01) 
# and Quercus boyntonii (QUBO; GSNAP4 alignment with Quercus robur reference) NextRAD samples

# In addition to procesing SNP datasets, this script reads in the QUAC and QUBO microsatellite
# genind files as well, to compare results between marker types

# Resampling arrays have 3 dimensions: rows are sample numbers, columns are allele frequncy categories,
# and Z dimension is replicates. nrow of the array is the number of individuals-1, because the sample
# function doesn't work for vectors of length 1

library(adegenet)
library(RColorBrewer)
library(scales)

# %%%% FUNCTIONS %%%% ----
# Read in relevant functions required for resampling analyses
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)
source("exSituRepresentation/functions_exSituRepresentation.R")

# %%%% QUAC %%%% ----
# ---- MSATs ----
# READ IN GENIND FILE (GCC_QUAC_ZAIN repo; QUAC_wK_garden_wild_clean.gen) ----
genpop.filePath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/Garden_Wild/"
setwd(genpop.filePath)
QUAC.MSAT.genind <- read.genepop("QUAC_wK_garden_wild_clean.gen", ncode = 3)
# Correct popNames: pop1 is Garden, pop2 is Wild
pop(QUAC.MSAT.genind) <- gsub("pop1", "garden", pop(QUAC.MSAT.genind))
pop(QUAC.MSAT.genind) <- gsub("pop2", "wild", pop(QUAC.MSAT.genind))
# Number of loci
nLoc(QUAC.MSAT.genind)

# CREATE WILD SAMPLE MATRIX ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUAC.MSAT.wild <- seq(from=length(which(pop(QUAC.MSAT.genind)=="garden"))+1, to=nInd(QUAC.MSAT.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUAC.MSAT.wildMat <- 
  QUAC.MSAT.genind@tab[QUAC.MSAT.wild,which(colSums(QUAC.MSAT.genind@tab[QUAC.MSAT.wild,], na.rm = TRUE) > 0)]

# CREATE SAMPLING RESULTS ARRAY ----
# Specify number of replicates
num_reps <- 5 
# Build microsatellite sampling array
samplingResults_QUAC.MSAT <- 
  sapply(1:num_reps, function(a) exSitu_Resample(wildMat=QUAC.MSAT.wildMat), simplify = "array")

# CALCULATE MEANS AND PLOT ----

# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min_95_QUAC.MSAT <- min(which(apply(samplingResults_QUAC.MSAT[,1,],1,mean) > 95)); min_95_QUAC.MSAT

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUAC.MSAT[,1,], 1, mean)
total_sd <- apply(samplingResults_QUAC.MSAT[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUAC.MSAT[,2,], 1, mean)
v.com_sd <- apply(samplingResults_QUAC.MSAT[,2,], 1, sd)

com_means <- apply(samplingResults_QUAC.MSAT[,3,], 1, mean)
com_sd <- apply(samplingResults_QUAC.MSAT[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUAC.MSAT[,4,], 1, mean)
lowfr_sd <- apply(samplingResults_QUAC.MSAT[,4,], 1, sd)

rare_means <- apply(samplingResults_QUAC.MSAT[,5,], 1, mean)
rare_sd <- apply(samplingResults_QUAC.MSAT[,5,], 1, sd)

# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4))
# Plots all sets of points onto single graph, as well as 95% threshold line
plotColors <- c("red","firebrick","darkorange3","coral","deeppink4")
plotColors[2:5] <- alpha(plotColors[2:5], 0.5)
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture",
     main="Microsatellites")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend(x=136.1294, y=60, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.3)
# Lines for 95% threshold
abline(h=95, col="black", lty=3)
abline(v=min_95_QUAC.MSAT, col="black")

# ---- SNPs ----
# READ IN GENIND FILE (QUAC DNFA; R0, min-maf=0; first SNP per locus; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops//"
setwd(genpop.filePath)
QUAC.SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUAC.SNP.genind)

# CREATE WILD SAMPLE MATRIX AND ALLELE FREQUENCY VECTOR ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUAC.wild <- seq(from=length(which(pop(QUAC.SNP.genind)=="garden"))+1, to=nInd(QUAC.SNP.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUAC.SNP.wildMat <- 
  QUAC.SNP.genind@tab[QUAC.wild,which(colSums(QUAC.SNP.genind@tab[QUAC.wild,], na.rm = TRUE) > 0)]

# CREATE SAMPLING RESULTS ARRAY ----
# Specify number of replicates
num_reps <- 5 
# Build SNP sampling array
samplingResults_QUAC.SNP <- 
  sapply(1:num_reps, function(a) exSitu_Resample(wildMat=QUAC.SNP.wildMat), simplify = "array")

# CALCULATE MEANS AND PLOT ----
# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min_95_QUAC.SNP <- min(which(apply(samplingResults_QUAC.SNP[,1,],1,mean) > 95)); min_95_QUAC.SNP

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUAC.SNP[,1,], 1, mean)
total_sd <- apply(samplingResults_QUAC.SNP[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUAC.SNP[,2,], 1, mean)
v.com_sd <- apply(samplingResults_QUAC.SNP[,2,], 1, sd)

com_means <- apply(samplingResults_QUAC.SNP[,3,], 1, mean)
com_sd <- apply(samplingResults_QUAC.SNP[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUAC.SNP[,4,], 1, mean)
lowfr_sd <- apply(samplingResults_QUAC.SNP[,4,], 1, sd)

rare_means <- apply(samplingResults_QUAC.SNP[,5,], 1, mean)
rare_sd <- apply(samplingResults_QUAC.SNP[,5,], 1, sd)
# Plots all sets of points onto single graph, as well as 95% threshold line
plotColors <- c("red","firebrick","darkorange3","coral","deeppink4")
plotColors[2:5] <- alpha(plotColors[2:5], 0.5)
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture",
     main="SNPs")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend(x=83, y=60, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.3)
# Lines for 95% threshold
abline(h=95, col="black", lty=3)
abline(v=min_95_QUAC.SNP, col="black")

# %%%% QUBO %%%% ----
# ---- MSATs ----
# READ IN GENIND FILE (Southeast Oaks repo; Qb_total.gen file) ----
genpop.filePath <- 
  "~/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
setwd(genpop.filePath)
QUBO.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_total.gen"), ncode=3)
# Correct popNames: last population (IMLS4_MP1_IMLS336_C05) is Garden; rest are Wild
pop(QUBO.MSAT.genind) <- gsub("IMLS4_MP1_IMLS336_C05", "garden", pop(QUBO.MSAT.genind))
# Number of loci
nLoc(QUBO.MSAT.genind)

# CREATE WILD SAMPLE MATRIX ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUBO.MSAT.wild <- which(pop(QUBO.MSAT.genind) != "garden")
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUBO.MSAT.wildMat <- 
  QUBO.MSAT.genind@tab[QUBO.MSAT.wild,which(colSums(QUBO.MSAT.genind@tab[QUBO.MSAT.wild,], na.rm = TRUE) > 0)]

# CREATE SAMPLING RESULTS ARRAY ----
# Specify number of replicates
num_reps <- 5 
# Build microsatellite sampling array
samplingResults_QUBO.MSAT <- 
  sapply(1:num_reps, function(a) exSitu_Resample(wildMat=QUBO.MSAT.wildMat), simplify = "array")

# CALCULATE MEANS AND PLOT ----

# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min_95_QUBO.MSAT <- min(which(apply(samplingResults_QUBO.MSAT[,1,],1,mean) > 95)); min_95_QUBO.MSAT

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUBO.MSAT[,1,], 1, mean)
total_sd <- apply(samplingResults_QUBO.MSAT[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUBO.MSAT[,2,], 1, mean)
v.com_sd <- apply(samplingResults_QUBO.MSAT[,2,], 1, sd)

com_means <- apply(samplingResults_QUBO.MSAT[,3,], 1, mean)
com_sd <- apply(samplingResults_QUBO.MSAT[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUBO.MSAT[,4,], 1, mean)
lowfr_sd <- apply(samplingResults_QUBO.MSAT[,4,], 1, sd)

rare_means <- apply(samplingResults_QUBO.MSAT[,5,], 1, mean)
rare_sd <- apply(samplingResults_QUBO.MSAT[,5,], 1, sd)

# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1), oma=rep(0.3,4))
# Plots all sets of points onto single graph, as well as 95% threshold line
# plotColors <- brewer.pal(n=5, name="Dark2")
plotColors <- c("red","firebrick","darkorange3","coral","deeppink4")
plotColors[2:5] <- alpha(plotColors[2:5], 0.5)
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture",
     main="Microsatellites")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend(x=200, y=65, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.3)
# Lines for 95% threshold
abline(h=95, col="black", lty=3)
abline(v=min_95_QUBO.MSAT, col="black")

# ---- SNPs ----
# READ IN GENIND FILE (QUBO GSNAP4 alignment; R0, min-maf=0; first SNP per locus; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R0_NOMAF.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# Number of loci
nLoc(QUBO.R0_NOMAF.genind)

# CREATE WILD SAMPLE MATRIX ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUBO.wild <- seq(from=length(which(pop(QUBO.R0_NOMAF.genind)=="garden"))+1, to=nInd(QUBO.R0_NOMAF.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUBO.wildMat <- 
  QUBO.R0_NOMAF.genind@tab[QUBO.wild,which(colSums(QUBO.R0_NOMAF.genind@tab[QUBO.wild,], na.rm = TRUE) > 0)]

# CREATE SAMPLING RESULTS ARRAY ----
num_reps <- 5 
# Build SNP sampling array
samplingResults_QUBO.SNP <- 
  sapply(1:num_reps, function(a) exSitu_Resample(wildMat=QUBO.SNP.wildMat), simplify = "array")

# CALCULATE MEANS AND PLOT ----
# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min_95_QUBO.SNP <- min(which(apply(samplingResults_QUBO.SNP[,1,],1,mean) > 95)); min_95_QUBO.SNP

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUBO.SNP[,1,], 1, mean)
total_sd <- apply(samplingResults_QUBO.SNP[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUBO.SNP[,2,], 1, mean)
v.com_sd <- apply(samplingResults_QUBO.SNP[,2,], 1, sd)

com_means <- apply(samplingResults_QUBO.SNP[,3,], 1, mean)
com_sd <- apply(samplingResults_QUBO.SNP[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUBO.SNP[,4,], 1, mean)
lowfr_sd <- apply(samplingResults_QUBO.SNP[,4,], 1, sd)

rare_means <- apply(samplingResults_QUBO.SNP[,5,], 1, mean)
rare_sd <- apply(samplingResults_QUBO.SNP[,5,], 1, sd)
# Plots all sets of points onto single graph, as well as 95% threshold line
# plotColors <- brewer.pal(n=5, name="Dark2")
plotColors <- c("red","firebrick","darkorange3","coral","deeppink4")
plotColors[2:5] <- alpha(plotColors[2:5], 0.5)
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture",
     main="SNPs")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend(x=80, y=65, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.3)
# Lines for 95% threshold
abline(h=95, col="black", lty=3)
abline(v=min_95_QUBO.SNP, col="black")
