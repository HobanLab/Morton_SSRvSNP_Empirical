# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% RESAMPLING WITH MICROSATELLITES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script generates Resampling arrays, and plots their results, 
# for Quercus acerifolia (QUAC; optimized Stacks de novo assembly, m 7, M/n 4, gt-alpha 0.01) 
# and Quercus boyntonii (QUBO; GSNAP4 alignment with Quercus robur reference) NextRAD samples

# In addition to procesing SNP datasets, this script reads in the QUAC and QUBO microsatellite
# genind files as well, to compare results between marker types

library(adegenet)
library(RColorBrewer)
library(scales)

# %%%% FUNCTIONS %%%% ----
# Function for reporting capture rates, using a sample matrix and a vector of allele frequencies
get.allele.cat.NEW <- function(freq.vector, sample.mat){
  # Total alleles
  # Determine how many alleles in the sample (i.e. greater than 0) are found in the frequency vector 
  total <- length(which(names(which(freq.vector > 0)) %in% names(which(colSums(sample.mat, na.rm = TRUE) > 0))))/length(which(freq.vector > 0))*100
  # Very common alleles (greater than 10%)
  v_com <- length(which(names(which(freq.vector > 10)) %in% names(which(colSums(sample.mat, na.rm = TRUE) > 0))))/length(which(freq.vector > 10))*100
  # Common alleles (greater than 5%)
  com <- length(which(names(which(freq.vector > 5)) %in% names(which(colSums(sample.mat, na.rm = TRUE) > 0))))/length(which(freq.vector > 5))*100
  # Low frequency alleles (between 1% and 10%)
  low_freq <- length(which(names(which(freq.vector < 10 & freq.vector > 1)) %in% names(which(colSums(sample.mat, na.rm = TRUE) > 0))))/length(which(freq.vector < 10 & freq.vector > 1))*100
  # Rare alleles (less than 1%)
  rare <- length(which(names(which(freq.vector < 1 & freq.vector > 0)) %in% names(which(colSums(sample.mat, na.rm = TRUE) > 0))))/length(which(freq.vector < 1 & freq.vector > 0))*100
  # Concatenate values to a vector, and return
  return(c(total,v_com,com,low_freq,rare))
}

# %%%% QUAC %%%% ----
# ---- MSATs ----
# READ IN GENIND FILE (QUAC_insitu_exsitu repo; QUAC_garden_wild_clean.gen) ----
genpop.filePath <- 
  "~/Documents/peripheralProjects/QUAC_insitu_exsitu/QUAC_data_files/QUAC_adegenet_files/Garden_Wild/"
setwd(genpop.filePath)
QUAC.MSAT.genind <- read.genepop("QUAC_garden_wild_clean.gen", quiet = TRUE, ncode = 3)
# Correct popNames: pop1 is Garden, pop2 is Wild
pop(QUAC.MSAT.genind) <- gsub("pop1", "garden", pop(QUAC.MSAT.genind))
pop(QUAC.MSAT.genind) <- gsub("pop2", "wild", pop(QUAC.MSAT.genind))
pop(QUAC.MSAT.genind)
# Number of loci
nLoc(QUAC.MSAT.genind)

# CREATE WILD SAMPLE MATRIX AND ALLELE FREQUENCY VECTOR ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUAC.MSAT.wild <- seq(from=length(which(pop(QUAC.MSAT.genind)=="garden"))+1, to=nInd(QUAC.MSAT.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUAC.MSAT.wildMat <- 
  QUAC.MSAT.genind@tab[QUAC.MSAT.wild,which(colSums(QUAC.MSAT.genind@tab[QUAC.MSAT.wild,], na.rm = TRUE) > 0)]
# Generate wild allele frequency vector
QUAC.MSAT_wildFreqs <- colSums(QUAC.MSAT.wildMat, na.rm = TRUE)/(length(QUAC.MSAT.wild)*2)*100

# CREATE SAMPLING RESULTS ARRAY ----
# Build a 3D array, with rows being sample numbers and columns being allele frequency categories
# 3rd dimension is replicates. nrow is number of individuals-1 
# (because sample doesn't work for vectors of length 1)
num_reps <- 50 # In Hoban et al. 2020 (where QUBO is considered), num_reps = 75,000
list_allele_cat <- c("tot","v_com","com","low_freq","rare")
samplingResults_QUAC.MSAT <- array(dim=c(nrow(QUAC.MSAT.wildMat)-1,length(list_allele_cat),num_reps))
colnames(samplingResults_QUAC.MSAT) <- list_allele_cat

# For each replicate (which is the third dimension, in the samplingResults array)...
for(i in 1:num_reps){
  # loop through sampling from 2 to the maximum number of wild individuals
  for(j in 2:nrow(QUAC.MSAT.wildMat)){
    # Create a sample of the wild allele matrix, of "j" size
    samp <- QUAC.MSAT.wildMat[sample(nrow(QUAC.MSAT.wildMat), size=j, replace = FALSE),]
    # Calculate how many alleles of each category that sample captures,
    # and place those percentages into the row of the samplingResults array
    samplingResults_QUAC.MSAT[j-1,,i] <- get.allele.cat.NEW(QUAC.MSAT_wildFreqs,samp)
  }
}
str(samplingResults_QUAC.MSAT)

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
# READ IN GENIND FILE (QUAC DNFA; R0, min-maf=0; 1 SNP/locus; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_TwoPops/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUAC.R0_NOMAF.genind)

# CREATE WILD SAMPLE MATRIX AND ALLELE FREQUENCY VECTOR ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUAC.wild <- seq(from=length(which(pop(QUAC.R0_NOMAF.genind)=="garden"))+1, to=nInd(QUAC.R0_NOMAF.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUAC.wildMat <- 
  QUAC.R0_NOMAF.genind@tab[QUAC.wild,which(colSums(QUAC.R0_NOMAF.genind@tab[QUAC.wild,], na.rm = TRUE) > 0)]
# Generate wild allele frequency vector
QUAC_wildFreqs <- colSums(QUAC.wildMat, na.rm = TRUE)/(length(QUAC.wild)*2)*100

# CREATE SAMPLING RESULTS ARRAY ----
# Build a 3D array, with rows being sample numbers and columns being allele frequency categories
# 3rd dimension is replicates. nrow is number of individuals-1 
# (because sample doesn't work for vectors of length 1)
num_reps <- 5 # In Hoban et al. 2020 (where QUBO is considered), num_reps = 75,000
list_allele_cat <- c("tot","v_com","com","low_freq","rare")
samplingResults_QUAC <- array(dim=c(nrow(QUAC.wildMat)-1,length(list_allele_cat),num_reps))
colnames(samplingResults_QUAC) <- list_allele_cat

# For each replicate (which is the third dimension, in the samplingResults array)...
for(i in 1:num_reps){
  # loop through sampling from 2 to the maximum number of wild individuals
  for(j in 2:nrow(QUAC.wildMat)){
    # Create a sample of the wild allele matrix, of "j" size
    samp <- QUAC.wildMat[sample(nrow(QUAC.wildMat), size=j, replace = FALSE),]
    # Calculate how many alleles of each category that sample captures,
    # and place those percentages into the row of the samplingResults array
    samplingResults_QUAC[j-1,,i] <- get.allele.cat.NEW(QUAC_wildFreqs,samp)
  }
}
str(samplingResults_QUAC)

# CALCULATE MEANS AND PLOT ----
# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min_95_QUAC.SNP <- min(which(apply(samplingResults_QUAC[,1,],1,mean) > 95)); min_95_QUAC.SNP

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUAC[,1,], 1, mean)
total_sd <- apply(samplingResults_QUAC[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUAC[,2,], 1, mean)
v.com_sd <- apply(samplingResults_QUAC[,2,], 1, sd)

com_means <- apply(samplingResults_QUAC[,3,], 1, mean)
com_sd <- apply(samplingResults_QUAC[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUAC[,4,], 1, mean)
lowfr_sd <- apply(samplingResults_QUAC[,4,], 1, sd)

rare_means <- apply(samplingResults_QUAC[,5,], 1, mean)
rare_sd <- apply(samplingResults_QUAC[,5,], 1, sd)
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
QUBO.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_total.gen"), quiet = TRUE, ncode=3)
# Correct popNames: last population (IMLS4_MP1_IMLS336_C05) is Garden; rest are Wild
pop(QUBO.MSAT.genind) <- gsub("IMLS4_MP1_IMLS336_C05", "garden", pop(QUBO.MSAT.genind))
# Number of loci
nLoc(QUBO.MSAT.genind)

# CREATE WILD SAMPLE MATRIX AND ALLELE FREQUENCY VECTOR ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUBO.MSAT.wild <- which(pop(QUBO.MSAT.genind) != "garden")
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUBO.MSAT.wildMat <- 
  QUBO.MSAT.genind@tab[QUBO.MSAT.wild,which(colSums(QUBO.MSAT.genind@tab[QUBO.MSAT.wild,], na.rm = TRUE) > 0)]
# Generate wild allele frequency vector
QUBO.MSAT_wildFreqs <- colSums(QUBO.MSAT.wildMat, na.rm = TRUE)/(length(QUBO.MSAT.wild)*2)*100

# CREATE SAMPLING RESULTS ARRAY ----
# Build a 3D array, with rows being sample numbers and columns being allele frequency categories
# 3rd dimension is replicates. nrow is number of individuals-1 
# (because sample doesn't work for vectors of length 1)
num_reps <- 50 # In Hoban et al. 2020 (where QUBO is considered), num_reps = 75,000
list_allele_cat <- c("tot","v_com","com","low_freq","rare")
samplingResults_QUBO.MSAT <- array(dim=c(nrow(QUBO.MSAT.wildMat)-1,length(list_allele_cat),num_reps))
colnames(samplingResults_QUBO.MSAT) <- list_allele_cat

# For each replicate (which is the third dimension, in the samplingResults array)...
for(i in 1:num_reps){
  # loop through sampling from 2 to the maximum number of wild individuals
  for(j in 2:nrow(QUBO.MSAT.wildMat)){
    # Create a sample of the wild allele matrix, of "j" size
    samp <- QUBO.MSAT.wildMat[sample(nrow(QUBO.MSAT.wildMat), size=j, replace = FALSE),]
    # Calculate how many alleles of each category that sample captures,
    # and place those percentages into the row of the samplingResults array
    samplingResults_QUBO.MSAT[j-1,,i] <- get.allele.cat.NEW(QUBO.MSAT_wildFreqs,samp)
  }
}
str(samplingResults_QUBO.MSAT)

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
# READ IN GENIND FILE (QUBO GSNAP4 alignment; R0, min-maf=0; 1 SNP/locus; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_TwoPops/"
setwd(genpop.filePath)
QUBO.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# Number of loci
nLoc(QUBO.R0_NOMAF.genind)

# CREATE WILD SAMPLE MATRIX AND ALLELE FREQUENCY VECTOR ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUBO.wild <- seq(from=length(which(pop(QUBO.R0_NOMAF.genind)=="garden"))+1, to=nInd(QUBO.R0_NOMAF.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUBO.wildMat <- 
  QUBO.R0_NOMAF.genind@tab[QUBO.wild,which(colSums(QUBO.R0_NOMAF.genind@tab[QUBO.wild,], na.rm = TRUE) > 0)]
# Generate wild allele frequency vector
QUBO_wildFreqs <- colSums(QUBO.wildMat, na.rm = TRUE)/(length(QUBO.wild)*2)*100

# CREATE SAMPLING RESULTS ARRAY ----
# Build a 3D array, with rows being sample numbers and columns being allele frequency categories
# 3rd dimension is replicates. nrow is number of individuals-1 
# (because sample doesn't work for vectors of length 1)
num_reps <- 5 # In Hoban et al. 2020 (where QUBO is considered), num_reps = 75,000
list_allele_cat <- c("tot","v_com","com","low_freq","rare")
samplingResults_QUBO <- array(dim=c(nrow(QUBO.wildMat)-1,length(list_allele_cat),num_reps))
colnames(samplingResults_QUBO) <- list_allele_cat

# For each replicate (which is the third dimension, in the samplingResults array)...
for(i in 1:num_reps){
  # loop through sampling from 2 to the maximum number of wild individuals
  for(j in 2:nrow(QUBO.wildMat)){
    # Create a sample of the wild allele matrix, of "j" size
    samp <- QUBO.wildMat[sample(nrow(QUBO.wildMat), size=j, replace = FALSE),]
    # Calculate how many alleles of each category that sample captures,
    # and place those percentages into the row of the samplingResults array
    samplingResults_QUBO[j-1,,i] <- get.allele.cat.NEW(QUBO_wildFreqs,samp)
  }
}
str(samplingResults_QUBO)

# CALCULATE MEANS AND PLOT ----
# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min_95_QUBO.SNP <- min(which(apply(samplingResults_QUBO[,1,],1,mean) > 95)); min_95_QUBO.SNP

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUBO[,1,], 1, mean)
total_sd <- apply(samplingResults_QUBO[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUBO[,2,], 1, mean)
v.com_sd <- apply(samplingResults_QUBO[,2,], 1, sd)

com_means <- apply(samplingResults_QUBO[,3,], 1, mean)
com_sd <- apply(samplingResults_QUBO[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUBO[,4,], 1, mean)
lowfr_sd <- apply(samplingResults_QUBO[,4,], 1, sd)

rare_means <- apply(samplingResults_QUBO[,5,], 1, mean)
rare_sd <- apply(samplingResults_QUBO[,5,], 1, sd)
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
