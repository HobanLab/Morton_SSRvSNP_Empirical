# %%%%%%%%%%%%%%%%%%
# %%% RESAMPLING %%%
# %%%%%%%%%%%%%%%%%%

# This script generates Resampling arrays, and plots their results, 
# for Quercus acerifolia (QUAC; optimized Stacks de novo assembly, m 7, M/n 4, gt-alpha 0.01) 
# and Quercus boyntonii (QUBO; GSNAP4 alignment with Quercus robur reference) NextRAD samples

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
QUAC.wild.mat <- 
  QUAC.R0_NOMAF.genind@tab[QUAC.wild,which(colSums(QUAC.R0_NOMAF.genind@tab[QUAC.wild,], na.rm = TRUE) > 0)]
# Generate wild allele frequency vector
QUAC_wildFreqs <- colSums(QUAC.wild.mat, na.rm = TRUE)/(length(QUAC.wild)*2)*100

# CREATE SAMPLING RESULTS ARRAY ----
# Build a 3D array, with rows being sample numbers and columns being allele frequency categories
# 3rd dimension is replicates. nrow is number of individuals-1 
# (because sample doesn't work for vectors of length 1)
num_reps <- 5 # In Hoban et al. 2020 (where QUBO is considered), num_reps = 75,000
list_allele_cat <- c("tot","v_com","com","low_freq","rare")
samplingResults_QUAC <- array(dim=c(nrow(QUAC.wild.mat)-1,length(list_allele_cat),num_reps))
colnames(samplingResults_QUAC) <- list_allele_cat

# For each replicate (which is the third dimension, in the samplingResults array)...
for(i in 1:num_reps){
  # loop through sampling from 2 to the maximum number of wild individuals
  for(j in 2:nrow(QUAC.wild.mat)){
    # Create a sample of the wild allele matrix, of "j" size
    samp <- QUAC.wild.mat[sample(nrow(QUAC.wild.mat), size=j, replace = FALSE),]
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
min_95_QUAC <- min(which(apply(samplingResults_QUAC[,1,],1,mean) > 95)); min_95_QUAC

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
# Pick colors, with transparency for values other than Total
plotColors <- brewer.pal(n=5, name="Dark2")
# plotColors <- c("red","firebrick","darkorange3","coral","deeppink4")
plotColors[2:5] <- alpha(plotColors[2:5], 0.5)
# Plots all sets of points onto single graph, as well as 95% threshold line
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture",
     main="QUAC (R0, NOMAF) Resampling")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend(x=83, y=20.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.4)
# Lines for 95% threshold
abline(h=95, col="black", lty=3)
abline(v=min_95_QUAC, col="black")

# %%%% QUBO %%%% ----
# READ IN GENIND FILE (QUBO DNFA; R0, min-maf=0; 1 SNP/locus; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_TwoPops/"
setwd(genpop.filePath)
QUBO.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# Number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUBO.R0_NOMAF.genind)

# CREATE WILD SAMPLE MATRIX AND ALLELE FREQUENCY VECTOR ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUBO.wild <- seq(from=length(which(pop(QUBO.R0_NOMAF.genind)=="garden"))+1, to=nInd(QUBO.R0_NOMAF.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUBO.wild.mat <- 
  QUBO.R0_NOMAF.genind@tab[QUBO.wild,which(colSums(QUBO.R0_NOMAF.genind@tab[QUBO.wild,], na.rm = TRUE) > 0)]
# Generate wild allele frequency vector
QUBO_wildFreqs <- colSums(QUBO.wild.mat, na.rm = TRUE)/(length(QUBO.wild)*2)*100

# CREATE SAMPLING RESULTS ARRAY ----
# Build a 3D array, with rows being sample numbers and columns being allele frequency categories
# 3rd dimension is replicates. nrow is number of individuals-1 
# (because sample doesn't work for vectors of length 1)
num_reps <- 5 # In Hoban et al. 2020 (where QUBO is considered), num_reps = 75,000
list_allele_cat <- c("tot","v_com","com","low_freq","rare")
samplingResults_QUBO <- array(dim=c(nrow(QUBO.wild.mat)-1,length(list_allele_cat),num_reps))
colnames(samplingResults_QUBO) <- list_allele_cat

# For each replicate (which is the third dimension, in the samplingResults array)...
for(i in 1:num_reps){
  # loop through sampling from 2 to the maximum number of wild individuals
  for(j in 2:nrow(QUBO.wild.mat)){
    # Create a sample of the wild allele matrix, of "j" size
    samp <- QUBO.wild.mat[sample(nrow(QUBO.wild.mat), size=j, replace = FALSE),]
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
min_95_QUBO <- min(which(apply(samplingResults_QUBO[,1,],1,mean) > 95)); min_95_QUBO

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
# Pick colors, with transparency for values other than Total
plotColors <- brewer.pal(n=5, name="Dark2")
# plotColors <- c("red","firebrick","darkorange3","coral","deeppink4")
plotColors[2:5] <- alpha(plotColors[2:5], 0.5)
# Plots all sets of points onto single graph, as well as 95% threshold line
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture",
     main="QUBO (R0, NOMAF) Resampling")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend(x=83, y=20.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.4)
abline(h=95, col="black", lty=3)
abline(v=min_95_QUBO, col="black")
