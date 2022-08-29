# %%%%%%%%%%%%%%%%%%
# %%% RESAMPLING %%%
# %%%%%%%%%%%%%%%%%%

# This script generates Resampling arrays, and plots their results, 
# for Quercus acerifolia (QUAC; optimized Stacks de novo assembly, m 7, M/n 4, gt-alpha 0.01) 
# and Quercus boyntonii (QUBO; GSNAP4 alignment with Quercus robur reference) NextRAD samples

# Resampling arrays have 3 dimensions: rows are sample numbers, columns are allele frequncy categories,
# and Z dimension is replicates. nrow of the array is the number of individuals-1, because the sample
# function doesn't work for vectors of length 1

library(adegenet)
library(RColorBrewer)
library(scales)
library(parallel)

# %%%% FUNCTIONS AND PARALLELIZATION %%%% ----
# Read in relevant functions required for resampling analyses
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)
source("exSituRepresentation/functions_exSituRepresentation.R")
# Set up relevant cores
num_cores <- detectCores() - 1 ; cl <- makeCluster(num_cores)

# %%%% QUAC %%%% ----
# ---- FIRST SNP PER LOCUS ----
# READ IN GENIND FILE (QUAC DNFA; R0, min-maf=0; first SNP per locus; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.1SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.1SNP.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUAC.1SNP.genind)

# CREATE WILD SAMPLE MATRIX ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUAC.wild <- seq(from=length(which(pop(QUAC.1SNP.genind)=="garden"))+1, to=nInd(QUAC.1SNP.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUAC.1SNP.wildMat <- 
  QUAC.1SNP.genind@tab[QUAC.wild,which(colSums(QUAC.1SNP.genind@tab[QUAC.wild,], na.rm = TRUE) > 0)]

# CREATE SAMPLING RESULTS ARRAY ----
# Specify number of replicates
num_reps <- 50 
# Export relevant functions and variables
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "QUAC.1SNP.wildMat"))

# Run resampling in parallel, to generate an array
samplingResults_QUAC.1SNP <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(wildMat=QUAC.1SNP.wildMat), simplify = "array")
# Close cores
# stopCluster(cl)
# Examine sampling array
str(samplingResults_QUAC.1SNP)

# CALCULATE MEANS AND PLOT ----
# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min95_QUAC.1SNP <- min(which(apply(samplingResults_QUAC.1SNP[,1,],1,mean) > 95)); min95_QUAC.1SNP

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUAC.1SNP[,1,], 1, mean)
total_sd <- apply(samplingResults_QUAC.1SNP[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUAC.1SNP[,2,], 1, mean)
v.com_sd <- apply(samplingResults_QUAC.1SNP[,2,], 1, sd)

com_means <- apply(samplingResults_QUAC.1SNP[,3,], 1, mean)
com_sd <- apply(samplingResults_QUAC.1SNP[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUAC.1SNP[,4,], 1, mean)
lowfr_sd <- apply(samplingResults_QUAC.1SNP[,4,], 1, sd)

rare_means <- apply(samplingResults_QUAC.1SNP[,5,], 1, mean)
rare_sd <- apply(samplingResults_QUAC.1SNP[,5,], 1, sd)
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

# ---- MICROHAPLOTYPE SNPS ----
# READ IN GENIND FILE (QUAC DNFA; R0, min-maf=0; SNPs filtered microhaplotype-wise; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUAC.HapSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.HapSNPs.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUAC.HapSNPs.genind)

# CREATE WILD SAMPLE MATRIX ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUAC.wild <- seq(from=length(which(pop(QUAC.HapSNPs.genind)=="garden"))+1, to=nInd(QUAC.HapSNPs.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUAC.HapSNPs.wildMat <- 
  QUAC.HapSNPs.genind@tab[QUAC.wild,which(colSums(QUAC.HapSNPs.genind@tab[QUAC.wild,], na.rm = TRUE) > 0)]

# CREATE SAMPLING RESULTS ARRAY ----
# Specify number of replicates
num_reps <- 10 
# Export relevant functions and variables
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "QUAC.HapSNPs.wildMat"))

# Run resampling in parallel, to generate an array
samplingResults_QUAC.HapSNPs <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(wildMat=QUAC.HapSNPs.wildMat), simplify = "array")
# Close cores
# stopCluster(cl)
# Examine sampling array
str(samplingResults_QUAC.HapSNPs)

# CALCULATE MEANS AND PLOT ----
# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min95_QUAC.HapSNPs <- min(which(apply(samplingResults_QUAC.HapSNPs[,1,],1,mean) > 95)); min95_QUAC.HapSNPs

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUAC.HapSNPs[,1,], 1, mean)
total_sd <- apply(samplingResults_QUAC.HapSNPs[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUAC.HapSNPs[,2,], 1, mean)
v.com_sd <- apply(samplingResults_QUAC.HapSNPs[,2,], 1, sd)

com_means <- apply(samplingResults_QUAC.HapSNPs[,3,], 1, mean)
com_sd <- apply(samplingResults_QUAC.HapSNPs[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUAC.HapSNPs[,4,], 1, mean)
lowfr_sd <- apply(samplingResults_QUAC.HapSNPs[,4,], 1, sd)

rare_means <- apply(samplingResults_QUAC.HapSNPs[,5,], 1, mean)
rare_sd <- apply(samplingResults_QUAC.HapSNPs[,5,], 1, sd)
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
# ---- FIRST SNP PER LOCUS ----
# READ IN GENIND FILE (QUBO DNFA; R0, min-maf=0; first SNP per locus; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.1SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.1SNP.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# Number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUBO.1SNP.genind)

# CREATE WILD SAMPLE MATRIX ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUBO.wild <- seq(from=length(which(pop(QUBO.1SNP.genind)=="garden"))+1, to=nInd(QUBO.1SNP.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUBO.1SNP.wildMat <- 
  QUBO.1SNP.genind@tab[QUBO.wild,which(colSums(QUBO.1SNP.genind@tab[QUBO.wild,], na.rm = TRUE) > 0)]

# CREATE SAMPLING RESULTS ARRAY ----
# Specify number of replicates
num_reps <- 50 
# Export relevant functions and variables
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "QUBO.1SNP.wildMat"))

# Run resampling in parallel, to generate an array
samplingResults_QUBO.1SNP <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(wildMat=QUBO.1SNP.wildMat), simplify = "array")
# Close cores
# stopCluster(cl)
# Examine sampling array
str(samplingResults_QUBO.1SNP)

# CALCULATE MEANS AND PLOT ----
# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min95_QUBO.1SNP <- min(which(apply(samplingResults_QUBO.1SNP[,1,],1,mean) > 95)); min95_QUBO.1SNP

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUBO.1SNP[,1,], 1, mean)
total_sd <- apply(samplingResults_QUBO.1SNP[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUBO.1SNP[,2,], 1, mean)
v.com_sd <- apply(samplingResults_QUBO.1SNP[,2,], 1, sd)

com_means <- apply(samplingResults_QUBO.1SNP[,3,], 1, mean)
com_sd <- apply(samplingResults_QUBO.1SNP[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUBO.1SNP[,4,], 1, mean)
lowfr_sd <- apply(samplingResults_QUBO.1SNP[,4,], 1, sd)

rare_means <- apply(samplingResults_QUBO.1SNP[,5,], 1, mean)
rare_sd <- apply(samplingResults_QUBO.1SNP[,5,], 1, sd)
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

# ---- MICROHAPLOTYPES ----
# READ IN GENIND FILE (QUBO DNFA; R0, min-maf=0; SNPs filtered microhaplotype-wise; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUBO.HapSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.HapSNPs.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# Number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUBO.HapSNPs.genind)

# CREATE WILD SAMPLE MATRIX ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUBO.wild <- seq(from=length(which(pop(QUBO.HapSNPs.genind)=="garden"))+1, to=nInd(QUBO.HapSNPs.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUBO.HapSNPs.wildMat <- 
  QUBO.HapSNPs.genind@tab[QUBO.wild,which(colSums(QUBO.HapSNPs.genind@tab[QUBO.wild,], na.rm = TRUE) > 0)]

# CREATE SAMPLING RESULTS ARRAY ----
# Specify number of replicates
num_reps <- 50 
# Export relevant functions and variables
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "QUBO.HapSNPs.wildMat"))

# Run resampling in parallel, to generate an array
samplingResults_QUBO.HapSNPs <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(wildMat=QUBO.HapSNPs.wildMat), simplify = "array")
# Close cores
stopCluster(cl)
# Examine sampling array
str(samplingResults_QUBO.HapSNPs)

# CALCULATE MEANS AND PLOT ----
# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min95_QUBO.HapSNPs <- min(which(apply(samplingResults_QUBO.HapSNPs[,1,],1,mean) > 95)); min95_QUBO.HapSNPs

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUBO.HapSNPs[,1,], 1, mean)
total_sd <- apply(samplingResults_QUBO.HapSNPs[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUBO.HapSNPs[,2,], 1, mean)
v.com_sd <- apply(samplingResults_QUBO.HapSNPs[,2,], 1, sd)

com_means <- apply(samplingResults_QUBO.HapSNPs[,3,], 1, mean)
com_sd <- apply(samplingResults_QUBO.HapSNPs[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUBO.HapSNPs[,4,], 1, mean)
lowfr_sd <- apply(samplingResults_QUBO.HapSNPs[,4,], 1, sd)

rare_means <- apply(samplingResults_QUBO.HapSNPs[,5,], 1, mean)
rare_sd <- apply(samplingResults_QUBO.HapSNPs[,5,], 1, sd)
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
