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
# READ IN GENIND FILE (QUAC DNFA; R0, min-maf=0; first SNP per locus; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.R0_NOMAF.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUAC.R0_NOMAF.genind)

# CREATE WILD SAMPLE MATRIX ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUAC.wild <- seq(from=length(which(pop(QUAC.R0_NOMAF.genind)=="garden"))+1, to=nInd(QUAC.R0_NOMAF.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUAC.wildMat <- 
  QUAC.R0_NOMAF.genind@tab[QUAC.wild,which(colSums(QUAC.R0_NOMAF.genind@tab[QUAC.wild,], na.rm = TRUE) > 0)]

# CREATE SAMPLING RESULTS ARRAY ----
# Specify number of replicates
num_reps <- 50 
# Export relevant functions and variables
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "QUAC.wildMat"))

# Run resampling in parallel, to generate an array
samplingResults_QUAC <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(wildMat=QUAC.wildMat), simplify = "array")
# Close cores
# stopCluster(cl)
# Examine sampling array
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
# READ IN GENIND FILE (QUBO DNFA; R0, min-maf=0; first SNP per locus; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R0_NOMAF.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# Number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUBO.R0_NOMAF.genind)

# CREATE WILD SAMPLE MATRIX ----
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUBO.wild <- seq(from=length(which(pop(QUBO.R0_NOMAF.genind)=="garden"))+1, to=nInd(QUBO.R0_NOMAF.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUBO.wildMat <- 
  QUBO.R0_NOMAF.genind@tab[QUBO.wild,which(colSums(QUBO.R0_NOMAF.genind@tab[QUBO.wild,], na.rm = TRUE) > 0)]

# CREATE SAMPLING RESULTS ARRAY ----
# Specify number of replicates
num_reps <- 50 
# Export relevant functions and variables
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "QUBO.wildMat"))

# Run resampling in parallel, to generate an array
samplingResults_QUBO <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(wildMat=QUBO.wildMat), simplify = "array")
# Close cores
stopCluster(cl)
# Examine sampling array
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
