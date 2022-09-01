# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% RESAMPLING: MICROHAPLOTYPE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses the ex situ representation functions to conduct resampling
# for Quercus boyntonii (QUBO; GSNAP4 alignment with Quercus robur reference) NextRAD samples.
# (Quercus acerifolia de novo haplotype dataset is too large to be analyzed)

# Two QUBO datasets are analyzed: one with the first SNP written for each locus, and one in which
# SNPs are filtered haplotype-wise (unshared SNPs are pruned, and filters applied for entire SNP locus, 
# rather than by each variant site). Neither dataset uses a filter for missing data (R0). 
# Given the size of haplotype dataests, it's recommend this script run in the background (i.e. not
# in the RStudio GUI), and results are exported to arrays.

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
# CODE COMMENTED OUT: MICROHAPLOTYPE DATASET TOO LARGE TO BE PROCESSED ON HOBAN LAB SERVER
# genpop.filePath <- 
#   "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_HapSNPs_2Pops/"
# setwd(genpop.filePath)
# QUAC.HapSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# # Correct popNames
# pop(QUAC.HapSNPs.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# 
# # CREATE SAMPLING RESULTS ARRAY ----
# # Specify number of replicates
# num_reps <- 25
# # Export relevant functions and variables
# clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
#                               "num_reps", "QUAC.HapSNPs.genind"))
# 
# # Run resampling in parallel, to generate an array
# samplingResults_QUAC.HapSNPs <- 
#   parSapply(cl, 1:num_reps, 
#             function(a) exSitu_Resample(gen.obj=QUAC.HapSNPs.genind), simplify = "array")
# 
# # CALCULATE MEANS AND EXPORT ----
# # Average results across replicates (slices) of the sampling array, to determine
# # the minimum number of samples required to capture 95% wild genetic diversity
# # (We average samplingResults[,1,], since this column contains the total genetic diversity)
# min95_QUAC.HapSNPs <- min(which(apply(samplingResults_QUAC.HapSNPs[,1,],1,mean) > 95)) 
# print(min95_QUAC.HapSNPs)
# 
# # Export resampling array, for later analysis
# saveRDS(samplingResults_QUAC.HapSNPs, 
#         file=paste0(SSRvSNP.wd,"exSituRespresentation/Resampling/QUAC_HapSNPs_resamplingArray.Rdata"))


# %%%% QUBO %%%% ----
# ---- FIRST SNP PER LOCUS ----
# READ IN GENIND FILE (QUBO DNFA; R0, min-maf=0; first SNP per locus; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.1SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.1SNP.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# CREATE SAMPLING RESULTS ARRAY ----
# Specify number of replicates
num_reps <- 50 
# Export relevant functions and variables
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "QUBO.1SNP.genind"))

# Run resampling in parallel, to generate an array
samplingResults_QUBO.1SNP <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(gen.obj=QUBO.1SNP.genind), simplify = "array")

# CALCULATE MEANS AND EXPORT ----
# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min_95_QUBO.1SNP <- min(which(apply(samplingResults_QUBO.1SNP[,1,],1,mean) > 95))
print(min_95_QUBO.1SNP)
# Measuring the standard deviation of the minimum sampling size
min95_SD_QUBO.1SNP <- apply(samplingResults_QUBO.1SNP[,1,],1,sd)[min_95_QUBO.1SNP]
print(min95_SD_QUBO.1SNP)
# Export resampling array, for later analysis
saveRDS(samplingResults_QUBO.1SNP, 
        file=paste0(SSRvSNP.wd,"exSituRepresentation/Resampling/QUBO_1SNP_resamplingArray.Rdata"))

# PLOTTING ----
# Read in resampling array
samplingResults_QUBO.1SNP <- readRDS(file=paste0(SSRvSNP.wd,"exSituRepresentation/Resampling/QUBO_1SNP_resamplingArray.Rdata"))
# Set graphing parameters to plot to charts vertically (for 1 SNP per locus and microhaplotypes)
par(mfcol=c(2,1), oma=rep(0.2,4))
# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUBO.1SNP[,1,], 1, mean); print(total_means)
total_sd <- apply(samplingResults_QUBO.1SNP[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUBO.1SNP[,2,], 1, mean); print(v.com_means)
v.com_sd <- apply(samplingResults_QUBO.1SNP[,2,], 1, sd)

com_means <- apply(samplingResults_QUBO.1SNP[,3,], 1, mean); print(com_means)
com_sd <- apply(samplingResults_QUBO.1SNP[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUBO.1SNP[,4,], 1, mean); print(lowfr_means)
lowfr_sd <- apply(samplingResults_QUBO.1SNP[,4,], 1, sd)

rare_means <- apply(samplingResults_QUBO.1SNP[,5,], 1, mean); print(rare_means)
rare_sd <- apply(samplingResults_QUBO.1SNP[,5,], 1, sd)
# Pick colors, with transparency for values other than Total
plotColors <- c("red","red4","darkorange3","coral","purple")
plotColors[2:5] <- alpha(plotColors[2:5], 0.5)
# Plots all sets of points onto single graph, as well as 95% threshold line
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture",
     main="QUBO (R0, NOMAF, First SNP per locus) Resampling (50 replicates)")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.15)
# Lines for 95% threshold
abline(h=95, col="black", lty=3)
abline(v=min95_QUBO.1SNP, col="black")

# ---- MICROHAPLOTYPES ----
# READ IN GENIND FILE (QUBO DNFA; R0, min-maf=0; SNPs filtered microhaplotype-wise; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUBO.HapSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.HapSNPs.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# CREATE SAMPLING RESULTS ARRAY ----
# Specify number of replicates
num_reps <- 50 
# Export relevant functions and variables
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "QUBO.HapSNPs.genind"))

# Run resampling in parallel, to generate an array
samplingResults_QUBO.HapSNPs <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(gen.obj = QUBO.HapSNPs.genind), simplify = "array")
# Close cores
stopCluster(cl)

# CALCULATE MEANS AND EXPORT ----
# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min_95_QUBO.HapSNPs <- min(which(apply(samplingResults_QUBO.HapSNPs[,1,],1,mean) > 95)) 
print(min_95_QUBO.HapSNPs)
# Measuring the standard deviation of the minimum sampling size
min_95SD_QUBO.HapSNPs <- apply(samplingResults_QUBO.HapSNPs[,1,],1,sd)[min_95_QUBO.HapSNPs]
print(min_95SD_QUBO.HapSNPs)
# Export resampling array, for later analysis
saveRDS(samplingResults_QUBO.HapSNPs, 
        file=paste0(SSRvSNP.wd,"exSituRepresentation/Resampling/QUBO_HapSNPs_resamplingArray.Rdata"))

# PLOTTING ----
# Read in resampling array
samplingResults_QUBO.HapSNPs <- readRDS(file=paste0(SSRvSNP.wd,"exSituRepresentation/Resampling/QUBO_HapSNPs_resamplingArray.Rdata"))
# Calculate means and standard deviations, for each capture rate category
total_means <- apply(samplingResults_QUBO.HapSNPs[,1,], 1, mean); print(total_means)
total_sd <- apply(samplingResults_QUBO.HapSNPs[,1,], 1, sd)

v.com_means <- apply(samplingResults_QUBO.HapSNPs[,2,], 1, mean); print(v.com_means)
v.com_sd <- apply(samplingResults_QUBO.HapSNPs[,2,], 1, sd)

com_means <- apply(samplingResults_QUBO.HapSNPs[,3,], 1, mean); print(com_means)
com_sd <- apply(samplingResults_QUBO.HapSNPs[,3,], 1, sd)

lowfr_means <- apply(samplingResults_QUBO.HapSNPs[,4,], 1, mean); print(lowfr_means)
lowfr_sd <- apply(samplingResults_QUBO.HapSNPs[,4,], 1, sd)

rare_means <- apply(samplingResults_QUBO.HapSNPs[,5,], 1, mean); print(rare_means)
rare_sd <- apply(samplingResults_QUBO.HapSNPs[,5,], 1, sd)
# Pick colors, with transparency for values other than Total
plotColors <- c("red","red4","darkorange3","coral","purple")
plotColors[2:5] <- alpha(plotColors[2:5], 0.5)
# Plots all sets of points onto single graph, as well as 95% threshold line
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture",
     main="QUBO (R0, NOMAF, SNPs Filtered Haplotype-wise) Resampling (50 replicates)")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.15)
# Lines for 95% threshold
abline(h=95, col="black", lty=3)
abline(v=min95_QUBO.HapSNPs, col="black")
