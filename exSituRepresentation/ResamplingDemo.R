# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% RESAMPLING DEMONSTRATION %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in the 3 functions declared in the functions_exSituRepresentation.R script,
# and demonstrates the outputs of each function using both a microsatellite genind object
# and a SNP genind object (both QUAC)

library(adegenet)
library(RColorBrewer)
library(scales)
library(parallel)

# %%%% FUNCTIONS %%%% ----
# Read in relevant functions required for resampling analyses
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)
source("exSituRepresentation/functions_exSituRepresentation.R")

# %%%% MICROSATELLITES %%%% ----
# ---- PROCESS GENIND ----
# Read in genind file (GCC_QUAC_ZAIN repo; QUAC_wK_garden_wild_clean.gen)
QUAC.MSAT.genpop.filePath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/Garden_Wild/"
setwd(QUAC.MSAT.genpop.filePath)
QUAC.MSAT.genind <- read.genepop("QUAC_wK_garden_wild_clean.gen", quiet = TRUE, ncode = 3)
pop(QUAC.MSAT.genind)
# Correct popNames: pop1 is Garden, pop2 is Wild
pop(QUAC.MSAT.genind) <- gsub("pop1", "garden", pop(QUAC.MSAT.genind))
pop(QUAC.MSAT.genind) <- gsub("pop2", "wild", pop(QUAC.MSAT.genind))
# Create a matrix of strictly wild samples
QUAC.MSAT.wildMat <- QUAC.MSAT.genind@tab[which(pop(QUAC.MSAT.genind) == "wild"),]
nrow(QUAC.MSAT.wildMat)
# Create a vector of wild allele frequencies
QUAC.MSAT.freqVector <- colSums(QUAC.MSAT.wildMat, na.rm = TRUE)/(nrow(QUAC.MSAT.wildMat)*2)*100
hist(QUAC.MSAT.freqVector)

# ---- DEMONSTRATE FUNCTIONS ----
# getAlleleCategories: examines the representation of a matrix of individuals, based on a vector of allele frequencies
# The allele frequency vector represents the frequencies over the ENTIRE wild population
# All samples
getAlleleCategories(QUAC.MSAT.freqVector, QUAC.MSAT.wildMat)
# 75% of samples
getAlleleCategories(QUAC.MSAT.freqVector, QUAC.MSAT.wildMat[1:129,])
# 50% of samples
getAlleleCategories(QUAC.MSAT.freqVector, QUAC.MSAT.wildMat[1:86,])
# 25% of samples
getAlleleCategories(QUAC.MSAT.freqVector, QUAC.MSAT.wildMat[1:43,])

# exSitu_Sample: wrapper for getAlleleCategories. When given the ENTIRE sample matrix, first calculates
# the allele frequency vector, then samples the matrix based on the num_samples argument
exSitu_Sample(QUAC.MSAT.wildMat, numSamples = 43)
exSitu_Sample(QUAC.MSAT.wildMat, numSamples = 129)

# exSitu_Resample: wrapper for exSitu_Sample. Iterates exSitu_Sample for every numSamples value from 
# 2 to nrow(sampleMatrix), using sapply
exSitu_Resample(QUAC.MSAT.wildMat)

# ---- PARALLELIZATION ----
# exSitu_Resample generates a matrix. To build an array, use sapply (or parSapply, to do so in parallel)

# WITHOUT PARALLELIZATION
num_reps <- 10
QUAC.MSAT.samplingResults <- 
  sapply(1:num_reps, function(a) exSitu_Resample(wildMat=QUAC.MSAT.wildMat), simplify = "array")

str(QUAC.MSAT.samplingResults)
QUAC.MSAT.samplingResults[,,1]

# WITH PARALLELIZATION
# Set up cores, and declare replicates variable
num_cores <- detectCores() - 1 ; cl <- makeCluster(num_cores)
num_reps <- 5000
# Export relevant functions and variables for parallelization
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "QUAC.MSAT.wildMat"))
# Build array using parSapply
QUAC.MSAT.samplingResults <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(wildMat=QUAC.MSAT.wildMat), simplify = "array")

str(QUAC.MSAT.samplingResults)
QUAC.MSAT.samplingResults[,,1]
# Average minimum sampling size to capture 95% of total allelic diversity
min_95_QUAC.MSAT <- min(which(apply(QUAC.MSAT.samplingResults[,1,],1,mean) > 95)); min_95_QUAC.MSAT
# Close cores
stopCluster(cl)

# %%%% SNPS %%%% ----
# ---- PROCESS GENIND ----
# Read in genind file (QUAC DNFA; R0, min-maf=0; 1 SNP/locus; 2 populations, garden and wild)
QUAC.SNP.genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(QUAC.SNP.genpop.filePath)
QUAC.SNP.genind <- read.genepop(paste0(QUAC.SNP.genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Create a matrix of strictly wild samples
QUAC.SNP.wildMat <- QUAC.SNP.genind@tab[which(pop(QUAC.SNP.genind) == "wild"),]
nrow(QUAC.SNP.wildMat)
# Create a vector of wild allele frequencies
QUAC.SNP.freqVector <- colSums(QUAC.SNP.wildMat, na.rm = TRUE)/(nrow(QUAC.SNP.wildMat)*2)*100
hist(QUAC.SNP.freqVector)

# ---- DEMONSTRATE FUNCTIONS ----
# getAlleleCategories: examines the representation of a matrix of individuals, based on a vector of allele frequencies
# The allele frequency vector represents the frequencies over the ENTIRE wild population
# All samples
getAlleleCategories(QUAC.SNP.freqVector, QUAC.SNP.wildMat)
# 75% of samples
getAlleleCategories(QUAC.SNP.freqVector, QUAC.SNP.wildMat[1:72,])
# 50% of samples
getAlleleCategories(QUAC.SNP.freqVector, QUAC.SNP.wildMat[1:48,])
# 25% of samples
getAlleleCategories(QUAC.SNP.freqVector, QUAC.SNP.wildMat[1:24,])

# exSitu_Sample: wrapper for getAlleleCategories. When given the ENTIRE sample matrix, first calculates
# the allele frequency vector, then samples the matrix based on the num_samples argument
exSitu_Sample(QUAC.SNP.wildMat, numSamples = 24)
exSitu_Sample(QUAC.SNP.wildMat, numSamples = 72)

# exSitu_Resample: wrapper for exSitu_Sample. Iterates exSitu_Sample for every numSamples value from 
# 2 to nrow(sampleMatrix), using sapply
exSitu_Resample(QUAC.SNP.wildMat)

# ---- PARALLELIZATION ----
# exSitu_Resample generates a matrix. To build an array, use sapply (or parSapply, to do so in parallel)

# WITHOUT PARALLELIZATION
num_reps <- 10
QUAC.SNP.samplingResults <- 
  sapply(1:num_reps, function(a) exSitu_Resample(wildMat=QUAC.SNP.wildMat), simplify = "array")

str(QUAC.SNP.samplingResults)
QUAC.SNP.samplingResults[,,1]

# WITH PARALLELIZATION
# Set up cores, and declare replicates variable
num_cores <- detectCores() - 1 ; cl <- makeCluster(num_cores)
num_reps <- 1000
# Export relevant functions and variables for parallelization
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "QUAC.SNP.wildMat"))
# Build array using parSapply
QUAC.SNP.samplingResults <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(wildMat=QUAC.SNP.wildMat), simplify = "array")

str(QUAC.SNP.samplingResults)
QUAC.SNP.samplingResults[,,1]
# Average minimum sampling size to capture 95% of total allelic diversity
min_95_QUAC.SNP <- min(which(apply(QUAC.SNP.samplingResults[,1,],1,mean) > 95)); min_95_QUAC.SNP
# Close cores
stopCluster(cl)
