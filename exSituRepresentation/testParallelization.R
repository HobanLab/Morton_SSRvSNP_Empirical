# Test parallelization script

library(adegenet)
library(parallel)

# Generate a list of matrices of sampling results using parallelized sapply (parSapply)

# Source required functions
setwd("~/Documents/SSRvSNP/Code/exSituRepresentation/")
source("functions_exSituRepresentation.R")
# Set up relevant cores
num_cores <- detectCores() - 1 ; cl <- makeCluster(num_cores)

# Read in QUAC genind object, and create sample matrix
# Read in genind file (QUAC DNFA; R0, min-maf=0; first SNP/locus; 2 populations, garden and wild) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.R0_NOMAF.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Make a vector corresponding to wild individuals in the genind file (to extract them)
# (Wild individuals are those that don't have a population named "garden")
QUAC.wild <- seq(from=length(which(pop(QUAC.R0_NOMAF.genind)=="garden"))+1, to=nInd(QUAC.R0_NOMAF.genind))
# Create a matrix of ONLY wild individuals with present alleles, from original genind object
QUAC.SNP.wildMat <- 
  QUAC.R0_NOMAF.genind@tab[QUAC.wild,which(colSums(QUAC.R0_NOMAF.genind@tab[QUAC.wild,], na.rm = TRUE) > 0)]
print("Genind object read into environment!")

# Read in the adegenet library on the clusters (may not be necessary)
clusterEvalQ(cl, library(adegenet))
# Set number of replicates
num_reps <- 1000 ; print(paste0("Number of replicates: ", num_reps))
print(paste0("Number of cores: ", num_cores))
# Export relevant functions and variables
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "QUAC.SNP.wildMat"))

# Print start time
print(paste0("Resamping start time: ",Sys.time()))

# Run resampling in parallel, to generate an array
samplingResults_QUAC.SNP <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(wildMat=QUAC.SNP.wildMat), simplify = "array")
# Close cores
stopCluster(cl)

# Print stop time
print(paste0("Resamping stop time: ",Sys.time()))

# Calculate and report the minimum number of samples
min_95_QUAC.SNP <- min(which(apply(samplingResults_QUAC.SNP[,1,],1,mean) > 95)) 
print(paste0("Minimum sampling size: ",min_95_QUAC.SNP))

# Export the resampling results array to R data file, in the main code repository
setwd("~/Documents/SSRvSNP/Code/exSituRepresentation/")
saveRDS(samplingResults_QUAC.SNP, file = "QUAC.SNP.samplingArray.Rdata")
