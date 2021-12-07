# Extracting de novo assembly metrics from text files generated using stacks-dist-extract
# These are single values (atomics) that are appended to rows of a matrix
# Each row of the matrix represents an assembly parameter set; columns are assembly metrics

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% BUILD ASSEMBLY MATRIX OBJECT %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assemblyFilepath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/QUAC/analysis/assemblyMetrics.Rdata"
# This code is commented out because it only needs to be run once
# Create matrix of appropriate length
# assemblyMetrics.mat <- matrix(ncol = 5)
# colnames(assemblyMetrics.mat) <- c("depth_of_cov","assembled_loci","polymorphic_loci","number_of_snps","pcr_duplication_rate")
# # Save empty matrix to file, to be read in later on
# saveRDS(assemblyMetrics.mat, file = assemblyFilepath)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% EXTRACT/CALCULATE ASSEMBLY METRICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Extract coverage values for each sample and calculate mean
depth.of.cov <- mean(read.table("metric-depth_of_cov", header=FALSE, skip=2)[,1])
# Read in total assembled loci
assembled.loci <- as.numeric(unique(read.table("metric-assembled_loci", header=TRUE, skip=3)))
# Read in total polymorphic loci (take sum of all loci that aren't monomorphic)
polymorphic.loci <- sum(read.table("metric-polymorphic_loci", header=FALSE)[,1])
# Read in total number of SNPs
number.of.snps <- as.numeric(unique(read.table("metric-number_of_SNPs", header=TRUE, skip=3)))
# Read in PCR duplication rate
pcr.duplication.rate <- mean((read.table("metric-pcr_duplication_rate", header=TRUE, skip=2))[,1])
# Combine values into a vector
values <- c(depth.of.cov, assembled.loci, polymorphic.loci, number.of.snps, pcr.duplication.rate)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% ADD ASSEMBLY METRICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Read in existing matrix
assemblyMetrics.mat <- readRDS(assemblyFilepath)

# If the first matrix value is NA, write values to first row. Otherwise, append the values
if(is.na(assemblyMetrics.mat[1,1])){
  assemblyMetrics.mat[1,] <- values
} else{
  assemblyMetrics.mat <- rbind(assemblyMetrics.mat, values)
  # If the matrix has been filled (80 rows), give the rows names according to the assembly parameters
  if(nrow(assemblyMetrics.mat)==80){
    rownames(assemblyMetrics.mat) <- read.csv("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/QUAC/output/QUAC_assemblies", header=FALSE)[,1]
  }
}

# Save the assembly metrics matrix object
saveRDS(assemblyMetrics.mat, file = assemblyFilepath)
