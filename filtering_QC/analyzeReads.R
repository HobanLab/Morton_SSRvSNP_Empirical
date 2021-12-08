# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% ANALYZE READ COUNTS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Extract total and per sample read values from file
# This file was generated using stacks-dist-extract process_radtags.log per_barcode_raw_read_counts

# Extract 2nd column: total number of reads per sample
totalReads <- (read.table("metric-retainedReads_allsamples", header=TRUE)[,2])

# Extract 3rd column: retained reads per sample
retainedReads <- (read.table("metric-retainedReads_allsamples", header=TRUE)[,3])
# Calculate each sample's retained read proportion
retainedSampleProportion <- (retainedReads/totalReads)*100

# Calculate the number of retained reads in entire dataset
allReads <- sum(retainedReads)
# Calculate each sample's contribution to the entire dataset
totalSampleProportion <- (retainedReads/allReads)*100

# Combine vectors
read.mat <- cbind(read.table("metric-retainedReads_allsamples", header=TRUE)[,1:3], 
                  retainedSampleProportion, totalSampleProportion)

# Data exploration
# Retained reads
mean(retainedReads); sd(retainedReads)
read.mat[which.min(read.mat[,3]),]; read.mat[which.max(read.mat[,3]),]
# Retained read proportion: mean, standard deviation, min/max
mean(retainedSampleProportion); sd(retainedSampleProportion)
read.mat[which.min(read.mat[,4]),]; read.mat[which.max(read.mat[,4]),]
# Total sample proportion: mean, standard deviation, min/max
mean(totalSampleProportion); sd(totalSampleProportion)
read.mat[which.min(read.mat[,5]),]; read.mat[which.max(read.mat[,5]),]
