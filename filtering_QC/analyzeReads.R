# %%%%%%%%%%%%%%%%%%%%%
# %%% ANALYZE READS %%%
# %%%%%%%%%%%%%%%%%%%%%
# Extract total and per sample read values from file
# This file was generated using stacks-dist-extract process_radtags.log per_barcode_raw_read_counts

# Extract 2nd column: total number of reads per sample
totalReads <- (read.table("metric-retainedReads_allsamples", header=TRUE)[,2])

# Extract 3rd column: retained reads per sample
retainedReads <- (read.table("metric-retainedReads_allsamples", header=TRUE)[,3])

# Calculate each sample's retained read proportion
retainedSampleProportion <- (retainedReads/totalReads)*100
# Data exploration
print(c("Minimum retained read proportion", min(retainedSampleProportion), "Maximum retained read proportion", max(retainedSampleProportion)))

# Calculate the number of retained reads in entire dataset
allReads <- sum(retainedReads)

# Calculate each sample's contribution to the entire dataset
totalSampleProportion <- (retainedReads/allReads)*100
# Data exploration
print(c("Minimum read contribution", min(totalSampleProportion), "Maximum read contribution", max(totalSampleProportion)))

# Combine vectors
read.mat <- cbind(read.table("metric-retainedReads_allsamples", header=TRUE)[,1:3], 
                  retainedSampleProportion, totalSampleProportion)
