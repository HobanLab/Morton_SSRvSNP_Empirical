# Analyzing output files from stacks-dist-extract
# This is to summarize the results of the final QUAC de novo assembly, QUBO_m7_M5_n5_ga0.01

# %%% Coverage values: generated during ustacks %%%
distTable_coverage <- read.table("metric-depth_of_cov", header=TRUE, skip=1, sep ="\t", check.names = FALSE)
# Mean and sd coverage
mean(distTable_coverage$`depth of cov`); sd(distTable_coverage$`depth of cov`)
# Min and max coverage samples
distTable_coverage$sample[which.min(distTable_coverage$`depth of cov`)];min(distTable_coverage$`depth of cov`)
distTable_coverage$sample[which.max(distTable_coverage$`depth of cov`)];max(distTable_coverage$`depth of cov`)
# Mean % reads incorporated
mean(distTable_coverage$`% reads incorporated`)

# %%% Weighted coverage values: generated during gstacks %%%
distTable_wCoverage <- read.table("metric-weighted_cov", header=TRUE, skip=2, sep ="\t")
# Mean weighted coverage and PCR duplication rate
mean(distTable_wCoverage$mean_cov); sd(distTable_wCoverage$mean_cov)
mean(distTable_wCoverage$pcr_dupl_rate); sd(distTable_wCoverage$pcr_dupl_rate)
# Min and max weighted coverage samples
distTable_wCoverage$sample[which.min(distTable_wCoverage$mean_cov)];min(distTable_wCoverage$mean_cov)
distTable_wCoverage$sample[which.max(distTable_wCoverage$mean_cov)];max(distTable_wCoverage$mean_cov)

# %%% Assembled loci: generated during populations %%%
distTable_assembledLoci <- read.table("metric-assembled_loci", header=TRUE, skip=1, sep ="\t")
# Number of assembled loci: 5,410
unique(distTable_assembledLoci$n_loci)
# Mean frequency of missing loci
mean(distTable_assembledLoci$frequency_missing)
# Min and max frequency missing loci
distTable_assembledLoci$sample[which.min(distTable_assembledLoci$frequency_missing)];min(distTable_assembledLoci$frequency_missing)
distTable_assembledLoci$sample[which.max(distTable_assembledLoci$frequency_missing)];max(distTable_assembledLoci$frequency_missing)

# %%% Polymorphic loci: generated during populations %%%
distTable_polymorphicLoci <- read.table("metric-polymorphic_loci", header=TRUE, skip=1, sep ="\t")

# %%% Number of SNPs: generated during populations %%%
distTable_numberOfSNPs <- read.table("metric-number_of_SNPs", header=TRUE, skip=1, sep ="\t")
# Number of polymorphic sites: 79,291
unique(distTable_numberOfSNPs$n_sites)
# Mean frequency of missing SNPs
mean(distTable_numberOfSNPs$frequency_missing)
# Min and max frequency missing SNPs
distTable_numberOfSNPs$sample[which.min(distTable_numberOfSNPs$frequency_missing)];min(distTable_numberOfSNPs$frequency_missing)
distTable_numberOfSNPs$sample[which.max(distTable_numberOfSNPs$frequency_missing)];max(distTable_numberOfSNPs$frequency_missing)
