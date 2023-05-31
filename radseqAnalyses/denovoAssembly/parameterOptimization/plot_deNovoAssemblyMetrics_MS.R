# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% PLOTTING DE NOVO ASSEMBLY METRICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script plots the de novo assembly metrics for both Quercus acerifolia (QUAC) and Q. boyntonii (QUBO).
# It is functionally identical to the plot_QUAC_assemblyMetrics.R and plot_QUBO_assemblyMetrics.R scripts
# (in the respective QUAV/QUBO folders), but generates PDF images for the manuscript.

# De novo assemblies were generated using the Stacks denovo_map.pl script, and the parameter values specified 
# in the denovo.params file (which are the m, M, n, and gt-alpha parameter values). For each of the 80 parameter
# combinations, a denovo assembly of a subset of samples from each species was utilized, and the following metrics
# were calculated with respect to each de novo assembly parameter:
#   1. Depth of coverage
#   2. Number of assembled loci 
#   3. Number of polymorphic loci 
#   4. Number of SNPs
# These 4 metrics were generated using the stacks-dist-extract function (see the QU(AC/BO)_extractMetrics.sh scripts).
# R80 loci (loci shared across at least 80% of the samples) were used for these analyses.

# We sought the lowest possible parameter values which maximized these 4 metric values. We also considered the number
# of duplicate loci (difference in loci seen between duplicate samples), but that analysis isn't included in this script.

# Specify path to the directory (on the lab server), where plots (PDFs) will be saved
imageOutDir <- "/home/akoontz/Documents/SSRvSNP/Documentation/Images/MolEcol_202305_Images/"

# %%%% QUAC %%%% ----
# %%%% READ IN PARAMETERS AND ASSEMBLY METRICS ----
# De novo assembly parameter values
QUAC.assemblyParameters <- read.csv2("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/denovo.params", header=FALSE, 
                                     sep = ",", col.names = c("m","M","n","gt-alpha"))
# Assembly metrics
QUAC.assemblyFilepath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/QUAC/analysis/assemblyMetrics.Rdata"
QUAC.assemblyMetrics.mat <- readRDS(QUAC.assemblyFilepath)
# Combine parameters and assembly metrics into single data.frame
QUAC.assembly.mat <- cbind(QUAC.assemblyParameters, QUAC.assemblyMetrics.mat)

# %%%% PLOTTING ----
# Parameter: m ----
# Call pdf command, to save resampling plots to disk
pdf(file = paste0(imageOutDir, "supplement/Fig_S2.pdf"), width = 9, height = 7.5)
# Set plotting parameters: two rows, two columns (4 graphs per window)
par(mfcol=c(2,2), mar = c(3,2,4,2), oma = c(1, 3, 1, .5), mgp = c(2, 0.6, 0))
# Coverage vs. m
boxplot(depth_of_cov~m,data=QUAC.assembly.mat, xlab="m", ylab="Depth of coverage")
title(main="Depth of Coverage across m", line=0.5, cex.main=0.8)
# Assembled loci vs. m
boxplot(assembled_loci~m,data=QUAC.assembly.mat, xlab="m", ylab="Assembled loci")
title(main="Assembled loci (R80) across m", line=0.5, cex.main=0.8)
# Polymorphic loci vs. m
boxplot(polymorphic_loci~m,data=QUAC.assembly.mat, xlab="m", ylab="Polymorphic loci")
title(main="Polymorphic loci (R80) across m", line=0.5, cex.main=0.8)
# Number of SNPs vs. m
boxplot(number_of_snps~m,data=QUAC.assembly.mat, xlab="m", ylab="Number of SNPs")
title(main="Number of SNPs (total, R80 loci) across m", line=0.5, cex.main=0.8)
# Title
mtext(side=3, line=24.5, at=-3.5, adj=0, cex=1.2, "Q. acerifolia: assembly metrics across m values")
# Turn off plotting device, to save plot
dev.off()

# Parameter: M/n ----
# Call pdf command, to save resampling plots to disk
pdf(file = paste0(imageOutDir, "supplement/Fig_S3.pdf"), width = 9, height = 7.5)
# Set plotting parameters: two rows, two columns (4 graphs per window)
par(mfcol=c(2,2), mar = c(3,2,4,2), oma = c(1, 3, 1, .5), mgp = c(2, 0.6, 0))
# Coverage vs. M/n
boxplot(depth_of_cov~M,data=QUAC.assembly.mat, xlab="M/n", ylab="Depth of coverage")
title(main="Depth of Coverage across M/n", line=0.5, cex.main=0.8)
# Assembled loci vs. M/n
boxplot(assembled_loci~M,data=QUAC.assembly.mat, xlab="M/n", ylab="Assembled loci")
title(main="Assembled loci across M/n", line=0.5, cex.main=0.8)
# Polymorphic loci vs. M/n
boxplot(polymorphic_loci~M,data=QUAC.assembly.mat, xlab="M/n", ylab="Polymorphic loci")
title(main="Polymorphic loci across M/n", line=0.5, cex.main=0.8)
# Number of SNPs vs. M/n
boxplot(number_of_snps~M,data=QUAC.assembly.mat, xlab="M/n", ylab="Number of SNPs")
title(main="Number of SNPs (total) across M/n", line=0.5, cex.main=0.8)
# Title
mtext(side=3, line=24.5, at=-6, adj=0, cex=1.2, "Q. acerifolia: assembly metrics across M/n values")
# Turn off plotting device, to save plot
dev.off()

# Parameter: gt-alpha ----
# Call pdf command, to save resampling plots to disk
pdf(file = paste0(imageOutDir, "supplement/Fig_S4.pdf"), width = 9, height = 7.5)
# Set plotting parameters: two rows, two columns (4 graphs per window)
par(mfcol=c(2,2), mar = c(3,2,4,2), oma = c(1, 3, 1, .5), mgp = c(2, 0.6, 0))
# Coverage vs. gt-alpha
boxplot(depth_of_cov~gt.alpha,data=QUAC.assembly.mat, xlab="gt-alpha", ylab="Depth of coverage")
title(main="Depth of Coverage across gt-alpha", line=0.5, cex.main=0.8)
# Assembled loci vs. gt-alpha
boxplot(assembled_loci~gt.alpha,data=QUAC.assembly.mat, xlab="gt-alpha", ylab="Assembled loci")
title(main="Assembled loci across gt-alpha", line=0.5, cex.main=0.8)
# Polymorphic loci vs. gt-alpha
boxplot(polymorphic_loci~gt.alpha,data=QUAC.assembly.mat, xlab="gt-alpha", ylab="Polymorphic loci")
title(main="Polymorphic loci across gt-alpha", line=0.5, cex.main=0.8)
# Number of SNPs vs. gt-alpha
boxplot(number_of_snps~gt.alpha,data=QUAC.assembly.mat, xlab="gt-alpha", ylab="Number of SNPs")
title(main="Number of SNPs (total) across gt-alpha", line=0.5, cex.main=0.8)
# Title
mtext(side=3, line=24.5, at=-1.3, adj=0, cex=1.2, "Q. acerifolia: assembly metrics across gt-alpha values")
# Turn off plotting device, to save plot
dev.off()

# %%%% QUBO %%%% ----
# %%%% READ IN PARAMETERS AND ASSEMBLY METRICS ----
# De novo assembly parameter values
QUBO.assemblyParameters <- read.csv2("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/denovo.params", header=FALSE, 
                                     sep = ",", col.names = c("m","M","n","gt-alpha"))
# Assembly metrics
QUBO.assemblyFilepath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/QUBO/analysis/QUBO_assemblyMetrics.Rdata"
QUBO.assemblyMetrics.mat <- readRDS(QUBO.assemblyFilepath)
# Combine parameters and assembly metrics into single data.frame
QUBO.assembly.mat <- cbind(QUBO.assemblyParameters, QUBO.assemblyMetrics.mat)

# %%%% PLOTTING ----
# Parameter: m ----
# Call pdf command, to save resampling plots to disk
pdf(file = paste0(imageOutDir, "supplement/Fig_S5.pdf"), width = 9, height = 7.5)
# Set plotting parameters: two rows, two columns (4 graphs per window)
par(mfcol=c(2,2), mar = c(3,2,4,2), oma = c(1, 3, 1, .5), mgp = c(2, 0.6, 0))
# Coverage vs. m
boxplot(depth_of_cov~m,data=QUBO.assembly.mat, xlab="m", ylab="Depth of coverage")
title(main="Depth of Coverage across m", line=0.5, cex.main=0.8)
# Assembled loci vs. m
boxplot(assembled_loci~m,data=QUBO.assembly.mat, xlab="m", ylab="Assembled loci")
title(main="Assembled loci (R80) across m", line=0.5, cex.main=0.8)
# Polymorphic loci vs. m
boxplot(polymorphic_loci~m,data=QUBO.assembly.mat, xlab="m", ylab="Polymorphic loci")
title(main="Polymorphic loci (R80) across m", line=0.5, cex.main=0.8)
# Number of SNPs vs. m
boxplot(number_of_snps~m,data=QUBO.assembly.mat, xlab="m", ylab="Number of SNPs")
title(main="Number of SNPs (total, R80 loci) across m", line=0.5, cex.main=0.8)
# Title
mtext(side=3, line=24.5, at=-3.3, adj=0, cex=1.2, "Q. boyntonii: assembly metrics across m values")
# Turn off plotting device, to save plot
dev.off()

# Parameter: M/n ----
# Call pdf command, to save resampling plots to disk
pdf(file = paste0(imageOutDir, "supplement/Fig_S6.pdf"), width = 9, height = 7.5)
# Set plotting parameters: two rows, two columns (4 graphs per window)
par(mfcol=c(2,2), mar = c(3,2,4,2), oma = c(1, 3, 1, .5), mgp = c(2, 0.6, 0))
# Coverage vs. M/n
boxplot(depth_of_cov~M,data=QUBO.assembly.mat, xlab="M/n", ylab="Depth of coverage")
title(main="Depth of Coverage across M/n", line=0.5, cex.main=0.8)
# Assembled loci vs. M/n
boxplot(assembled_loci~M,data=QUBO.assembly.mat, xlab="M/n", ylab="Assembled loci")
title(main="Assembled loci across M/n", line=0.5, cex.main=0.8)
# Polymorphic loci vs. M/n
boxplot(polymorphic_loci~M,data=QUBO.assembly.mat, xlab="M/n", ylab="Polymorphic loci")
title(main="Polymorphic loci across M/n", line=0.5, cex.main=0.8)
# Number of SNPs vs. M/n
boxplot(number_of_snps~M,data=QUBO.assembly.mat, xlab="M/n", ylab="Number of SNPs")
title(main="Number of SNPs (total) across M/n", line=0.5, cex.main=0.8)
# Title
mtext(side=3, line=24.5, at=-6, adj=0, cex=1.2, "Q. boyntonii: assembly metrics across M/n values")
# Turn off plotting device, to save plot
dev.off()

# Parameter: gt-alpha ----
# Call pdf command, to save resampling plots to disk
pdf(file = paste0(imageOutDir, "supplement/Fig_S7.pdf"), width = 9, height = 7.5)
# Set plotting parameters: two rows, two columns (4 graphs per window)
par(mfcol=c(2,2), mar = c(3,2,4,2), oma = c(1, 3, 1, .5), mgp = c(2, 0.6, 0))
# Coverage vs. gt-alpha
boxplot(depth_of_cov~gt.alpha,data=QUBO.assembly.mat, xlab="gt-alpha", ylab="Depth of coverage")
title(main="Depth of Coverage across gt-alpha", line=0.5, cex.main=0.8)
# Assembled loci vs. gt-alpha
boxplot(assembled_loci~gt.alpha,data=QUBO.assembly.mat, xlab="gt-alpha", ylab="Assembled loci")
title(main="Assembled loci across gt-alpha", line=0.5, cex.main=0.8)
# Polymorphic loci vs. gt-alpha
boxplot(polymorphic_loci~gt.alpha,data=QUBO.assembly.mat, xlab="gt-alpha", ylab="Polymorphic loci")
title(main="Polymorphic loci across gt-alpha", line=0.5, cex.main=0.8)
# Number of SNPs vs. gt-alpha
boxplot(number_of_snps~gt.alpha,data=QUBO.assembly.mat, xlab="gt-alpha", ylab="Number of SNPs")
title(main="Number of SNPs (total) across gt-alpha", line=0.5, cex.main=0.8)
# Title
mtext(side=3, line=24.5, at=-1.3, adj=0, cex=1.2, "Q. boyntonii: assembly metrics across gt-alpha values")
# Turn off plotting device, to save plot
dev.off()
