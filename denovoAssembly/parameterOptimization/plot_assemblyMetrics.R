# Script for plotting assembly metrics values against assembly parameters, to optimize de novo assembly

# QUAC assembly metrics
assemblyFilepath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/QUAC/analysis/assemblyMetrics.Rdata"
assemblyMetrics.mat <- readRDS(assemblyFilepath)

# QUAC parameters
assemblyParameters <- read.csv2("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/denovo.params", header=FALSE, 
                                sep = ",", col.names = c("m","M","n","gt-alpha"))

# Combine parameters and assembly metrics into single data.frame
assembly.mat <- cbind(assemblyParameters, assemblyMetrics.mat)

# PLOTTING
subtitle <- "R80: loci shared across 80% of samples"

# Coverage vs. m
boxplot(depth_of_cov~m,data=assembly.mat, xlab="m", ylab="Depth of coverage")
title(main="QUAC: Depth of Coverage across m", line=3)
mtext(side=3, line=1.5, at=2, adj=0, cex=1, subtitle)

# Coverage vs. gt-alpha
boxplot(depth_of_cov~gt.alpha,data=assembly.mat, xlab="gt-alpha", ylab="Depth of coverage")
title(main="QUAC: Depth of Coverage across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1.1, adj=0, cex=1, subtitle)

# Assembled loci vs. M/n
boxplot(assembled_loci~M,data=assembly.mat, xlab="M/n", ylab="Assembled loci")
title(main="QUAC: Assembled loci across M/n", line=3)
mtext(side=3, line=1.5, at=3, adj=0, cex=1, subtitle)

# Assembled loci vs. gt-alpha
boxplot(assembled_loci~gt.alpha,data=assembly.mat, xlab="gt-alpha", ylab="Assembled loci")
title(main="QUAC: Assembled loci across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1.1, adj=0, cex=1, subtitle)

# Polymorphic loci vs. M/n
boxplot(polymorphic_loci~M,data=assembly.mat, xlab="M/n", ylab="Polymorphic loci")
title(main="QUAC: Polymorphic loci across M/n", line=3)
mtext(side=3, line=1.5, at=3, adj=0, cex=1, subtitle)

# Polymorphic loci vs. gt-alpha
boxplot(polymorphic_loci~gt.alpha,data=assembly.mat, xlab="gt-alpha", ylab="Polymorphic loci")
title(main="QUAC: Polymorphic loci across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1.1, adj=0, cex=1, subtitle)

# Number of SNPs vs. M/n
boxplot(number_of_snps~M,data=assembly.mat, xlab="M/n", ylab="Number of SNPs")
title(main="QUAC: Number of SNPs (total) across M/n", line=3)
mtext(side=3, line=1.5, at=3, adj=0, cex=1, subtitle)

# Number of SNPs vs. gt-alpha
boxplot(number_of_snps~gt.alpha,data=assembly.mat, xlab="gt-alpha", ylab="Number of SNPs")
title(main="QUAC: Number of SNPs (total) across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1.1, adj=0, cex=1, subtitle)
