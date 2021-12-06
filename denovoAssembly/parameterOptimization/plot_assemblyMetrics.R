# Script for plotting assembly metrics values against assembly parameters, to optimize de novo assembly

# QUAC assembly metrics
quac.assemblyFilepath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/QUAC/analysis/assemblyMetrics.Rdata"
quac.assemblyMetrics.mat <- readRDS(quac.assemblyFilepath)

# QUAC parameters
quac.assemblyParameters <- read.csv2("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/denovo.params", header=FALSE, 
                                sep = ",", col.names = c("m","M","n","gt-alpha"))

# Combine parameters and assembly metrics into single data.frame
quac.assembly.mat <- cbind(quac.assemblyParameters, quac.assemblyMetrics.mat)

# PLOTTING
subtitle <- "R80: loci shared across 80% of samples"

# ----QUAC COVERAGE----
# Coverage vs. m
boxplot(depth_of_cov~m,data=quac.assembly.mat, xlab="m", ylab="Depth of coverage")
title(main="QUAC: Depth of Coverage across m", line=3)
mtext(side=3, line=1.5, at=2, adj=0, cex=1, subtitle)

# ----QUAC ASSEMBLED LOCI----
# Assembled loci vs. m
boxplot(assembled_loci~m,data=quac.assembly.mat, xlab="m", ylab="Assembled loci")
title(main="QUAC: Assembled loci across m", line=3)
mtext(side=3, line=1.5, at=2, adj=0, cex=1, subtitle)

# Assembled loci vs. M/n
boxplot(assembled_loci~M,data=quac.assembly.mat, xlab="M/n", ylab="Assembled loci")
title(main="QUAC: Assembled loci across M/n", line=3)
mtext(side=3, line=1.5, at=3, adj=0, cex=1, subtitle)

# ----QUAC POLYMORPHIC LOCI----
# Polymorphic loci vs. m
boxplot(polymorphic_loci~m,data=quac.assembly.mat, xlab="m", ylab="Polymorphic loci")
title(main="QUAC: Polymorphic loci across m", line=3)
mtext(side=3, line=1.5, at=2, adj=0, cex=1, subtitle)

# Polymorphic loci vs. M/n
boxplot(polymorphic_loci~M,data=quac.assembly.mat, xlab="M/n", ylab="Polymorphic loci")
title(main="QUAC: Polymorphic loci across M/n", line=3)
mtext(side=3, line=1.5, at=3, adj=0, cex=1, subtitle)

# Polymorphic loci vs. gt-alpha
boxplot(polymorphic_loci~gt.alpha,data=quac.assembly.mat, xlab="gt-alpha", ylab="Polymorphic loci")
title(main="QUAC: Polymorphic loci across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1.1, adj=0, cex=1, subtitle)

# ----QUAC NUMBER OF SNPS----
# Number of SNPs vs. m
boxplot(number_of_snps~m,data=quac.assembly.mat, xlab="m", ylab="Number of SNPs")
title(main="QUAC: Number of SNPs (total) across m", line=3)
mtext(side=3, line=1.5, at=2, adj=0, cex=1, subtitle)

# Number of SNPs vs. M/n
boxplot(number_of_snps~M,data=quac.assembly.mat, xlab="M/n", ylab="Number of SNPs")
title(main="QUAC: Number of SNPs (total) across M/n", line=3)
mtext(side=3, line=1.5, at=3, adj=0, cex=1, subtitle)

# Number of SNPs vs. gt-alpha
boxplot(number_of_snps~gt.alpha,data=quac.assembly.mat, xlab="gt-alpha", ylab="Number of SNPs")
title(main="QUAC: Number of SNPs (total) across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1.1, adj=0, cex=1, subtitle)

# DEPRECATED (UNUSEFUL) PLOTS----
# Coverage vs. gt-alpha
boxplot(depth_of_cov~gt.alpha,data=quac.assembly.mat, xlab="gt-alpha", ylab="Depth of coverage")
title(main="QUAC: Depth of Coverage across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1.1, adj=0, cex=1, subtitle)

# Assembled loci vs. gt-alpha
boxplot(assembled_loci~gt.alpha,data=quac.assembly.mat, xlab="gt-alpha", ylab="Assembled loci")
title(main="QUAC: Assembled loci across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1.1, adj=0, cex=1, subtitle)