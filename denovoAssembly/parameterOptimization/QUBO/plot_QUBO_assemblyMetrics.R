# Script for plotting assembly metrics values against assembly parameters, to optimize de novo assembly

# Commands below are used to construct a matrix of assembly parameters and output metrics
# QUBO assembly metrics
qubo.assemblyFilepath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/QUBO/analysis/QUBO_assemblyMetrics.Rdata"
qubo.assemblyMetrics.mat <- readRDS(qubo.assemblyFilepath)

# QUBO parameters
qubo.assemblyParameters <- read.csv2("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_paramOpt/denovo.params", header=FALSE, 
                                sep = ",", col.names = c("m","M","n","gt-alpha"))

# Combine parameters and assembly metrics into single data.frame
qubo.assembly.mat <- cbind(qubo.assemblyParameters, qubo.assemblyMetrics.mat)

# PLOTTING
subtitle <- "R80: loci shared across 80% of samples"

# ----QUBO m----
# Coverage vs. m
boxplot(depth_of_cov~m,data=qubo.assembly.mat, xlab="m", ylab="Depth of coverage")
title(main="Depth of Coverage across m", line=3)
mtext(side=3, line=1.5, at=1.2, adj=0, cex=0.75, subtitle)

# Weighted coverage vs. m
boxplot(weighted_cov~m,data=qubo.assembly.mat, xlab="m", ylab="Weighted coverage")
title(main="Weighted coverage across m", line=3)
mtext(side=3, line=1.5, at=1.2, adj=0, cex=0.75, subtitle)

# Assembled loci vs. m
boxplot(assembled_loci~m,data=qubo.assembly.mat, xlab="m", ylab="Assembled loci")
title(main="Assembled loci across m", line=3)
mtext(side=3, line=1.5, at=1.2, adj=0, cex=0.75, subtitle)

# Polymorphic loci vs. m
boxplot(polymorphic_loci~m,data=qubo.assembly.mat, xlab="m", ylab="Polymorphic loci")
title(main="Polymorphic loci across m", line=3)
mtext(side=3, line=1.5, at=1.2, adj=0, cex=0.75, subtitle)

# Number of SNPs vs. m
boxplot(number_of_snps~m,data=qubo.assembly.mat, xlab="m", ylab="Number of SNPs")
title(main="Number of SNPs (total) across m", line=3)
mtext(side=3, line=1.5, at=1.2, adj=0, cex=0.75, subtitle)

# ----QUBO M/n----
# Coverage vs. M/n
boxplot(depth_of_cov~M,data=qubo.assembly.mat, xlab="M/n", ylab="Depth of coverage")
title(main="Depth of Coverage across M/n", line=3)
mtext(side=3, line=1.5, at=1.2, adj=0, cex=0.75, subtitle)

# Weighted coverage vs. M/n
boxplot(weighted_cov~M,data=qubo.assembly.mat, xlab="M/n", ylab="Weighted coverage")
title(main="Weighted coverage across M/n", line=3)
mtext(side=3, line=1.5, at=1.2, adj=0, cex=0.75, subtitle)

# Assembled loci vs. M/n
boxplot(assembled_loci~M,data=qubo.assembly.mat, xlab="M/n", ylab="Assembled loci")
title(main="Assembled loci across M/n", line=3)
mtext(side=3, line=1.5, at=1.2, adj=0, cex=0.75, subtitle)

# Polymorphic loci vs. M/n
boxplot(polymorphic_loci~M,data=qubo.assembly.mat, xlab="M/n", ylab="Polymorphic loci")
title(main="Polymorphic loci across M/n", line=3)
mtext(side=3, line=1.5, at=1.2, adj=0, cex=0.75, subtitle)

# Number of SNPs vs. M/n
boxplot(number_of_snps~M,data=qubo.assembly.mat, xlab="M/n", ylab="Number of SNPs")
title(main="Number of SNPs (total) across M/n", line=3)
mtext(side=3, line=1.5, at=1.2, adj=0, cex=0.75, subtitle)

# ----QUBO gt-alpha----
# Coverage vs. gt-alpha
boxplot(depth_of_cov~gt.alpha,data=qubo.assembly.mat, xlab="gt-alpha", ylab="Depth of coverage")
title(main="Depth of Coverage across gt-alpha", line=3)
mtext(side=3, line=1.5, at=0.75, adj=0, cex=0.75, subtitle)

# Weighted coverage vs. gt-alpha
boxplot(weighted_cov~gt.alpha,data=qubo.assembly.mat, xlab="gt-alpha", ylab="Weighted coverage")
title(main="Weighted coverage across gt-alpha", line=3)
mtext(side=3, line=1.5, at=0.75, adj=0, cex=0.75, subtitle)

# Assembled loci vs. gt-alpha
boxplot(assembled_loci~gt.alpha,data=qubo.assembly.mat, xlab="gt-alpha", ylab="Assembled loci")
title(main="Assembled loci across gt-alpha", line=3)
mtext(side=3, line=1.5, at=0.75, adj=0, cex=0.75, subtitle)

# Polymorphic loci vs. gt-alpha
boxplot(polymorphic_loci~gt.alpha,data=qubo.assembly.mat, xlab="gt-alpha", ylab="Polymorphic loci")
title(main="Polymorphic loci across gt-alpha", line=3)
mtext(side=3, line=1.5, at=0.75, adj=0, cex=0.75, subtitle)

# Number of SNPs vs. gt-alpha
boxplot(number_of_snps~gt.alpha,data=qubo.assembly.mat, xlab="gt-alpha", ylab="Number of SNPs")
title(main="Number of SNPs (total) across gt-alpha", line=3)
mtext(side=3, line=1.5, at=0.75, adj=0, cex=0.75, subtitle)


# Two rows, two columns: summary charts for each assembly parameter
# Capture original graphing parameters
def.par <- par()
par(mfcol=c(2,2), mar = c(5,2,4,2), oma = c(2, 4, .5, .5), mgp = c(2, 0.6, 0))

# ----QUBO m----
# Coverage vs. m
boxplot(depth_of_cov~m,data=qubo.assembly.mat, xlab="m", ylab="Depth of coverage")
title(main="Depth of Coverage across m", line=3)
mtext(side=3, line=1.5, at=1.8, adj=0, cex=0.75, subtitle)

# Assembled loci vs. m
boxplot(assembled_loci~m,data=qubo.assembly.mat, xlab="m", ylab="Assembled loci")
title(main="Assembled loci across m", line=3)
mtext(side=3, line=1.5, at=1.8, adj=0, cex=0.75, subtitle)

# Polymorphic loci vs. m
boxplot(polymorphic_loci~m,data=qubo.assembly.mat, xlab="m", ylab="Polymorphic loci")
title(main="Polymorphic loci across m", line=3)
mtext(side=3, line=1.5, at=1.8, adj=0, cex=0.75, subtitle)

# Number of SNPs vs. m
boxplot(number_of_snps~m,data=qubo.assembly.mat, xlab="m", ylab="Number of SNPs")
title(main="Number of SNPs (total) across m", line=3)
mtext(side=3, line=1.5, at=1.8, adj=0, cex=0.75, subtitle)

# ----QUBO M/n----
# Coverage vs. M/n
boxplot(depth_of_cov~M,data=qubo.assembly.mat, xlab="M/n", ylab="Depth of coverage")
title(main="Depth of Coverage across M/n", line=3)
mtext(side=3, line=1.5, at=1.8, adj=0, cex=0.75, subtitle)

# Assembled loci vs. M/n
boxplot(assembled_loci~M,data=qubo.assembly.mat, xlab="M/n", ylab="Assembled loci")
title(main="Assembled loci across M/n", line=3)
mtext(side=3, line=1.5, at=1.8, adj=0, cex=0.75, subtitle)

# Polymorphic loci vs. M/n
boxplot(polymorphic_loci~M,data=qubo.assembly.mat, xlab="M/n", ylab="Polymorphic loci")
title(main="Polymorphic loci across M/n", line=3)
mtext(side=3, line=1.5, at=1.8, adj=0, cex=0.75, subtitle)

# Number of SNPs vs. M/n
boxplot(number_of_snps~M,data=qubo.assembly.mat, xlab="M/n", ylab="Number of SNPs")
title(main="Number of SNPs (total) across M/n", line=3)
mtext(side=3, line=1.5, at=1.8, adj=0, cex=0.75, subtitle)

# ----QUBO gt-alpha----
# Coverage vs. gt-alpha
boxplot(depth_of_cov~gt.alpha,data=qubo.assembly.mat, xlab="gt-alpha", ylab="Depth of coverage")
title(main="Depth of Coverage across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1, adj=0, cex=0.75, subtitle)

# Assembled loci vs. gt-alpha
boxplot(assembled_loci~gt.alpha,data=qubo.assembly.mat, xlab="gt-alpha", ylab="Assembled loci")
title(main="Assembled loci across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1, adj=0, cex=0.75, subtitle)

# Polymorphic loci vs. gt-alpha
boxplot(polymorphic_loci~gt.alpha,data=qubo.assembly.mat, xlab="gt-alpha", ylab="Polymorphic loci")
title(main="Polymorphic loci across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1, adj=0, cex=0.75, subtitle)

# Number of SNPs vs. gt-alpha
boxplot(number_of_snps~gt.alpha,data=qubo.assembly.mat, xlab="gt-alpha", ylab="Number of SNPs")
title(main="Number of SNPs (total) across gt-alpha", line=3)
mtext(side=3, line=1.5, at=1, adj=0, cex=0.75, subtitle)
