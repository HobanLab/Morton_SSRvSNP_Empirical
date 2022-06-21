# %%%%%%%%%%%%%%%%%%%
# %%% MEASURE FST %%%
# %%%%%%%%%%%%%%%%%%%

# This script takes Fst values calculated in Stacks and plots them in a table using heatmap
# For both QUAC (optimized de novo assembly) and QUBO (aligned to the Q. robur reference genome),
# Fst maps are generated using a single (random) SNP per locus, and using all SNPs per locus

# Set plotting window to stack 2 Fst charts vertically
par(mfcol=c(2,1), oma=c(1,1,1,1))
par(mfcol=c(2,1), oma=rep(0.1,4))

# %%%% QUAC %%%% ----
# SINGLE (RANDOM) SNP PER LOCUS
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF/")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUAC.fst.mat <- as.matrix(read.table("populations.fst_summary.tsv", 
                                     header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUAC.fst.mat <- rbind(QUAC.fst.mat, rep(NA,5))
# Update row names
rownames(QUAC.fst.mat) <- colnames(QUAC.fst.mat)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUAC.fst.mat), y=1:nrow(QUAC.fst.mat), z=t(QUAC.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="QUAC Fst Values: Single (random) SNP/locus")
# Add boundary lines
grid(nx=ncol(QUAC.fst.mat), ny=nrow(QUAC.fst.mat), col="black", lty=1)
# Include sample names, to understand the context of genetic distances
axis(1, 1:ncol(QUAC.fst.mat), colnames(QUAC.fst.mat), cex.axis=1.2, tick=FALSE)
text(1, c(1:5), labels=rownames(QUAC.fst.mat), cex=1.2)
for(x in 1:ncol(QUAC.fst.mat)){
  for(y in 1:nrow(QUAC.fst.mat)){
    text(x, y, QUAC.fst.mat[y,x], cex=1.5)
  }
}

# ALL SNPS PER LOCUS
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF_AllSNPs/")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUAC.fst.mat <- as.matrix(read.table("populations.fst_summary.tsv", 
                                     header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUAC.fst.mat <- rbind(QUAC.fst.mat, rep(NA,5))
# Update row names
rownames(QUAC.fst.mat) <- colnames(QUAC.fst.mat)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUAC.fst.mat), y=1:nrow(QUAC.fst.mat), z=t(QUAC.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="QUAC Fst Values: All SNPs/locus")
# Add boundary lines
grid(nx=ncol(QUAC.fst.mat), ny=nrow(QUAC.fst.mat), col="black", lty=1)
# Include sample names, to understand the context of genetic distances
axis(1, 1:ncol(QUAC.fst.mat), colnames(QUAC.fst.mat), cex.axis=1.2, tick=FALSE)
text(1, c(1:5), labels=rownames(QUAC.fst.mat), cex=1.2)
for(x in 1:ncol(QUAC.fst.mat)){
  for(y in 1:nrow(QUAC.fst.mat)){
    text(x, y, QUAC.fst.mat[y,x], cex=1.5)
  }
}

# %%%% QUBO %%%% ----
# SINGLE (RANDOM) SNP PER LOCUS
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0_NOMAF/")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUBO.fst.mat <- as.matrix(read.table("populations.fst_summary.tsv", 
                                     header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUBO.fst.mat <- rbind(QUBO.fst.mat, rep(NA,10))
# Update row names
rownames(QUBO.fst.mat) <- colnames(QUBO.fst.mat)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUBO.fst.mat), y=1:nrow(QUBO.fst.mat), z=t(QUBO.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="QUBO Fst Values: Single (random) SNP/locus")
# Add boundary lines
grid(nx=ncol(QUBO.fst.mat), ny=nrow(QUBO.fst.mat), col="black", lty=1)
# Include sample names, to understand the context of genetic distances
axis(1, 1:ncol(QUBO.fst.mat), colnames(QUBO.fst.mat), cex.axis=1.2, tick=FALSE)
text(1, c(1:11), labels=rownames(QUBO.fst.mat), cex=1.2)
for(x in 1:ncol(QUBO.fst.mat)){
  for(y in 1:nrow(QUBO.fst.mat)){
    text(x, y, QUBO.fst.mat[y,x], cex=1.5)
  }
}

# ALL SNPS PER LOCUS
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0_NOMAF_AllSNPs/")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUBO.fst.mat <- as.matrix(read.table("populations.fst_summary.tsv", 
                                     header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUBO.fst.mat <- rbind(QUBO.fst.mat, rep(NA,10))
# Update row names
rownames(QUBO.fst.mat) <- colnames(QUBO.fst.mat)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUBO.fst.mat), y=1:nrow(QUBO.fst.mat), z=t(QUBO.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="QUBO Fst Values: All SNPs/locus")
# Add boundary lines
grid(nx=ncol(QUBO.fst.mat), ny=nrow(QUBO.fst.mat), col="black", lty=1)
# Include sample names, to understand the context of genetic distances
axis(1, 1:ncol(QUBO.fst.mat), colnames(QUBO.fst.mat), cex.axis=1.2, tick=FALSE)
text(1, c(1:11), labels=rownames(QUBO.fst.mat), cex=1.2)
for(x in 1:ncol(QUBO.fst.mat)){
  for(y in 1:nrow(QUBO.fst.mat)){
    text(x, y, QUBO.fst.mat[y,x], cex=1.5)
  }
}
