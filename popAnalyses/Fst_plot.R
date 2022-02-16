# Fst plotting/heatmap
# Script for plotting Fst charts generated for de novo, reference, and hybrid Stacks datasets
# These Fst values are derived for wild populations only, and are pulled from the populations.fst_summary.tsv Stacks file

# Set plotting window to stack 3 Fst charts vertically
par(mfcol=c(3,1), oma=c(1,1,1,1))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% Quercus acerifolia (QUAC) %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# De novo final assembly (DNFA)----
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUAC.fst.mat <- as.matrix(read.table("output/populations_wild/populations.fst_summary.tsv", header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUAC.fst.mat <- rbind(QUAC.fst.mat, rep(NA,5))
# Update row and column names
rownames(QUAC.fst.mat) <- colnames(QUAC.fst.mat) <- unique(QUAC_wild.popnames)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUAC.fst.mat), y=1:nrow(QUAC.fst.mat), z=t(QUAC.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="QUAC Fst Values: De novo (6,850 loci)")
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

# Reference dataset (filtered, GSNAP aligned; REF)----
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUAC.fst.mat <- as.matrix(read.table("GSNAP/output/populations_wild/populations.fst_summary.tsv", header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUAC.fst.mat <- rbind(QUAC.fst.mat, rep(NA,5))
# Update row and column names
rownames(QUAC.fst.mat) <- colnames(QUAC.fst.mat) <- unique(QUAC_wild.popnames)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUAC.fst.mat), y=1:nrow(QUAC.fst.mat), z=t(QUAC.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="Reference (17,844 loci)")
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

# Hybrid dataset (optimized de novo consensus loci, BWA aligned; HYB)----
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_denovoConsensus/QUAC/")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUAC.fst.mat <- as.matrix(read.table("output/populations_wild/populations.fst_summary.tsv", header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUAC.fst.mat <- rbind(QUAC.fst.mat, rep(NA,5))
# Update row and column names
rownames(QUAC.fst.mat) <- colnames(QUAC.fst.mat) <- unique(QUAC_wild.popnames)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUAC.fst.mat), y=1:nrow(QUAC.fst.mat), z=t(QUAC.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="Hybrid (1,701 loci)")
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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% Quercus boyntonii (QUBO) %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# QUBO de novo final assembly (DNFA)----
# Feed actual wild population names (rather than wild1, wild2, etc.)
QUBO_wild.popnames <- c("Oakbrook","Worldsong","Irondale","Blue Trail Kings Chair", "EBSCO Parking lot",
                     "Peavine Falls","Wattsville","Moss Rock Preserve","EBSCO Ridge","Hinds Road","Pop11")

# De novo final assembly (DNFA)----
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUBO.fst.mat <- as.matrix(read.table("output/populations_wild/populations.fst_summary.tsv", header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUBO.fst.mat <- rbind(QUBO.fst.mat, rep(NA,11))
# Update row and column names
rownames(QUBO.fst.mat) <- colnames(QUBO.fst.mat) <- unique(QUBO_wild.popnames)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUBO.fst.mat), y=1:nrow(QUBO.fst.mat), z=t(QUBO.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="QUBO Fst Values: De novo (5,621 loci)")
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

# Reference dataset (filtered, GSNAP aligned; REF)----
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUBO.fst.mat <- as.matrix(read.table("GSNAP/output/populations_wild/populations.fst_summary.tsv", header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUBO.fst.mat <- rbind(QUBO.fst.mat, rep(NA,11))
# Update row and column names
rownames(QUBO.fst.mat) <- colnames(QUBO.fst.mat) <- unique(QUBO_wild.popnames)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUBO.fst.mat), y=1:nrow(QUBO.fst.mat), z=t(QUBO.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="Reference (40,998 loci)")
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

# Hybrid dataset (optimized de novo consensus loci, BWA aligned; HYB)----
setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_denovoConsensus/QUBO/")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUBO.fst.mat <- as.matrix(read.table("output/populations_wild/populations.fst_summary.tsv", header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUBO.fst.mat <- rbind(QUBO.fst.mat, rep(NA,11))
# Update row and column names
rownames(QUBO.fst.mat) <- colnames(QUBO.fst.mat) <- unique(QUBO_wild.popnames)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUBO.fst.mat), y=1:nrow(QUBO.fst.mat), z=t(QUBO.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="Hybrid (2,160 loci)")
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
