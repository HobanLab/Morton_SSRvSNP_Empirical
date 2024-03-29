# Fst plotting/heatmap

# Specify sample names for plotting
# These have to be truncated, in order for the plot to be legible
fst.names <- c(rep("QUAC_AA1",2),rep("QUAC_SLM1",2),"QUAC_SLM2","QUAC_KSM","QUBO_ABG","QUBO_TMA",
               "QUBO_BTK", rep("QUBO_EBS1",2),"QUBO_EBS2")

# Specify which Fst matrix to pull values from: subset (1,000), or all (47,266) loci
# This output file is generated by the populations module in Stacks, when specifying the parameter "--fstats"
# Subset
#setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/PreliminarySubset/output/populations_1000/")
plt.title <- "Preliminary Subset Fst Values: Subset (1,000 loci)"
# All
#setwd("/RAID1/IMLS_GCCO/Analysis/Stacks/PreliminarySubset/output/populations_all/")
plt.title <- "Preliminary Subset Fst Values: Subset (All loci)"


# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
fst.mat <- as.matrix(read.table("populations.fst_summary.tsv", header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
fst.mat <- rbind(fst.mat, rep(NA,12))
# Update row and column names
rownames(fst.mat) <- colnames(fst.mat) <- fst.names


# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(fst.mat), y=1:nrow(fst.mat), z=t(fst.mat), axes=FALSE, xlab="", ylab="", 
      main=plt.title)
# Add boundary lines
grid(nx=ncol(fst.mat), ny=nrow(fst.mat), col="black", lty=1)
# Include sample names, to understand the context of genetic distances
axis(1, 1:ncol(fst.mat), colnames(fst.mat), cex.axis=0.65, tick=FALSE)
axis(2, 1:nrow(fst.mat), rownames(fst.mat), cex.axis=0.65, gap.axis=0.1, tick=FALSE, pos=1.2)
for(x in 1:ncol(fst.mat)){
  for(y in 1:nrow(fst.mat)){
    text(x, y, fst.mat[y,x])
  }
}
