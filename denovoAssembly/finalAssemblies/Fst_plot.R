# Fst plotting/heatmap

# Working directory set to /RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/

# QUAC de novo final assembly (DNFA)----
# Feed actual wild population names (rather than wild1, wild2, etc.)
QUAC_DNFA.names <- c("garden","Porter Mt.","Magazine Mt.","Pryor Mt.","Sugarloaf Mtns.",
                     "Kessler Mt.")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUAC.fst.mat <- as.matrix(read.table("QUAC/output/populations_Sum/QUAC_populations.fst_summary.tsv", header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUAC.fst.mat <- rbind(QUAC.fst.mat, rep(NA,6))
# Update row and column names
rownames(QUAC.fst.mat) <- colnames(QUAC.fst.mat) <- unique(QUAC_DNFA.names)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUAC.fst.mat), y=1:nrow(QUAC.fst.mat), z=t(QUAC.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="QUAC de novo: Fst Values")
# Add boundary lines
grid(nx=ncol(QUAC.fst.mat), ny=nrow(QUAC.fst.mat), col="black", lty=1)
# Include sample names, to understand the context of genetic distances
axis(1, 1:ncol(QUAC.fst.mat), colnames(QUAC.fst.mat), cex.axis=0.65, tick=FALSE)
axis(2, 1:nrow(QUAC.fst.mat), rownames(QUAC.fst.mat), cex.axis=0.65, gap.axis=0.1, tick=FALSE, pos=1.2)
for(x in 1:ncol(QUAC.fst.mat)){
  for(y in 1:nrow(QUAC.fst.mat)){
    text(x, y, QUAC.fst.mat[y,x])
  }
}

# QUBO de novo final assembly (DNFA)----
# Feed actual wild population names (rather than wild1, wild2, etc.)
QUBO_DNFA.names <- c("garden","Oakbrook","Worldsong","Irondale","Blue Trail Kings Chair", "EBSCO Parking lot",
                     "Peavine Falls","Wattsville","Moss Rock Preserve","EBSCO Ridge","Hinds Road","Pop11")

# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUBO.fst.mat <- as.matrix(read.table("QUBO/output/populations_Sum/QUBO_populations.fst_summary.tsv", header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUBO.fst.mat <- rbind(QUBO.fst.mat, rep(NA,12))
# Update row and column names
rownames(QUBO.fst.mat) <- colnames(QUBO.fst.mat) <- unique(QUBO_DNFA.names)

# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUBO.fst.mat), y=1:nrow(QUBO.fst.mat), z=t(QUBO.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="QUBO de novo: Fst Values")
# Add boundary lines
grid(nx=ncol(QUBO.fst.mat), ny=nrow(QUBO.fst.mat), col="black", lty=1)
# Include sample names, to understand the context of genetic distances
axis(1, 1:ncol(QUBO.fst.mat), colnames(QUBO.fst.mat), cex.axis=0.65, tick=FALSE)
axis(2, 1:nrow(QUBO.fst.mat), rownames(QUBO.fst.mat), cex.axis=0.65, gap.axis=0.1, tick=FALSE, pos=1.2)
for(x in 1:ncol(QUBO.fst.mat)){
  for(y in 1:nrow(QUBO.fst.mat)){
    text(x, y, QUBO.fst.mat[y,x])
  }
}
