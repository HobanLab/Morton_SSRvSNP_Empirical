# Extracting de novo assembly metrics from text files generated using stacks-dist-extract

# Extract coverage values for each sample and calculate mean
depth.of.cov <- mean(read.table("depth_of_cov", header=FALSE, skip=2)[,1])
# Read in total assembled loci
assembled.loci <- unique(read.table("assembled_loci", header=TRUE, skip=3))
# Read in total polymorphic loci (take sum of all loci that aren't monomorphic)
polymorphic.loci <- sum(read.table("polymorphic_loci", header=FALSE)[,1])
# Read in total number of SNPs
number.of.snps <- unique(read.table("number_of_SNPs", header=TRUE, skip=3))

