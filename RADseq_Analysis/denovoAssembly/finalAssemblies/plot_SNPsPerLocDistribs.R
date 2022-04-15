# Read SNPs per locus prefilters
SNPsPerLocPre_table <- read.table("metric-SNPS_perLoc_Pre", header=TRUE, skip=1, sep ="\t", check.names = FALSE)

# Barplot: currently not working
barplot(as.matrix(SNPsPerLocPre_table), beside = TRUE, xlab="Number of SNPs", ylab="Number of loci", main="Distribution of number of SNPs per locus")
