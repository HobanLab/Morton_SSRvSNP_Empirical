# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% POPULATION GENETIC STATISTICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script calculates popgulation genetic statistics (heterozygosity, allelic richness, Fst)
# for the two study species of the SSRvSNP study: Quercus acerifolia (QUAC) and Q. boyntonii (QUBO)

library(adegenet)
library(pegas)
library(hierfstat)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ---- QUERCUS ACERIFOLIA ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%

# READ IN/PROCESS GEN FILE----
genpop.filePath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_sum/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
QUAC.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Read in the Stacks popmap values, and use these to replace @pop values 
# (using the pops accessor; this is necessary because original pop names are incorrect)
pop(QUAC.genind) <- factor(read.table("../../../QUAC_popmap2", header=FALSE)[,2])

# HARDY-WEINBERG EQUILIBRIUM----
hw.test(QUAC.genind) # Stalls out--this is likely too many loci to process

# HETEROZYGOSITY----
Hs(QUAC.genind)

# Barplot for expected heterozygosity, SNP markers
barplot(Hs(QUAC.genind), beside = TRUE, 
        ylim = c(0,0.06), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf", 
                  "Kessler"), 
        main = "QUAC Heterozygosity: SNPs", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0)

# ALLELIC RICHNESS----
# Values per population
apply(allelic.richness(QUAC.genind)$Ar, 2, mean, na.rm=TRUE)

# FST----
# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUAC.fst.mat <- as.matrix(read.table("populations.fst_summary.tsv", header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUAC.fst.mat <- rbind(QUAC.fst.mat, rep(NA,6))
# Update row and column names
# rownames(QUAC.fst.mat) <- colnames(QUAC.fst.mat) <- levels(unique(pop(QUAC.genind)))
QUAC_popNames <- c("Garden", "Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf", 
                   "Kessler")
rownames(QUAC.fst.mat) <- colnames(QUAC.fst.mat) <- QUAC_popNames
# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUAC.fst.mat), y=1:nrow(QUAC.fst.mat), z=t(QUAC.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="QUAC Fst Values: SNPs")
# Add boundary lines
grid(nx=ncol(QUAC.fst.mat), ny=nrow(QUAC.fst.mat), col="black", lty=1)
# Include sample names, to understand the context of genetic distances
axis(1, 1:ncol(QUAC.fst.mat), colnames(QUAC.fst.mat), cex.axis=1.2, tick=FALSE)
text(1, c(1:6), labels=rownames(QUAC.fst.mat), cex=1.2)
for(x in 1:ncol(QUAC.fst.mat)){
  for(y in 1:nrow(QUAC.fst.mat)){
    text(x, y, QUAC.fst.mat[y,x], cex=1.5)
  }
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# ---- QUERCUS BOYNTONII ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# READ IN/PROCESS GEN FILE----
genpop.filePath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP/output/populations_sum/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
QUBO.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Read in the Stacks popmap values, and use these to replace @pop values 
# (using the pops accessor; this is necessary because original pop names are incorrect)
pop(QUBO.genind) <- factor(read.table("../../../QUBO_popmap2", header=FALSE)[,2])

# HARDY-WEINBERG EQUILIBRIUM----
hw.test(QUBO.genind) # Stalls out--this is likely too many loci to process

# HETEROZYGOSITY----
Hs(QUBO.genind)

# Barplot for expected heterozygosity, SNP markers
barplot(Hs(QUBO.genind), beside = TRUE, 
        ylim = c(0,0.06), col = c("darkseagreen1", rep("darkgreen", 11)),
        names = c("Garden", "Oakbrook", "Worldsong", "Irondale", "BTKC", 
                  "EBSCO_PL", "Peavine", "Wattsville", "MossRock", "EBSCO_Ridge", "HindsRoad", "Pop11"), 
        main = "QUBO Heterozygosity: SNPs", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0)

# ALLELIC RICHNESS----
# Values per population
apply(allelic.richness(QUBO.genind)$Ar, 2, mean, na.rm=TRUE)

# FST----
# Stacks exports an Fst table when the populations command is given the --fstats argument. Read this in
QUBO.fst.mat <- as.matrix(read.table("populations.fst_summary.tsv", header=TRUE, row.names=1, sep = "\t"))
# Add a row at the bottom to make matrix symmetrical (nrow=ncol)
QUBO.fst.mat <- rbind(QUBO.fst.mat, rep(NA,12))
# Update row and column names
# rownames(QUBO.fst.mat) <- colnames(QUBO.fst.mat) <- levels(unique(pop(QUBO.genind)))
QUBO_popNames <- c("Garden", "Oakbrook", "Worldsong", "Irondale", "BTKC", 
                   "EBSCO_PL", "Peavine", "Wattsville", "MossRock", "EBSCO_Ridge", "HindsRoad", "Pop11") 
rownames(QUBO.fst.mat) <- colnames(QUBO.fst.mat) <- QUBO_popNames
# Use image command to plot a heatmap
# First two arguments specify the boundaries of the heatmap; z provides actual values
# z is transposed in order to plot numeral values later on
image(x=1:ncol(QUBO.fst.mat), y=1:nrow(QUBO.fst.mat), z=t(QUBO.fst.mat), axes=FALSE, xlab="", ylab="", 
      main="QUBO Fst Values: SNPs")
# Add boundary lines
grid(nx=ncol(QUBO.fst.mat), ny=nrow(QUBO.fst.mat), col="black", lty=1)
# Include sample names, to understand the context of genetic distances
axis(1, 1:ncol(QUBO.fst.mat), colnames(QUBO.fst.mat), cex.axis=1.1, tick=FALSE)
text(1, c(1:12), labels=rownames(QUBO.fst.mat), cex=1.1)
for(x in 1:ncol(QUBO.fst.mat)){
  for(y in 1:nrow(QUBO.fst.mat)){
    text(x, y, QUBO.fst.mat[y,x], cex=1.1)
  }
}
