# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% POPULATION GENETIC STATISTICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script calculates the following population genetic statistics:
# heterozygosity, allele counts (i.e. number of variant sites across all loci), and allelic richness
# It does this for both QUAC (optimized de novo assembly) and QUBO 
# (aligned to the Q. robur reference genome) NextRAD datasets

# Metrics are calculated 1) between all wild populations, and 2) between garden individuals
# and a collection of all wild individuals (generated using the repool function)

library(adegenet)
library(hierfstat)

# %%%% QUAC %%%% ----

# READ IN GENIND FILE (QUAC DNFA; R0, min-maf=0; 1 SNP/locus)
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# Separate QUAC populations
QUAC.R0.NOMAF.pops <- seppop(QUAC.R0_NOMAF.genind)
# Repool wild populations into a single genind. Then bind this to garden
# This code takes an extremely long time to run--I stopped it after 7 minutes. Might have to do this in Stacks
# QUAC.wild <- repool(QUAC.R0.NOMAF.G_W$porterMt, QUAC.R0.NOMAF.G_W$magazineMt, QUAC.R0.NOMAF.G_W$pryorMt,
#                     QUAC.R0.NOMAF.G_W$sugarloaf_midlandPeak, QUAC.R0.NOMAF.G_W$kessler_shaleBarrenRidge)
# Rename QUAC.wild, so that all individuals have the pop name "wild"
# QUAC.R0.NOMAF.G_W <- repool(QUAC.R0.NOMAF.pops$garden, QUAC.wild)
# HETEROZYGOSITY
QUAC_HZ <- Hs(QUAC.R0_NOMAF.genind)

# Barplot for expected heterozygosity, SNP markers
barplot(QUAC_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf", 
                  "Kessler"), 
        main = "QUAC Heterozygosity: SNPs", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# ALLELE COUNTS, ALLELIC RICHNESS
# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUAC.R0_NOMAF.genind)
# Allelic richness: values per population
apply(allelic.richness(QUAC.R0_NOMAF.genind)$Ar, 2, mean, na.rm=TRUE)

# %%%% QUBO %%%% ----
# READ IN GENIND FILE (QUBO DNFA; R0, min-maf=0; 1 SNP/locus) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF/"
setwd(genpop.filePath)
QUBO.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# Number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUBO.R0_NOMAF.genind)

# HETEROZYGOSITY
QUBO_HZ <- Hs(QUBO.R0_NOMAF.genind)

# Barplot for expected heterozygosity, SNP markers
barplot(QUBO_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 11)),
        names = c("Garden", "Oakbrook", "Worldsong", "Irondale", "BTKC", 
                  "EBSCO_PL", "Peavine", "Wattsville", "MossRock", "EBSCO_Ridge", "HindsRoad", "Pop11"), 
        main = "QUBO Heterozygosity: SNPs", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# ALLELE COUNTS, ALLELIC RICHNESS
# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUBO.R0_NOMAF.genind)
# Allelic richness: values per population
apply(allelic.richness(QUBO.R0_NOMAF.genind)$Ar, 2, mean, na.rm=TRUE)
