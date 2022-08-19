# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% POPULATION GENETIC STATISTICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script calculates the following population genetic statistics:
# heterozygosity, allele counts (i.e. number of variant sites across all loci), and allelic richness
# It does this for both QUAC (optimized de novo assembly) and QUBO 
# (aligned to the Q. robur reference genome) NextRAD datasets

# TO DO: Make this script calculate heterozygostiy, allele counts, and allelic richness in the
# following scenarios
# 1. MSAT: all samples
# 2. MSAT: subsampled down to NextRAD samples (for garden and wild--so, you need to update your current code)
# 3. SNP: all samples (garden and wild)
# 4. SNP: wild samples
# 5. SNP: subsampled down to MSAT samples (for garden and wild--so, you need to update your current code)

# Metrics are calculated 1) between garden individuals and a collection of all wild individuals ("wild"), 
# and 2) between all wild individuals. These groupings are achieved using the Stacks populations module

library(adegenet)
library(hierfstat)

# %%%% QUAC %%%% ----
# ---- MSATS ----
# Read in genind file (GCC_QUAC_ZAIN repo; QUAC_wK_garden_wild_clean.gen)
QUAC.MSAT.genpop.filePath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/Garden_Wild/"
setwd(QUAC.MSAT.genpop.filePath)
QUAC.MSAT.genind <- read.genepop("QUAC_wK_garden_wild_clean.gen", quiet = TRUE, ncode = 3)
# Correct popNames: pop1 is Garden, pop2 is Wild
pop(QUAC.MSAT.genind) <- gsub("pop1", "garden", pop(QUAC.MSAT.genind))
pop(QUAC.MSAT.genind) <- gsub("pop2", "wild", pop(QUAC.MSAT.genind))

# Heterozygosity
QUAC.MSAT_HZ <- Hs(QUAC.MSAT.genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.MSAT_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: SNPs", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Allele counts, allelic richness
ncol(QUAC.MSAT.genind@tab)
# Allelic richness: values per population
apply(allelic.richness(QUAC.MSAT.genind)$Ar, 2, mean, na.rm=TRUE)

# SUBSET MSATS ----
# Correct sample names: read in tissue database names from GCC_QUAC_ZAIN repository, and rename genind object rows
QUAC.MSAT.sampleNames_filepath <- "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Data_Frames/QUAC_allpop_clean_df.csv"
QUAC.MSAT.sampleNames <- unlist(read.csv2(QUAC.MSAT.sampleNames_filepath, header = TRUE, sep=",")[1])
rownames(QUAC.MSAT.genind@tab) <- QUAC.MSAT.sampleNames
# Create a matrix of strictly wild samples

# Heterozygosity
QUAC.MSAT_HZ <- Hs(QUAC.MSAT.genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.MSAT_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: SNPs", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Allele counts, allelic richness
ncol(QUAC.MSAT.genind@tab)
# Allelic richness: values per population
apply(allelic.richness(QUAC.MSAT.genind)$Ar, 2, mean, na.rm=TRUE)




# ---- SNPS ----
# GARDEN AND WILD ----
# Read in genind file: Optimized de novo assembly; R0, NOMAF, first SNP/locus, 2 populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.R0_NOMAF.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])

# Heterozygosity
QUAC.SNP_HZ <- Hs(QUAC.R0_NOMAF.genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: SNPs", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Allele counts, allelic richness
# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUAC.R0_NOMAF.genind)
# Allelic richness: values per population
apply(allelic.richness(QUAC.R0_NOMAF.genind)$Ar, 2, mean, na.rm=TRUE)

# WILD ONLY ----
# Read in genind file: Optimized de novo assembly; R0, NOMAF, first SNP/locus, only wild populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF_1SNP/"
setwd(genpop.filePath)
QUAC.Wild_R0.NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.Wild_R0.NOMAF.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# Heterozygosity
QUAC.SNP.Wild_HZ <- Hs(QUAC.Wild_R0.NOMAF.genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP.Wild_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf", 
                  "Kessler"), 
        main = "QUAC Heterozygosity: SNPs", 
        xlab = "Wild population", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# ALLELE COUNTS, ALLELIC RICHNESS
# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUAC.Wild_R0.NOMAF.genind)
# Allelic richness: values per population
apply(allelic.richness(QUAC.Wild_R0.NOMAF.genind)$Ar, 2, mean, na.rm=TRUE)

# ---- MSATS & SNPS, SUBSET ----


# %%%% QUBO %%%% ----
# ---- GARDEN AND WILD ----
# MSATS ----


# SNPS ----
# Read in genind file: GSNAP4 alignment with Quercus robur genome; R0, NOMAF, first SNP/locus, 2 populations 
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO_R0.NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO_R0.NOMAF.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# Heterozygosity
QUBO_HZ <- Hs(QUBO_R0.NOMAF.genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Heterozygosity: SNPs", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Allele counts, allelic richness
# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUBO_R0.NOMAF.genind)
# Allelic richness: values per population
apply(allelic.richness(QUBO_R0.NOMAF.genind)$Ar, 2, mean, na.rm=TRUE)

# ---- WILD ONLY ----
# Read in genind file: GSNAP4 alignment with Quercus robur genome; R0, NOMAF, first SNP/locus, only wild populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF_1SNP/"
setwd(genpop.filePath)
QUBO.Wild_R0.NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.Wild_R0.NOMAF.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# Heterozygosity
QUBO.Wild_HZ <- Hs(QUBO.Wild_R0.NOMAF.genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.Wild_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf", 
                  "Kessler"), 
        main = "QUAC Heterozygosity: SNPs", 
        xlab = "Wild source", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# ALLELE COUNTS, ALLELIC RICHNESS
# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUBO.Wild_R0.NOMAF.genind)
# Allelic richness: values per population
apply(allelic.richness(QUBO.Wild_R0.NOMAF.genind)$Ar, 2, mean, na.rm=TRUE)
