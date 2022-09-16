# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% POPULATION GENETIC STATISTICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script calculates the following population genetic statistics:
# heterozygosity, allele counts (i.e. number of variant sites across all loci), and allelic richness
# It does this for both QUAC (optimized de novo assembly) and QUBO 
# (aligned to the Q. robur reference genome) NextRAD datasets

# The script performs these calculations with different sets of samples and loci, as outlined below:
# 1. Microsatellites: all samples, garden and wild populations
# 2. SNPs: all samples, garden and wild populations, R0 filter (no filter for missing data)
# 3. SNPs: all samples, garden and wild populations, R80 filter (loci present in 80% of samples), 1st SNP/locus
# 4. SNPs: all samples, garden and wild populations, R80 filter (loci present in 80% of samples), haplotype-wise SNPs
# 5. SNPs: all samples, wild populations only, R0 filter (no filter for missing data)
# 6. SNPs: all samples, wild populations only, R80 filter (loci present in 80% of samples)
# 7. Microsatellites: subset samples (only those shared with SNP dataset), garden and wild populations
# 8. SNPs: subset samples (only those shared with SNP dataset), garden and wild populations, R0 filter
# 9. SNPs: subset samples (only those shared with SNP dataset), garden and wild populations, R80 filter

# Each of these scenarios is explored for both species

library(adegenet)
library(hierfstat)
# Set working directory
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)

# %%%% QUAC %%%% ----
# ---- MSATS ----
# %%% GARDEN AND WILD ----
# Read in genind file (GCC_QUAC_ZAIN repo; QUAC_wK_garden_wild_clean.gen)
QUAC.MSAT.genpop.filePath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/Garden_Wild/"
setwd(QUAC.MSAT.genpop.filePath)
QUAC.MSAT_genind <- read.genepop("QUAC_wK_garden_wild_clean.gen", ncode = 3)
# Correct popNames: pop1 is Garden, pop2 is Wild
pop(QUAC.MSAT_genind) <- gsub("pop1", "garden", pop(QUAC.MSAT_genind))
pop(QUAC.MSAT_genind) <- gsub("pop2", "wild", pop(QUAC.MSAT_genind))

# Heterozygosity
QUAC.MSAT_HZ <- Hs(QUAC.MSAT_genind); print(QUAC.MSAT_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.MSAT_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: MSATs (Complete)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUAC.MSAT_genind); ncol(QUAC.MSAT_genind@tab)
# Allelic richness: values per population
QUAC.MSAT_AR <- apply(allelic.richness(QUAC.MSAT_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUAC.MSAT_AR)

# ---- SNPS ----
# %%% GARDEN AND WILD ----
# R0 ----
# Read in genind file: Optimized de novo assembly; R0, NOMAF, first SNP/locus, 2 populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.R0_genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.R0_genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])

# Heterozygosity
QUAC.SNP.R0_HZ <- Hs(QUAC.SNP.R0_genind); print(QUAC.SNP.R0_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP.R0_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: SNPs (Complete), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUAC.SNP.R0_genind); ncol(QUAC.SNP.R0_genind@tab)
# Allelic richness: values per population
QUAC.SNP.R0_AR <- apply(allelic.richness(QUAC.SNP.R0_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUAC.SNP.R0_AR)

# R80 ----
# 1st SNP/locus ----
# Read in genind file: Optimized de novo assembly; R80, NOMAF, first SNP/locus, 2 populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.R80_genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.R80_genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])

# Heterozygosity
QUAC.SNP.R80_HZ <- Hs(QUAC.SNP.R80_genind); print(QUAC.SNP.R80_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP.R80_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: SNPs (Complete), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUAC.SNP.R80_genind); ncol(QUAC.SNP.R80_genind@tab)
# Allelic richness: values per population
QUAC.SNP.R80_AR <- apply(allelic.richness(QUAC.SNP.R80_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUAC.SNP.R80_AR)

# Microhaplotypes ----
# Read in genind file: Optimized de novo assembly; R80, NOMAF, haplotype-wise SNPs, 2 populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.R80.HapSNP_genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.R80.HapSNP_genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])

# Heterozygosity
QUAC.SNP.R80.HapSNP_HZ <- Hs(QUAC.SNP.R80.HapSNP_genind); print(QUAC.SNP.R80.HapSNP_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP.R80.HapSNP_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: SNPs (Complete), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUAC.SNP.R80.HapSNP_genind); ncol(QUAC.SNP.R80.HapSNP_genind@tab)
# Allelic richness: values per population
QUAC.SNP.R80.HapSNP_AR <- apply(allelic.richness(QUAC.SNP.R80.HapSNP_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUAC.SNP.R80.HapSNP_AR)

# %%% WILD ONLY ----
# R0 ----
# Read in genind file: Optimized de novo assembly; R0, NOMAF, first SNP/locus, only wild populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF_1SNP/"
setwd(genpop.filePath)
QUAC.Wild.SNP.R0_genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.Wild.SNP.R0_genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# Heterozygosity
QUAC.SNP.Wild.R0_HZ <- Hs(QUAC.Wild.SNP.R0_genind); print(QUAC.SNP.Wild.R0_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP.Wild.R0_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf", 
                  "Kessler"), 
        main = "QUAC Wild Heterozygosity: SNPs (Complete), R0", 
        xlab = "Wild population", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUAC.Wild.SNP.R0_genind); ncol(QUAC.Wild.SNP.R0_genind@tab)
# Allelic richness: values per population
QUAC.SNP.Wild.R0_AR <- apply(allelic.richness(QUAC.Wild.SNP.R0_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUAC.SNP.Wild.R0_AR)

# R80 ----
# Read in genind file: Optimized de novo assembly; R80, NOMAF, first SNP/locus, only wild populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R80_NOMAF_1SNP/"
setwd(genpop.filePath)
QUAC.Wild.SNP.R80_genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.Wild.SNP.R80_genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# Heterozygosity
QUAC.SNP.Wild.R80_HZ <- Hs(QUAC.Wild.SNP.R80_genind); print(QUAC.SNP.Wild.R80_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP.Wild.R80_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf", 
                  "Kessler"), 
        main = "QUAC Wild Heterozygosity: SNPs (Complete), R80", 
        xlab = "Wild population", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUAC.Wild.SNP.R80_genind); ncol(QUAC.Wild.SNP.R80_genind@tab)
# Allelic richness: values per population
QUAC.SNP.Wild.R80_AR <- apply(allelic.richness(QUAC.Wild.SNP.R80_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUAC.SNP.Wild.R80_AR)

# ---- SUBSET ----
# MSAT: read in Tissue database names from GCC_QUAC_ZAIN repository, and rename MSAT genind matrix
QUAC.MSAT.tissueNames_filepath <- "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Data_Frames/QUAC_allpop_clean_df.csv"
QUAC.MSAT.tissueNames <- unlist(read.csv2(QUAC.MSAT.tissueNames_filepath, header = TRUE, sep=",")[1])
rownames(QUAC.MSAT_genind@tab) <- QUAC.MSAT.tissueNames

# SNP: read in Tissue database names, and rename SNP genind matrices (both R0 and R80)
# This file was created (by Austin K.), and can be found on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
QUAC.SNP.tissueNames_filepath <- paste0(SSRvSNP.wd,"exSituRepresentation/Resampling/QUAC_SNP_TissueNames.csv")
QUAC.SNP.tissueNames <- unlist(read.csv2(QUAC.SNP.tissueNames_filepath, header = TRUE, sep = ",")[3])
rownames(QUAC.SNP.R0_genind@tab) <- rownames(QUAC.SNP.R80_genind@tab) <- QUAC.SNP.tissueNames

# Subset SNP sample names by those that are also seen within the MSAT samples 
QUAC_sharedSamples <- sort(QUAC.SNP.tissueNames[which(QUAC.SNP.tissueNames %in% QUAC.MSAT.tissueNames)])
# Subset MSAT and SNP genind matrices to strictly shared samples, dropping now absent alleles
QUAC.MSAT.subset_genind <- QUAC.MSAT_genind[QUAC_sharedSamples,, drop=TRUE]
QUAC.SNP.R0.subset_genind <- QUAC.SNP.R0_genind[QUAC_sharedSamples,, drop=TRUE]
QUAC.SNP.R80.subset_genind <- QUAC.SNP.R80_genind[QUAC_sharedSamples,, drop=TRUE]

# Heterozygosity
QUAC.MSAT.subset_HZ <- Hs(QUAC.MSAT.subset_genind); print(QUAC.MSAT.subset_HZ)
QUAC.SNP.R0.subset_HZ <- Hs(QUAC.SNP.R0.subset_genind); print(QUAC.SNP.R0.subset_HZ)
QUAC.SNP.R80.subset_HZ <- Hs(QUAC.SNP.R80.subset_genind); print(QUAC.SNP.R80.subset_HZ)
# Barplot for expected heterozygosity, MSAT subset
barplot(QUAC.MSAT.subset_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: MSATs (Subset)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Barplot for expected heterozygosity, SNP R0 subset
barplot(QUAC.SNP.R0.subset_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: SNPs (Subset), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Barplot for expected heterozygosity, SNP R80 subset
barplot(QUAC.SNP.R80.subset_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: SNPs (Subset), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUAC.MSAT.subset_genind); ncol(QUAC.MSAT.subset_genind@tab)
nLoc(QUAC.SNP.R0.subset_genind); ncol(QUAC.SNP.R0.subset_genind@tab)
nLoc(QUAC.SNP.R80.subset_genind); ncol(QUAC.SNP.R80.subset_genind@tab)
# Allelic richness: values per population
QUAC.MSAT.subset_AR <- apply(allelic.richness(QUAC.MSAT.subset_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUAC.MSAT.subset_AR)
QUAC.SNP.R0.subset_AR <- apply(allelic.richness(QUAC.SNP.R0.subset_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUAC.SNP.R0.subset_AR)
QUAC.SNP.R80.subset_AR <- apply(allelic.richness(QUAC.SNP.R80.subset_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUAC.SNP.R80.subset_AR)

# %%%% QUBO %%%% ----
# ---- MSATS ----
# %%% GARDEN AND WILD ----
# Read in genind file (Southeast Oaks repo; genetic_data/Qb_total.gen file)
genpop.filePath <- 
  "~/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
setwd(genpop.filePath)
QUBO.MSAT_genind <- read.genepop(paste0(genpop.filePath,"Qb_total.gen"), ncode=3)
# Correct popNames: last population (IMLS4_MP1_IMLS336_C05) is garden; rest (9) are wild
levels(QUBO.MSAT_genind@pop) <- c(rep("wild",9), "garden") 

# Heterozygosity
QUBO.MSAT_HZ <- Hs(QUBO.MSAT_genind); print(QUBO.MSAT_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.MSAT_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: MSATs (Complete)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUBO.MSAT_genind); ncol(QUBO.MSAT_genind@tab)
# Allelic richness: values per population
QUBO.MSAT_AR <- apply(allelic.richness(QUBO.MSAT_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUBO.MSAT_AR)

# ---- SNPS ----
# %%% GARDEN AND WILD ----
# R0 ----
# Read in genind file: GSNAP4 alignment with Quercus robur genome; R0, NOMAF, first SNP/locus, 2 populations 
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.R0_genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.R0_genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# Heterozygosity
QUBO.SNP.R0_HZ <- Hs(QUBO.SNP.R0_genind); print(QUBO.SNP.R0_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.SNP.R0_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs (Complete), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUBO.SNP.R0_genind); ncol(QUBO.SNP.R0_genind@tab)
# Allelic richness: values per population
QUBO.SNP.R0_AR <- apply(allelic.richness(QUBO.SNP.R0_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUBO.SNP.R0_AR)

# R80 ----
# 1st SNP/locus ----
# Read in genind file: GSNAP4 alignment with Quercus robur genome; R80, NOMAF, first SNP/locus, 2 populations 
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.R80_genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.R80_genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# Heterozygosity
QUBO.SNP.R80_HZ <- Hs(QUBO.SNP.R80_genind); print(QUBO.SNP.R80_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.SNP.R80_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs (Complete), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUBO.SNP.R80_genind); ncol(QUBO.SNP.R80_genind@tab)
# Allelic richness: values per population
QUBO.SNP.R80_AR <- apply(allelic.richness(QUBO.SNP.R80_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUBO.SNP.R80_AR)

# Microhaplotypes ----
# Read in genind file: GSNAP4 alignment with Quercus robur genome; R80, NOMAF, haplotype-wise SNPs, 2 populations 
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.R80.HapSNP_genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.R80.HapSNP_genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# Heterozygosity
QUBO.SNP.R80.HapSNP_HZ <- Hs(QUBO.SNP.R80.HapSNP_genind); print(QUBO.SNP.R80.HapSNP_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.SNP.R80.HapSNP_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs (Complete), R80 (Haplotype-wise SNPs)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUBO.SNP.R80.HapSNP_genind); ncol(QUBO.SNP.R80.HapSNP_genind@tab)
# Allelic richness: values per population
QUBO.SNP.R80.HapSNP_AR <- apply(allelic.richness(QUBO.SNP.R80.HapSNP_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUBO.SNP.R80.HapSNP_AR)

# %%% WILD ONLY ----
# R0 ----
# Read in genind file: GSNAP4 alignment with Quercus robur genome; R0, NOMAF, first SNP/locus, only wild populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0_NOMAF_1SNP/"
setwd(genpop.filePath)
QUBO.Wild.SNP.R0_genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.Wild.SNP.R0_genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# Heterozygosity
QUBO.Wild.SNP.R0_HZ <- Hs(QUBO.Wild.SNP.R0_genind); print(QUBO.Wild.SNP.R0_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.Wild.SNP.R0_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf", 
                  "Kessler"), 
        main = "QUBO Wild Heterozygosity: SNPs (Complete), R0", 
        xlab = "Wild source", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUBO.Wild.SNP.R0_genind); ncol(QUBO.Wild.SNP.R0_genind@tab)
# Allelic richness: values per population
QUBO.Wild.SNP.R0_AR <- apply(allelic.richness(QUBO.Wild.SNP.R0_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUBO.Wild.SNP.R0_AR)

# R80 ----
# Read in genind file: GSNAP4 alignment with Quercus robur genome; R80, NOMAF, first SNP/locus, only wild populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80_NOMAF_1SNP/"
setwd(genpop.filePath)
QUBO.Wild.SNP.R80_genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.Wild.SNP.R80_genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# Heterozygosity
QUBO.Wild.SNP.R80_HZ <- Hs(QUBO.Wild.SNP.R80_genind); print(QUBO.Wild.SNP.R80_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.Wild.SNP.R80_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf", 
                  "Kessler"), 
        main = "QUBO Wild Heterozygosity: SNPs (Complete), R80", 
        xlab = "Wild source", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUBO.Wild.SNP.R80_genind); ncol(QUBO.Wild.SNP.R80_genind@tab)
# Allelic richness: values per population
QUBO.Wild.SNP.R80_AR <- apply(allelic.richness(QUBO.Wild.SNP.R80_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUBO.Wild.SNP.R80_AR)

# ---- SUBSET ----
# MSAT: Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.MSAT_genind@tab), 
                                       function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.MSAT_genind@tab) <- QUBO.MSAT.sampleNames

# SNP: Remove QUBO_W_ headers from sample names
QUBO.SNP.sampleNames <- gsub("QUBO_G_",replacement = "", row.names(QUBO.SNP.R0_genind@tab))
QUBO.SNP.sampleNames <- gsub("QUBO_W_",replacement = "", QUBO.SNP.sampleNames)
# Replace SH-Q names in SNP list with IMLS names
# These were determined by Austin K., and are outlined on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
# Only 1 of the 11 SH_Q garden samples has an IMLS sample name (SHQ2177); others are unshared 
QUBO.SNP.sampleNames <- gsub("SH_Q2177",replacement = "IMLS338", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2178",replacement = "IMLS312", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2179",replacement = "IMLS062", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2180",replacement = "IMLS051", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2181",replacement = "IMLS011", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2182",replacement = "IMLS144", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2183",replacement = "IMLS170", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2184",replacement = "IMLS005", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2186",replacement = "IMLS017", QUBO.SNP.sampleNames)
# Rename sample matrices (both R0 and R80)
rownames(QUBO.SNP.R0_genind@tab) <- rownames(QUBO.SNP.R80_genind@tab) <- QUBO.SNP.sampleNames

# Subset SNP sample names by those that are also seen within the MSAT samples
QUBO_sharedSamples <- sort(QUBO.SNP.sampleNames[which(QUBO.SNP.sampleNames %in% QUBO.MSAT.sampleNames)])
# Subset MSAT and SNP wild matrix objects to strictly shared samples
QUBO.MSAT.subset_genind <- QUBO.MSAT_genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.R0.subset_genind <- QUBO.SNP.R0_genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.R80.subset_genind <- QUBO.SNP.R80_genind[QUBO_sharedSamples,, drop=TRUE]

# Heterozygosity
QUBO.MSAT.subset_HZ <- Hs(QUBO.MSAT.subset_genind); print(QUBO.MSAT.subset_HZ)
QUBO.SNP.R0.subset_HZ <- Hs(QUBO.SNP.R0.subset_genind); print(QUBO.SNP.R0.subset_HZ)
QUBO.SNP.R80.subset_HZ <- Hs(QUBO.SNP.R80.subset_genind); print(QUBO.SNP.R80.subset_HZ)
# Barplot for expected heterozygosity, MSAT subset
barplot(QUBO.MSAT.subset_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: MSATs (Subset)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Barplot for expected heterozygosity, SNP R0 subset
barplot(QUBO.SNP.R0.subset_HZ, beside = TRUE, 
        ylim = c(0,0.), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs (Subset), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Barplot for expected heterozygosity, SNP R80 subset
barplot(QUBO.SNP.R80.subset_HZ, beside = TRUE, 
        ylim = c(0,0.), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs (Subset), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
nLoc(QUBO.MSAT.subset_genind); ncol(QUBO.MSAT.subset_genind@tab)
nLoc(QUBO.SNP.R0.subset_genind); ncol(QUBO.SNP.R0.subset_genind@tab)
nLoc(QUBO.SNP.R80.subset_genind); ncol(QUBO.SNP.R80.subset_genind@tab)
# Allelic richness: values per population
QUBO.MSAT.subset_AR <- apply(allelic.richness(QUBO.MSAT.subset_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUBO.MSAT.subset_AR)
QUBO.SNP.R0.subset_AR <- apply(allelic.richness(QUBO.SNP.R0.subset_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUBO.SNP.R0.subset_AR)
QUBO.SNP.R80.subset_AR <- apply(allelic.richness(QUBO.SNP.R80.subset_genind)$Ar, 2, mean, na.rm=TRUE)
print(QUBO.SNP.R80.subset_AR)
