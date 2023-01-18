# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% POPULATION GENETIC STATISTICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script calculates the following population genetic statistics:
# heterozygosity, allele counts (i.e. number of variant sites across all loci), and allelic richness
# It does this for both QUAC (optimized de novo assembly) and QUBO 
# (aligned to the Q. robur reference genome) NextRAD datasets

# The script performs these calculations with different sets of samples and loci
# For each species, analyses are divided by using garden and wild samples and using only wild samples, as outlined below:
# GARDEN AND WILD SAMPLES
# 1. Microsatellites (Complete)
# 2. SNPs, De novo, R0 (Complete)
# 3. SNPs, De novo, R80 (Complete)
# 4. SNPs, De novo, R80 (Complete)
# 5. SNPs, Reference, R0, 1st SNP/locus (Complete)
# 6. SNPs, Reference, R80, microhaplotype-wise SNPs (Complete)
# 7. Microsatellites (Subset)
# 8. SNPs, De novo, R0 (Subset)
# 9. SNPs, De novo, R80 (Subset)
# 10. SNPs, Reference, R0 (Subset)
# 11. SNPs, Reference, R80 (Subset)
# WILD SAMPLES
# 1. QUAC: SNPs, De novo, R0 (Complete)
# 2. QUAC: SNPs, De novo, R80 (Complete)
# 3. QUBO: SNPs, Reference, R0 (Complete)
# 4. QUBO: SNPs, Reference, R80 (Complete)
# Each of the garden and wild scenarios is explored for both species. 
# However, for QUBO, microhaplotype SNPs are analyzed for de novo datasets, rather than reference datasets
# Complete: all samples originally included in analysis, for specified marker type
# Subset: only samples shared between MSAT and SNP analyses

library(adegenet)
library(hierfstat)
# Set working directory
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)

# %%%% QUAC %%%% ----
# %%% GARDEN AND WILD ----
# ---- MSATS ----
# Read in genind file (GCC_QUAC_ZAIN repo; QUAC_wK_garden_wild_clean.gen)
QUAC.MSAT_filePath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/"
QUAC.MSAT_genind <- read.genepop(paste0(QUAC.MSAT_filePath, "QUAC_wK_allpop_clean.gen"), ncode = 3)
# Correct popNames: samples with popName pattern QAc-G- are garden 
levels(QUAC.MSAT_genind@pop)[grep(pattern = "QAc-G-", levels(QUAC.MSAT_genind@pop))] <- 
  rep("garden", length(grep(pattern = "QAc-G-", levels(QUAC.MSAT_genind@pop))))
# Correct popNames: samples with popName pattern QAc-W- are wild
levels(QUAC.MSAT_genind@pop)[grep(pattern = "QAc-W-", levels(QUAC.MSAT_genind@pop))] <- 
  rep("wild", length(grep(pattern = "QAc-W-", levels(QUAC.MSAT_genind@pop))))

# Heterozygosity
QUAC.MSAT_HZ <- Hs(QUAC.MSAT_genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.MSAT_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: MSATs (Complete)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUAC.MSAT_LC <- nLoc(QUAC.MSAT_genind); QUAC.MSAT_AC <- ncol(QUAC.MSAT_genind@tab)
# Allelic richness: values per population
QUAC.MSAT_AR <- apply(allelic.richness(QUAC.MSAT_genind)$Ar, 2, mean, na.rm=TRUE)
# Print pop. gen. results
print("%%% QUAC MSAT, GARDEN AND WILD (COMPLETE) %%%")
print(paste0("Heterozygosity: ", QUAC.MSAT_HZ))
print(paste0("Loci count: ", QUAC.MSAT_LC))
print(paste0("Allele count: ", QUAC.MSAT_AC))
print(paste0("Allelic richness: ", QUAC.MSAT_AR))

# ---- SNPS ----
# DE NOVO ----
# R0 ----
# Read in genind file: Optimized de novo assembly; R0, NOMAF, first SNP/locus, 2 populations
QUAC.SNP.DN.R0_filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops/"
QUAC.SNP.DN.R0_genind <- read.genepop(paste0(QUAC.SNP.DN.R0_filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R0_genind) <- 
  factor(read.table(paste0(QUAC.SNP.DN.R0_filePath,"QUAC_popmap_GardenWild"), header=FALSE)[,2])

# Heterozygosity
QUAC.SNP.DN.R0_HZ <- Hs(QUAC.SNP.DN.R0_genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP.DN.R0_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: SNPs, De novo (Complete), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUAC.SNP.DN.R0_LC <- nLoc(QUAC.SNP.DN.R0_genind); QUAC.SNP.DN.R0_AC <- ncol(QUAC.SNP.DN.R0_genind@tab)
# Allelic richness: values per population
QUAC.SNP.DN.R0_AR <- apply(allelic.richness(QUAC.SNP.DN.R0_genind)$Ar, 2, mean, na.rm=TRUE)
# Print pop. gen. results
print("%%% QUAC SNP (DE NOVO, R0), GARDEN AND WILD (COMPLETE) %%%")
print(paste0("Heterozygosity: ", QUAC.SNP.DN.R0_HZ))
print(paste0("Loci count: ", QUAC.SNP.DN.R0_LC))
print(paste0("Allele count: ", QUAC.SNP.DN.R0_AC))
print(paste0("Allelic richness: ", QUAC.SNP.DN.R0_AR))

# R80 ----
# 1st SNP/locus ----
# Read in genind file: Optimized de novo assembly; R80, NOMAF, first SNP/locus, 2 populations
QUAC.SNP.DN.R80_filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops/"
QUAC.SNP.DN.R80_genind <- read.genepop(paste0(QUAC.SNP.DN.R80_filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R80_genind) <- factor(read.table(paste0(QUAC.SNP.DN.R80_filePath,"QUAC_popmap_GardenWild"), header=FALSE)[,2])

# Heterozygosity
QUAC.SNP.DN.R80_HZ <- Hs(QUAC.SNP.DN.R80_genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP.DN.R80_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: SNPs, De novo (Complete), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUAC.SNP.DN.R80_LC <- nLoc(QUAC.SNP.DN.R80_genind); QUAC.SNP.DN.R80_AC <- ncol(QUAC.SNP.DN.R80_genind@tab)
# Allelic richness: values per population
QUAC.SNP.DN.R80_AR <- apply(allelic.richness(QUAC.SNP.DN.R80_genind)$Ar, 2, mean, na.rm=TRUE)
# Print pop. gen. results
print("%%% QUAC SNP (DE NOVO, R80), GARDEN AND WILD (COMPLETE) %%%")
print(paste0("Heterozygosity: ", QUAC.SNP.DN.R80_HZ))
print(paste0("Loci count: ", QUAC.SNP.DN.R80_LC))
print(paste0("Allele count: ", QUAC.SNP.DN.R80_AC))
print(paste0("Allelic richness: ", QUAC.SNP.DN.R80_AR))

# # Microhaplotypes ----
# # Read in genind file: Optimized de novo assembly; R80, NOMAF, haplotype-wise SNPs, 2 populations
# QUAC.SNP.DN.R80.HapSNP_filePath <-
#   "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_HapSNPs_2Pops/"
# QUAC.SNP.DN.R80.HapSNP_genind <- read.genepop(paste0(QUAC.SNP.DN.R80.HapSNP_filePath,"populations.snps.gen"))
# # Correct popNames
# pop(QUAC.SNP.DN.R80.HapSNP_genind) <-
#   factor(read.table(paste0(QUAC.SNP.DN.R80.HapSNP_filePath,"QUAC_popmap_GardenWild"), header=FALSE)[,2])
# 
# # Heterozygosity
# QUAC.SNP.DN.R80.HapSNP_HZ <- Hs(QUAC.SNP.DN.R80.HapSNP_genind)
# # Barplot for expected heterozygosity, SNP markers
# barplot(QUAC.SNP.DN.R80.HapSNP_HZ, beside = TRUE,
#         ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
#         names = c("Garden", "Wild"),
#         main = "QUAC Garden and Wild Heterozygosity: SNPs, De novo, (Complete), R80 (Haplotype-wise SNPs)",
#         xlab = "Population Type", ylab = "Expected Heterozygosity")
# abline(h = 0, lwd=2)
# 
# # Loci count, allele count
# QUAC.SNP.DN.R80.HapSNP_LC <- nLoc(QUAC.SNP.DN.R80.HapSNP_genind) 
# QUAC.SNP.DN.R80.HapSNP_AC <- ncol(QUAC.SNP.DN.R80.HapSNP_genind@tab)
# # Allelic richness: values per population
# QUAC.SNP.DN.R80.HapSNP_AR <- apply(allelic.richness(QUAC.SNP.DN.R80.HapSNP_genind)$Ar, 2, mean, na.rm=TRUE)
# # Print pop. gen. results
# print("%%% QUAC SNP (DE NOVO, HAPLOTYPE-WISE), GARDEN AND WILD (COMPLETE) %%%")
# print(paste0("Heterozygosity: ", QUAC.SNP.DN.R80.HapSNP_HZ))
# print(paste0("Loci count: ", QUAC.SNP.DN.R80.HapSNP_LC))
# print(paste0("Allele count: ", QUAC.SNP.DN.R80.HapSNP_AC))
# print(paste0("Allelic richness: ", QUAC.SNP.DN.R80.HapSNP_AR))

# REFERENCE ----
# R0 ----
# Read in genind file: GSNAP4 alignment with Quercus rubra genome; R0, NOMAF, first SNP/locus, 2 populations
QUAC.SNP.REF.R0_filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R0_NOMAF_1SNP_2Pops/"
QUAC.SNP.REF.R0_genind <- read.genepop(paste0(QUAC.SNP.REF.R0_filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R0_genind) <- 
  factor(read.table(paste0(QUAC.SNP.REF.R0_filePath,"QUAC_popmap_GardenWild"), header=FALSE)[,2])

# Heterozygosity
QUAC.SNP.REF.R0_HZ <- Hs(QUAC.SNP.REF.R0_genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP.REF.R0_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: SNPs, Reference (Complete), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUAC.SNP.REF.R0_LC <- nLoc(QUAC.SNP.REF.R0_genind); QUAC.SNP.REF.R0_AC <- ncol(QUAC.SNP.REF.R0_genind@tab)
# Allelic richness: values per population
QUAC.SNP.REF.R0_AR <- apply(allelic.richness(QUAC.SNP.REF.R0_genind)$Ar, 2, mean, na.rm=TRUE)
# Print pop. gen. results
print("%%% QUAC SNP (REFERENCE, R0), GARDEN AND WILD (COMPLETE) %%%")
print(paste0("Heterozygosity: ", QUAC.SNP.REF.R0_HZ))
print(paste0("Loci count: ", QUAC.SNP.REF.R0_LC))
print(paste0("Allele count: ", QUAC.SNP.REF.R0_AC))
print(paste0("Allelic richness: ", QUAC.SNP.REF.R0_AR))

# R80 ----
# Read in genind file: GSNAP4 alignment with Quercus rubra genome; R80, NOMAF, first SNP/locus, 2 populations
QUAC.SNP.REF.R80_filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R80_NOMAF_1SNP_2Pops/"
QUAC.SNP.REF.R80_genind <- read.genepop(paste0(QUAC.SNP.REF.R80_filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R80_genind) <- factor(read.table(paste0(QUAC.SNP.REF.R80_filePath,"QUAC_popmap_GardenWild"), header=FALSE)[,2])

# Heterozygosity
QUAC.SNP.REF.R80_HZ <- Hs(QUAC.SNP.REF.R80_genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP.REF.R80_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: SNPs, Reference (Complete), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUAC.SNP.REF.R80_LC <- nLoc(QUAC.SNP.REF.R80_genind); QUAC.SNP.REF.R80_AC <- ncol(QUAC.SNP.REF.R80_genind@tab)
# Allelic richness: values per population
QUAC.SNP.REF.R80_AR <- apply(allelic.richness(QUAC.SNP.REF.R80_genind)$Ar, 2, mean, na.rm=TRUE)
# Print pop. gen. results
print("%%% QUAC SNP (REFERENCE, R80), GARDEN AND WILD (COMPLETE) %%%")
print(paste0("Heterozygosity: ", QUAC.SNP.REF.R80_HZ))
print(paste0("Loci count: ", QUAC.SNP.REF.R80_LC))
print(paste0("Allele count: ", QUAC.SNP.REF.R80_AC))
print(paste0("Allelic richness: ", QUAC.SNP.REF.R80_AR))

# ---- SUBSET ----
# MSAT: read in Tissue database names from GCC_QUAC_ZAIN repository, and rename MSAT genind matrix
QUAC.MSAT.tissueNames_filepath <- "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Data_Frames/QUAC_allpop_clean_df.csv"
QUAC.MSAT.tissueNames <- unlist(read.csv2(QUAC.MSAT.tissueNames_filepath, header = TRUE, sep=",")[1])
rownames(QUAC.MSAT_genind@tab) <- QUAC.MSAT.tissueNames

# SNP: read in Tissue database names, and rename SNP genind matrices (both R0 and R80)
# This file was created (by Austin K.), and can be found on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
QUAC.SNP.tissueNames_filepath <- paste0(SSRvSNP.wd,"exSituRepresentation/Resampling/QUAC_SNP_TissueNames.csv")
QUAC.SNP.tissueNames <- unlist(read.csv2(QUAC.SNP.tissueNames_filepath, header = TRUE, sep = ",")[3])
rownames(QUAC.SNP.DN.R0_genind@tab) <- rownames(QUAC.SNP.DN.R80_genind@tab) <- QUAC.SNP.tissueNames
rownames(QUAC.SNP.REF.R0_genind@tab) <- rownames(QUAC.SNP.REF.R80_genind@tab) <- QUAC.SNP.tissueNames

# Subset SNP sample names by those that are also seen within the MSAT samples 
QUAC_sharedSamples <- sort(QUAC.SNP.tissueNames[which(QUAC.SNP.tissueNames %in% QUAC.MSAT.tissueNames)])
# Subset MSAT and SNP genind matrices to strictly shared samples, dropping now absent alleles
QUAC.MSAT.subset_genind <- QUAC.MSAT_genind[QUAC_sharedSamples,, drop=TRUE]
# De novo
QUAC.SNP.DN.R0.subset_genind <- QUAC.SNP.DN.R0_genind[QUAC_sharedSamples,, drop=TRUE]
QUAC.SNP.DN.R80.subset_genind <- QUAC.SNP.DN.R80_genind[QUAC_sharedSamples,, drop=TRUE]
# Reference
QUAC.SNP.REF.R0.subset_genind <- QUAC.SNP.REF.R0_genind[QUAC_sharedSamples,, drop=TRUE]
QUAC.SNP.REF.R80.subset_genind <- QUAC.SNP.REF.R80_genind[QUAC_sharedSamples,, drop=TRUE]

# Heterozygosity
QUAC.MSAT.subset_HZ <- Hs(QUAC.MSAT.subset_genind)
# De novo
QUAC.SNP.DN.R0.subset_HZ <- Hs(QUAC.SNP.DN.R0.subset_genind)
QUAC.SNP.DN.R80.subset_HZ <- Hs(QUAC.SNP.DN.R80.subset_genind)
# Reference
QUAC.SNP.REF.R0.subset_HZ <- Hs(QUAC.SNP.REF.R0.subset_genind)
QUAC.SNP.REF.R80.subset_HZ <- Hs(QUAC.SNP.REF.R80.subset_genind)
# Barplot for expected heterozygosity, MSAT subset
barplot(QUAC.MSAT.subset_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: MSATs (Subset)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# De novo
# Barplot for expected heterozygosity, SNP R0 subset
barplot(QUAC.SNP.DN.R0.subset_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: SNPs, De novo (Subset), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Barplot for expected heterozygosity, SNP R80 subset
barplot(QUAC.SNP.DN.R80.subset_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: SNPs, De novo (Subset), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Reference
# Barplot for expected heterozygosity, SNP R0 subset
barplot(QUAC.SNP.REF.R0.subset_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: SNPs, Reference (Subset), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Barplot for expected heterozygosity, SNP R80 subset
barplot(QUAC.SNP.REF.R80.subset_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Garden and Wild Heterozygosity: SNPs, Reference (Subset), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUAC.MSAT.subset_LC <- nLoc(QUAC.MSAT.subset_genind) 
QUAC.MSAT.subset_AC <- ncol(QUAC.MSAT.subset_genind@tab)
# De novo
QUAC.SNP.DN.R0.subset_LC <- nLoc(QUAC.SNP.DN.R0.subset_genind) 
QUAC.SNP.DN.R0.subset_AC <- ncol(QUAC.SNP.DN.R0.subset_genind@tab)
QUAC.SNP.DN.R80.subset_LC <- nLoc(QUAC.SNP.DN.R80.subset_genind) 
QUAC.SNP.DN.R80.subset_AC <- ncol(QUAC.SNP.DN.R80.subset_genind@tab)
# Reference
QUAC.SNP.REF.R0.subset_LC <- nLoc(QUAC.SNP.REF.R0.subset_genind) 
QUAC.SNP.REF.R0.subset_AC <- ncol(QUAC.SNP.REF.R0.subset_genind@tab)
QUAC.SNP.REF.R80.subset_LC <- nLoc(QUAC.SNP.REF.R80.subset_genind) 
QUAC.SNP.REF.R80.subset_AC <- ncol(QUAC.SNP.REF.R80.subset_genind@tab)

# Allelic richness: values per population
QUAC.MSAT.subset_AR <- apply(allelic.richness(QUAC.MSAT.subset_genind)$Ar, 2, mean, na.rm=TRUE)
# De novo
QUAC.SNP.DN.R0.subset_AR <- apply(allelic.richness(QUAC.SNP.DN.R0.subset_genind)$Ar, 2, mean, na.rm=TRUE)
QUAC.SNP.DN.R80.subset_AR <- apply(allelic.richness(QUAC.SNP.DN.R80.subset_genind)$Ar, 2, mean, na.rm=TRUE)
# Reference
QUAC.SNP.REF.R0.subset_AR <- apply(allelic.richness(QUAC.SNP.REF.R0.subset_genind)$Ar, 2, mean, na.rm=TRUE)
QUAC.SNP.REF.R80.subset_AR <- apply(allelic.richness(QUAC.SNP.REF.R80.subset_genind)$Ar, 2, mean, na.rm=TRUE)

# Print pop. gen. results
# MSAT
print("%%% QUAC MSAT, GARDEN AND WILD (SUBSET) %%%")
print(paste0("Heterozygosity: ", QUAC.MSAT.subset_HZ))
print(paste0("Loci count: ", QUAC.MSAT.subset_LC))
print(paste0("Allele count: ", QUAC.MSAT.subset_AC))
print(paste0("Allelic richness: ", QUAC.MSAT.subset_AR))
# SNP, De novo
# R0
print("%%% QUAC SNP (DE NOVO, R0), GARDEN AND WILD (SUBSET) %%%")
print(paste0("Heterozygosity: ", QUAC.SNP.DN.R0.subset_HZ))
print(paste0("Loci count: ", QUAC.SNP.DN.R0.subset_LC))
print(paste0("Allele count: ", QUAC.SNP.DN.R0.subset_AC))
print(paste0("Allelic richness: ", QUAC.SNP.DN.R0.subset_AR))
# R80
print("%%% QUAC SNP (DE NOVO, R80), GARDEN AND WILD (SUBSET) %%%")
print(paste0("Heterozygosity: ", QUAC.SNP.DN.R80.subset_HZ))
print(paste0("Loci count: ", QUAC.SNP.DN.R80.subset_LC))
print(paste0("Allele count: ", QUAC.SNP.DN.R80.subset_AC))
print(paste0("Allelic richness: ", QUAC.SNP.DN.R80.subset_AR))
# SNP, Reference
# R0
print("%%% QUAC SNP (REFERENCE, R0), GARDEN AND WILD (SUBSET) %%%")
print(paste0("Heterozygosity: ", QUAC.SNP.REF.R0.subset_HZ))
print(paste0("Loci count: ", QUAC.SNP.REF.R0.subset_LC))
print(paste0("Allele count: ", QUAC.SNP.REF.R0.subset_AC))
print(paste0("Allelic richness: ", QUAC.SNP.REF.R0.subset_AR))
# R80
print("%%% QUAC SNP (REFERENCE, R80), GARDEN AND WILD (SUBSET) %%%")
print(paste0("Heterozygosity: ", QUAC.SNP.REF.R80.subset_HZ))
print(paste0("Loci count: ", QUAC.SNP.REF.R80.subset_LC))
print(paste0("Allele count: ", QUAC.SNP.REF.R80.subset_AC))
print(paste0("Allelic richness: ", QUAC.SNP.REF.R80.subset_AR))

# %%% WILD ONLY ----
# SNP, DE NOVO ----
# R0 ----
# # Read in genind file: Optimized de novo assembly; R0, NOMAF, first SNP/locus, only wild populations
# QUAC.Wild.SNP.DN.R0_filePath <-
#   "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF_1SNP/"
# QUAC.Wild.SNP.DN.R0_genind <- read.genepop(paste0(QUAC.Wild.SNP.DN.R0_filePath,"populations.snps.gen"))
# # Correct popNames
# pop(QUAC.Wild.SNP.DN.R0_genind) <- factor(read.table(paste0(QUAC.Wild.SNP.DN.R0_filePath,"QUAC_popmap_wild"), header=FALSE)[,2])
# 
# # Heterozygosity
# QUAC.Wild.SNP.DN.R0_HZ <- Hs(QUAC.Wild.SNP.DN.R0_genind)
# # Barplot for expected heterozygosity, SNP markers
# barplot(QUAC.Wild.SNP.DN.R0_HZ, beside = TRUE,
#         ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
#         names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf",
#                   "Kessler"),
#         main = "QUAC Wild Heterozygosity: SNPs, De novo (Complete), R0",
#         xlab = "Wild population", ylab = "Expected Heterozygosity")
# abline(h = 0, lwd=2)
# 
# # Loci count, allele count
# QUAC.Wild.SNP.DN.R0_LC <- nLoc(QUAC.Wild.SNP.DN.R0_genind)
# QUAC.Wild.SNP.DN.R0_AC <- ncol(QUAC.Wild.SNP.DN.R0_genind@tab)
# # Allelic richness: values per population
# QUAC.Wild.SNP.DN.R0_AR <- apply(allelic.richness(QUAC.Wild.SNP.DN.R0_genind)$Ar, 2, mean, na.rm=TRUE)
# # Print pop. gen. results
# print("%%% QUAC SNP (DE NOVO, R0), WILD (COMPLETE) %%%")
# print(paste0("Heterozygosity: ", QUAC.Wild.SNP.DN.R0_HZ))
# print(paste0("Loci count: ", QUAC.Wild.SNP.DN.R0_LC))
# print(paste0("Allele count: ", QUAC.Wild.SNP.DN.R0_AC))
# print(paste0("Allelic richness: ", QUAC.Wild.SNP.DN.R0_AR))
# 
# # R80 ----
# # Read in genind file: Optimized de novo assembly; R80, NOMAF, first SNP/locus, only wild populations
# QUAC.Wild.SNP.DN.R80_filePath <-
#   "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R80_NOMAF_1SNP/"
# QUAC.Wild.SNP.DN.R80_genind <- read.genepop(paste0(QUAC.Wild.SNP.DN.R80_filePath,"populations.snps.gen"))
# # Correct popNames
# pop(QUAC.Wild.SNP.DN.R80_genind) <- factor(read.table(paste0(QUAC.Wild.SNP.DN.R80_filePath,"QUAC_popmap_wild"), header=FALSE)[,2])
# 
# # Heterozygosity
# QUAC.Wild.SNP.DN.R80_HZ <- Hs(QUAC.Wild.SNP.DN.R80_genind)
# # Barplot for expected heterozygosity, SNP markers
# barplot(QUAC.Wild.SNP.DN.R80_HZ, beside = TRUE,
#         ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
#         names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf",
#                   "Kessler"),
#         main = "QUAC Wild Heterozygosity: SNPs, De novo (Complete), R80",
#         xlab = "Wild population", ylab = "Expected Heterozygosity")
# abline(h = 0, lwd=2)
# 
# # Loci count, allele count
# QUAC.Wild.SNP.DN.R80_LC <- nLoc(QUAC.Wild.SNP.DN.R80_genind)
# QUAC.Wild.SNP.DN.R80_AC <- ncol(QUAC.Wild.SNP.DN.R80_genind@tab)
# # Allelic richness: values per population
# QUAC.Wild.SNP.DN.R80_AR <- apply(allelic.richness(QUAC.Wild.SNP.DN.R80_genind)$Ar, 2, mean, na.rm=TRUE)
# # Print pop. gen. results
# print("%%% QUAC SNP (DE NOVO, R80), WILD (COMPLETE) %%%")
# print(paste0("Heterozygosity: ", QUAC.Wild.SNP.DN.R80_HZ))
# print(paste0("Loci count: ", QUAC.Wild.SNP.DN.R80_LC))
# print(paste0("Allele count: ", QUAC.Wild.SNP.DN.R80_AC))
# print(paste0("Allelic richness: ", QUAC.Wild.SNP.DN.R80_AR))

# %%%% QUBO %%%% ----
# %%% GARDEN AND WILD ----
# ---- MSATS ----
# Read in genind file (Southeast Oaks repo; genetic_data/Qb_total.gen file)
QUBO.MSAT_filePath <- 
  "~/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
QUBO.MSAT_genind <- read.genepop(paste0(QUBO.MSAT_filePath,"Qb_total.gen"), ncode=3)
# Correct popNames: last population (IMLS4_MP1_IMLS336_C05) is garden; rest (9) are wild
levels(QUBO.MSAT_genind@pop) <- c(rep("wild",9), "garden") 

# Heterozygosity
QUBO.MSAT_HZ <- Hs(QUBO.MSAT_genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.MSAT_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: MSATs (Complete)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUBO.MSAT_LC <- nLoc(QUBO.MSAT_genind); QUBO.MSAT_AC <- ncol(QUBO.MSAT_genind@tab)
# Allelic richness: values per population
QUBO.MSAT_AR <- apply(allelic.richness(QUBO.MSAT_genind)$Ar, 2, mean, na.rm=TRUE)
# Print pop. gen. results
print("%%% QUBO MSAT, GARDEN AND WILD (COMPLETE) %%%")
print(paste0("Heterozygosity: ", QUBO.MSAT_HZ))
print(paste0("Loci count: ", QUBO.MSAT_LC))
print(paste0("Allele count: ", QUBO.MSAT_AC))
print(paste0("Allelic richness: ", QUBO.MSAT_AR))

# ---- SNPS ----
# DE NOVO ----
# R0 ----
# Read in genind file: optimized de novo assembly; R0, NOMAF, first SNP/locus, 2 populations 
QUBO.SNP.DN.R0_filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R0_NOMAF_1SNP_2Pops/"
QUBO.SNP.DN.R0_genind <- read.genepop(paste0(QUBO.SNP.DN.R0_filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R0_genind) <- 
  factor(read.table(paste0(QUBO.SNP.DN.R0_filePath, "QUBO_popmap_GardenWild"), header=FALSE)[,2])

# Heterozygosity
QUBO.SNP.DN.R0_HZ <- Hs(QUBO.SNP.DN.R0_genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.SNP.DN.R0_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs, De novo (Complete), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUBO.SNP.DN.R0_LC <- nLoc(QUBO.SNP.DN.R0_genind); QUBO.SNP.DN.R0_AC <- ncol(QUBO.SNP.DN.R0_genind@tab)
# Allelic richness: values per population
QUBO.SNP.DN.R0_AR <- apply(allelic.richness(QUBO.SNP.DN.R0_genind)$Ar, 2, mean, na.rm=TRUE)
# Print pop. gen. results
print("%%% QUBO SNP (DE NOVO, R0), GARDEN AND WILD (COMPLETE) %%%")
print(paste0("Heterozygosity: ", QUBO.SNP.DN.R0_HZ))
print(paste0("Loci count: ", QUBO.SNP.DN.R0_LC))
print(paste0("Allele count: ", QUBO.SNP.DN.R0_AC))
print(paste0("Allelic richness: ", QUBO.SNP.DN.R0_AR))

# R80 ----
# Read in genind file: optimized de novo assembly; R80, NOMAF, first SNP/locus, 2 populations 
QUBO.SNP.DN.R80_filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R80_NOMAF_1SNP_2Pops/"
QUBO.SNP.DN.R80_genind <- read.genepop(paste0(QUBO.SNP.DN.R80_filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R80_genind) <- 
  factor(read.table(paste0(QUBO.SNP.DN.R80_filePath,"QUBO_popmap_GardenWild"), header=FALSE)[,2])

# Heterozygosity
QUBO.SNP.DN.R80_HZ <- Hs(QUBO.SNP.DN.R80_genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.SNP.DN.R80_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs, De novo (Complete), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUBO.SNP.DN.R80_LC <- nLoc(QUBO.SNP.DN.R80_genind); QUBO.SNP.DN.R80_AC <- ncol(QUBO.SNP.DN.R80_genind@tab)
# Allelic richness: values per population
QUBO.SNP.DN.R80_AR <- apply(allelic.richness(QUBO.SNP.DN.R80_genind)$Ar, 2, mean, na.rm=TRUE)
# Print pop. gen. results
print("%%% QUBO SNP (DE NOVO, R80), GARDEN AND WILD (COMPLETE) %%%")
print(paste0("Heterozygosity: ", QUBO.SNP.DN.R80_HZ))
print(paste0("Loci count: ", QUBO.SNP.DN.R80_LC))
print(paste0("Allele count: ", QUBO.SNP.DN.R80_AC))
print(paste0("Allelic richness: ", QUBO.SNP.DN.R80_AR))

# REFERENCE ----
# R0 ----
# Read in genind file: GSNAP4 alignment with Quercus robur genome; R0, NOMAF, first SNP/locus, 2 populations 
QUBO.SNP.REF.R0_filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
QUBO.SNP.REF.R0_genind <- read.genepop(paste0(QUBO.SNP.REF.R0_filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R0_genind) <- 
  factor(read.table(paste0(QUBO.SNP.REF.R0_filePath,"QUBO_popmap_GardenWild"), header=FALSE)[,2])

# Heterozygosity
QUBO.SNP.REF.R0_HZ <- Hs(QUBO.SNP.REF.R0_genind)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.SNP.REF.R0_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs, Reference (Complete), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUBO.SNP.REF.R0_LC <- nLoc(QUBO.SNP.REF.R0_genind); QUBO.SNP.REF.R0_AC <- ncol(QUBO.SNP.REF.R0_genind@tab)
# Allelic richness: values per population
QUBO.SNP.REF.R0_AR <- apply(allelic.richness(QUBO.SNP.REF.R0_genind)$Ar, 2, mean, na.rm=TRUE)
# Print pop. gen. results
print("%%% QUBO SNP (REFERENCE, R0), GARDEN AND WILD (COMPLETE) %%%")
print(paste0("Heterozygosity: ", QUBO.SNP.REF.R0_HZ))
print(paste0("Loci count: ", QUBO.SNP.REF.R0_LC))
print(paste0("Allele count: ", QUBO.SNP.REF.R0_AC))
print(paste0("Allelic richness: ", QUBO.SNP.REF.R0_AR))

# R80 ----
# 1st SNP/locus ----
# Read in genind file: GSNAP4 alignment with Quercus robur genome; R80, NOMAF, first SNP/locus, 2 populations 
QUBO.SNP.REF.R80_filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_1SNP_2Pops/"
QUBO.SNP.REF.R80_genind <- read.genepop(paste0(QUBO.SNP.REF.R80_filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R80_genind) <- 
  factor(read.table(paste0(QUBO.SNP.REF.R80_filePath,"QUBO_popmap_GardenWild"), header=FALSE)[,2])

# Heterozygosity
QUBO.SNP.REF.R80_HZ <- Hs(QUBO.SNP.REF.R80_genind); print(QUBO.SNP.REF.R80_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.SNP.REF.R80_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs, Reference (Complete), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUBO.SNP.REF.R80_LC <- nLoc(QUBO.SNP.REF.R80_genind); QUBO.SNP.REF.R80_AC <- ncol(QUBO.SNP.REF.R80_genind@tab)
# Allelic richness: values per population
QUBO.SNP.REF.R80_AR <- apply(allelic.richness(QUBO.SNP.REF.R80_genind)$Ar, 2, mean, na.rm=TRUE)
# Print pop. gen. results
print("%%% QUBO SNP (REFERENCE, R80), GARDEN AND WILD (COMPLETE) %%%")
print(paste0("Heterozygosity: ", QUBO.SNP.REF.R80_HZ))
print(paste0("Loci count: ", QUBO.SNP.REF.R80_LC))
print(paste0("Allele count: ", QUBO.SNP.REF.R80_AC))
print(paste0("Allelic richness: ", QUBO.SNP.REF.R80_AR))

# # Microhaplotypes ----
# # Read in genind file: GSNAP4 alignment with Quercus robur genome; R80, NOMAF, haplotype-wise SNPs, 2 populations
# QUBO.SNP.REF.R80.HapSNP_filePath <-
#   "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_HapSNPs_2Pops/"
# QUBO.SNP.REF.R80.HapSNP_genind <- read.genepop(paste0(QUBO.SNP.REF.R80.HapSNP_filePath,"populations.snps.gen"))
# # Correct popNames
# pop(QUBO.SNP.REF.R80.HapSNP_genind) <-
# factor(read.table(paste0(QUBO.SNP.REF.R80.HapSNP_filePath,"QUBO_popmap_GardenWild"), header=FALSE)[,2])
# 
# # Heterozygosity
# QUBO.SNP.REF.R80.HapSNP_HZ <- Hs(QUBO.SNP.REF.R80.HapSNP_genind)
# # Barplot for expected heterozygosity, SNP markers
# barplot(QUBO.SNP.REF.R80.HapSNP_HZ, beside = TRUE,
#         ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
#         names = c("Garden", "Wild"),
#         main = "QUBO Garden and Wild Heterozygosity: SNPs, Reference (Complete), R80 (Haplotype-wise SNPs)",
#         xlab = "Population Type", ylab = "Expected Heterozygosity")
# abline(h = 0, lwd=2)
# 
# # Loci count, allele count
# QUBO.SNP.REF.R80.HapSNP_LC <- nLoc(QUBO.SNP.REF.R80.HapSNP_genind) 
# QUBO.SNP.REF.R80.HapSNP_AC <- ncol(QUBO.SNP.REF.R80.HapSNP_genind@tab)
# # Allelic richness: values per population
# QUBO.SNP.REF.R80.HapSNP_AR <- apply(allelic.richness(QUBO.SNP.REF.R80.HapSNP_genind)$Ar, 2, mean, na.rm=TRUE)
# # Print pop. gen. results
# print("%%% QUBO SNP (REFERENCE, R80), GARDEN AND WILD (COMPLETE) %%%")
# print(paste0("Heterozygosity: ", QUBO.SNP.REF.R80.HapSNP_HZ))
# print(paste0("Loci count: ", QUBO.SNP.REF.R80.HapSNP_LC))
# print(paste0("Allele count: ", QUBO.SNP.REF.R80.HapSNP_AC))
# print(paste0("Allelic richness: ", QUBO.SNP.REF.R80.HapSNP_AR))

# ---- SUBSET ----
# MSAT: Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.MSAT_genind@tab), 
                                       function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.MSAT_genind@tab) <- QUBO.MSAT.sampleNames

# SNP: Remove QUBO_W_ headers from sample names
QUBO.SNP.sampleNames <- gsub("QUBO_G_",replacement = "", row.names(QUBO.SNP.DN.R0_genind@tab))
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
rownames(QUBO.SNP.DN.R0_genind@tab) <- rownames(QUBO.SNP.DN.R80_genind@tab) <- QUBO.SNP.sampleNames
rownames(QUBO.SNP.REF.R0_genind@tab) <- rownames(QUBO.SNP.REF.R80_genind@tab) <- QUBO.SNP.sampleNames

# Subset SNP sample names by those that are also seen within the MSAT samples
QUBO_sharedSamples <- sort(QUBO.SNP.sampleNames[which(QUBO.SNP.sampleNames %in% QUBO.MSAT.sampleNames)])
# Subset MSAT and SNP wild matrix objects to strictly shared samples
QUBO.MSAT.subset_genind <- QUBO.MSAT_genind[QUBO_sharedSamples,, drop=TRUE]
# De novo
QUBO.SNP.DN.R0.subset_genind <- QUBO.SNP.DN.R0_genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.DN.R80.subset_genind <- QUBO.SNP.DN.R80_genind[QUBO_sharedSamples,, drop=TRUE]
# Reference
QUBO.SNP.REF.R0.subset_genind <- QUBO.SNP.REF.R0_genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.REF.R80.subset_genind <- QUBO.SNP.REF.R80_genind[QUBO_sharedSamples,, drop=TRUE]

# Heterozygosity
QUBO.MSAT.subset_HZ <- Hs(QUBO.MSAT.subset_genind)
# De novo
QUBO.SNP.DN.R0.subset_HZ <- Hs(QUBO.SNP.DN.R0.subset_genind)
QUBO.SNP.DN.R80.subset_HZ <- Hs(QUBO.SNP.DN.R80.subset_genind)
# Reference
QUBO.SNP.REF.R0.subset_HZ <- Hs(QUBO.SNP.REF.R0.subset_genind)
QUBO.SNP.REF.R80.subset_HZ <- Hs(QUBO.SNP.REF.R80.subset_genind)
# Barplot for expected heterozygosity, MSAT subset
barplot(QUBO.MSAT.subset_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: MSATs (Subset)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# De novo
# Barplot for expected heterozygosity, SNP R0 subset
barplot(QUBO.SNP.DN.R0.subset_HZ, beside = TRUE, 
        ylim = c(0,0.), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs, De novo (Subset), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Barplot for expected heterozygosity, SNP R80 subset
barplot(QUBO.SNP.DN.R80.subset_HZ, beside = TRUE, 
        ylim = c(0,0.), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs, De novo (Subset), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Reference
# Barplot for expected heterozygosity, SNP R0 subset
barplot(QUBO.SNP.REF.R0.subset_HZ, beside = TRUE, 
        ylim = c(0,0.), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs, Reference (Subset), R0", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Barplot for expected heterozygosity, SNP R80 subset
barplot(QUBO.SNP.REF.R80.subset_HZ, beside = TRUE, 
        ylim = c(0,0.), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Garden and Wild Heterozygosity: SNPs, Reference (Subset), R80", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Loci count, allele count
QUBO.MSAT.subset_LC <- nLoc(QUBO.MSAT.subset_genind) 
QUBO.MSAT.subset_AC <- ncol(QUBO.MSAT.subset_genind@tab)
# De novo
QUBO.SNP.DN.R0.subset_LC <- nLoc(QUBO.SNP.DN.R0.subset_genind) 
QUBO.SNP.DN.R0.subset_AC <- ncol(QUBO.SNP.DN.R0.subset_genind@tab)
QUBO.SNP.DN.R80.subset_LC <- nLoc(QUBO.SNP.DN.R80.subset_genind) 
QUBO.SNP.DN.R80.subset_AC <- ncol(QUBO.SNP.DN.R80.subset_genind@tab)
# Reference
QUBO.SNP.REF.R0.subset_LC <- nLoc(QUBO.SNP.REF.R0.subset_genind) 
QUBO.SNP.REF.R0.subset_AC <- ncol(QUBO.SNP.REF.R0.subset_genind@tab)
QUBO.SNP.REF.R80.subset_LC <- nLoc(QUBO.SNP.REF.R80.subset_genind) 
QUBO.SNP.REF.R80.subset_AC <- ncol(QUBO.SNP.REF.R80.subset_genind@tab)

# Allelic richness: values per population
QUBO.MSAT.subset_AR <- apply(allelic.richness(QUBO.MSAT.subset_genind)$Ar, 2, mean, na.rm=TRUE)
# De novo
QUBO.SNP.DN.R0.subset_AR <- apply(allelic.richness(QUBO.SNP.DN.R0.subset_genind)$Ar, 2, mean, na.rm=TRUE)
QUBO.SNP.DN.R80.subset_AR <- apply(allelic.richness(QUBO.SNP.DN.R80.subset_genind)$Ar, 2, mean, na.rm=TRUE)
# Reference
QUBO.SNP.REF.R0.subset_AR <- apply(allelic.richness(QUBO.SNP.REF.R0.subset_genind)$Ar, 2, mean, na.rm=TRUE)
QUBO.SNP.REF.R80.subset_AR <- apply(allelic.richness(QUBO.SNP.REF.R80.subset_genind)$Ar, 2, mean, na.rm=TRUE)

# Print pop. gen. results
# MSAT
print("%%% QUBO MSAT, GARDEN AND WILD (SUBSET) %%%")
print(paste0("Heterozygosity: ", QUBO.MSAT.subset_HZ))
print(paste0("Loci count: ", QUBO.MSAT.subset_LC))
print(paste0("Allele count: ", QUBO.MSAT.subset_AC))
print(paste0("Allelic richness: ", QUBO.MSAT.subset_AR))
# SNP, De novo
# R0
print("%%% QUBO SNP (DE NOVO, R0), GARDEN AND WILD (SUBSET) %%%")
print(paste0("Heterozygosity: ", QUBO.SNP.DN.R0.subset_HZ))
print(paste0("Loci count: ", QUBO.SNP.DN.R0.subset_LC))
print(paste0("Allele count: ", QUBO.SNP.DN.R0.subset_AC))
print(paste0("Allelic richness: ", QUBO.SNP.DN.R0.subset_AR))
# R80
print("%%% QUBO SNP (DE NOVO, R80), GARDEN AND WILD (SUBSET) %%%")
print(paste0("Heterozygosity: ", QUBO.SNP.DN.R80.subset_HZ))
print(paste0("Loci count: ", QUBO.SNP.DN.R80.subset_LC))
print(paste0("Allele count: ", QUBO.SNP.DN.R80.subset_AC))
print(paste0("Allelic richness: ", QUBO.SNP.DN.R80.subset_AR))
# SNP, Reference
# R0
print("%%% QUBO SNP (REFERENCE, R0), GARDEN AND WILD (SUBSET) %%%")
print(paste0("Heterozygosity: ", QUBO.SNP.REF.R0.subset_HZ))
print(paste0("Loci count: ", QUBO.SNP.REF.R0.subset_LC))
print(paste0("Allele count: ", QUBO.SNP.REF.R0.subset_AC))
print(paste0("Allelic richness: ", QUBO.SNP.REF.R0.subset_AR))
# R80
print("%%% QUBO SNP (REFERENCE, R80), GARDEN AND WILD (SUBSET) %%%")
print(paste0("Heterozygosity: ", QUBO.SNP.REF.R80.subset_HZ))
print(paste0("Loci count: ", QUBO.SNP.REF.R80.subset_LC))
print(paste0("Allele count: ", QUBO.SNP.REF.R80.subset_AC))
print(paste0("Allelic richness: ", QUBO.SNP.REF.R80.subset_AR))

# %%% WILD ONLY ----
# SNP, REFERENCE ----
# R0 ----
# # Read in genind file: GSNAP4 alignment with Quercus robur genome; R0, NOMAF, first SNP/locus, only wild populations
# QUBO.Wild.SNP.DN.R0_filePath <-
#   "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0_NOMAF_1SNP/"
# QUBO.Wild.SNP.DN.R0_genind <- read.genepop(paste0(QUBO.Wild.SNP.DN.R0_filePath,"populations.snps.gen"))
# # Correct popNames
# pop(QUBO.Wild.SNP.DN.R0_genind) <- factor(read.table(paste0(QUBO.Wild.SNP.DN.R0_filePath, "QUBO_popmap_wild"), header=FALSE)[,2])
# 
# # Heterozygosity
# QUBO.Wild.SNP.DN.R0_HZ <- Hs(QUBO.Wild.SNP.DN.R0_genind)
# # Barplot for expected heterozygosity, SNP markers
# barplot(QUBO.Wild.SNP.DN.R0_HZ, beside = TRUE,
#         ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
#         names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf",
#                   "Kessler"),
#         main = "QUBO Wild Heterozygosity: SNPs (Complete), R0",
#         xlab = "Wild source", ylab = "Expected Heterozygosity")
# abline(h = 0, lwd=2)
# 
# # Loci count, allele count
# QUBO.Wild.SNP.DN.R0_LC <- nLoc(QUBO.Wild.SNP.DN.R0_genind)
# QUBO.Wild.SNP.DN.R0_AC <- ncol(QUBO.Wild.SNP.DN.R0_genind@tab)
# # Allelic richness: values per population
# QUBO.Wild.SNP.DN.R0_AR <- apply(allelic.richness(QUBO.Wild.SNP.DN.R0_genind)$Ar, 2, mean, na.rm=TRUE)
# # Print pop. gen. results
# print("%%% QUBO SNP (DE NOVO, R0), WILD (COMPLETE) %%%")
# print(paste0("Heterozygosity: ", QUBO.Wild.SNP.DN.R0_HZ))
# print(paste0("Loci count: ", QUBO.Wild.SNP.DN.R0_LC))
# print(paste0("Allele count: ", QUBO.Wild.SNP.DN.R0_AC))
# print(paste0("Allelic richness: ", QUBO.Wild.SNP.DN.R0_AR))
# 
# # R80 ----
# # Read in genind file: GSNAP4 alignment with Quercus robur genome; R80, NOMAF, first SNP/locus, only wild populations
# QUBO.Wild.SNP.DN.R80_filePath <-
#   "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80_NOMAF_1SNP/"
# QUBO.Wild.SNP.DN.R80_genind <- read.genepop(paste0(QUBO.Wild.SNP.DN.R80_filePath,"populations.snps.gen"))
# # Correct popNames
# pop(QUBO.Wild.SNP.DN.R80_genind) <- factor(read.table(paste0(QUBO.Wild.SNP.DN.R80_filePath,"QUBO_popmap_wild"), header=FALSE)[,2])
# 
# # Heterozygosity
# QUBO.Wild.SNP.DN.R80_HZ <- Hs(QUBO.Wild.SNP.DN.R80_genind)
# # Barplot for expected heterozygosity, SNP markers
# barplot(QUBO.Wild.SNP.DN.R80_HZ, beside = TRUE,
#         ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
#         names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf",
#                   "Kessler"),
#         main = "QUBO Wild Heterozygosity: SNPs (Complete), R80",
#         xlab = "Wild source", ylab = "Expected Heterozygosity")
# abline(h = 0, lwd=2)
# 
# # Loci count, allele count
# QUBO.Wild.SNP.DN.R80_LC <- nLoc(QUBO.Wild.SNP.DN.R80_genind)
# QUBO.Wild.SNP.DN.R80_AC <- ncol(QUBO.Wild.SNP.DN.R80_genind@tab)
# # Allelic richness: values per population
# QUBO.Wild.SNP.DN.R80_AR <- apply(allelic.richness(QUBO.Wild.SNP.DN.R80_genind)$Ar, 2, mean, na.rm=TRUE)
# # Print pop. gen. results
# print("%%% QUBO SNP (DE NOVO, R80), WILD (COMPLETE) %%%")
# print(paste0("Heterozygosity: ", QUBO.Wild.SNP.DN.R80_HZ))
# print(paste0("Loci count: ", QUBO.Wild.SNP.DN.R80_LC))
# print(paste0("Allele count: ", QUBO.Wild.SNP.DN.R80_AC))
# print(paste0("Allelic richness: ", QUBO.Wild.SNP.DN.R80_AR))
