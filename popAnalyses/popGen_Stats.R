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
# Set working directory
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)

# %%%% QUAC %%%% ----
# ---- MSATS ----
# GARDEN AND WILD ----
# Read in genind file (GCC_QUAC_ZAIN repo; QUAC_wK_garden_wild_clean.gen)
QUAC.MSAT.genpop.filePath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/Garden_Wild/"
setwd(QUAC.MSAT.genpop.filePath)
QUAC.MSAT.genind <- read.genepop("QUAC_wK_garden_wild_clean.gen", ncode = 3)
# Correct popNames: pop1 is Garden, pop2 is Wild
pop(QUAC.MSAT.genind) <- gsub("pop1", "garden", pop(QUAC.MSAT.genind))
pop(QUAC.MSAT.genind) <- gsub("pop2", "wild", pop(QUAC.MSAT.genind))

# Heterozygosity
QUAC.MSAT_HZ <- Hs(QUAC.MSAT.genind); print(QUAC.MSAT_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.MSAT_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: MSATs (Complete)", 
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
QUAC.SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])

# Heterozygosity
QUAC.SNP_HZ <- Hs(QUAC.SNP.genind); print(QUAC.SNP_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP_HZ, beside = TRUE, 
        ylim = c(0,0.3), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: SNPs (Complete)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUAC.SNP.genind)
# Allelic richness: values per population
apply(allelic.richness(QUAC.SNP.genind)$Ar, 2, mean, na.rm=TRUE)

# WILD ONLY ----
# Read in genind file: Optimized de novo assembly; R0, NOMAF, first SNP/locus, only wild populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF_1SNP/"
setwd(genpop.filePath)
QUAC.Wild.SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.Wild.SNP.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# Heterozygosity
QUAC.SNP.Wild_HZ <- Hs(QUAC.Wild.SNP.genind); print(QUAC.SNP.Wild_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUAC.SNP.Wild_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf", 
                  "Kessler"), 
        main = "QUAC Heterozygosity: SNPs", 
        xlab = "Wild population", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUAC.Wild.SNP.genind)
# Allelic richness: values per population
apply(allelic.richness(QUAC.Wild.SNP.genind)$Ar, 2, mean, na.rm=TRUE)

# ---- SUBSET ----
# MSAT: read in tissue database names from GCC_QUAC_ZAIN repository, and rename genind object rows
QUAC.MSAT.tissueNames_filepath <- "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Data_Frames/QUAC_allpop_clean_df.csv"
QUAC.MSAT.tissueNames <- unlist(read.csv2(QUAC.MSAT.tissueNames_filepath, header = TRUE, sep=",")[1])
rownames(QUAC.MSAT.genind@tab) <- QUAC.MSAT.tissueNames

# SNP: SNP Tissue sample names, and rename SNP genind matrix
# This file was created (by Austin K.), and can be found on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
QUAC.SNP.tissueNames_filepath <- paste0(SSRvSNP.wd,"exSituRepresentation/QUAC_SNP_TissueNames.csv")
QUAC.SNP.tissueNames <- unlist(read.csv2(QUAC.SNP.tissueNames_filepath, header = TRUE, sep = ",")[3])
rownames(QUAC.SNP.genind@tab) <- QUAC.SNP.tissueNames

# Subset SNP sample names by those that are also seen within the MSAT samples 
QUAC_sharedSamples <- sort(QUAC.SNP.tissueNames[which(QUAC.SNP.tissueNames %in% QUAC.MSAT.tissueNames)])
# Subset MSAT and SNP genind matrices to strictly shared samples, dropping now absent alleles
QUAC.SNP.genind.subset <- QUAC.SNP.genind[QUAC_sharedSamples,, drop=TRUE]
QUAC.MSAT.genind.subset <- QUAC.MSAT.genind[QUAC_sharedSamples,, drop=TRUE]

# Heterozygosity
QUAC.SNP.subset_HZ <- Hs(QUAC.SNP.genind.subset); print(QUAC.SNP.subset_HZ)
QUAC.MSAT.subset_HZ <- Hs(QUAC.MSAT.genind.subset); print(QUAC.MSAT.subset_HZ)
# Barplot for expected heterozygosity, SNP subset
barplot(QUAC.SNP.subset_HZ, beside = TRUE, 
        ylim = c(0,0.), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: SNPs (Subset)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Barplot for expected heterozygosity, MSAT subset
barplot(QUAC.MSAT.subset_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUAC Heterozygosity: MSATs (Subset)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUAC.SNP.genind.subset)
nLoc(QUAC.MSAT.genind.subset)
# Allelic richness: values per population
apply(allelic.richness(QUAC.SNP.genind.subset)$Ar, 2, mean, na.rm=TRUE)
apply(allelic.richness(QUAC.MSAT.genind.subset)$Ar, 2, mean, na.rm=TRUE)

# %%%% QUBO %%%% ----
# ---- MSATS ----
# GARDEN AND WILD ----
# Read in genind file (Southeast Oaks repo; genetic_data/Qb_total.gen file)
genpop.filePath <- 
  "~/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
setwd(genpop.filePath)
QUBO.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_total.gen"), quiet = TRUE, ncode=3)
# Correct popNames: last population (IMLS4_MP1_IMLS336_C05) is Garden; rest are Wild
pop(QUBO.MSAT.genind) <- gsub("IMLS4_MP1_IMLS336_C05", "garden", pop(QUBO.MSAT.genind))

# Heterozygosity
QUBO.MSAT_HZ <- Hs(QUBO.MSAT.genind); print(QUBO.MSAT_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.MSAT_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Heterozygosity: MSATs (Complete)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Allele counts, allelic richness
ncol(QUBO.MSAT.genind@tab)
# Allelic richness: values per population
apply(allelic.richness(QUBO.MSAT.genind)$Ar, 2, mean, na.rm=TRUE)

# ---- SNPS ----
# GARDEN AND WILD ----
# Read in genind file: GSNAP4 alignment with Quercus robur genome; R0, NOMAF, first SNP/locus, 2 populations 
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# Heterozygosity
QUBO.SNP_HZ <- Hs(QUBO.SNP.genind); print(QUBO.SNP_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.SNP_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Heterozygosity: SNPs (Complete)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUBO.SNP.genind)
# Allelic richness: values per population
apply(allelic.richness(QUBO.SNP.genind)$Ar, 2, mean, na.rm=TRUE)

# WILD ONLY ----
# Read in genind file: GSNAP4 alignment with Quercus robur genome; R0, NOMAF, first SNP/locus, only wild populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0_NOMAF_1SNP/"
setwd(genpop.filePath)
QUBO.Wild.SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.Wild.SNP.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# Heterozygosity
QUBO.Wild.SNP_HZ <- Hs(QUBO.Wild.SNP.genind); print(QUBO.Wild.SNP_HZ)
# Barplot for expected heterozygosity, SNP markers
barplot(QUBO.Wild.SNP_HZ, beside = TRUE, 
        ylim = c(0,0.75), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Porter Mt.", "Magazine Mt.", "Pryor Mt.", "Sugarloaf", 
                  "Kessler"), 
        main = "QUAC Heterozygosity: SNPs", 
        xlab = "Wild source", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUBO.Wild.SNP.genind)
# Allelic richness: values per population
apply(allelic.richness(QUBO.Wild.SNP.genind)$Ar, 2, mean, na.rm=TRUE)

# ---- SUBSET ----
# MSAT: Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.MSAT.genind), function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.MSAT.genind) <- QUBO.MSAT.sampleNames

# SNP: Remove QUBO_W_ headers from sample names
QUBO.SNP.sampleNames <- gsub("QUBO_G_",replacement = "", row.names(QUBO.SNP.genind@tab))
QUBO.SNP.sampleNames <- gsub("QUBO_W_",replacement = "", QUBO.SNP.sampleNames)
# Replace SH-Q names in SNP list with IMLS names
# These were determined by Austin K., and are outlined on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
# Only 1 of the 11 SH_Q garden samples has an IMLS sample name (SHQ2177)
QUBO.SNP.sampleNames <- gsub("SH_Q2177",replacement = "IMLS338", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2178",replacement = "IMLS312", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2179",replacement = "IMLS062", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2180",replacement = "IMLS051", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2181",replacement = "IMLS011", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2182",replacement = "IMLS144", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2183",replacement = "IMLS170", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2184",replacement = "IMLS005", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2186",replacement = "IMLS017", QUBO.SNP.sampleNames)
# Rename sample matrix
rownames(QUBO.SNP.genind@tab) <- QUBO.SNP.sampleNames

# ---- GENERATE DATASET OF SHARED SAMPLES ----
# Subset SNP sample names by those that are also seen within the MSAT samples
QUBO_sharedSamples <- sort(QUBO.SNP.sampleNames[which(QUBO.SNP.sampleNames %in% QUBO.MSAT.sampleNames)])
# Subset MSAT and SNP wild matrix objects to strictly shared samples
QUBO.SNP.genind.subset <- QUBO.SNP.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.MSAT.genind.subset <- QUBO.MSAT.genind[QUBO_sharedSamples,, drop=TRUE]

# Heterozygosity
QUBO.SNP.subset_HZ <- Hs(QUBO.SNP.genind.subset); print(QUBO.SNP.subset_HZ)
QUBO.MSAT.subset_HZ <- Hs(QUBO.MSAT.genind.subset); print(QUBO.MSAT.subset_HZ)
# Barplot for expected heterozygosity, SNP subset
barplot(QUBO.SNP.subset_HZ, beside = TRUE, 
        ylim = c(0,0.), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Heterozygosity: SNPs (Subset)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)
# Barplot for expected heterozygosity, MSAT subset
barplot(QUBO.MSAT.subset_HZ, beside = TRUE, 
        ylim = c(0,0.7), col = c("darkseagreen1", rep("darkgreen", 5)),
        names = c("Garden", "Wild"), 
        main = "QUBO Heterozygosity: MSATs (Subset)", 
        xlab = "Population Type", ylab = "Expected Heterozygosity")
abline(h = 0, lwd=2)

# Allele counts: number of variant sites (number of loci, and 1 SNP/locus)
nLoc(QUBO.SNP.genind.subset)
nLoc(QUBO.MSAT.genind.subset)
# Allelic richness: values per population
apply(allelic.richness(QUBO.SNP.genind.subset)$Ar, 2, mean, na.rm=TRUE)
apply(allelic.richness(QUBO.MSAT.genind.subset)$Ar, 2, mean, na.rm=TRUE)
