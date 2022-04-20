# %%%%%%%%%%%%%%%%%%%%%%%%
# %%% EX SITU SAMPLING %%%
# %%%%%%%%%%%%%%%%%%%%%%%%

# This script analyzes ex situ conservation metrics (allelic capture) and runs resampling analyses
# for the two study species of the SSRvSNP study: Quercus acerifolia (QUAC) and Q. boyntonii (QUBO)

library(adegenet)
library(hierfstat)
library(RColorBrewer)

# %%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUERCUS ACERIFOLIA %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%

# READ IN AND PROCESS GENEPOP FILE----
# Set the filepath variable below to specify which populations directory to source the genpop file from
# (Note that Stacks will write this files with the suffix ".genepop", but it needs to be ".gen")
# Using the QUAC Reference dataset (aligned using GSNAP), and the "summary" populations 
# This dataset uses an -R 80 parameter, meaning loci had to be present in 80% of all samples to be retained
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_sum/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
QUAC.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Read in the Stacks popmap values, and use these to replace @pop values 
# (using the pops accessor; this is necessary because original pop names are incorrect)
pop(QUAC.genind) <- factor(read.table("../../../QUAC_popmap2", header=FALSE)[,2])
# Genind/genpop tab slot contains a matrix of allele counts
# Total loci in assembly (after Stacks populations -R80 filtering) is 14,033. But, this includes monomorphic loci
# Genpop file only includes polymorphic loci, of which there are 6,361
nLoc(QUAC.genind)
# Each locus contains two alleles, which leads to 6,361 * 2 = 12,722 columns of the sample x allele matrix
ncol(QUAC.genind@tab)
# Create vectors corresponding to sample numbers of garden and wild individuals using seq
QUAC.garden <- seq_len(length(which(pop(QUAC.genind)=="garden")))
QUAC.wild <- seq(from=length(which(pop(QUAC.genind)=="garden"))+1, to=nInd(QUAC.genind))
QUAC.wild.N <- length(QUAC.wild)
QUAC.garden.N <- length(QUAC.garden)
# rownames(QUAC.genind@tab[QUAC.garden,]) # Demonstration
# rownames(QUAC.genind@tab[QUAC.wild,]) # Demonstration

# ALLELE CATEGORIES----
# Categorize wild alleles, then determine how many of each allele category gardens are capturing

# ALL ALLELES----
# Wild allele frequency vector
QUAC_wildFreqs <- colSums(QUAC.genind@tab[QUAC.wild,], na.rm = TRUE)/(QUAC.wild.N*2)*100

# Total
length(which(names(which(QUAC_wildFreqs > 0)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs > 0))*100
# Very common
length(which(names(which(QUAC_wildFreqs > 10)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs > 10))*100
# Common
length(which(names(which(QUAC_wildFreqs > 5)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs > 5))*100
# Low frequency
length(which(names(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1))*100
# Rare
length(which(names(which(QUAC_wildFreqs < 1)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs < 1))*100
# Alternate rare calculation, corresponding to QUHA approach 
# length(which(names(which(QUAC_wildFreqs < 1 & QUAC_wildFreqs > min(QUAC_wildFreqs))) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs < 1 & QUAC_wildFreqs > min(QUAC_wildFreqs)))*100

# Somehow, without using names, we get the same values when analyzing the entire dataset (but not subsets...)
# # Total
# length(which(which(QUAC_wildFreqs > 0) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)))/length(which(QUAC_wildFreqs > 0))*100
# # Very common
# length(which(which(QUAC_wildFreqs > 10) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)))/length(which(QUAC_wildFreqs > 10))*100
# # Common
# length(which(which(QUAC_wildFreqs > 5) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)))/length(which(QUAC_wildFreqs > 5))*100
# # Low frequency
# length(which(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)))/length(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1))*100
# # Rare
# length(which(which(QUAC_wildFreqs < 1) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)))/length(which(QUAC_wildFreqs < 1))*100

# 100% Rare allele capture exploration----
# Vector of wild alleles, along with their frequencies
QUAC_rareWildAlleles <- QUAC_wildFreqs[which(QUAC_wildFreqs < 1 & QUAC_wildFreqs > 0)]
# Number of rare alleles (with and without absent alleles)
length(QUAC_rareWildAlleles) ; length(which(QUAC_rareWildAlleles > 0))

# How many rare wild alleles captured?
length(which(names(which(QUAC_wildFreqs < 1)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0))))

# Names of the alleles that are absent (0) in the garden
names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) == 0))

# Matching the names of the rare wild alleles and the names of the garden alleles that are absent
unique(match(names(QUAC_rareWildAlleles), names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) == 0))))
unique(names(QUAC_rareWildAlleles) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) == 0)))
# None of the alleles that are missing in the garden samples match the names of rare wild alleles

# Building a list of garden+wild sample combinations, for samples containing rare wild alleles
garden_wild_rareMatches <- list()
for(i in 1:length(QUAC_rareWildAlleles)){
  garden_wild_rareMatches[[i]] <- which(QUAC.genind@tab[,names(QUAC_rareWildAlleles)[i]] != 0)
}
garden_wild_rareMatches
unique(sapply(garden_wild_rareMatches, length))

# Checking out which alleles are rare in the gardens
QUAC_gardenFreqs <- colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE)/(QUAC.garden.N*2)*100
length(which(QUAC_gardenFreqs < 1)) # 1,439 rare garden alleles
# Garden rare alleles
QUAC_rareGardenAlleles <- QUAC_gardenFreqs[which(QUAC_gardenFreqs < 1 & QUAC_gardenFreqs > 0)]
# Number of rare alleles (with and without absent alleles); very similar numbers between garden and wild samples
length(QUAC_rareGardenAlleles) ; length(which(QUAC_rareGardenAlleles > 0))

# Building a list of garden+wild sample combinations, for samples containing rare garden alleles
wild_garden_rareMatches <- list()
for(i in 1:length(QUAC_rareGardenAlleles)){
  wild_garden_rareMatches[[i]] <- which(QUAC.genind@tab[,names(QUAC_rareGardenAlleles)[i]] != 0)
}
wild_garden_rareMatches
unique(sapply(wild_garden_rareMatches, length))

# REMOVING WILD SINGLETONS/DOUBLETONS----
# How many wild alleles only show up once or twice, i.e. colSums = 0, 1, or 2?
length(which(colSums(QUAC.genind@tab[QUAC.wild,], na.rm = TRUE) <= 2)) #2,673 alleles 
length(which(QUAC_wildFreqs <1)) #1,445 alleles have frequencies less than 1
# Removing singletons/doubletons suggests that we will remove all instances of rare alleles (and some low frequency alleles)
# Demo with allele 7593_61.04
QUAC.genind@tab[QUAC.wild,"7593_61.04"]
unique(QUAC.genind@tab[QUAC.wild,"7593_61.04"])
QUAC_wildFreqs["7593_61.04"]

# Subset wild allele frequency vector to only contain alleles with colSums greater than 3
QUAC_wildFreqs <- colSums(QUAC.genind@tab[QUAC.wild,which(colSums(QUAC.genind@tab[QUAC.wild,], na.rm = TRUE) >= 3)], na.rm = TRUE)/(QUAC.wild.N*2)*100

# What percentage of total wild alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs > 0)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,names(QUAC_wildFreqs)], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs > 0))*100
# What percentage of the very common wild alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs > 10)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,names(QUAC_wildFreqs)], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs > 10))*100
# What percentage of the common wild alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs > 5)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,names(QUAC_wildFreqs)], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs > 5))*100
# What percentage of the low frequency alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,names(QUAC_wildFreqs)], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1))*100
# What percentage of the rare alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs < 1)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,names(QUAC_wildFreqs)], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs < 1))*100

# DEMO: ALLELES 11-30----
# Number of garden and wild individuals with alleles 11-30
colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE)
colSums(QUAC.genind@tab[QUAC.wild,11:30], na.rm = TRUE)

# Wild allele frequency vector
QUAC_wildFreqs <- colSums(QUAC.genind@tab[QUAC.wild,11:30], na.rm = TRUE)/(QUAC.wild.N*2)*100

QUAC_wildFreqs <- QUAC_wildFreqs[which(QUAC_wildFreqs > 0)]

# Check out vectors: wild frequencies, and numbers of garden alleles present
QUAC_wildFreqs
colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE)

# Buildup: Very common wild alleles
names(which(QUAC_wildFreqs > 10))
# Which garden alleles are present
names(which(colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE) > 0))
# Are the very common wild alleles are seen in the (present) garden alleles
names(which(QUAC_wildFreqs > 10)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE) > 0))
# Which very common wild alleles are seen in the (present) garden alleles
which(names(which(QUAC_wildFreqs > 10)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE) > 0)))
# How many very common wild alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs > 10)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE) > 0))))
# What percentage of the very common wild alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs > 10)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs > 10))*100

# What percentage of total wild alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs > 0)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs > 0))*100
# What percentage of the very common wild alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs > 10)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs > 10))*100
# What percentage of the common wild alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs > 5)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs > 5))*100
# What percentage of the low frequency alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1))*100
# What percentage of the rare alleles are seen in the (present) garden alleles
length(which(names(which(QUAC_wildFreqs < 1)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,11:30], na.rm = TRUE) > 0))))/length(which(QUAC_wildFreqs < 1))*100

# OLD----
# # %%% DEMONSTRATION OF DIFFERENCES IN APPROACHES WITH AND WITHOUT NAMES %%%
# # Note the different lengths below
# length(QUAC_wildFreqs)
# length(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0))
# # If vectors are of different lengths, can you match their indices using %in% ? 
# # No: demonstrations below
# 
# # Wild alleles 9:20, and garden colSums 10:20
# QUAC_wildFreqs[9:20]
# which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]
# 
# # Without names
# # Total: should be 75%, getting 16.6%
# length(which(which(QUAC_wildFreqs[9:20] > 0) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]))/length(which(QUAC_wildFreqs[9:20] > 0))*100
# # Very common: should be 83.3%, getting 16.6%
# length(which(which(QUAC_wildFreqs[9:20] > 10) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]))/length(which(QUAC_wildFreqs[9:20] > 10))*100
# # Common: should be 83/3%, getting 16.6%
# length(which(which(QUAC_wildFreqs[9:20] > 5) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]))/length(which(QUAC_wildFreqs[9:20] > 5))*100
# # Low frequency: should be 50%, getting 25%
# length(which(which(QUAC_wildFreqs[9:20] < 10 & QUAC_wildFreqs[9:20] > 1) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]))/length(which(QUAC_wildFreqs[9:20] < 10 & QUAC_wildFreqs[9:20] > 1))*100
# # Rare: should be 100%, getting 0%
# length(which(which(QUAC_wildFreqs[9:20] < 1) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]))/length(which(QUAC_wildFreqs[9:20] < 1))*100
# 
# # Using names
# # Total
# length(which(names(which(QUAC_wildFreqs[9:20] > 0)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20])))/length(which(QUAC_wildFreqs[9:20] > 0))*100
# # Very common
# length(which(names(which(QUAC_wildFreqs[9:20] > 10)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20])))/length(which(QUAC_wildFreqs[9:20] > 10))*100
# # Common
# length(which(names(which(QUAC_wildFreqs[9:20] > 5)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20])))/length(which(QUAC_wildFreqs[9:20] > 5))*100
# # Low frequency
# length(which(names(which(QUAC_wildFreqs[9:20] < 10 & QUAC_wildFreqs[9:20] > 1)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20])))/length(which(QUAC_wildFreqs[9:20] < 10 & QUAC_wildFreqs[9:20] > 1))*100
# # Rare
# length(which(names(which(QUAC_wildFreqs[9:20] < 1)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20])))/length(which(QUAC_wildFreqs[9:20] < 1))*100
# 
# # Wild alleles 10:21, and garden colSums 10:20
# QUAC_wildFreqs[10:21]
# which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]
# 
# # Without names
# # Total: should be 83.3%, getting 16.7%
# length(which(which(QUAC_wildFreqs[10:21] > 0) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]))/length(which(QUAC_wildFreqs[10:21] > 0))*100
# # Very common: should be 100%, getting 16.7%
# length(which(which(QUAC_wildFreqs[10:21] > 10) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]))/length(which(QUAC_wildFreqs[10:21] > 10))*100
# # Common: should be 100%, getting 16.7%
# length(which(which(QUAC_wildFreqs[10:21] > 5) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]))/length(which(QUAC_wildFreqs[10:21] > 5))*100
# # Low frequency: should be 50%, getting 0%
# length(which(which(QUAC_wildFreqs[10:21] < 10 & QUAC_wildFreqs[10:21] > 1) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]))/length(which(QUAC_wildFreqs[10:21] < 10 & QUAC_wildFreqs[10:21] > 1))*100
# # Rare: should be 100%, getting 50%
# length(which(which(QUAC_wildFreqs[10:21] < 1) %in% which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20]))/length(which(QUAC_wildFreqs[10:21] < 1))*100
# 
# # Using names
# # Total
# length(which(names(which(QUAC_wildFreqs[10:21] > 0)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20])))/length(which(QUAC_wildFreqs[10:21] > 0))*100
# # Very common: should be 100%, getting 16.7%
# length(which(names(which(QUAC_wildFreqs[10:21] > 10)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20])))/length(which(QUAC_wildFreqs[10:21] > 10))*100
# # Common: should be 100%, getting 16.7%
# length(which(names(which(QUAC_wildFreqs[10:21] > 5)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20])))/length(which(QUAC_wildFreqs[10:21] > 5))*100
# # Low frequency: should be 50%, getting 0%
# length(which(names(which(QUAC_wildFreqs[10:21] < 10 & QUAC_wildFreqs[10:21] > 1)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20])))/length(which(QUAC_wildFreqs[10:21] < 10 & QUAC_wildFreqs[10:21] > 1))*100
# # Rare: should be 100%, getting 50%
# length(which(names(which(QUAC_wildFreqs[10:21] < 1)) %in% names(which(colSums(QUAC.genind@tab[QUAC.garden,], na.rm = TRUE) > 0)[10:20])))/length(which(QUAC_wildFreqs[10:21] < 1))*100
#
# # %%% DEMONSTRATION OF OLD APPROACH WITHOUT NAMES, USING ALLELES 11-20 %%%
# # Number of garden and wild individuals with alleles 11-20
# colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) 
# colSums(QUAC.genind@tab[QUAC.wild,11:20], na.rm = TRUE)
# 
# # Wild allele frequency vector
# QUAC_wildFreqs <- colSums(QUAC.genind@tab[QUAC.wild,11:20], na.rm = TRUE)/(QUAC.wild.N*2)*100
# 
# # Check out vectors: wild frequencies, and numbers of garden alleles present
# QUAC_wildFreqs
# colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE)
# 
# # Buildup: Very common wild alleles
# which(QUAC_wildFreqs > 10)
# # Which garden alleles are present
# which(colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) > 0)
# # Are the very common wild alleles are seen in the (present) garden alleles
# which(QUAC_wildFreqs > 10) %in% which(colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) > 0)
# # Which very common wild alleles are seen in the (present) garden alleles
# which(which(QUAC_wildFreqs > 10) %in% which(colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) > 0))
# # How many very common wild alleles are seen in the (present) garden alleles
# length(which(which(QUAC_wildFreqs > 10) %in% which(colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) > 0)))
# # What percentage of the very common wild alleles are seen in the (present) garden alleles
# length(which(which(QUAC_wildFreqs > 10) %in% which(colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) > 0)))/length(which(QUAC_wildFreqs > 10))*100
# 
# # What percentage of total wild alleles are seen in the (present) garden alleles
# length(which(which(QUAC_wildFreqs > 0) %in% which(colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) > 0)))/length(which(QUAC_wildFreqs > 0))*100
# # What percentage of the very common wild alleles are seen in the (present) garden alleles (again)
# length(which(which(QUAC_wildFreqs > 10) %in% which(colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) > 0)))/length(which(QUAC_wildFreqs > 10))*100
# # What percentage of the common wild alleles are seen in the (present) garden alleles
# length(which(which(QUAC_wildFreqs > 5) %in% which(colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) > 0)))/length(which(QUAC_wildFreqs > 5))*100
# # Which low frequency alleles are seen (present) in the garden alleles
# which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1) %in% which(colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) > 0)
# # What percentage of the low frequency alleles are seen in the (present) garden alleles
# length(which(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1) %in% which(colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) > 0)))/length(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1))*100
# # What percentage of the rare alleles are seen in the (present) garden alleles
# length(which(which(QUAC_wildFreqs < 1) %in% which(colSums(QUAC.genind@tab[QUAC.garden,11:20], na.rm = TRUE) > 0)))/length(which(QUAC_wildFreqs < 1))*100
# 
# # %%% DEMONSTRATIONS OF DIFFERENT ALLELES BEING RETAING BASED ON DIFFERENT SEPARATION TECHNIQUES %%%
# # Subsetting wild samples using seppop/repool
# QUAC.genind.wild <- seppop(QUAC.genind)[2:6]
# # (Using repool will subset the resulting genind object to only include loci present in the combined populations)
# QUAC.genind.wild <- repool(QUAC.genind.wild$porterMt,QUAC.genind.wild$magazineMt,QUAC.genind.wild$pryorMt,
#                            QUAC.genind.wild$sugarloaf_midlandPeak,QUAC.genind.wild$kessler_shaleBarrenRidge)
# 
# # Subsetting wild samples using repool_new function from Sean
# repool_new <- function(genind_obj,vect_pops){
#   genind_obj_sep <- seppop(genind_obj)
#   genind_obj_merge <- genind_obj_sep[[vect_pops[1]]]
#   for (i in 1:(length(vect_pops)-1)) genind_obj_merge<-repool(genind_obj_merge,genind_obj_sep[[vect_pops[i+1]]])
#   genind_obj_merge
# }
# QUAC.genind.wild.NEW <- repool_new(QUAC.genind,2:6)
# 
# # NOTE: there are differences in allele matrices depending on how genind/genpop objects are constructed
# ncol(QUAC.genind@tab)
# ncol(QUAC.genpop@tab)
# # Without na.rm argument, colSums don't match across genind and genpop objects
# identical(colSums(QUAC.genind@tab),colSums(QUAC.genpop@tab))
# identical(colSums(QUAC.genind@tab, na.rm = TRUE),colSums(QUAC.genpop@tab, na.rm = TRUE))
# # Subset objects have NAs removed
# ncol(QUAC.genind.wild@tab)
# ncol(QUAC.genind.wild.NEW@tab)
# ncol(QUAC.genpop.wild@tab)

# RESAMPLING----
# Create a matrix of ONLY wild individuals with present alleles
QUAC.wild.mat <- QUAC.genind@tab[QUAC.wild,which(colSums(QUAC.genind@tab[QUAC.wild,], na.rm = TRUE) > 0)]

# Processing duplicates
# Only one duplicate, in QUAC Wild data set
which(rownames(QUAC.wild.mat) == "QUAC_W_DUP_SH_Q2121")
QUAC.wild.mat[65:66,1:7] # QUAC_W_SH_Q2121 and QUAC_W_DUP_SH_Q2121
length(which(QUAC.wild.mat[65,] != QUAC.wild.mat[66,]))
# 874/12,159 alleles (7.19%) are different between these duplicate samples
# For now, just removing duplicate sample
QUAC.wild.mat <- QUAC.wild.mat[-66,]

# Function for measuring wild allelic capture of samples
# The 2 arguments are: vector of wild allele frequencies, and sample matrix (of the wild genind object, nrow=2:96)
get.allele.cat <- function(freq.vector, sample.mat){
  # Total alleles
  # Determine how many alleles in the sample (i.e. greater than 0) are found in the frequency vector 
  total <- length(which(names(which(freq.vector > 0)) %in% names(which(colSums(sample.mat, na.rm = TRUE) > 0))))/length(which(freq.vector > 0))*100
  # Very common alleles (greater than 10%)
  v_com <- length(which(names(which(freq.vector > 10)) %in% names(which(colSums(sample.mat, na.rm = TRUE) > 0))))/length(which(freq.vector > 10))*100
  # Common alleles (greater than 5%)
  com <- length(which(names(which(freq.vector > 5)) %in% names(which(colSums(sample.mat, na.rm = TRUE) > 0))))/length(which(freq.vector > 5))*100
  # Low frequency alleles (between 1% and 10%)
  low_freq <- length(which(names(which(freq.vector < 10 & freq.vector > 1)) %in% names(which(colSums(sample.mat, na.rm = TRUE) > 0))))/length(which(freq.vector < 10 & freq.vector > 1))*100
  # Rare alleles (less than 1%)
  rare <- length(which(names(which(freq.vector < 1)) %in% names(which(colSums(sample.mat, na.rm = TRUE) > 0))))/length(which(freq.vector < 1))*100
  # Concatentate values to a vector, and return
  return(c(total,v_com,com,low_freq,rare))
}

# %%% ENTIRE DATASET %%%
# Build a 3D array, with rows being sample numbers and columns being allele frequency categories
# 3rd dimension is replicates. nrow is number of individuals-1 because sample doesn't work for vectors of length 1
num_reps <- 5 # In Hoban et al. 2020 (where QUBO is considered), num_reps = 75,000
list_allele_cat <- c("tot","v_com","com","low_freq","rare")
samplingResults <- array(dim=c(nrow(QUAC.wild.mat)-1,length(list_allele_cat),num_reps))
colnames(samplingResults) <- list_allele_cat
# Generate wild allele frequency vector
QUAC_wildFreqs <- colSums(QUAC.genind@tab[QUAC.wild,], na.rm = TRUE)/(QUAC.wild.N*2)*100
# Remove garden specific alleles.
# Failing to do this will prevent rare capture rates from reaching 100%, even when all wild individuals are resampled
QUAC_wildFreqs <- QUAC_wildFreqs[which(QUAC_wildFreqs > 0)]

# For each replicate (which is the third dimension, in the samplingResults array)...
for(i in 1:num_reps){
  # prog.bar(i,num_reps) <-- Use browser() to figure out how this fxn is working. Then comment it out!
  # loop through sampling from 2 to the maximum number of wild individuals
  for(j in 2:nrow(QUAC.wild.mat)){
    # Create a sample of the wild allele matrix, of "j" size
    samp <- QUAC.wild.mat[sample(nrow(QUAC.wild.mat), size=j, replace = FALSE),]
    # Calculate how many alleles of each category that sample captures,
    # and place those percentages into the row of the samplingResults array
    samplingResults[j-1,,i] <- get.allele.cat(QUAC_wildFreqs,samp)
  }
}
samplingResults[,,1]
samplingResults[,,2]
samplingResults[,,3]
samplingResults[,,4]
samplingResults[,,5]

# OUTDATED----
# Function for measuring wild allelic capture of samples
# The 2 arguments are: vector of wild allele frequencies, and sample matrix (of the wild genind object, nrow=2:96)
get.allele.cat.OLD <- function(freq.vector, sample.mat){
  # Total alleles
  # Determine how many alleles in the sample (i.e. greater than 0) are found in the frequency vector
  total <- length(which(which(freq.vector > 0) %in% which(colSums(sample.mat, na.rm = TRUE) > 0)))/length(which(freq.vector > 0))*100
  # Very common alleles (greater than 10%)
  v_com <- length(which(which(freq.vector > 10) %in% which(colSums(sample.mat, na.rm = TRUE) > 0)))/length(which(freq.vector > 10))*100
  # Common alleles (greater than 5%)
  com <- length(which(which(freq.vector > 5) %in% which(colSums(sample.mat, na.rm = TRUE) > 0)))/length(which(freq.vector > 5))*100
  # Low frequency alleles (between 1% and 10%)
  low_freq <- length(which(which(freq.vector < 10 & freq.vector > 1) %in% which(colSums(sample.mat, na.rm = TRUE) > 0)))/length(which(freq.vector < 10 & freq.vector > 1))*100
  # Rare alleles (less than 1%)
  rare <- length(which(which(freq.vector < 1) %in% which(colSums(sample.mat, na.rm = TRUE) > 0)))/length(which(freq.vector < 1))*100
  # Concatentate values to a vector, and return
  return(c(total,v_com,com,low_freq,rare))
}

# Build a 3D array, with rows being sample numbers and columns being allele frequency categories
# 3rd dimension is replicates. nrow is number of individuals-1 because sample doesn't work for vectors of length 1
num_reps <- 5 # In Hoban et al. 2020 (where QUBO is considered), num_reps = 75,000
list_allele_cat <- c("tot","v_com","com","low_freq","rare")
samplingResults.OLD <- array(dim=c(nrow(QUAC.wild.mat)-1,length(list_allele_cat),num_reps))
colnames(samplingResults.OLD) <- list_allele_cat
# Generate wild allele frequency vector
QUAC_wildFreqs <- colSums(QUAC.genind@tab[QUAC.wild,], na.rm = TRUE)/(QUAC.wild.N*2)*100
# Remove garden specific alleles.
# Failing to do this will prevent rare capture rates from reaching 100%, even when all wild individuals are resampled
QUAC_wildFreqs <- QUAC_wildFreqs[which(QUAC_wildFreqs > 0)]

# For each replicate (which is the third dimension, in the samplingResults array)...
for(i in 1:num_reps){
  # loop through sampling from 2 to the maximum number of wild individuals
  for(j in 2:nrow(QUAC.wild.mat)){
    # Create a sample of the wild allele matrix, of "j" size
    samp <- QUAC.wild.mat[sample(nrow(QUAC.wild.mat), size=j, replace = FALSE),]
    # Calculate how many alleles of each category that sample captures,
    # and place those percentages into the row of the samplingResults array
    samplingResults.OLD[j-1,,i] <- get.allele.cat.OLD(QUAC_wildFreqs,samp)
  }
}
samplingResults.OLD[,,1]
samplingResults.OLD[,,2]
samplingResults.OLD[,,3]
samplingResults.OLD[,,4]
samplingResults.OLD[,,5]

# # Old function for measuring wild allelic capture of samples
# # The two arguments are the vector of frequencies of wild alleles, and the sample of the wild genind object
# get.allele.cat <- function(freq.vector, sample.mat){
#   # Total alleles
#   # Determine how many alleles in the sample (i.e. greater than 0) are found in the frequency vector 
#   total <- length(which(colSums(sample.mat, na.rm = TRUE) > 0) %in% freq.vector)/length(freq.vector)*100
#   # Very common alleles (greater than 10%)
#   v_com <- length(which(which(colSums(sample.mat, na.rm = TRUE) > 0) %in% which(freq.vector > 10)))/length(which(freq.vector > 10))*100
#   # Common alleles (greater than 5%)
#   com <- length(which(which(colSums(sample.mat, na.rm = TRUE) > 0) %in% which(freq.vector > 5)))/length(which(freq.vector > 5))*100
#   # Low frequency alleles (between 1% and 10%)
#   low_freq <- length(which(which(colSums(sample.mat, na.rm = TRUE) > 0) %in% which(freq.vector < 10 & freq.vector > 1)))/length(which(freq.vector < 10 & freq.vector > 1))*100
#   # Rare alleles (less than 1%)
#   rare <- length(which(which(colSums(sample.mat, na.rm = TRUE) > 0) %in% which(freq.vector < 1)))/length(which(freq.vector < 1))*100
#   # Concatentate values to a vector, and return
#   return(c(total,v_com,com,low_freq,rare))
# }

# RECREATING ABOVE RESAMPLING CODE USING APPLY FUNCTIONS----

# First, trying to recreate the results of the inner loop (j), which samples the wild allele matrix
# and then calculates the captured allele frequencies of that sample
samp <- QUAC.wild.mat[sample(nrow(QUAC.wild.mat), size=2, replace = FALSE),]

# 2 approaches: use a lambda function within some kind of apply function,
# or create a new function and then try to apply that over the wild allele matrix

# Lambda function
# Works for 2, but not multiple values
NEW.samp <- apply(QUAC.wild.mat,2, function(x) sample(x, size=2, replace = FALSE))

NEW.samp <- apply(QUAC.wild.mat,2, function(x) sample(x, size=seq(from=2, to=nrow(x)), replace = FALSE))
# The 2 calls below still only use 2 as the sample size
NEW.samp <- apply(QUAC.wild.mat,2, function(x) sample(x, size=2:96, replace = FALSE))

NEW.samp <- apply(QUAC.wild.mat,2, function(x) sample(x, size=seq(from=2, to=nrow(QUAC.wild.mat)), replace = FALSE))

# New function
resample.test <- function(wild.matrix, sample.number, allele.frequencies){
  # Create a sample of the wild allele matrix, of "j" size
  # wild.sample <- wild.matrix[sample(nrow(wild.matrix), size=sample.number, replace = FALSE)]
  wild.sample <- wild.matrix[sample(wild.matrix, size=sample.number, replace = FALSE),]
  # Calculate how many alleles of each category that sample captures,
  # and output those percentages
  output <- get.allele.cat(allele.frequencies, wild.sample)
}

# Single example
test <- resample.test(QUAC.wild.mat, sample.number = 2, allele.frequencies = QUAC.wildAlleleFreqs)
# Again, this doesn't work, because I'm trying to pass a vector to the size argument of 'sample'...
test.resamplingMat <- apply(QUAC.wild.mat, 2, resample.test, sample.number=seq(from=2, to=nrow(QUAC.wild.mat)), allele.frequencies=QUAC.wildAlleleFreqs)
# Attempts using replicate function
replicate(96, sample(QUAC.wild.mat, size = 2:nrow(QUAC.wild.mat), replace = FALSE), simplify=FALSE)
replicate(96, sample(QUAC.wild.mat, replace = FALSE), simplify=FALSE)
yo <- replicate(96, QUAC.wild.mat[sample(nrow(QUAC.wild.mat), size=2:nrow(QUAC.wild.mat), replace = FALSE),], simplify=FALSE)
# Variable samples per matrix, but not ordered...size argument is itself called from sample?
yo <- replicate(96, QUAC.wild.mat[sample(nrow(QUAC.wild.mat), sample(2:nrow(QUAC.wild.mat), 1L), replace = FALSE),], simplify=FALSE)
# 2 samples per matrix
yo <- replicate(96, QUAC.wild.mat[sample(nrow(QUAC.wild.mat), 2:nrow(QUAC.wild.mat), replace = FALSE),], simplify=FALSE)
# 96 samples per matrix
yo <- replicate(96, QUAC.wild.mat[sample(nrow(QUAC.wild.mat), 96, replace = FALSE),], simplify=FALSE)

# (tapply/lambda function example)
input <- 1:10
grouping <- rep(letters[1:2], 5)
tapply(input, grouping, function(x) sd(x)/sqrt(length(x)))

# RESAMPLING: PLOTTING----
total_means <- apply(samplingResults[,1,], 1, mean)
total_sd <- apply(samplingResults[,1,], 1, sd)

v.com_means <- apply(samplingResults[,2,], 1, mean)
v.com_sd <- apply(samplingResults[,2,], 1, sd)

com_means <- apply(samplingResults[,3,], 1, mean)
com_sd <- apply(samplingResults[,3,], 1, sd)

lowfr_means <- apply(samplingResults[,4,], 1, mean)
lowfr_sd <- apply(samplingResults[,4,], 1, sd)

rare_means <- apply(samplingResults[,5,], 1, mean)
rare_sd <- apply(samplingResults[,5,], 1, sd)

plotColors <- brewer.pal(n=5, name="Dark2")
# plotColors <- c("red","firebrick","darkorange3","coral","deeppink4")
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend("bottomright", inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.2)

# Zoomed in 
plot(v.com_means, ylim=c(0,110), xlim=c(0,5), col=plotColors[1], pch=16, main="QUAC: Resampling Curve (10 replicates)",
     xlab="Number of samples", ylab="Mean percent captured")
points(com_means, col=plotColors[2], pch=16)
points(lowfr_means, col=plotColors[3], pch=16)
points(rare_means, col=plotColors[4], pch=16)
legend("bottomright", inset = 0.05, legend = c("Very common: >10%","Common: >5%","Low frequency: 1--10%", "Rare: <1%"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2)

# How many samples reach 95% genetic capture threshold?
# Extract just the results for all alleles
samplingResults[,1,]
# Below line averages across reps (supposedly? Seems like MARGIN should be 2, not 1...)
apply(samplingResults[,1,],1,mean)
# 60 samples captures 95% of total genetic diversity


# %%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUERCUS BOYNTONII %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%

# READ IN AND PROCESS GENEPOP FILE----
# Set the filepath variable below to specify which populations directory to source the genpop file from
# (Note that Stacks will write this files with the suffix ".genepop", but it needs to be ".gen")
# Using the QUBO Reference dataset (aligned using GSNAP), and the "summary" populations 
# This dataset uses an -R 80 parameter, meaning loci had to be present in 80% of all samples to be retained
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP/output/populations_sum/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
QUBO.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Read in the Stacks popmap values, and use these to replace @pop values 
# (using the pops accessor; this is necessary because original pop names are incorrect)
pop(QUBO.genind) <- factor(read.table("../../../QUBO_popmap2", header=FALSE)[,2])
# Create vectors corresponding to sample numbers of garden and wild individuals using seq
QUBO.garden <- seq_len(length(which(pop(QUBO.genind)=="garden")))
QUBO.wild <- seq(from=length(which(pop(QUBO.genind)=="garden"))+1, to=nInd(QUBO.genind))
QUBO.wild.N <- length(QUBO.wild)

# ALLELE CATEGORIES----
# Categorize wild alleles, then determine how many of each allele category gardens are capturing

# ALL ALLELES----
# Wild allele frequency vector
QUBO_wildFreqs <- colSums(QUBO.genind@tab[QUBO.wild,], na.rm = TRUE)/(QUBO.wild.N*2)*100

# Total
length(which(names(which(QUBO_wildFreqs > 0)) %in% names(which(colSums(QUBO.genind@tab[QUBO.garden,], na.rm = TRUE) > 0))))/length(which(QUBO_wildFreqs > 0))*100
# Very common
length(which(names(which(QUBO_wildFreqs > 10)) %in% names(which(colSums(QUBO.genind@tab[QUBO.garden,], na.rm = TRUE) > 0))))/length(which(QUBO_wildFreqs > 10))*100
# Common
length(which(names(which(QUBO_wildFreqs > 5)) %in% names(which(colSums(QUBO.genind@tab[QUBO.garden,], na.rm = TRUE) > 0))))/length(which(QUBO_wildFreqs > 5))*100
# Low frequency
length(which(names(which(QUBO_wildFreqs < 10 & QUBO_wildFreqs > 1)) %in% names(which(colSums(QUBO.genind@tab[QUBO.garden,], na.rm = TRUE) > 0))))/length(which(QUBO_wildFreqs < 10 & QUBO_wildFreqs > 1))*100
# Rare
length(which(names(which(QUBO_wildFreqs < 1)) %in% names(which(colSums(QUBO.genind@tab[QUBO.garden,], na.rm = TRUE) > 0))))/length(which(QUBO_wildFreqs < 1))*100

# REMOVING WILD SINGLETONS/DOUBLETONS----
# Subset wild allele frequency vector to only contain alleles with colSums greater than 3
QUBO_wildFreqs <- colSums(QUBO.genind@tab[QUBO.wild,which(colSums(QUBO.genind@tab[QUBO.wild,], na.rm = TRUE) >= 3)], na.rm = TRUE)/(QUBO.wild.N*2)*100

# What percentage of total wild alleles are seen in the (present) garden alleles
length(which(names(which(QUBO_wildFreqs > 0)) %in% names(which(colSums(QUBO.genind@tab[QUBO.garden,names(QUBO_wildFreqs)], na.rm = TRUE) > 0))))/length(which(QUBO_wildFreqs > 0))*100
# What percentage of the very common wild alleles are seen in the (present) garden alleles
length(which(names(which(QUBO_wildFreqs > 10)) %in% names(which(colSums(QUBO.genind@tab[QUBO.garden,names(QUBO_wildFreqs)], na.rm = TRUE) > 0))))/length(which(QUBO_wildFreqs > 10))*100
# What percentage of the common wild alleles are seen in the (present) garden alleles
length(which(names(which(QUBO_wildFreqs > 5)) %in% names(which(colSums(QUBO.genind@tab[QUBO.garden,names(QUBO_wildFreqs)], na.rm = TRUE) > 0))))/length(which(QUBO_wildFreqs > 5))*100
# What percentage of the low frequency alleles are seen in the (present) garden alleles
length(which(names(which(QUBO_wildFreqs < 10 & QUBO_wildFreqs > 1)) %in% names(which(colSums(QUBO.genind@tab[QUBO.garden,names(QUBO_wildFreqs)], na.rm = TRUE) > 0))))/length(which(QUBO_wildFreqs < 10 & QUBO_wildFreqs > 1))*100
# What percentage of the rare alleles are seen in the (present) garden alleles
length(which(names(which(QUBO_wildFreqs < 1)) %in% names(which(colSums(QUBO.genind@tab[QUBO.garden,names(QUBO_wildFreqs)], na.rm = TRUE) > 0))))/length(which(QUBO_wildFreqs < 1))*100
