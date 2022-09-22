# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% EX SITU REPRESENTATION RATES: SNP %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses the functions declared in the functions_exSituRepresentation.R file to generate ex situ representation values 
# for Quercus acerifolia (QUAC; optimized Stacks de novo assembly, m 7, M/n 4, gt-alpha 0.01) 
# and Quercus boyntonii (QUBO; GSNAP4 alignment with Quercus robur reference) samples

# Ex situ representation rates (i.e. how well do gardens individuals contain, or represent, wild alleles)
# are impacted by methodology ("approach") and by the filters utilized
# in the Stacks populations module, which can filter the loci or the SNPs included in an analysis.

# %%% METHODOLOGY %%%
# We explored 2 different approaches, each of which had 2 possible alternatives

# -SINGLE GENIND (TOGETHER): a single genind file containing both garden and wild samples was 
# analyzed (separating samples out by rows). This used the reportAllelicRepresentation_Together* functions
# -TWO GENINDS (SEPARATE): two genind files, one for garden, and one for wild, were compared 
# by measuring the allele names across each (which were found to be consistent). 
# This used the reportAllelicRepresentation_Separate* functions
# -_Partial: the Together/Separate_Partial functions only look for the beginning section of
# an allele name (which was found to correspond to the locus number, in the locus catalog
# built by gstacks). 
# -(Complete): functions that do not include "_Partial" look for whole alleles names
# (i.e. LocusNumber_SNPposition.SNP)

# %%% FILTERS %%%
# Filters are described below, and code sections are broken out according to the filters used

# -R: corresponds the percentage of individuals a locus needs to be present in, 
# to be included in the analysis (0, 80%, or 100%)
# -NOMAF: means "no minor allele frequency": no threshold for a nucleotide's presence. 
# If an analysis doesn't have NOMAF, then the minor allele frequency filter is 1%
# -AllSNPs: means every SNP is written for every locus.
# -1SNP: means only the first SNP is written for every locus
# -H: SNPs are filtered haplotype-wise (unshared SNPs are pruned to reduce haplotype-wise missing data).
# Otherwise: only a single (random) SNP is written per locus
# -TwoPops: means all wild samples are grouped into a single population (called "wild") in the Stacks 
# popmap file. Otherwise, Stacks popmaps assign each wild sample to its source population

# TO DO: add a check to make sure garden individuals are in a population named "garden"

library(adegenet)

# %%%% FUNCTIONS %%%% ----
# Read in relevant functions required for resampling analyses
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)
source("exSituRepresentation/functions_exSituRepresentation.R")

# %%%% QUAC %%%% ----
# ---- SINGLE GENIND (TOGETHER) ----
# R0 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0/"
setwd(genpop.filePath)
QUAC.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R0 Representation rates
reportAllelicRepresentation_Together(QUAC.R0.genind)
reportAllelicRepresentation_Together_Partial(QUAC.R0.genind)
nLoc(QUAC.R0.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80/"
setwd(genpop.filePath)
QUAC.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R80 Representation rates
reportAllelicRepresentation_Together(QUAC.R80.genind)
reportAllelicRepresentation_Together_Partial(QUAC.R80.genind)
nLoc(QUAC.R80.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R100/"
setwd(genpop.filePath)
QUAC.R100.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R100 Representation rates
reportAllelicRepresentation_Together(QUAC.R100.genind)
reportAllelicRepresentation_Together_Partial(QUAC.R100.genind)
nLoc(QUAC.R100.genind)

# ----TWO GENINDS (SEPARATE) ----
# R0 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R0/"
setwd(genpop.filePath)
QUAC.R0.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0/"
setwd(genpop.filePath)
QUAC.R0.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R0 Representation rates
reportAllelicRepresentation_Separate(QUAC.R0.garden.genind, QUAC.R0.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUAC.R0.garden.genind, QUAC.R0.wild.genind)
nLoc(QUAC.R0.garden.genind)
nLoc(QUAC.R0.wild.genind)

# R80 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R80/"
setwd(genpop.filePath)
QUAC.R80.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R80/"
setwd(genpop.filePath)
QUAC.R80.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R80 Representation rates
reportAllelicRepresentation_Separate(QUAC.R80.garden.genind, QUAC.R80.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUAC.R80.garden.genind, QUAC.R80.wild.genind)
nLoc(QUAC.R80.garden.genind)
nLoc(QUAC.R80.wild.genind)

# R100 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R100/"
setwd(genpop.filePath)
QUAC.R100.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R100/"
setwd(genpop.filePath)
QUAC.R100.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R100 Representation rates
reportAllelicRepresentation_Separate(QUAC.R100.garden.genind, QUAC.R100.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUAC.R100.garden.genind, QUAC.R100.wild.genind)
nLoc(QUAC.R100.garden.genind)
nLoc(QUAC.R100.wild.genind)

# ---- NO MINOR ALLELE FREQUENCY (NOMAF) ----
# ---- SINGLE GENIND (TOGETHER) ----
# R0 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.R0_NOMAF.genind)
reportAllelicRepresentation_Together_Partial(QUAC.R0_NOMAF.genind)
nLoc(QUAC.R0_NOMAF.genind); ncol(QUAC.R0_NOMAF.genind@tab)

# R0, ALL SNPS, MULTIPLE WILD POPULATIONS----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUAC.R0_NOMAF_AllSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF_AllSNPs.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.R0_NOMAF_AllSNPs.genind)
reportAllelicRepresentation_Together_Partial(QUAC.R0_NOMAF_AllSNPs.genind)

# R0, FIRST SNP, MULTIPLE WILD POPULATIONS ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP/"
setwd(genpop.filePath)
QUAC.R0_NOMAF_1SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF_1SNP.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.R0_NOMAF_1SNP.genind)
reportAllelicRepresentation_Together_Partial(QUAC.R0_NOMAF_1SNP.genind)

# R0, TWO POPULATIONS----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_TwoPops/"
setwd(genpop.filePath)
QUAC.R0_NOMAF_TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF_TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.R0_NOMAF_TwoPops.genind)
reportAllelicRepresentation_Together_Partial(QUAC.R0_NOMAF_TwoPops.genind)
nLoc(QUAC.R0_NOMAF_TwoPops.genind)

# R0, TWO POPULATIONS, FIRST SNP ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.R0_NOMAF_1SNP.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.R0_NOMAF_1SNP.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.R0_NOMAF_1SNP.TwoPops.genind)
nLoc(QUAC.R0_NOMAF_1SNP.TwoPops.genind)
# Exploration of total and wild allele frequency proportions
# Functions come from Simulated repository (https://github.com/akoontz11/Morton_SSRvSNP_Simulations/blob/main/RScripts/functions_SSRvSNP_Sim.R)
getWildAlleleFreqProportions(QUAC.R0_NOMAF_1SNP.TwoPops.genind)
getTotalAlleleFreqProportions(QUAC.R0_NOMAF_1SNP.TwoPops.genind)

# R0, TWO POPULATIONS, ALL SNPS ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_AllSNPs_2Pops/"
setwd(genpop.filePath)
QUAC.R0_NOMAF_AllSNPs.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.R0_NOMAF_AllSNPs.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.R0_NOMAF_AllSNPs.TwoPops.genind)
nLoc(QUAC.R0_NOMAF_AllSNPs.TwoPops.genind)

# R0, TWO POPULATIONS, HAPLOTYPE-WISE SNP FILTER ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUAC.R0_NOMAF_H.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.R0_NOMAF_H.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.R0_NOMAF_H.TwoPops.genind)
nLoc(QUAC.R0_NOMAF_H.TwoPops.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF/"
setwd(genpop.filePath)
QUAC.R80_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80_NOMAF.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.R80_NOMAF.genind)
reportAllelicRepresentation_Together_Partial(QUAC.R80_NOMAF.genind)
nLoc(QUAC.R80_NOMAF.genind)

# R80, ALL SNPS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUAC.R80_NOMAF_AllSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80_NOMAF_AllSNPs.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.R80_NOMAF_AllSNPs.genind)
reportAllelicRepresentation_Together_Partial(QUAC.R80_NOMAF_AllSNPs.genind)
nLoc(QUAC.R80_NOMAF_AllSNPs.genind)

# R80, TWO POPULATIONS, FIRST SNP ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.R80_NOMAF_1SNP.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.R80_NOMAF_1SNP.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.R80_NOMAF_1SNP.TwoPops.genind)
nLoc(QUAC.R80_NOMAF_1SNP.TwoPops.genind)

# R80, TWO POPULATIONS, ALL SNPS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_AllSNPs_2Pops/"
setwd(genpop.filePath)
QUAC.R80_NOMAF_AllSNPs.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.R80_NOMAF_AllSNPs.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.R80_NOMAF_AllSNPs.TwoPops.genind)
nLoc(QUAC.R80_NOMAF_AllSNPs.TwoPops.genind)

# R80, TWO POPULATIONS, HAPLOTYPE-WISE SNP FILTER ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUAC.R80_NOMAF_H.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.R80_NOMAF_H.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.R80_NOMAF_H.TwoPops.genind)
nLoc(QUAC.R80_NOMAF_H.TwoPops.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R100_NOMAF/"
setwd(genpop.filePath)
QUAC.R100_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100_NOMAF.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R100_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.R100_NOMAF.genind)
reportAllelicRepresentation_Together_Partial(QUAC.R100_NOMAF.genind)
nLoc(QUAC.R100_NOMAF.genind)

# R100, ALL SNPS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R100_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUAC.R100_NOMAF_AllSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100_NOMAF_AllSNPs.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R100_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.R100_NOMAF_AllSNPs.genind)
reportAllelicRepresentation_Together_Partial(QUAC.R100_NOMAF_AllSNPs.genind)
nLoc(QUAC.R100_NOMAF_AllSNPs.genind)

# ----TWO GENINDS (SEPARATE) ----
# R0 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R0_NOMAF/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R0 Representation rates
reportAllelicRepresentation_Separate(QUAC.R0_NOMAF.garden.genind, QUAC.R0_NOMAF.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUAC.R0_NOMAF.garden.genind, QUAC.R0_NOMAF.wild.genind)
nLoc(QUAC.R0_NOMAF.garden.genind)
nLoc(QUAC.R0_NOMAF.wild.genind)

# R0, ALL SNPS ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R0_NOMAF_AllSNPs//"
setwd(genpop.filePath)
QUAC.R0_NOMAF.AllSNPs.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF.AllSNPs.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.AllSNPs.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF.AllSNPs.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R0 Representation rates
reportAllelicRepresentation_Separate(QUAC.R0_NOMAF.AllSNPs.garden.genind, QUAC.R0_NOMAF.AllSNPs.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUAC.R0_NOMAF.AllSNPs.garden.genind, QUAC.R0_NOMAF.AllSNPs.wild.genind)
nLoc(QUAC.R0_NOMAF.AllSNPs.garden.genind)
nLoc(QUAC.R0_NOMAF.AllSNPs.wild.genind)

# R80 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R80_NOMAF/"
setwd(genpop.filePath)
QUAC.R80_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80_NOMAF.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R80_NOMAF/"
setwd(genpop.filePath)
QUAC.R80_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80_NOMAF.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R80 Representation rates
reportAllelicRepresentation_Separate(QUAC.R80_NOMAF.garden.genind, QUAC.R80_NOMAF.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUAC.R80_NOMAF.garden.genind, QUAC.R80_NOMAF.wild.genind)
nLoc(QUAC.R80_NOMAF.garden.genind)
nLoc(QUAC.R80_NOMAF.wild.genind)

# R100 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R100_NOMAF/"
setwd(genpop.filePath)
QUAC.R100_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100_NOMAF.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R100_NOMAF/"
setwd(genpop.filePath)
QUAC.R100_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100_NOMAF.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R100 Representation rates
reportAllelicRepresentation_Separate(QUAC.R100_NOMAF.garden.genind, QUAC.R100_NOMAF.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUAC.R100_NOMAF.garden.genind, QUAC.R100_NOMAF.wild.genind)
nLoc(QUAC.R100_NOMAF.garden.genind)
nLoc(QUAC.R100_NOMAF.wild.genind)

# %%%% QUBO %%%% ----
# ---- SINGLE GENIND (TOGETHER) ----
# R0 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0/"
setwd(genpop.filePath)
QUBO.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R0 Representation rates
reportAllelicRepresentation_Together(QUBO.R0.genind)
reportAllelicRepresentation_Together_Partial(QUBO.R0.genind)
nLoc(QUBO.R0.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80/"
setwd(genpop.filePath)
QUBO.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R80 Representation rates
reportAllelicRepresentation_Together(QUBO.R80.genind)
reportAllelicRepresentation_Together_Partial(QUBO.R80.genind)
nLoc(QUBO.R80.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R100/"
setwd(genpop.filePath)
QUBO.R100.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R100 Representation rates
reportAllelicRepresentation_Together(QUBO.R100.genind)
reportAllelicRepresentation_Together_Partial(QUBO.R100.genind)
nLoc(QUBO.R100.genind)

# ----TWO GENINDS (SEPARATE) ----
# R0 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R0/"
setwd(genpop.filePath)
QUBO.R0.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0/"
setwd(genpop.filePath)
QUBO.R0.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R0 Representation rates
reportAllelicRepresentation_Separate(QUBO.R0.garden.genind, QUBO.R0.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUBO.R0.garden.genind, QUBO.R0.wild.genind)
nLoc(QUBO.R0.garden.genind)
nLoc(QUBO.R0.wild.genind)

# R80 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R80/"
setwd(genpop.filePath)
QUBO.R80.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80/"
setwd(genpop.filePath)
QUBO.R80.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R80 Representation rates
reportAllelicRepresentation_Separate(QUBO.R80.garden.genind, QUBO.R80.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUBO.R80.garden.genind, QUBO.R80.wild.genind)
nLoc(QUBO.R80.garden.genind)
nLoc(QUBO.R80.wild.genind)

# R100 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R100/"
setwd(genpop.filePath)
QUBO.R100.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R100/"
setwd(genpop.filePath)
QUBO.R100.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R100 Representation rates
reportAllelicRepresentation_Separate(QUBO.R100.garden.genind, QUBO.R100.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUBO.R100.garden.genind, QUBO.R100.wild.genind)
nLoc(QUBO.R100.garden.genind)
nLoc(QUBO.R100.wild.genind)

# ---- NO MINOR ALLELE FREQUENCY (NOMAF) ----
# ---- SINGLE GENIND (TOGETHER) ----
# R0 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF/"
setwd(genpop.filePath)
QUBO.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF.genind)
reportAllelicRepresentation_Together_Partial(QUBO.R0_NOMAF.genind)
nLoc(QUBO.R0_NOMAF.genind)

# R0, ALL SNPS----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_AllSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF_AllSNPs.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF_AllSNPs.genind)
reportAllelicRepresentation_Together_Partial(QUBO.R0_NOMAF_AllSNPs.genind)

# R0, FIRST SNP----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_1SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF_1SNP.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF_1SNP.genind)
reportAllelicRepresentation_Together_Partial(QUBO.R0_NOMAF_1SNP.genind)
nLoc(QUBO.R0_NOMAF_1SNP.genind)

# R0, TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_TwoPops/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF_TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF_TwoPops.genind)
reportAllelicRepresentation_Together_Partial(QUBO.R0_NOMAF_TwoPops.genind)
nLoc(QUBO.R0_NOMAF_TwoPops.genind)

# R0, TWO POPULATIONS, FIRST SNP ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_1SNP.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R0_NOMAF_1SNP.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF_1SNP.TwoPops.genind)
nLoc(QUBO.R0_NOMAF_1SNP.TwoPops.genind)

# Exploration of total and wild allele frequency proportions
# Functions come from Simulated repository (https://github.com/akoontz11/Morton_SSRvSNP_Simulations/blob/main/RScripts/functions_SSRvSNP_Sim.R)
getWildAlleleFreqProportions(QUBO.R0_NOMAF_1SNP.TwoPops.genind)
getTotalAlleleFreqProportions(QUBO.R0_NOMAF_1SNP.TwoPops.genind)

# R0, TWO POPULATIONS, ALL SNPS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_AllSNPs_2Pops/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_1SNP.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R0_NOMAF_1SNP.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF_1SNP.TwoPops.genind)
nLoc(QUBO.R0_NOMAF_1SNP.TwoPops.genind)

# R0, TWO POPULATIONS, HAPLOTYPE-WISE SNP FILTER ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_H.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R0_NOMAF_H.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF_H.TwoPops.genind)
nLoc(QUBO.R0_NOMAF_H.TwoPops.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF/"
setwd(genpop.filePath)
QUBO.R80_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80_NOMAF.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R80_NOMAF.genind)
reportAllelicRepresentation_Together_Partial(QUBO.R80_NOMAF.genind)
nLoc(QUBO.R80_NOMAF.genind)

# R80, ALL SNPS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUBO.R80_NOMAF_AllSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R80_NOMAF_AllSNPs.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R80_NOMAF_AllSNPs.genind)
reportAllelicRepresentation_Together_Partial(QUBO.R80_NOMAF_AllSNPs.genind)
nLoc(QUBO.R80_NOMAF_AllSNPs.genind)

# R80, TWO POPULATIONS, FIRST SNP ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.R80_NOMAF_1SNP.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R80_NOMAF_1SNP.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R80_NOMAF_1SNP.TwoPops.genind)
nLoc(QUBO.R80_NOMAF_1SNP.TwoPops.genind)

# R80, TWO POPULATIONS, ALL SNPS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_AllSNPs_2Pops/"
setwd(genpop.filePath)
QUBO.R80_NOMAF_AllSNPs.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R80_NOMAF_AllSNPs.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R80_NOMAF_AllSNPs.TwoPops.genind)
nLoc(QUBO.R80_NOMAF_1SNP.TwoPops.genind)

# R80, TWO POPULATIONS, HAPLOTYPE-WISE SNP FILTER ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUBO.R80_NOMAF_H.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R80_NOMAF_H.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R80_NOMAF_H.TwoPops.genind)
nLoc(QUBO.R80_NOMAF_H.TwoPops.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R100_NOMAF/"
setwd(genpop.filePath)
QUBO.R100_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100_NOMAF.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R100_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R100_NOMAF.genind)
reportAllelicRepresentation_Together_Partial(QUBO.R100_NOMAF.genind)
nLoc(QUBO.R100_NOMAF.genind)

# R100, ALL SNPS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R100_NOMAF/"
setwd(genpop.filePath)
QUBO.R100_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100_NOMAF.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R100_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R100_NOMAF.genind)
reportAllelicRepresentation_Together_Partial(QUBO.R100_NOMAF.genind)

# ---- TWO GENINDS (SEPARATE) ----
# R0 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R0_NOMAF/"
setwd(genpop.filePath)
QUBO.R0_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0_NOMAF/"
setwd(genpop.filePath)
QUBO.R0_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R0_NOMAF Representation rates
reportAllelicRepresentation_Separate(QUBO.R0_NOMAF.garden.genind, QUBO.R0_NOMAF.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUBO.R0_NOMAF.garden.genind, QUBO.R0_NOMAF.wild.genind)
nLoc(QUBO.R0_NOMAF.garden.genind)
nLoc(QUBO.R0_NOMAF.wild.genind)

# R0, ALL SNPS----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R0_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_AllSNPs.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF_AllSNPs.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_AllSNPs.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF_AllSNPs.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R0_NOMAF_AllSNPs Representation rates
reportAllelicRepresentation_Separate(QUBO.R0_NOMAF_AllSNPs.garden.genind, QUBO.R0_NOMAF_AllSNPs.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUBO.R0_NOMAF_AllSNPs.garden.genind, QUBO.R0_NOMAF_AllSNPs.wild.genind)
nLoc(QUBO.R0_NOMAF_AllSNPs.garden.genind)
nLoc(QUBO.R0_NOMAF_AllSNPs.wild.genind)

# R80 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R80_NOMAF/"
setwd(genpop.filePath)
QUBO.R80_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80_NOMAF.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80_NOMAF/"
setwd(genpop.filePath)
QUBO.R80_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80_NOMAF.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R80_NOMAF Representation rates
reportAllelicRepresentation_Separate(QUBO.R80_NOMAF.garden.genind, QUBO.R80_NOMAF.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUBO.R80_NOMAF.garden.genind, QUBO.R80_NOMAF.wild.genind)
nLoc(QUBO.R80_NOMAF.garden.genind)
nLoc(QUBO.R80_NOMAF.wild.genind)

# R100 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R100_NOMAF/"
setwd(genpop.filePath)
QUBO.R100_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100_NOMAF.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R100_NOMAF/"
setwd(genpop.filePath)
QUBO.R100_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100_NOMAF.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R100_NOMAF Representation rates
reportAllelicRepresentation_Separate(QUBO.R100_NOMAF.garden.genind, QUBO.R100_NOMAF.wild.genind)
reportAllelicRepresentation_Separate_Partial(QUBO.R100_NOMAF.garden.genind, QUBO.R100_NOMAF.wild.genind)
nLoc(QUBO.R100_NOMAF.garden.genind)
nLoc(QUBO.R100_NOMAF.wild.genind)
