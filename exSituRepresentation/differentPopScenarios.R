# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% EX SITU CAPTURE RATES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses the functions below to generate ex situ capture values 
# for Quercus acerifolia (QUAC; optimized Stacks de novo assembly, m 7, M/n 4, gt-alpha 0.01) 
# and Quercus boyntonii (QUBO; GSNAP4 reference alignment with Quercus robur reference) NextRAD samples

# Ex situ capture rates (i.e. how well/completely do gardens individuals contain, or capture, wild alleles)
# are impacted by methodology ("approach") and by the filters utilized
# in the Stacks populations module, which can filter the loci or the SNPs included in an analysis.

# %%% METHODOLOGY %%%
# We explored 2 different approaches, each of which had two possible alternatives

# -SINGLE GENIND (TOGETHER): a single genind file containing both garden and wild samples was 
# analyzed (separating samples out by rows). This used the reportAllelicCapture_Together* functions
# -TWO GENINDS (SEPARATE): two genind files, one for garden, and one for wild, were compared 
# by measuring the allele names across each (which were found to be consistent). 
# This used the reportAllelicCapture_Separate* functions
# -Partial: the Together/Separate_Partial functions only look for the beginning section of
# an allele name (which was found to correspond to the locus number, in the locus catalog
# built by gstacks). Functions that do not include "Partial" look for whole alleles names
# (i.e. LocusNumber_SNPposition.SNP)

# %%% FILTERS %%%
# Filters are described below, and code sections are broken out according to the filters used

# -R: corresponds the percentage of individuals a locus needs to be present in, 
# to be included in the analysis (0, 80%, or 100%)
# -NOMAF: means "no minor allele frequency": no threshold for a nucleotide's presence. 
# If an analysis doesn't have NOMAF, then the minor allele frequency filter is 1%
# -AllSNPs: means every SNP is written for every locus. 
# Otherwise, only a single (random) SNP is written per locus

library(adegenet)

# %%%% FUNCTIONS %%%% ----
# Function for reporting allelic capture rates, using a single genind object
reportAllelicCapture_Together <- function(gen.obj){
  # Generate numberical vectors corresponding to garden and wild rows, for later calculations
  garden.Rows <- seq_len(length(which(pop(gen.obj)=="garden")))
  wild.Rows <- seq(from=length(which(pop(gen.obj)=="garden"))+1, to=nInd(gen.obj))
  garden.N <- length(garden.Rows)
  wild.N <- length(wild.Rows)
  # Build the wild allele frequency vector
  wildFreqs <- colSums(gen.obj@tab[wild.Rows,], na.rm = TRUE)/(wild.N*2)*100
  # Calculate capture rates
  # Total
  total <- length(which(names(which(wildFreqs > 0)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 0))*100
  # Very common
  veryCommon <- length(which(names(which(wildFreqs > 10)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 10))*100
  # Common
  common <- length(which(names(which(wildFreqs > 5)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 5))*100
  # Low frequency
  lowFrequency <- length(which(names(which(wildFreqs < 10 & wildFreqs > 1)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs < 10 & wildFreqs > 1))*100
  # Rare
  rare <- length(which(names(which(wildFreqs < 1 & wildFreqs > 0)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs < 1 & wildFreqs > 0))*100
  # Build list of rates
  captureRates <- c(total,veryCommon,common,lowFrequency,rare)
  names(captureRates) <- c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%)","Rare (<1%)")
  # Print capture rates and return
  return(captureRates)
}

# Function for reporting allelic capture rates, using a single genind object
# This version of the function strips the characters following the underscore, from the allele names
# Therefore, only whole RAD loci are being compared across garden and wild individuals
reportAllelicCapture_Together_Partial <- function(gen.obj){
  # Generate numberical vectors corresponding to garden and wild rows, for later calculations
  garden.Rows <- seq_len(length(which(pop(gen.obj)=="garden")))
  wild.Rows <- seq(from=length(which(pop(gen.obj)=="garden"))+1, to=nInd(gen.obj))
  garden.N <- length(garden.Rows)
  wild.N <- length(wild.Rows)
  # Rename wild frequencies, dropping the portion of allele names following the underscore
  colnames(gen.obj@tab) <- gsub(pattern = "_[0-9]{1,4}.[0-9]{1,2}", replacement = "", colnames(gen.obj@tab))
  # Build the wild allele frequency vector
  wildFreqs <- colSums(gen.obj@tab[wild.Rows,], na.rm = TRUE)/(wild.N*2)*100
  # Calculate capture rates
  # Total
  total <- length(which(unique(names(which(wildFreqs > 0))) %in% unique(names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0)))))/length(unique(names(which(wildFreqs > 0))))*100 
  # Very common
  veryCommon <- length(which(unique(names(which(wildFreqs > 10))) %in% unique(names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0)))))/length(unique(names(which(wildFreqs > 10))))*100
  # Common
  common <- length(which(unique(names(which(wildFreqs > 5))) %in% unique(names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0)))))/length(unique(names(which(wildFreqs > 5))))*100
  # Low frequency
  lowFrequency <- length(which(unique(names(which(wildFreqs < 10 & wildFreqs > 1))) %in% unique(names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0)))))/length(unique(names(which(wildFreqs < 10 & wildFreqs > 1))))*100
  # Rare
  rare <- length(which(unique(names(which(wildFreqs < 1 & wildFreqs > 0))) %in% unique(names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0)))))/length(unique(names(which(wildFreqs < 1 & wildFreqs > 0))))*100
  # Build list of rates
  captureRates <- c(total,veryCommon,common,lowFrequency,rare)
  names(captureRates) <- c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%)","Rare (<1%)")
  # Print capture rates and return
  return(captureRates)
}

# Function for reporting allelic capture rates, using two genind objects (one for garden, one for wild)
reportAllelicCapture_Separate <- function(gen.obj.garden, gen.obj.wild){
  # Build the wild allele frequency vector
  wildFreqs <- colSums(gen.obj.wild@tab, na.rm = TRUE)/(nInd(gen.obj.wild)*2)*100
  # Calculate capture rates
  # Total
  total <- length(which(names(which(wildFreqs > 0)) %in% colnames(gen.obj.garden@tab)))/length(which(wildFreqs > 0))*100
  # Very common
  veryCommon <- length(which(names(which(wildFreqs > 10)) %in% colnames(gen.obj.garden@tab)))/length(which(wildFreqs > 10))*100
  # Common
  common <- length(which(names(which(wildFreqs > 5)) %in% colnames(gen.obj.garden@tab)))/length(which(wildFreqs > 5))*100
  # Low frequency
  lowFrequency <- length(which(names(which(wildFreqs < 10 & wildFreqs > 1)) %in% colnames(gen.obj.garden@tab)))/length(which(wildFreqs < 10 & wildFreqs > 1))*100
  # Rare
  rare <- length(which(names(which(wildFreqs < 1 & wildFreqs > 0)) %in% colnames(gen.obj.garden@tab)))/length(which(wildFreqs < 1 & wildFreqs > 0))*100
  # Build list of rates
  captureRates <- c(total,veryCommon,common,lowFrequency,rare)
  names(captureRates) <- c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%)","Rare (<1%)")
  # Print capture rates and return
  return(captureRates)
}

# Function for reporting allelic capture rates, using two genind objects (one for garden, one for wild)
# This version of the function strips the characters following the underscore, from the allele names
# Therefore, only whole RAD loci are being compared across garden and wild individuals
reportAllelicCapture_Separate_Partial <- function(gen.obj.garden, gen.obj.wild){
  # Build the wild allele frequency vector
  wildFreqs <- colSums(gen.obj.wild@tab, na.rm = TRUE)/(nInd(gen.obj.wild)*2)*100
  # Rename wild frequencies, dropping the portion of allele names following the underscore
  names(wildFreqs) <- gsub(pattern = "_[0-9]{1,4}.[0-9]{1,2}", replacement = "", names(wildFreqs))
  # Rename garden colnames, dropping the portion of allele names following the underscore
  colnames(gen.obj.garden@tab) <- gsub(pattern = "_[0-9]{1,4}.[0-9]{1,2}", replacement = "", colnames(gen.obj.garden@tab))
  # Calculate capture rates
  # Total
  total <- length(which(unique(names(which(wildFreqs > 0))) %in% unique(colnames(gen.obj.garden@tab))))/length(unique(names(which(wildFreqs > 0))))*100
  # Very common
  veryCommon <- length(which(unique(names(which(wildFreqs > 10))) %in% unique(colnames(gen.obj.garden@tab))))/length(unique(names(which(wildFreqs > 10))))*100
  # Common
  common <- length(which(unique(names(which(wildFreqs > 5))) %in% unique(colnames(gen.obj.garden@tab))))/length(unique(names(which(wildFreqs > 5))))*100
  # Low frequency
  lowFrequency <- length(which(unique(names(which(wildFreqs < 10 & wildFreqs > 1))) %in% unique(colnames(gen.obj.garden@tab))))/length(unique(names(which(wildFreqs < 10 & wildFreqs > 1))))*100
  # Rare
  rare <- length(which(unique(names(which(wildFreqs < 1 & wildFreqs > 0))) %in% unique(colnames(gen.obj.garden@tab))))/length(unique(names(which(wildFreqs < 1 & wildFreqs > 0))))*100
  # Build list of rates
  captureRates <- c(total,veryCommon,common,lowFrequency,rare)
  names(captureRates) <- c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%)","Rare (<1%)")
  # Print capture rates and return
  return(captureRates)
}

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
# R0 Capture rates
reportAllelicCapture_Together(QUAC.R0.genind)
reportAllelicCapture_Together_Partial(QUAC.R0.genind)
nLoc(QUAC.R0.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80/"
setwd(genpop.filePath)
QUAC.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R80 Capture rates
reportAllelicCapture_Together(QUAC.R80.genind)
reportAllelicCapture_Together_Partial(QUAC.R80.genind)
nLoc(QUAC.R80.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R100/"
setwd(genpop.filePath)
QUAC.R100.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R100 Capture rates
reportAllelicCapture_Together(QUAC.R100.genind)
reportAllelicCapture_Together_Partial(QUAC.R100.genind)
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

# R0 Capture rates
reportAllelicCapture_Separate(QUAC.R0.garden.genind, QUAC.R0.wild.genind)
reportAllelicCapture_Separate_Partial(QUAC.R0.garden.genind, QUAC.R0.wild.genind)
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

# R80 Capture rates
reportAllelicCapture_Separate(QUAC.R80.garden.genind, QUAC.R80.wild.genind)
reportAllelicCapture_Separate_Partial(QUAC.R80.garden.genind, QUAC.R80.wild.genind)
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

# R100 Capture rates
reportAllelicCapture_Separate(QUAC.R100.garden.genind, QUAC.R100.wild.genind)
reportAllelicCapture_Separate_Partial(QUAC.R100.garden.genind, QUAC.R100.wild.genind)
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
# R0_NOMAF Capture rates
reportAllelicCapture_Together(QUAC.R0_NOMAF.genind)
reportAllelicCapture_Together_Partial(QUAC.R0_NOMAF.genind)
nLoc(QUAC.R0_NOMAF.genind)

# R0, ALL SNPS----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUAC.R0_NOMAF_AllSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF_AllSNPs.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R0_NOMAF Capture rates
reportAllelicCapture_Together(QUAC.R0_NOMAF_AllSNPs.genind)
reportAllelicCapture_Together_Partial(QUAC.R0_NOMAF_AllSNPs.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF/"
setwd(genpop.filePath)
QUAC.R80_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80_NOMAF.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R80_NOMAF Capture rates
reportAllelicCapture_Together(QUAC.R80_NOMAF.genind)
reportAllelicCapture_Together_Partial(QUAC.R80_NOMAF.genind)
nLoc(QUAC.R80_NOMAF.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R100_NOMAF/"
setwd(genpop.filePath)
QUAC.R100_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100_NOMAF.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R100_NOMAF Capture rates
reportAllelicCapture_Together(QUAC.R100_NOMAF.genind)
reportAllelicCapture_Together_Partial(QUAC.R100_NOMAF.genind)
nLoc(QUAC.R100_NOMAF.genind)

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

# R0 Capture rates
reportAllelicCapture_Separate(QUAC.R0_NOMAF.garden.genind, QUAC.R0_NOMAF.wild.genind)
reportAllelicCapture_Separate_Partial(QUAC.R0_NOMAF.garden.genind, QUAC.R0_NOMAF.wild.genind)
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

# R0 Capture rates
reportAllelicCapture_Separate(QUAC.R0_NOMAF.AllSNPs.garden.genind, QUAC.R0_NOMAF.AllSNPs.wild.genind)
reportAllelicCapture_Separate_Partial(QUAC.R0_NOMAF.AllSNPs.garden.genind, QUAC.R0_NOMAF.AllSNPs.wild.genind)
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

# R80 Capture rates
reportAllelicCapture_Separate(QUAC.R80_NOMAF.garden.genind, QUAC.R80_NOMAF.wild.genind)
reportAllelicCapture_Separate_Partial(QUAC.R80_NOMAF.garden.genind, QUAC.R80_NOMAF.wild.genind)
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

# R100 Capture rates
reportAllelicCapture_Separate(QUAC.R100_NOMAF.garden.genind, QUAC.R100_NOMAF.wild.genind)
reportAllelicCapture_Separate_Partial(QUAC.R100_NOMAF.garden.genind, QUAC.R100_NOMAF.wild.genind)
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
# R0 Capture rates
reportAllelicCapture_Together(QUBO.R0.genind)
reportAllelicCapture_Together_Partial(QUBO.R0.genind)
nLoc(QUBO.R0.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80/"
setwd(genpop.filePath)
QUBO.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R80 Capture rates
reportAllelicCapture_Together(QUBO.R80.genind)
reportAllelicCapture_Together_Partial(QUBO.R80.genind)
nLoc(QUBO.R80.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R100/"
setwd(genpop.filePath)
QUBO.R100.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R100 Capture rates
reportAllelicCapture_Together(QUBO.R100.genind)
reportAllelicCapture_Together_Partial(QUBO.R100.genind)
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

# R0 Capture rates
reportAllelicCapture_Separate(QUBO.R0.garden.genind, QUBO.R0.wild.genind)
reportAllelicCapture_Separate_Partial(QUBO.R0.garden.genind, QUBO.R0.wild.genind)
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

# R80 Capture rates
reportAllelicCapture_Separate(QUBO.R80.garden.genind, QUBO.R80.wild.genind)
reportAllelicCapture_Separate_Partial(QUBO.R80.garden.genind, QUBO.R80.wild.genind)
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

# R100 Capture rates
reportAllelicCapture_Separate(QUBO.R100.garden.genind, QUBO.R100.wild.genind)
reportAllelicCapture_Separate_Partial(QUBO.R100.garden.genind, QUBO.R100.wild.genind)
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
# R0_NOMAF Capture rates
reportAllelicCapture_Together(QUBO.R0_NOMAF.genind)
reportAllelicCapture_Together_Partial(QUBO.R0_NOMAF.genind)
nLoc(QUBO.R0_NOMAF.genind)

# R0, ALL SNPS----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_AllSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF_AllSNPs.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R0_NOMAF Capture rates
reportAllelicCapture_Together(QUBO.R0_NOMAF_AllSNPs.genind)
reportAllelicCapture_Together_Partial(QUBO.R0_NOMAF_AllSNPs.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF/"
setwd(genpop.filePath)
QUBO.R80_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80_NOMAF.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R80_NOMAF Capture rates
reportAllelicCapture_Together(QUBO.R80_NOMAF.genind)
reportAllelicCapture_Together_Partial(QUBO.R80_NOMAF.genind)
nLoc(QUBO.R80_NOMAF.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R100_NOMAF/"
setwd(genpop.filePath)
QUBO.R100_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100_NOMAF.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R100_NOMAF Capture rates
reportAllelicCapture_Together(QUBO.R100_NOMAF.genind)
reportAllelicCapture_Together_Partial(QUBO.R100_NOMAF.genind)
nLoc(QUBO.R100_NOMAF.genind)

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

# R0_NOMAF Capture rates
reportAllelicCapture_Separate(QUBO.R0_NOMAF.garden.genind, QUBO.R0_NOMAF.wild.genind)
reportAllelicCapture_Separate_Partial(QUBO.R0_NOMAF.garden.genind, QUBO.R0_NOMAF.wild.genind)
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

# R0_NOMAF_AllSNPs Capture rates
reportAllelicCapture_Separate(QUBO.R0_NOMAF_AllSNPs.garden.genind, QUBO.R0_NOMAF_AllSNPs.wild.genind)
reportAllelicCapture_Separate_Partial(QUBO.R0_NOMAF_AllSNPs.garden.genind, QUBO.R0_NOMAF_AllSNPs.wild.genind)
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

# R80_NOMAF Capture rates
reportAllelicCapture_Separate(QUBO.R80_NOMAF.garden.genind, QUBO.R80_NOMAF.wild.genind)
reportAllelicCapture_Separate_Partial(QUBO.R80_NOMAF.garden.genind, QUBO.R80_NOMAF.wild.genind)
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

# R100_NOMAF Capture rates
reportAllelicCapture_Separate(QUBO.R100_NOMAF.garden.genind, QUBO.R100_NOMAF.wild.genind)
reportAllelicCapture_Separate_Partial(QUBO.R100_NOMAF.garden.genind, QUBO.R100_NOMAF.wild.genind)
nLoc(QUBO.R100_NOMAF.garden.genind)
nLoc(QUBO.R100_NOMAF.wild.genind)
