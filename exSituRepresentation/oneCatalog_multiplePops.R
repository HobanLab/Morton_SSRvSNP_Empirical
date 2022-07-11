# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% ONE CATALOG, MULTIPLE POPULATIONS MODULES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Comparison of wild and garden alleles separated using populations
# In this example, the loci catalog is constructed using both wild and garden individuals
# These individuals are separated out (via the populations module), and then the resulting 
# genind objects are compared, to see what alleles overlap

# This script is exploratory: it was used to form what's been termed the "Separate"
# approach for measuring ex situ capture. For more formal, organized ex situ capture results
# (with both the "Separate" and "Together" approaches, and with the impacts of different filters)
# check out the exSituCaptureRates.R script

# Filepaths refer to filtered reference dataset, aligned using GSNAP, with -R 80 filtering
library(adegenet)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% REFERENCE -- FILTERED, GSNAP %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# QUAC----
# Wild
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_wild/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
wild.QUAC.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(wild.QUAC.genind) <- factor(read.table("../../../QUAC_popmap_wild", header=FALSE)[,2])
nLoc(wild.QUAC.genind); ncol(wild.QUAC.genind@tab)

# Garden
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_garden/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
garden.QUAC.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(garden.QUAC.genind) <- factor(read.table("../../../QUAC_popmap_garden", header=FALSE)[,2])
nLoc(garden.QUAC.genind); ncol(garden.QUAC.genind@tab)

# Vector of TRUE/FALSE values, for whether wild allele names is seen in vector of garden allele names
colnames(wild.QUAC.genind@tab) %in% colnames(garden.QUAC.genind@tab)
length(colnames(wild.QUAC.genind@tab) %in% colnames(garden.QUAC.genind@tab)) 
length(which(colnames(wild.QUAC.genind@tab) %in% colnames(garden.QUAC.genind@tab))) # 1,780/15,972

# Wild allele frequency vector
QUAC_wildFreqs <- colSums(wild.QUAC.genind@tab, na.rm = TRUE)/(nInd(wild.QUAC.genind)*2)*100

# Total
length(which(names(which(QUAC_wildFreqs > 0)) %in% colnames(garden.QUAC.genind@tab)))/length(which(QUAC_wildFreqs > 0))*100
# Very common
length(which(names(which(QUAC_wildFreqs > 10)) %in% colnames(garden.QUAC.genind@tab)))/length(which(QUAC_wildFreqs > 10))*100
# Common
length(which(names(which(QUAC_wildFreqs > 5)) %in% colnames(garden.QUAC.genind@tab)))/length(which(QUAC_wildFreqs > 5))*100
# Low frequency
length(which(names(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1)) %in% colnames(garden.QUAC.genind@tab)))/length(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1))*100
# Rare
length(which(names(which(QUAC_wildFreqs < 1)) %in% colnames(garden.QUAC.genind@tab)))/length(which(QUAC_wildFreqs < 1))*100

# QUBO----
# Wild
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP/output/populations_wild/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
wild.QUBO.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(wild.QUBO.genind) <- factor(read.table("../../../QUBO_popmap_wild", header=FALSE)[,2])
nLoc(wild.QUBO.genind); ncol(wild.QUBO.genind@tab)

# Garden
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP/output/populations_garden/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
garden.QUBO.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(garden.QUBO.genind) <- factor(read.table("../../../QUBO_popmap_garden", header=FALSE)[,2])
nLoc(garden.QUBO.genind); ncol(garden.QUBO.genind@tab)

# Vector of TRUE/FALSE values, for whether wild allele names is seen in vector of garden allele names
colnames(wild.QUBO.genind@tab) %in% colnames(garden.QUBO.genind@tab)
length(colnames(wild.QUBO.genind@tab) %in% colnames(garden.QUBO.genind@tab)) 
length(which(colnames(wild.QUBO.genind@tab) %in% colnames(garden.QUBO.genind@tab))) # 2,598/33,070

# Wild allele frequency vector
QUBO_wildFreqs <- colSums(wild.QUBO.genind@tab, na.rm = TRUE)/(nInd(wild.QUBO.genind)*2)*100

# Total
length(which(names(which(QUBO_wildFreqs > 0)) %in% colnames(garden.QUBO.genind@tab)))/length(which(QUBO_wildFreqs > 0))*100
# Very common
length(which(names(which(QUBO_wildFreqs > 10)) %in% colnames(garden.QUBO.genind@tab)))/length(which(QUBO_wildFreqs > 10))*100
# Common
length(which(names(which(QUBO_wildFreqs > 5)) %in% colnames(garden.QUBO.genind@tab)))/length(which(QUBO_wildFreqs > 5))*100
# Low frequency
length(which(names(which(QUBO_wildFreqs < 10 & QUBO_wildFreqs > 1)) %in% colnames(garden.QUBO.genind@tab)))/length(which(QUBO_wildFreqs < 10 & QUBO_wildFreqs > 1))*100
# Rare
length(which(names(which(QUBO_wildFreqs < 1)) %in% colnames(garden.QUBO.genind@tab)))/length(which(QUBO_wildFreqs < 1))*100

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% DE NOVO -- FINAL ASSEMBLIES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# QUAC----
# Wild
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
wild.QUAC.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(wild.QUAC.genind) <- factor(read.table("../../QUAC_popmap_wild", header=FALSE)[,2])
nLoc(wild.QUAC.genind); ncol(wild.QUAC.genind@tab)

# Garden
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
garden.QUAC.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(garden.QUAC.genind) <- factor(read.table("../../QUAC_popmap_garden", header=FALSE)[,2])
nLoc(garden.QUAC.genind); ncol(garden.QUAC.genind@tab)

# Vector of TRUE/FALSE values, for whether wild allele names is seen in vector of garden allele names
colnames(wild.QUAC.genind@tab) %in% colnames(garden.QUAC.genind@tab)
length(colnames(wild.QUAC.genind@tab) %in% colnames(garden.QUAC.genind@tab)) 
length(which(colnames(wild.QUAC.genind@tab) %in% colnames(garden.QUAC.genind@tab))) # 964/12,098

# Wild allele frequency vector
QUAC_wildFreqs <- colSums(wild.QUAC.genind@tab, na.rm = TRUE)/(nInd(wild.QUAC.genind)*2)*100

# Total
length(which(names(which(QUAC_wildFreqs > 0)) %in% colnames(garden.QUAC.genind@tab)))/length(which(QUAC_wildFreqs > 0))*100
# Very common
length(which(names(which(QUAC_wildFreqs > 10)) %in% colnames(garden.QUAC.genind@tab)))/length(which(QUAC_wildFreqs > 10))*100
# Common
length(which(names(which(QUAC_wildFreqs > 5)) %in% colnames(garden.QUAC.genind@tab)))/length(which(QUAC_wildFreqs > 5))*100
# Low frequency
length(which(names(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1)) %in% colnames(garden.QUAC.genind@tab)))/length(which(QUAC_wildFreqs < 10 & QUAC_wildFreqs > 1))*100
# Rare
length(which(names(which(QUAC_wildFreqs < 1)) %in% colnames(garden.QUAC.genind@tab)))/length(which(QUAC_wildFreqs < 1))*100

# QUBO----
# Wild
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_wild/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
wild.QUBO.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(wild.QUBO.genind) <- factor(read.table("../../QUBO_popmap_wild", header=FALSE)[,2])
nLoc(wild.QUBO.genind); ncol(wild.QUBO.genind@tab)

# Garden
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_garden/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
garden.QUBO.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(garden.QUBO.genind) <- factor(read.table("../../QUBO_popmap_garden", header=FALSE)[,2])
nLoc(garden.QUBO.genind); ncol(garden.QUBO.genind@tab)

# Vector of TRUE/FALSE values, for whether wild allele names is seen in vector of garden allele names
colnames(wild.QUBO.genind@tab) %in% colnames(garden.QUBO.genind@tab)
length(colnames(wild.QUBO.genind@tab) %in% colnames(garden.QUBO.genind@tab)) 
length(which(colnames(wild.QUBO.genind@tab) %in% colnames(garden.QUBO.genind@tab))) # 818/10,122

# Wild allele frequency vector
QUBO_wildFreqs <- colSums(wild.QUBO.genind@tab, na.rm = TRUE)/(nInd(wild.QUBO.genind)*2)*100

# Total
length(which(names(which(QUBO_wildFreqs > 0)) %in% colnames(garden.QUBO.genind@tab)))/length(which(QUBO_wildFreqs > 0))*100
# Very common
length(which(names(which(QUBO_wildFreqs > 10)) %in% colnames(garden.QUBO.genind@tab)))/length(which(QUBO_wildFreqs > 10))*100
# Common
length(which(names(which(QUBO_wildFreqs > 5)) %in% colnames(garden.QUBO.genind@tab)))/length(which(QUBO_wildFreqs > 5))*100
# Low frequency
length(which(names(which(QUBO_wildFreqs < 10 & QUBO_wildFreqs > 1)) %in% colnames(garden.QUBO.genind@tab)))/length(which(QUBO_wildFreqs < 10 & QUBO_wildFreqs > 1))*100
# Rare
length(which(names(which(QUBO_wildFreqs < 1)) %in% colnames(garden.QUBO.genind@tab)))/length(which(QUBO_wildFreqs < 1))*100

# The values we're getting here are very low. However, these values don't reflect the ACTUAL
# capture of wild alleles, because the genotyping/SNP calling is still happening across ALL
# samples. What we need is to consider *just* the wild samples, genotype those, then *just* 
# the garden samples, genotype those, and then compare how many wild alleles are capture in
# the gardens

# The problem with this, though, is that it's going to lead to *extremely* low capture rates, 
# AND it's not even clear the loci will be comparable across the two datasets
