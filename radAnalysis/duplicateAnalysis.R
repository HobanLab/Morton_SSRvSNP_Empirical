# %%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% DUPLICATE ANALYSIS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%

# This script compares the number of loci and alleles present in genind objects
# generated with duplicates and without duplicates, in order to determine how many
# loci are generated as the result of random sequencing error.

# For both Quercus acerifolia (QUAC) and Quercus boyntonii (QUBO), 
# it reads in datasets from optimized de novo assemblies built using Stacks 
# (QUAC: m 7, M/n 5, gt-alpha 0.01; QUBO: m 7, M/n 5, gt-alpha 0.01)
# and datasets from reference alignments ("GSNAP4" alignment parameters; QUAC: Q. rubra genome;
# QUBO: Q. robur genome). Finally, for both de novo and reference datasets, it reads in loci
# unfiltered according to missing data (R0) and loci shared among at least 80% (R80) of all samples (garden and wild)

library(adegenet)

# %%%% QUAC %%%% ----
# ---- SNPS, DE NOVO ----
# R0 ----
# WITHOUT DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.DN.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R0.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])

# WITH DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops_wDup/"
setwd(genpop.filePath)
QUAC.SNP.DN.R0.wDup.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R0.wDup.genind) <- factor(read.table("QUAC_popmap_GardenWild_wDup", header=FALSE)[,2])

# Calculate differences in number of loci/alleles
nLoc(QUAC.SNP.DN.R0.wDup.genind) - nLoc(QUAC.SNP.DN.R0.genind)
ncol(QUAC.SNP.DN.R0.wDup.genind@tab) - ncol(QUAC.SNP.DN.R0.genind@tab)
# Calculate proportion of duplicate loci 
# (number of duplicate loci/number of loci in datasets with duplicates)
print(paste0("%%% QUAC: DN: R0: Duplicate loci percentage: ",
  ((nLoc(QUAC.SNP.DN.R0.wDup.genind) - nLoc(QUAC.SNP.DN.R0.genind))/nLoc(QUAC.SNP.DN.R0.wDup.genind))*100, "%"))

# R80 ----
# WITHOUT DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R80.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])

# WITH DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops_wDup/"
setwd(genpop.filePath)
QUAC.SNP.DN.R80.wDup.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R80.wDup.genind) <- factor(read.table("QUAC_popmap_GardenWild_wDup", header=FALSE)[,2])

# Calculate differences in number of loci/alleles
nLoc(QUAC.SNP.DN.R80.wDup.genind) - nLoc(QUAC.SNP.DN.R80.genind)
ncol(QUAC.SNP.DN.R80.wDup.genind@tab) - ncol(QUAC.SNP.DN.R80.genind@tab)
# Calculate proportion of duplicate loci 
# (number of duplicate loci/number of loci in datasets with duplicates)
print(paste0("%%% QUAC: DN: R80: Duplicate loci percentage: ",
             ((nLoc(QUAC.SNP.DN.R80.wDup.genind) - nLoc(QUAC.SNP.DN.R80.genind))/nLoc(QUAC.SNP.DN.R80.wDup.genind))*100, "%"))

# ---- SNPS, REFERENCE ----
# R0 ----
# WITHOUT DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.REF.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R0.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])

# WITH DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R0_NOMAF_1SNP_2Pops_wDup/"
setwd(genpop.filePath)
QUAC.SNP.REF.R0.wDup.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R0.wDup.genind) <- factor(read.table("QUAC_popmap_GardenWild_wDup", header=FALSE)[,2])

# Calculate differences in number of loci/alleles
nLoc(QUAC.SNP.REF.R0.wDup.genind) - nLoc(QUAC.SNP.REF.R0.genind)
ncol(QUAC.SNP.REF.R0.wDup.genind@tab) - ncol(QUAC.SNP.REF.R0.genind@tab)
# Calculate proportion of duplicate loci 
# (number of duplicate loci/number of loci in datasets with duplicates)
print(paste0("%%% QUAC: REF: R0: Duplicate loci percentage: ",
             ((nLoc(QUAC.SNP.REF.R0.wDup.genind) - nLoc(QUAC.SNP.REF.R0.genind))/nLoc(QUAC.SNP.REF.R0.wDup.genind))*100, "%"))

# R80 ----
# WITHOUT DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R80.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])

# WITH DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R80_NOMAF_1SNP_2Pops_wDup/"
setwd(genpop.filePath)
QUAC.SNP.REF.R80.wDup.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R80.wDup.genind) <- factor(read.table("QUAC_popmap_GardenWild_wDup", header=FALSE)[,2])

# Calculate differences in number of loci/alleles
nLoc(QUAC.SNP.REF.R80.wDup.genind) - nLoc(QUAC.SNP.REF.R80.genind)
ncol(QUAC.SNP.REF.R80.wDup.genind@tab) - ncol(QUAC.SNP.REF.R80.genind@tab)
# Calculate proportion of duplicate loci 
# (number of duplicate loci/number of loci in datasets with duplicates)
print(paste0("%%% QUAC: REF: R80: Duplicate loci percentage: ",
             ((nLoc(QUAC.SNP.REF.R80.wDup.genind) - nLoc(QUAC.SNP.REF.R80.genind))/nLoc(QUAC.SNP.REF.R80.wDup.genind))*100, "%"))

# %%%% QUBO %%%% ----
# ---- SNPS: DE NOVO ASSEMBLY ----
# R0 ----
# WITHOUT DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.DN.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R0.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# WITH DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R0_NOMAF_1SNP_2Pops_wDup/"
setwd(genpop.filePath)
QUBO.SNP.DN.R0.wDup.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R0.wDup.genind) <- factor(read.table("QUBO_popmap_GardenWild_wDup", header=FALSE)[,2])

# Calculate differences in number of loci/alleles
nLoc(QUBO.SNP.DN.R0.wDup.genind) - nLoc(QUBO.SNP.DN.R0.genind)
ncol(QUBO.SNP.DN.R0.wDup.genind@tab) - ncol(QUBO.SNP.DN.R0.genind@tab)
# Calculate proportion of duplicate loci 
# (number of duplicate loci/number of loci in datasets with duplicates)
print(paste0("%%% QUBO: DN: R0: Duplicate loci percentage: ",
             ((nLoc(QUBO.SNP.DN.R0.wDup.genind) - nLoc(QUBO.SNP.DN.R0.genind))/nLoc(QUBO.SNP.DN.R0.wDup.genind))*100, "%"))

# R80 ----
# WITHOUT DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R80.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# WITH DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R80_NOMAF_1SNP_2Pops_wDup/"
setwd(genpop.filePath)
QUBO.SNP.DN.R80.wDup.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R80.wDup.genind) <- factor(read.table("QUBO_popmap_GardenWild_wDup", header=FALSE)[,2])

# Calculate differences in number of loci/alleles
nLoc(QUBO.SNP.DN.R80.wDup.genind) - nLoc(QUBO.SNP.DN.R80.genind)
ncol(QUBO.SNP.DN.R80.wDup.genind@tab) - ncol(QUBO.SNP.DN.R80.genind@tab)
# Calculate proportion of duplicate loci 
# (number of duplicate loci/number of loci in datasets with duplicates)
print(paste0("%%% QUBO: DN: R80: Duplicate loci percentage: ",
             ((nLoc(QUBO.SNP.DN.R80.wDup.genind) - nLoc(QUBO.SNP.DN.R80.genind))/nLoc(QUBO.SNP.DN.R80.wDup.genind))*100, "%"))

# ---- SNPS: REFERENCE ----
# R0 ----
# WITHOUT DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.REF.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R0.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# WITH DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops_wDup/"
setwd(genpop.filePath)
QUBO.SNP.REF.R0.wDup.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R0.wDup.genind) <- factor(read.table("QUBO_popmap_GardenWild_wDup", header=FALSE)[,2])

# Calculate differences in number of loci/alleles
nLoc(QUBO.SNP.REF.R0.wDup.genind) - nLoc(QUBO.SNP.REF.R0.genind)
ncol(QUBO.SNP.REF.R0.wDup.genind@tab) - ncol(QUBO.SNP.REF.R0.genind@tab)
# Calculate proportion of duplicate loci 
# (number of duplicate loci/number of loci in datasets with duplicates)
print(paste0("%%% QUBO: REF: R0: Duplicate loci percentage: ",
             ((nLoc(QUBO.SNP.REF.R0.wDup.genind) - nLoc(QUBO.SNP.REF.R0.genind))/nLoc(QUBO.SNP.REF.R0.wDup.genind))*100, "%"))

# R80 ----
# WITHOUT DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R80.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# WITH DUPLICATES
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_1SNP_2Pops_wDup/"
setwd(genpop.filePath)
QUBO.SNP.REF.R80.wDup.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R80.wDup.genind) <- factor(read.table("QUBO_popmap_GardenWild_wDup", header=FALSE)[,2])

# Calculate differences in number of loci/alleles
nLoc(QUBO.SNP.REF.R80.wDup.genind) - nLoc(QUBO.SNP.REF.R80.genind)
ncol(QUBO.SNP.REF.R80.wDup.genind@tab) - ncol(QUBO.SNP.REF.R80.genind@tab)
# Calculate proportion of duplicate loci 
# (number of duplicate loci/number of loci in datasets with duplicates)
print(paste0("%%% QUBO: REF: R80: Duplicate loci percentage: ",
             ((nLoc(QUBO.SNP.REF.R80.wDup.genind) - nLoc(QUBO.SNP.REF.R80.genind))/nLoc(QUBO.SNP.REF.R80.wDup.genind))*100, "%"))
