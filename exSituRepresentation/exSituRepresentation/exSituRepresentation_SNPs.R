# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% EX SITU REPRESENTATION RATES: SNP %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses the functions declared in the functions_exSituRepresentation.R file to generate ex situ representation values 
# for Quercus acerifolia (QUAC, two datasets: optimized Stacks de novo assembly, 
# m 7, M/n 4, gt-alpha 0.01, and GSNAP4 alignment with Quercus rubra reference) 
# and Quercus boyntonii (QUBO; GSNAP4 alignment with Quercus robur reference) samples

# Many scenarios/methodologies/filters were explored when conducting ex situ representation analyses; to see 
# all of these, check out the exSituRepresentation_SNPs_Expanded.R script. The scenarios contained in this script
# are those we choose to include in our final Results, in one form or another. Some of the justifications for using
# these particular scenarios, out of all of those explores in the Expanded script, are offered below:

# -SINGLE GENIND OBJECT (TOGETHER): we decided that processing garden and wild samples together (creating a single genind)
# had a minimal impact on the SNP genotyping algorithm used by Stacks.
# -COMPLTE SNP LOCI: because we're interested in SNP values and positions in a RAD locus, we use complete locus names,
# rather than the _Partial analyses.
# NO MINOR ALLELE FREQUENCY (NOMAF): we found that including a minor allele frequency greatly impacted our rates of 
# gardens representing rare Wild alleles. Therefore, we turned of minor allele frequencies for the analyses we report.
# FIRST SNP PER LOCUS (1SNP): in order to maintain consistency across separate Stacks runs, we avoided writing a random SNP
# to each locus. Writing all SNPs to each locus led to computationally intensive datasets, and did not greatly impact 
# ex situ conservation (although some _AllSNPs scenarios are included here).

# Additional filters are described below, and code sections are broken out according to the filters used.
# %%% FILTERS %%%
# -R: corresponds the percentage of individuals a locus needs to be present in, 
# to be included in the analysis (0 or 80%)
# -AllSNPs: means every SNP is written for every locus.
# -1SNP: means only the first SNP is written for every locus
# -H: SNPs are filtered haplotype-wise (unshared SNPs are pruned to reduce haplotype-wise missing data).
# -TwoPops: means all wild samples are grouped into a single population (called "wild") in the Stacks 
# popmap file. Otherwise, Stacks popmaps assign each wild sample to its source population

library(adegenet)

# %%%% FUNCTIONS %%%% ----
# Read in relevant functions required for resampling analyses
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)
source("exSituRepresentation/functions_exSituRepresentation.R")

# %%%% QUAC %%%% ----
# %%%% DE NOVO %%%% ----
# %%%% R0 ----
# ---- ALL SNPS ----
# MULTIPLE WILD POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUAC.DN.R0_NOMAF_AllSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.DN.R0_NOMAF_AllSNPs.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.DN.R0_NOMAF_AllSNPs.genind)

# TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_AllSNPs_2Pops/"
setwd(genpop.filePath)
QUAC.DN.R0_NOMAF_AllSNPs.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.DN.R0_NOMAF_AllSNPs.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.DN.R0_NOMAF_AllSNPs.TwoPops.genind)

# ---- FIRST SNP ----
# MULTIPLE WILD POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP/"
setwd(genpop.filePath)
QUAC.DN.R0_NOMAF_1SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.DN.R0_NOMAF_1SNP.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.DN.R0_NOMAF_1SNP.genind)

# TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.DN.R0_NOMAF_1SNP.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.DN.R0_NOMAF_1SNP.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.DN.R0_NOMAF_1SNP.TwoPops.genind)
# Exploration of total and wild allele frequency proportions
# Functions come from Simulated repository 
# (https://github.com/akoontz11/Morton_SSRvSNP_Simulations/blob/main/RScripts/functions_SSRvSNP_Sim.R)
getWildAlleleFreqProportions(QUAC.DN.R0_NOMAF_1SNP.TwoPops.genind)
getTotalAlleleFreqProportions(QUAC.DN.R0_NOMAF_1SNP.TwoPops.genind)

# ---- HAPLOTYPE-WISE SNP FILTER ----
# TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUAC.DN.R0_NOMAF_H.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.DN.R0_NOMAF_H.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.DN.R0_NOMAF_H.TwoPops.genind)

# %%%% R80 ----
# ALL SNPS, TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_AllSNPs_2Pops/"
setwd(genpop.filePath)
QUAC.DN.R80_NOMAF_AllSNPs.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.DN.R80_NOMAF_AllSNPs.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.DN.R80_NOMAF_AllSNPs.TwoPops.genind)

# FIRST SNP, TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.DN.R80_NOMAF_1SNP.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.DN.R80_NOMAF_1SNP.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.DN.R80_NOMAF_1SNP.TwoPops.genind)
# Exploration of total and wild allele frequency proportions
# Functions come from Simulated repository 
# (https://github.com/akoontz11/Morton_SSRvSNP_Simulations/blob/main/RScripts/functions_SSRvSNP_Sim.R)
getWildAlleleFreqProportions(QUAC.DN.R80_NOMAF_1SNP.TwoPops.genind)
getTotalAlleleFreqProportions(QUAC.DN.R80_NOMAF_1SNP.TwoPops.genind)

# HAPLOTYPE-WISE SNP FILTER, TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUAC.DN.R80_NOMAF_H.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.DN.R80_NOMAF_H.TwoPops.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.DN.R80_NOMAF_H.TwoPops.genind)
# %%%% REFERENCE ALIGNMENT %%%% ----
# %%%% R0 ----
# FIRST SNP, TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.R.R0_NOMAF_AllSNPs.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.R.R0_NOMAF_AllSNPs.TwoPops.genind) <- factor(read.table("QUAC_popmap", header=FALSE)[,2])
# Representation rates
reportAllelicRepresentation_Together(QUAC.R.R0_NOMAF_AllSNPs.TwoPops.genind)

# %%%% R80 ----
# FIRST SNP, TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.R.R80_NOMAF_1SNP.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.R.R80_NOMAF_1SNP.TwoPops.genind) <- factor(read.table("QUAC_popmap", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUAC.R.R80_NOMAF_1SNP.TwoPops.genind)
1# Exploration of total and wild allele frequency proportions
# Functions come from Simulated repository 
# (https://github.com/akoontz11/Morton_SSRvSNP_Simulations/blob/main/RScripts/functions_SSRvSNP_Sim.R)
getWildAlleleFreqProportions(QUAC.R.R80_NOMAF_1SNP.TwoPops.genind)
getTotalAlleleFreqProportions(QUAC.R.R80_NOMAF_1SNP.TwoPops.genind)

# %%%% QUBO %%%% ----
# %%%% R0 ----
# ---- ALL SNPS ----
# MULTIPLE WILD POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_AllSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF_AllSNPs.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF_AllSNPs.genind)

# TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_AllSNPs_2Pops/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_1SNP.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R0_NOMAF_1SNP.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF_1SNP.TwoPops.genind)

# ---- FIRST SNP ----
# MULTIPLE WILD POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_1SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF_1SNP.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF_1SNP.genind)

# TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_1SNP.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R0_NOMAF_1SNP.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF_1SNP.TwoPops.genind)
# Exploration of total and wild allele frequency proportions
# Functions come from Simulated repository 
# (https://github.com/akoontz11/Morton_SSRvSNP_Simulations/blob/main/RScripts/functions_SSRvSNP_Sim.R)
getWildAlleleFreqProportions(QUBO.R0_NOMAF_1SNP.TwoPops.genind)
getTotalAlleleFreqProportions(QUBO.R0_NOMAF_1SNP.TwoPops.genind)

# ---- HAPLOTYPE-WISE SNP FILTER ----
# TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUBO.R0_NOMAF_H.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R0_NOMAF_H.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R0_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R0_NOMAF_H.TwoPops.genind)

# %%%% R80 ----
# ---- ALL SNPS ----
# MULTIPLE WILD POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_AllSNPs/"
setwd(genpop.filePath)
QUBO.R80_NOMAF_AllSNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R80_NOMAF_AllSNPs.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R80_NOMAF_AllSNPs.genind)

# TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_AllSNPs_2Pops/"
setwd(genpop.filePath)
QUBO.R80_NOMAF_AllSNPs.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R80_NOMAF_AllSNPs.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R80_NOMAF_AllSNPs.TwoPops.genind)

# ---- FIRST SNP ----
# TWO POPULATIONS ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.R80_NOMAF_1SNP.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R80_NOMAF_1SNP.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R80_NOMAF_1SNP.TwoPops.genind)
# Exploration of total and wild allele frequency proportions
# Functions come from Simulated repository 
# (https://github.com/akoontz11/Morton_SSRvSNP_Simulations/blob/main/RScripts/functions_SSRvSNP_Sim.R)
getWildAlleleFreqProportions(QUBO.R80_NOMAF_1SNP.TwoPops.genind)
getTotalAlleleFreqProportions(QUBO.R80_NOMAF_1SNP.TwoPops.genind)

# ---- HAPLOTYPE-WISE SNP FILTER ----
# TWO POPULATIONS
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_HapSNPs_2Pops/"
setwd(genpop.filePath)
QUBO.R80_NOMAF_H.TwoPops.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.R80_NOMAF_H.TwoPops.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# R80_NOMAF Representation rates
reportAllelicRepresentation_Together(QUBO.R80_NOMAF_H.TwoPops.genind)
