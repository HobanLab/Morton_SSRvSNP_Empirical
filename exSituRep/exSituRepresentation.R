# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% EX SITU REPRESENTATION RATES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses the functions declared in the functions_exSituRepresentation.R file to generate 
# ex situ representation values for Quercus acerifolia (QUAC) and Quercus boyntonii (QUBO). QUAC 
# datasets without Kessler individuals are analyzed.

# For both species, it reads in datasets from optimized de novo assemblies built using Stacks 
# (QUAC: m 7, M/n 5, gt-alpha 0.01; QUBO: m 7, M/n 5, gt-alpha 0.01)
# and datasets from reference alignments ("GSNAP4" alignment parameters; QUAC: Q. rubra genome;
# QUBO: Q. robur genome). Finally, for both de novo and reference datasets, it reads in loci
# unfiltered according to missing data (R0) and loci shared among at least 80% (R80) of all samples 
# (garden and wild)

# It does this for both the SNP (NextRAD) and MSAT datasets, using both complete datasets
# and partial datasets (samples subset to only those included within both studies)

library(adegenet)

# %%%% FUNCTIONS %%%% ----
# Read in relevant functions required for ex situ representation analyses
SSRvSNP.wd <- "/home/akoontz/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)
source("exSituRep/functions_exSituRepresentation.R")

# %%%% QUAC %%%% ----
# ---- MSATS: COMPLETE ----
QUAC.MSAT.filePath <- 
  "/home/akoontz/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/"
QUAC.MSAT.genind <- 
  read.genepop(paste0(QUAC.MSAT.filePath, "Adegenet_Files/QUAC_woK_allpop_clean.gen"), ncode = 3)
# Assign sample names: read in Tissue database names from GCC_QUAC_ZAIN repository
QUAC.MSAT.tissueNames <- 
  unlist(read.csv2(paste0(QUAC.MSAT.filePath, "Data_Frames/QUAC_woK_allpop_clean_df.csv"), header = TRUE, sep=",")[1])
rownames(QUAC.MSAT.genind@tab) <- QUAC.MSAT.tissueNames
# Correct popNames: samples with popName pattern QAc-G- are garden 
levels(QUAC.MSAT.genind@pop)[grep(pattern = "QAc-G-", levels(QUAC.MSAT.genind@pop))] <- 
  rep("garden", length(grep(pattern = "QAc-G-", levels(QUAC.MSAT.genind@pop))))
# Correct popNames: samples with popName pattern QAc-W- are wild
levels(QUAC.MSAT.genind@pop)[grep(pattern = "QAc-W-", levels(QUAC.MSAT.genind@pop))] <- 
  rep("wild", length(grep(pattern = "QAc-W-", levels(QUAC.MSAT.genind@pop))))
# Report representation of wild alleles in gardens
exSitu_Rep(QUAC.MSAT.genind)
# Exploration of total and wild allele frequency proportions.
getWildAlleleFreqProportions(QUAC.MSAT.genind)
getTotalAlleleFreqProportions(QUAC.MSAT.genind)

# ---- SNPS, DE NOVO: COMPLETE ----
# R0 ----
# Read in genind file: Optimized de novo assembly; R0, min-maf=0, first SNP/locus, 2 populations (garden and wild), no Kessler individuals
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops_NoK/"
QUAC.SNP.DN.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R0.genind) <- 
  factor(read.table(paste0(genpop.filePath, "QUAC_popmap_GardenWild_NoK"), header=FALSE)[,2])
# Report representation of wild alleles in gardens
exSitu_Rep(QUAC.SNP.DN.R0.genind)
# Exploration of total and wild allele frequency proportions
getWildAlleleFreqProportions(QUAC.SNP.DN.R0.genind)
getTotalAlleleFreqProportions(QUAC.SNP.DN.R0.genind)

# R80 ----
# Read in genind file: Optimized de novo assembly; R80, min-maf=0, first SNP/locus, 2 populations (garden and wild), no Kessler individuals
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops_NoK/"
QUAC.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath, "QUAC_popmap_GardenWild_NoK"), header=FALSE)[,2])
# Report representation of wild alleles in gardens
exSitu_Rep(QUAC.SNP.DN.R80.genind)
# Exploration of total and wild allele frequency proportions
getWildAlleleFreqProportions(QUAC.SNP.DN.R80.genind)
getTotalAlleleFreqProportions(QUAC.SNP.DN.R80.genind)

# ---- SNPS, REFERENCE: COMPLETE ----
# R0 ----
# Read in genind file: QUAC GSNAP4 alignment; R0, min-maf=0, first SNP/locus, 2 populations (garden and wild), no Kessler individuals
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R0_NOMAF_1SNP_2Pops_NoK/"
QUAC.SNP.REF.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R0.genind) <- 
  factor(read.table(paste0(genpop.filePath, "QUAC_popmap_GardenWild_NoK"), header=FALSE)[,2])
# Representation rates
exSitu_Rep(QUAC.SNP.REF.R0.genind)
# Exploration of total and wild allele frequency proportions
getWildAlleleFreqProportions(QUAC.SNP.REF.R0.genind)
getTotalAlleleFreqProportions(QUAC.SNP.REF.R0.genind)

# R80 ----
# Read in genind file: QUAC GSNAP4 alignment; R80, min-maf=0, first SNP/locus, 2 populations (garden and wild), no Kessler individuals
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R80_NOMAF_1SNP_2Pops_NoK/"
QUAC.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath,"QUAC_popmap_GardenWild_NoK"), header=FALSE)[,2])
# Representation rates
exSitu_Rep(QUAC.SNP.REF.R80.genind)
# Exploration of total and wild allele frequency proportions
getWildAlleleFreqProportions(QUAC.SNP.REF.R80.genind)
getTotalAlleleFreqProportions(QUAC.SNP.REF.R80.genind)

# ---- MSATS AND SNPS: SUBSET ----
# SNP: read in Tissue database names, and rename SNP genind matrix
# Austin Koontz created this file, and it lives on the Hoban Lab Drive 
# ("MSATcomparisons_TissueNames:QUAC_NextRAD_TissueDatabaseNames")
QUAC.SNP.tissueNames_filepath <- paste0(SSRvSNP.wd,"exSituRep/QUAC_SNP_TissueNames.csv")
QUAC.SNP.tissueNames <- unlist(read.csv2(QUAC.SNP.tissueNames_filepath, header = TRUE, sep = ",")[3])
rownames(QUAC.SNP.DN.R0.genind@tab) <- rownames(QUAC.SNP.REF.R0.genind@tab) <- QUAC.SNP.tissueNames
rownames(QUAC.SNP.DN.R80.genind@tab) <- rownames(QUAC.SNP.REF.R80.genind@tab) <- QUAC.SNP.tissueNames

# Subset SNP sample names by those that are also seen within the MSAT samples 
QUAC_sharedSamples <- sort(QUAC.SNP.tissueNames[which(QUAC.SNP.tissueNames %in% QUAC.MSAT.tissueNames)])
# Subset MSAT and SNP genind matrices to strictly shared samples, dropping now absent alleles
QUAC.MSAT_subset.genind <- QUAC.MSAT.genind[QUAC_sharedSamples,, drop=TRUE]
# R0
QUAC.SNP.DN.R0_subset.genind <- QUAC.SNP.DN.R0.genind[QUAC_sharedSamples,, drop=TRUE]
QUAC.SNP.REF.R0_subset.genind <- QUAC.SNP.REF.R0.genind[QUAC_sharedSamples,, drop=TRUE]
# R80
QUAC.SNP.DN.R80_subset.genind <- QUAC.SNP.DN.R80.genind[QUAC_sharedSamples,, drop=TRUE]
QUAC.SNP.REF.R80_subset.genind <- QUAC.SNP.REF.R80.genind[QUAC_sharedSamples,, drop=TRUE]

# Subset MSAT
# Wild allelic representation in gardens
exSitu_Rep(QUAC.MSAT_subset.genind)
# Allele frequency proportions
getWildAlleleFreqProportions(QUAC.MSAT_subset.genind)
getTotalAlleleFreqProportions(QUAC.MSAT_subset.genind)

# Subset SNP
# Wild allelic representation in gardens
# De novo
exSitu_Rep(QUAC.SNP.DN.R0_subset.genind)
exSitu_Rep(QUAC.SNP.DN.R80_subset.genind)
# Reference
exSitu_Rep(QUAC.SNP.REF.R0_subset.genind)
exSitu_Rep(QUAC.SNP.REF.R80_subset.genind)
# Allele frequency proportions
# De novo
getWildAlleleFreqProportions(QUAC.SNP.DN.R80_subset.genind)
getWildAlleleFreqProportions(QUAC.SNP.DN.R80_subset.genind)
getTotalAlleleFreqProportions(QUAC.SNP.DN.R80_subset.genind)
getTotalAlleleFreqProportions(QUAC.SNP.DN.R80_subset.genind)
# Reference
getWildAlleleFreqProportions(QUAC.SNP.REF.R80_subset.genind)
getWildAlleleFreqProportions(QUAC.SNP.REF.R80_subset.genind)
getTotalAlleleFreqProportions(QUAC.SNP.REF.R80_subset.genind)
getTotalAlleleFreqProportions(QUAC.SNP.REF.R80_subset.genind)

# %%%% QUBO %%%% ----
# ---- MSATS: COMPLETE ----
# Read in genind file (Southeast Oaks repo; genetic_data/Qb_total.gen file)
genpop.filePath <- 
  "/home/akoontz/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
QUBO.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_total.gen"), ncode=3)
# Correct popNames: last population (IMLS4_MP1_IMLS336_C05) is garden; rest (9) are wild
levels(QUBO.MSAT.genind@pop) <- c(rep("wild",9), "garden") 
# Report representation of wild alleles in gardens
exSitu_Rep(QUBO.MSAT.genind)
# Exploration of total and wild allele frequency proportions
getWildAlleleFreqProportions(QUBO.MSAT.genind)
getTotalAlleleFreqProportions(QUBO.MSAT.genind)

# ---- SNPS: DE NOVO ASSEMBLY: COMPLETE ----
# R0 ----
# Read in genind file: Optimized de novo assembly; R0, min-maf=0, first SNP/locus, 2 populations (garden and wild)
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R0_NOMAF_1SNP_2Pops/"
QUBO.SNP.DN.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R0.genind) <- 
  factor(read.table(paste0(genpop.filePath,"QUBO_popmap_GardenWild"), header=FALSE)[,2])
# Report representation of wild alleles in gardens
exSitu_Rep(QUBO.SNP.DN.R0.genind)
# Exploration of total and wild allele frequency proportions
getWildAlleleFreqProportions(QUBO.SNP.DN.R0.genind)
getTotalAlleleFreqProportions(QUBO.SNP.DN.R0.genind)

# R80 ----
# Read in genind file: Optimized de novo assembly; R80, min-maf=0, first SNP/locus, 2 populations (garden and wild)
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath,"QUBO_popmap_GardenWild"), header=FALSE)[,2])
# Report representation of wild alleles in gardens
exSitu_Rep(QUBO.SNP.DN.R80.genind)
# Exploration of total and wild allele frequency proportions
getWildAlleleFreqProportions(QUBO.SNP.DN.R80.genind)
getTotalAlleleFreqProportions(QUBO.SNP.DN.R80.genind)

# ---- SNPS: REFERENCE: COMPLETE ----
# R0 ----
# Read in genind file: QUBO GSNAP4 alignment; R0, min-maf=0, first SNP/locus, 2 populations (garden and wild)
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
QUBO.SNP.REF.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R0.genind) <- 
  factor(read.table(paste0(genpop.filePath,"QUBO_popmap_GardenWild"), header=FALSE)[,2])
# Report representation of wild alleles in gardens
exSitu_Rep(QUBO.SNP.REF.R0.genind)
# Exploration of total and wild allele frequency proportions
getWildAlleleFreqProportions(QUBO.SNP.REF.R0.genind)
getTotalAlleleFreqProportions(QUBO.SNP.REF.R0.genind)

# R80 ----
# Read in genind file: QUBO GSNAP4 alignment; R80, min-maf=0, first SNP/locus, 2 populations (garden and wild)
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_1SNP_2Pops/"
QUBO.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath,"QUBO_popmap_GardenWild"), header=FALSE)[,2])
# Report representation of wild alleles in gardens
exSitu_Rep(QUBO.SNP.REF.R80.genind)
# Exploration of total and wild allele frequency proportions
getWildAlleleFreqProportions(QUBO.SNP.REF.R80.genind)
getTotalAlleleFreqProportions(QUBO.SNP.REF.R80.genind)

# ---- MSATS AND SNPS: SUBSET ----
# MSAT: Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.MSAT.genind@tab), function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.MSAT.genind@tab) <- QUBO.MSAT.sampleNames

# SNP: Remove QUBO_W_ headers from sample names
QUBO.SNP.sampleNames <- gsub("QUBO_G_",replacement = "", row.names(QUBO.SNP.REF.R0.genind@tab))
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
# Rename sample matrices
rownames(QUBO.SNP.REF.R0.genind@tab) <- rownames(QUBO.SNP.DN.R0.genind@tab) <- QUBO.SNP.sampleNames
rownames(QUBO.SNP.REF.R80.genind@tab) <- rownames(QUBO.SNP.DN.R80.genind@tab) <- QUBO.SNP.sampleNames

# Subset SNP sample names by those that are also seen within the MSAT samples
QUBO_sharedSamples <- sort(QUBO.SNP.sampleNames[which(QUBO.SNP.sampleNames %in% QUBO.MSAT.sampleNames)])
# Subset MSAT and SNP wild matrix objects to strictly shared samples
QUBO.MSAT_subset.genind <- QUBO.MSAT.genind[QUBO_sharedSamples,, drop=TRUE]
# R0
QUBO.SNP.DN.R0_subset.genind <- QUBO.SNP.DN.R0.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.REF.R0_subset.genind <- QUBO.SNP.REF.R0.genind[QUBO_sharedSamples,, drop=TRUE]
# R80
QUBO.SNP.DN.R80_subset.genind <- QUBO.SNP.DN.R80.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.REF.R80_subset.genind <- QUBO.SNP.REF.R80.genind[QUBO_sharedSamples,, drop=TRUE]

# Subset MSAT 
# Wild allelic representation in gardens
exSitu_Rep(QUBO.MSAT_subset.genind)
# Allele frequency proportions
getWildAlleleFreqProportions(QUBO.MSAT_subset.genind)
getTotalAlleleFreqProportions(QUBO.MSAT_subset.genind)

# Subset SNP
# Wild allelic representation in gardens
# De novo
exSitu_Rep(QUBO.SNP.DN.R0_subset.genind)
exSitu_Rep(QUBO.SNP.DN.R80_subset.genind)
# Reference
exSitu_Rep(QUBO.SNP.REF.R0_subset.genind)
exSitu_Rep(QUBO.SNP.REF.R80_subset.genind)
# Allele frequency proportions
# De novo
getWildAlleleFreqProportions(QUBO.SNP.DN.R80_subset.genind)
getWildAlleleFreqProportions(QUBO.SNP.DN.R80_subset.genind)
getTotalAlleleFreqProportions(QUBO.SNP.DN.R80_subset.genind)
getTotalAlleleFreqProportions(QUBO.SNP.DN.R80_subset.genind)
# Reference
getWildAlleleFreqProportions(QUBO.SNP.REF.R80_subset.genind)
getWildAlleleFreqProportions(QUBO.SNP.REF.R80_subset.genind)
getTotalAlleleFreqProportions(QUBO.SNP.REF.R80_subset.genind)
getTotalAlleleFreqProportions(QUBO.SNP.REF.R80_subset.genind)
