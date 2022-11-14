# %%%%%%%%%%%%%%%%%%%%%%%%
# %%% FST CALCULATIONS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%

# This script generates Fst matrices from genind files generated for QUAC and QUBO, 
# and plots these matrices as heatmaps using the image command. 
# Pairwise Fst values are derived for wild populations only.

# Different methods for calculating Fst are offered, as outlined below
# 1. Fst_Nei_report: uses a genind file, from which a pairwise Fst is calculated using Nei's 1987 method
# 2. Fst_WC_report: uses a genind file, from which a pairwise Fst is calculated using Weir and Cockham's 1984 method
# 3. Fst_Stacks_report: reads in the populations.fst_summary.tsv file generated from the Stacks populations module (--fstats flag)

# For both species, this script reads in genind files from optimized de novo assemblies built using Stacks 
# (QUAC: m 7, M/n 5, gt-alpha 0.01; QUBO: m 7, M/n 5, gt-alpha 0.01)
# and genind files from reference alignments ("GSNAP4" alignment parameters; QUAC: Q. rubra genome;
# QUBO: Q. robur genome). Finally, for both de novo and refernce datasets, it reads in loci
# unfiltered according to missing data (R0) and loci shared among at least 80% (R80) of all samples (garden and wild)

# It does this for both the SNP (NextRAD) and corresponding MSAT datasets for both species. Both SNP and 
# MSAT datasets are pared down to just the samples shared between the two datasets, before calculating Fst

library(adegenet)
library(hierfstat)

# %%%% FUNCTIONS %%%% ----
# Set the working directory, in order to properly read in tissue and population names
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)
source("exSituRepresentation/functions_exSituRepresentation.R")

# Declare different functions for calculating Fst
# Function for plotting Fst values using a heatmap, given a matrix of values (and a title)
# The lower half of the fst_mat argument needs to have NA values, as does the diagonal 
Fst_plot <- function(fst_mat, title="Fst Values"){
  # Use image to plot a heatmap. First two arguments specify the boundaries of the heatmap; z provides actual values
  # z is transposed in order to plot numeral values later on
  image(x=1:ncol(fst_mat), y=1:nrow(fst_mat), z=t(fst_mat), axes=FALSE, xlab="", ylab="", 
        main=title)
  # Add boundary lines
  grid(nx=ncol(fst_mat), ny=nrow(fst_mat), col="black", lty=1)
  # Include sample names, to understand the context of genetic distances
  axis(1, 1:ncol(fst_mat), colnames(fst_mat), cex.axis=1.2, tick=FALSE)
  text(1, c(1:nrow(fst_mat)), labels=rownames(fst_mat), cex=1.2)
  # Go through matrix, and plot values on each cell
  for(x in 1:ncol(fst_mat)){
    for(y in 1:nrow(fst_mat)){
      text(x, y, fst_mat[y,x], cex=1.5)
    }
  }
}

# Function for calculating and plotting Fst values (Nei, 1987), from genind object (using hierfstat package)
Fst_Nei_report <- function(gen.obj, title){
  # Convert genind object to hierfstat format
  hierfstat.obj <- genind2hierfstat(gen.obj)
  # Calculate pairwise Fst (Nei, 1987)
  fst_mat <- pairwise.neifst(hierfstat.obj)
  # pairwise.neifst returns a full matrix, whereas we need only the upper half
  # The line below gets the matrix into the format needed for plotting (lower triangle values removed)
  fst_mat[lower.tri(fst_mat, diag = TRUE)] <- NA
  
  # Plot using the Fst_plot command
  Fst_plot(fst_mat, title=title)
  # Return matrix
  return(fst_mat)
}

# Function for calculating and plotting Fst values (Weir & Cockham, 1984), from genind object (using hierfstat package)
Fst_WC_report <- function(gen.obj, title){
  # Convert genind object to hierfstat format
  hierfstat.obj <- genind2hierfstat(gen.obj)
  # Calculate pairwise Fst (Weir & Cockham, 1984)
  fst_mat <- pairwise.WCfst(hierfstat.obj)
  # pairwise.neifst returns a full matrix, whereas we need only the upper half
  # The line below gets the matrix into the format needed for plotting (lower triangle values removed)
  fst_mat[lower.tri(fst_mat, diag = TRUE)] <- NA
  # Round the values in the Fst matrix to the 4th digit, to make it comparable to Nei values
  fst_mat <- round(fst_mat, digits = 4)
  
  # Plot using the Fst_plot command
  Fst_plot(fst_mat, title=title)
  # Return matrix
  return(fst_mat)
}

# Function for calculating and plotting Fst values, generated from Stacks populations module
# This function is mainly used for SNP R0 datasets, which are too large to be processed using other functions
Fst_Stacks_report <- function(filepath.fst_tab, title){
  # Read the populations.fst_summary.tsv file in the specified directory
  fst_mat <- as.matrix(read.table(paste0(filepath.fst_tab,"populations.fst_summary.tsv"), 
                                  header=TRUE, row.names=1, sep = "\t"))
  # Add a row at the bottom to make matrix symmetrical (nrow=ncol)
  fst_mat <- rbind(fst_mat, rep(NA, ncol(fst_mat)))
  # Update row names based on column names
  rownames(fst_mat) <- colnames(fst_mat)
  # Round the values in the Fst matrix to the 4th digit, to make it comparable to Nei values
  fst_mat <- round(fst_mat, digits = 4)
  
  # Plot using the Fst_plot command
  Fst_plot(fst_mat, title=title)
  # Return matrix
  return(fst_mat)
}

# %%%% QUAC %%%% ----
# ---- READ IN GENIND FILES ----
# MICROSATELLITE
QUAC.MSAT.genpop.filePath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/Garden_Wild/"
setwd(QUAC.MSAT.genpop.filePath)
QUAC.MSAT.genind <- read.genepop("QUAC_wK_garden_wild_clean.gen", ncode = 3)
# Correct popNames: pop1 is Garden, pop2 is Wild
pop(QUAC.MSAT.genind) <- gsub("pop1", "garden", pop(QUAC.MSAT.genind))
pop(QUAC.MSAT.genind) <- gsub("pop2", "wild", pop(QUAC.MSAT.genind))
# Subset to only wild individuals
QUAC.MSAT.genind <- QUAC.MSAT.genind[which(pop(QUAC.MSAT.genind)=="wild"),, drop=TRUE]

# SNP: DE NOVO, R0
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.DN.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R0.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Subset to only wild individuals
QUAC.SNP.DN.R0.genind <- QUAC.SNP.DN.R0.genind[which(pop(QUAC.SNP.DN.R0.genind)=="wild"),, drop=TRUE]

# SNP: DE NOVO, R80
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R80.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Subset to only wild individuals
QUAC.SNP.DN.R80.genind <- QUAC.SNP.DN.R80.genind[which(pop(QUAC.SNP.DN.R80.genind)=="wild"),, drop=TRUE]

# SNP: REFERENCE, R0
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.REF.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R0.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Subset to only wild individuals
QUAC.SNP.REF.R0.genind <- QUAC.SNP.REF.R0.genind[which(pop(QUAC.SNP.REF.R0.genind)=="wild"),, drop=TRUE]

# SNP: REFERENCE, R80
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R80.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Subset to only wild individuals
QUAC.SNP.REF.R80.genind <- QUAC.SNP.REF.R80.genind[which(pop(QUAC.SNP.REF.R80.genind)=="wild"),, drop=TRUE]

# ---- SUBSET GENIND FILES ----
# MSAT: read in Tissue database names from GCC_QUAC_ZAIN repository
QUAC.MSAT.tissueNames_filepath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Data_Frames/QUAC_allpop_clean_df.csv"
QUAC.MSAT.tissueNames <- unlist(read.csv2(QUAC.MSAT.tissueNames_filepath, header = TRUE, sep=",")[1])
# Remove garden samples, and rename MSAT matrix
QUAC.MSAT.tissueNames <- QUAC.MSAT.tissueNames[-grep(pattern = "QAc-G-", QUAC.MSAT.tissueNames)]
rownames(QUAC.MSAT.genind@tab) <- QUAC.MSAT.tissueNames

# SNP: read in Tissue database names
# This file was created (by Austin K.), and can be found on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
QUAC.SNP.tissueNames_filepath <- paste0(SSRvSNP.wd,"exSituRepresentation/Resampling/QUAC_SNP_TissueNames.csv")
QUAC.SNP.tissueNames <- unlist(read.csv2(QUAC.SNP.tissueNames_filepath, header = TRUE, sep = ",")[3])
# Remove garden samples, and rename SNP matrices
QUAC.SNP.tissueNames <- QUAC.SNP.tissueNames[-grep(pattern = "QAc-G-", QUAC.SNP.tissueNames)]
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

# ---- RENAME WILD POPULATIONS ----
# Read in wild population names from QUAC_SNP_TissuNames file (used above)
QUAC.wildPopNames <- unlist(read.csv2(QUAC.SNP.tissueNames_filepath, header = TRUE, sep = ",")[5])
# Remove populations corresponding to garden samples
QUAC.wildPopNames <- QUAC.wildPopNames[-(1:96)]

# Rename all subset genind objects (MSAT and SNP)
pop(QUAC.MSAT_subset.genind) <- pop(QUAC.SNP.DN.R0_subset.genind) <- pop(QUAC.SNP.DN.R80_subset.genind) <-
  pop(QUAC.SNP.REF.R0_subset.genind) <- pop(QUAC.SNP.REF.R80_subset.genind) <- QUAC.wildPopNames

# ---- CALCULATE FST VALUES ----
# Subset MSAT
Fst_Nei_report(QUAC.MSAT_subset.genind, title = "QUAC MSAT: Fst Values")

# Subset SNP
# De novo
# Fst_Nei_report(QUAC.SNP.DN.R0_subset.genind, title = "QUAC SNP De novo (R0): Fst Values")
Fst_Nei_report(QUAC.SNP.DN.R80_subset.genind, title = "QUAC SNP De novo (R80): Fst Values")
# Reference
# Fst_Nei_report(QUAC.SNP.REF.R0_subset.genind, title = "QUAC SNP De novo (R0): Fst Values")
Fst_Nei_report(QUAC.SNP.REF.R80_subset.genind, title = "QUAC SNP De novo (R80): Fst Values")

# %%%% QUBO %%%% ----
# ---- READ IN GENIND FILES ----
# MICROSATELLITE
genpop.filePath <- 
  "~/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
setwd(genpop.filePath)
QUBO.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_total.gen"), ncode=3)
# Correct popNames: last population (IMLS4_MP1_IMLS336_C05) is garden; rest (9) are wild
levels(QUBO.MSAT.genind@pop) <- c(rep("wild",9), "garden") 
# Subset to only wild individuals
QUBO.MSAT.genind <- QUBO.MSAT.genind[which(pop(QUBO.MSAT.genind)=="wild"),, drop=TRUE]

# SNP: DE NOVO, R0
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.DN.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R0.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# Subset to only wild individuals
QUBO.SNP.DN.R0.genind <- QUBO.SNP.DN.R0.genind[which(pop(QUBO.SNP.DN.R0.genind)=="wild"),, drop=TRUE]

# SNP: DE NOVO, R80
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R80.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# Subset to only wild individuals
QUBO.SNP.DN.R80.genind <- QUBO.SNP.DN.R80.genind[which(pop(QUBO.SNP.DN.R80.genind)=="wild"),, drop=TRUE]

# SNP: REFERENCE, R0
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.REF.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R0.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# Subset to only wild individuals
QUBO.SNP.REF.R0.genind <- QUBO.SNP.REF.R0.genind[which(pop(QUBO.SNP.REF.R0.genind)=="wild"),, drop=TRUE]

# SNP: REFERENCE, R80
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R80.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# Subset to only wild individuals
QUBO.SNP.REF.R80.genind <- QUBO.SNP.REF.R80.genind[which(pop(QUBO.SNP.REF.R80.genind)=="wild"),, drop=TRUE]

# ---- SUBSET GENIND FILES ----
# MSAT: Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.MSAT.genind@tab), function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.MSAT.genind@tab) <- QUBO.MSAT.sampleNames

# SNP: Remove QUBO_W_ headers from sample names
QUBO.SNP.sampleNames <- gsub("QUBO_W_",replacement = "", row.names(QUBO.SNP.REF.R0.genind@tab))
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
# MSAT
QUBO.MSAT_subset.genind <- QUBO.MSAT.genind[QUBO_sharedSamples,, drop=TRUE]
# SNP
# R0
QUBO.SNP.DN.R0_subset.genind <- QUBO.SNP.DN.R0.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.REF.R0_subset.genind <- QUBO.SNP.REF.R0.genind[QUBO_sharedSamples,, drop=TRUE]
# R80
QUBO.SNP.DN.R80_subset.genind <- QUBO.SNP.DN.R80.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.REF.R80_subset.genind <- QUBO.SNP.REF.R80.genind[QUBO_sharedSamples,, drop=TRUE]

# ---- RENAME WILD POPULATIONS ----
# Assign population names to subset QUBO wild samples
# This file was created (by Austin K.), and can be found on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
QUBO.wildPopNames_filepath <- paste0(SSRvSNP.wd,"popAnalyses/QUBO_sharedWildPops.csv")
QUBO.wildPopNames <- unlist(read.csv2(QUBO.wildPopNames_filepath, header = TRUE, sep = ",")[2])

# Rename all subset genind objects (MSAT and SNP)
pop(QUBO.MSAT_subset.genind) <- pop(QUBO.SNP.DN.R0_subset.genind) <- pop(QUBO.SNP.DN.R80_subset.genind) <-
  pop(QUBO.SNP.REF.R0_subset.genind) <- pop(QUBO.SNP.REF.R80_subset.genind) <- QUBO.wildPopNames

# ---- CALCULATE FST VALUES ----
# Subset MSAT
Fst_Nei_report(QUBO.MSAT_subset.genind, title = "QUBO MSAT: Fst Values")

# Subset SNP
# De novo
# Fst_Nei_report(QUBO.SNP.DN.R0_subset.genind, title = "QUBO SNP De novo (R0): Fst Values")
Fst_Nei_report(QUBO.SNP.DN.R80_subset.genind, title = "QUBO SNP De novo (R80): Fst Values")
# Reference
# Fst_Nei_report(QUBO.SNP.REF.R0_subset.genind, title = "QUBO SNP De novo (R0): Fst Values")
Fst_Nei_report(QUBO.SNP.REF.R80_subset.genind, title = "QUBO SNP De novo (R80): Fst Values")
