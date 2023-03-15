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
# For all of these, the mean pairwise Fst value is also calculated

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
  
  # Print mean pairwise Fst value
  print(paste0("%%% ", title, " %%%"))
  print(paste0("%%% Mean pairwise Fst: ", mean(fst_mat, na.rm=TRUE)))
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
  
  # Print mean pairwise Fst value
  print(paste0("%%% Mean pairwise Fst: ", mean(fst_mat, na.rm=TRUE)))
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
  
  # Print mean pairwise Fst value
  print(paste0("%%% Mean pairwise Fst: ", mean(fst_mat, na.rm=TRUE)))
  # Plot using the Fst_plot command
  Fst_plot(fst_mat, title=title)
  # Return matrix
  return(fst_mat)
}

# %%%% QUAC %%%% ----
# ---- READ IN GENIND FILES ----
# MICROSATELLITE ----
QUAC.MSAT.genpop.filePath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/"
QUAC.MSAT.genind <- read.genepop(paste0(QUAC.MSAT.genpop.filePath, "QUAC_woK_allpop_clean.gen"), ncode = 3)
# Specify filepath to GCC_QUAC_ZAIN dataframe, containing sample names and population names
QUAC.MSAT.dataframe_filepath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Data_Frames/QUAC_woK_allpop_clean_df.csv"
# Assign sample names: read in Tissue database names from GCC_QUAC_ZAIN repository
QUAC.MSAT.tissueNames <- unlist(read.csv2(QUAC.MSAT.dataframe_filepath, header = TRUE, sep=",")[1])
rownames(QUAC.MSAT.genind@tab) <- QUAC.MSAT.tissueNames
# Correct population names: read in a dataframe containing population values
pop(QUAC.MSAT.genind) <- unlist(read.csv2(QUAC.MSAT.dataframe_filepath, header = TRUE, sep=",")[2])
# Remove garden samples from tissueNames vector
QUAC.MSAT.tissueNames <- QUAC.MSAT.tissueNames[-grep(pattern = "QAc-G-", QUAC.MSAT.tissueNames)]
# Subset Complete MSAT genind object to just wild populations (last 5 populations)
QUAC.MSAT.genind <- QUAC.MSAT.genind[QUAC.MSAT.tissueNames,,drop=TRUE]

# SNP ----
# Read in Tissue database names, for subsetting genetic matrices later
# This file was created (by Austin K.), and can be found on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
QUAC.SNP.tissueNames_filepath <- paste0(SSRvSNP.wd,"exSituRepresentation/QUAC_SNP_TissueNames.csv")
QUAC.SNP.tissueNames <- unlist(read.csv2(QUAC.SNP.tissueNames_filepath, header = TRUE, sep = ",")[3])
# Remove garden samples from SNP tissue database names vector
QUAC.SNP.tissueNames <- QUAC.SNP.tissueNames[-grep(pattern = "QAc-G-", QUAC.SNP.tissueNames)]

# DE NOVO, R0
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF_1SNP_NoK/"
QUAC.SNP.DN.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R0.genind) <- factor(read.table(paste0(genpop.filePath, "QUAC_popmap_wild_NoK"), header=FALSE)[,2])
# Correct sample names
rownames(QUAC.SNP.DN.R0.genind@tab) <- QUAC.SNP.tissueNames

# SNP: DE NOVO, R80
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R80_NOMAF_1SNP_NoK/"
QUAC.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R80.genind) <- factor(read.table(paste0(genpop.filePath, "QUAC_popmap_wild_NoK"), header=FALSE)[,2])
# Correct sample names
rownames(QUAC.SNP.DN.R80.genind@tab) <- QUAC.SNP.tissueNames

# SNP: REFERENCE, R0
# Values cannot be calculated for original SNP REF R0 dataset, because Stacks population module
# is automatically killed (code 5806). Therefore, a set of 5,000 randomly selected "whitelisted" loci
# are chosen from a incomplete populations.sumstats.tsv file (from a different run, of the REF R0 dataset), 
# and these whitelisted loci are passed to the populations module for the successful run.
genpop.filePath <-
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_wild_R0_NOMAF_1SNP_NoK_WL/"
setwd(genpop.filePath)
QUAC.SNP.REF.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R0.genind) <- factor(read.table("QUAC_popmap_wild_NoK", header=FALSE)[,2])
# Correct sample names
rownames(QUAC.SNP.REF.R0.genind@tab) <- QUAC.SNP.tissueNames

# SNP: REFERENCE, R80
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_wild_R80_NOMAF_1SNP_NoK/"
QUAC.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R80.genind) <- factor(read.table(paste0(genpop.filePath, "QUAC_popmap_wild_NoK"), header=FALSE)[,2])
# Correct sample names
rownames(QUAC.SNP.REF.R80.genind@tab) <- QUAC.SNP.tissueNames

# ---- CALCULATE FST VALUES: COMPLETE ----
# MSAT
Fst_Nei_report(QUAC.MSAT.genind, title = "QUAC MSAT (Complete)")

# Subset SNP
# De novo
Fst_Nei_report(QUAC.SNP.DN.R0.genind, title = "QUAC SNP De novo (Complete, R0)")
Fst_Nei_report(QUAC.SNP.DN.R80.genind, title = "QUAC SNP De novo (Complete, R80)")
# Reference
Fst_Nei_report(QUAC.SNP.REF.R0.genind, title = "QUAC SNP Reference (Complete, R0)")
Fst_Nei_report(QUAC.SNP.REF.R80.genind, title = "QUAC SNP Reference (Complete, R80)")

# ---- SUBSET GENIND FILES ----
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

# ---- CALCULATE FST VALUES: SUBSET ----
# Subset MSAT
Fst_Nei_report(QUAC.MSAT_subset.genind, title = "QUAC MSAT (Subset)")

# Subset SNP
# De novo
Fst_Nei_report(QUAC.SNP.DN.R0_subset.genind, title = "QUAC SNP De novo (Subset, R0)")
Fst_Nei_report(QUAC.SNP.DN.R80_subset.genind, title = "QUAC SNP De novo (Subset, R80)")
# Reference
Fst_Nei_report(QUAC.SNP.REF.R0_subset.genind, title = "QUAC SNP Reference (Subset, R0)")
Fst_Nei_report(QUAC.SNP.REF.R80_subset.genind, title = "QUAC SNP Reference (Subset, R80)")

# %%%% QUBO %%%% ----
# ---- READ IN GENIND FILES ----
# MICROSATELLITE ----
# Read in the genind file from the SE oaks project repo. This genind contains only wild QUBO sampels
genpop.filePath <- 
  "~/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
QUBO.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_wild_w_ALL.gen"), ncode=3)
# Correct popNames: read in a dataframe which contains population names, in order of 
# samples in genind file. This document was created manually by Austin Koontz, and can be 
# found on the Hoban Lab Drive
QUBO.MSAT.dataframe_filepath <- 
  "~/Documents/SSRvSNP/Code/popAnalyses/QUBO_MSAT_Complete_WildPops.csv"
# Correct population names: read in a dataframe containing population values
pop(QUBO.MSAT.genind) <- unlist(read.csv2(QUBO.MSAT.dataframe_filepath, header = TRUE, sep=",")[2])
# Rename MSAT samples (for subsetting later on)
# Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.MSAT.genind@tab), function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.MSAT.genind@tab) <- QUBO.MSAT.sampleNames

# SNP ----
# DE NOVO, R0
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_wild_R0_NOMAF_1SNP_Ordered2/"
QUBO.SNP.DN.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R0.genind) <- factor(read.table(paste0(genpop.filePath, "QUBO_popmap_wild"), header=FALSE)[,2])

# SNP: DE NOVO, R80
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_wild_R80_NOMAF_1SNP_Ordered2/"
QUBO.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R80.genind) <- factor(read.table(paste0(genpop.filePath, "QUBO_popmap_wild"), header=FALSE)[,2])

# SNP: REFERENCE, R0
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0_NOMAF_1SNP_Ordered2/"
QUBO.SNP.REF.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R0.genind) <- factor(read.table(paste0(genpop.filePath, "QUBO_popmap_wild"), header=FALSE)[,2])

# SNP: REFERENCE, R80
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80_NOMAF_1SNP_Ordered2/"
QUBO.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R80.genind) <- factor(read.table(paste0(genpop.filePath, "QUBO_popmap_wild"), header=FALSE)[,2])

# Rename SNP samples (for subsetting later on)
# Remove QUBO_W_ headers from sample names
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
# Rename SNP sample matrices
rownames(QUBO.SNP.REF.R0.genind@tab) <- rownames(QUBO.SNP.DN.R0.genind@tab) <- QUBO.SNP.sampleNames
rownames(QUBO.SNP.REF.R80.genind@tab) <- rownames(QUBO.SNP.DN.R80.genind@tab) <- QUBO.SNP.sampleNames

# ---- CALCULATE FST VALUES: COMPLETE ----
# Complete MSAT
Fst_Nei_report(QUBO.MSAT.genind, title = "QUBO MSAT (Complete)")

# Complete SNP
# De novo
Fst_Nei_report(QUBO.SNP.DN.R0.genind, title = "QUBO SNP De novo (Complete, R0)")
Fst_Nei_report(QUBO.SNP.DN.R80.genind, title = "QUBO SNP De novo (Complete, R80)")
# Reference
Fst_Nei_report(QUBO.SNP.REF.R0.genind, title = "QUBO SNP Reference (Complete, R0)")
Fst_Nei_report(QUBO.SNP.REF.R80.genind, title = "QUBO SNP Reference (Complete, R80)")

# ---- SUBSET GENIND FILES ----
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

# ---- CALCULATE FST VALUES: SUBSET ----
# Subset MSAT
Fst_Nei_report(QUBO.MSAT_subset.genind, title = "QUBO MSAT (Subset)")

# Subset SNP
# De novo
Fst_Nei_report(QUBO.SNP.DN.R0_subset.genind, title = "QUBO SNP De novo (Subset, R0)")
Fst_Nei_report(QUBO.SNP.DN.R80_subset.genind, title = "QUBO SNP De novo (Subset, R80)")
# Reference
Fst_Nei_report(QUBO.SNP.REF.R0_subset.genind, title = "QUBO SNP Reference (Subset, R0)")
Fst_Nei_report(QUBO.SNP.REF.R80_subset.genind, title = "QUBO SNP Reference (Subset, R80)")
