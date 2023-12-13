# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% ALLELE FREQUENCY HISTOGRAMS: PLOTTING %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in microsatellite (MSAT) and SNP genind files for Quercus acerifolia (QUAC) 
# and Quercus boyntonii (QUBO). It then subsets those genind files (to generate identical samples 
# sets) and plots allele frequency histograms. These plots are included in the supplement for the 
# Evolutionary Applications manuscript submitted for the Empirical MSAT v. SNP project. For SNP datasets, 
# only R80 data is shown.

library(adegenet)
library(RColorBrewer)
library(scales)
library(parallel)

# %%%% FUNCTIONS AND VARIABLES %%%% ----
# Read in relevant functions required for resampling analyses
SSRvSNP.wd <- "/home/akoontz/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)
source("exSituRep/functions_exSituRepresentation.R")
# Declare custom function for plotting a single histogram of allele frequencies from a genind object
makeAlleleFreqHist <- function(gen.obj, title=gen.obj@other, yMax=nLoc(gen.obj)){
  # Make a vector of allele frequency values, from the genind object
  wildAlleleFreqs <- getWildFreqs(gen.obj)
  # Generate histogram
  hist(wildAlleleFreqs, main=title, ylim=c(0, yMax), xlab="", ylab="")
}
# Specify path to the directory (on the lab server), where resampling plots (PDFs) will be saved
imageOutDir <- "/home/akoontz/Documents/SSRvSNP/Documentation/Images/EvolApp_2023_Images/FINAL/supplement/"

# %%%% READ IN AND PROCESS GENIND FILES %%%% ----
# ---- QUAC ----
# MSATS ----
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

# SNPS, DE NOVO ----
# Read in genind file: Optimized de novo assembly; R80, min-maf=0, first SNP/locus, 2 populations (garden and wild), no Kessler individuals
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops_NoK/"
QUAC.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath, "QUAC_popmap_GardenWild_NoK"), header=FALSE)[,2])

# SNPS, REFERENCE ----
# Read in genind file: QUAC GSNAP4 alignment; R80, min-maf=0, first SNP/locus, 2 populations (garden and wild), no Kessler individuals
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_R80_NOMAF_1SNP_2Pops_NoK/"
QUAC.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath,"QUAC_popmap_GardenWild_NoK"), header=FALSE)[,2])

# SUBSET DATASETS ----
# SNP: read in Tissue database names, and rename SNP genind matrix
# Austin Koontz created this file, and it lives on the Hoban Lab Drive 
# ("MSATcomparisons_TissueNames:QUAC_NextRAD_TissueDatabaseNames")
QUAC.SNP.tissueNames_filepath <- paste0(SSRvSNP.wd,"exSituRep/QUAC_SNP_TissueNames.csv")
QUAC.SNP.tissueNames <- unlist(read.csv2(QUAC.SNP.tissueNames_filepath, header = TRUE, sep = ",")[3])
rownames(QUAC.SNP.DN.R80.genind@tab) <- rownames(QUAC.SNP.REF.R80.genind@tab) <- QUAC.SNP.tissueNames

# Subset SNP sample names by those that are also seen within the MSAT samples 
QUAC_sharedSamples <- sort(QUAC.SNP.tissueNames[which(QUAC.SNP.tissueNames %in% QUAC.MSAT.tissueNames)])
# Subset MSAT and SNP genind matrices to strictly shared samples, dropping now absent alleles
QUAC.MSAT_subset.genind <- QUAC.MSAT.genind[QUAC_sharedSamples,, drop=TRUE]
QUAC.SNP.DN.R80_subset.genind <- QUAC.SNP.DN.R80.genind[QUAC_sharedSamples,, drop=TRUE]
QUAC.SNP.REF.R80_subset.genind <- QUAC.SNP.REF.R80.genind[QUAC_sharedSamples,, drop=TRUE]

# ---- QUBO ----
# MSATS ----
# Read in genind file (Southeast Oaks repo; genetic_data/Qb_total.gen file)
genpop.filePath <- 
  "/home/akoontz/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
QUBO.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_total.gen"), ncode=3)
# Correct popNames: last population (IMLS4_MP1_IMLS336_C05) is garden; rest (9) are wild
levels(QUBO.MSAT.genind@pop) <- c(rep("wild",9), "garden") 

# SNPS: DE NOVO ASSEMBLY ----
# Read in genind file: Optimized de novo assembly; R80, min-maf=0, first SNP/locus, 2 populations (garden and wild)
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath,"QUBO_popmap_GardenWild"), header=FALSE)[,2])

# SNPS: REFERENCE ----
# Read in genind file: QUBO GSNAP4 alignment; R80, min-maf=0, first SNP/locus, 2 populations (garden and wild)
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_1SNP_2Pops/"
QUBO.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath,"QUBO_popmap_GardenWild"), header=FALSE)[,2])

# SUBSET DATASETS ----
# MSAT: Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.MSAT.genind@tab), function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.MSAT.genind@tab) <- QUBO.MSAT.sampleNames

# SNP: Remove QUBO_W_ headers from sample names
QUBO.SNP.sampleNames <- gsub("QUBO_G_",replacement = "", row.names(QUBO.SNP.REF.R80.genind@tab))
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
rownames(QUBO.SNP.REF.R80.genind@tab) <- rownames(QUBO.SNP.DN.R80.genind@tab) <- QUBO.SNP.sampleNames

# Subset SNP sample names by those that are also seen within the MSAT samples
QUBO_sharedSamples <- sort(QUBO.SNP.sampleNames[which(QUBO.SNP.sampleNames %in% QUBO.MSAT.sampleNames)])
# Subset MSAT and SNP wild matrix objects to strictly shared samples
QUBO.MSAT_subset.genind <- QUBO.MSAT.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.DN.R80_subset.genind <- QUBO.SNP.DN.R80.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.REF.R80_subset.genind <- QUBO.SNP.REF.R80.genind[QUBO_sharedSamples,, drop=TRUE]

# %%%% PLOTTING %%%% ----
# QUAC %%%% ----
# Call pdf command, to save resampling plots to disk
pdf(file = paste0(imageOutDir, "Fig_S8.pdf"), width = 9, height = 7.5)
# Set plotting window to stack 3 graphs vertically
par(mfcol=c(3,1), oma=rep(0.2,4), mar=c(2,4.5,4,1.5)+0.1)
# MSAT
makeAlleleFreqHist(QUAC.MSAT_subset.genind, title = "", yMax = 150)
mtext("Microsatellites", side=1, line=-7, cex=1.2)
# Plot title
title("Quercus acerifolia: Wild allele frequency distributions", line=1.5, cex.main=2)
# Update margins to reduce space above middle graph
par(mar=c(2,4.5,2,1.5)+0.1)
# SNP: DE NOVO ASSEMBLY (R80)
makeAlleleFreqHist(QUAC.SNP.DN.R80_subset.genind, title = "", yMax=6000)
mtext("SNPs: De novo (R80)", side=1, line=-7, cex=1.2)
# y-axis label
mtext(text="Number of alleles", side=2, line=3, cex=1, srt=90)
# Update margins to allow for more space below bottom graph
par(mar=c(4,4.5,2,1.5)+0.1)
# SNP: REFERENCE (R80)
makeAlleleFreqHist(QUAC.SNP.REF.R80_subset.genind, title = "", yMax=15000)
mtext("SNPs: Reference (R80)", side=1, line=-7, cex=1.2)
# x-axis label
mtext(text="Allele frequency category", side=1, line=2.5, cex=1)
# Turn off plotting device
dev.off()

# QUBO %%%% ----
# Call pdf command, to save resampling plots to disk
pdf(file = paste0(imageOutDir, "Fig_S9.pdf"), width = 9, height = 7.5)
# Set plotting window to stack 3 graphs vertically
par(mfcol=c(3,1), oma=rep(0.2,4), mar=c(2,4.5,4,1.5)+0.1)
# MSAT
makeAlleleFreqHist(QUBO.MSAT_subset.genind, title = "", yMax = 120)
mtext("Microsatellites", side=1, line=-7, cex=1.2)
# Plot title
title("Quercus boyntonii: Wild allele frequency distributions", line=1.5, cex.main=2)
# Update margins to reduce space above middle graph
par(mar=c(2,4.5,2,1.5)+0.1)
# SNP: DE NOVO ASSEMBLY (R80)
makeAlleleFreqHist(QUBO.SNP.DN.R80_subset.genind, title = "", yMax=5000)
mtext("SNPs: De novo (R80)", side=1, line=-7, cex=1.2)
# y-axis label
mtext(text="Number of alleles", side=2, line=3, cex=1, srt=90)
# Update margins to allow for more space below bottom graph
par(mar=c(4,4.5,2,1.5)+0.1)
# SNP: REFERENCE (R80)
makeAlleleFreqHist(QUBO.SNP.REF.R80_subset.genind, title = "", yMax=7000)
mtext("SNPs: Reference (R80)", side=1, line=-7, cex=1.2)
# x-axis label
mtext(text="Allele frequency category", side=1, line=2.5, cex=1)
# Turn off plotting device
dev.off()
