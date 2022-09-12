# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% EX SITU REPRESENTATION RATES: QUHA APPROACH %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses the approach from the Qharvardii_ex_situ analysis to generate ex situ representation values 
# for Quercus acerifolia (QUAC; optimized Stacks de novo assembly, m 7, M/n 4, gt-alpha 0.01) 
# and Quercus boyntonii (QUBO; GSNAP4 alignment with Quercus robur reference) samples

# In addition to analyzing the SNP (NextRAD) datasets, this script generates ex situ representation values for 
# corresponding QUAC/QUBO microsatellite (MSAT) datasets as well. It also subsets SNP and MSAT datasets to only include
# samples shared between the two, and measures ex situ representation in these subset datasets as well.

# Specifically, it uses the get.allele.cat function declared in the Fa_sample_funcs.R file, along with the structure
# used in the Qhav_ex_situ_code.R script, to calculate ex situ representation values. Both the get.allele.cat function
# and the Qhav_ex_situ_code.R script have been adapted to work with the QUAC and QUBO SNP and MSATdatasets 
# (single garden population, no n_to_drop values, etc.)

# The purpose of this script is to demonstrate the identical results between the two approaches used to measure
# ex situ representation, and thereby validate the results

library(adegenet)
# Set working directory for ex situ analyses
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)

# %%%% FUNCTIONS %%%% ----
# Adapted version of the original get.allele.cat function, without analysis by region or n_to_drop arguments
# Also, note that these categories do not take into account the minimum allele frequency. Therefore, total
# and rare allele calculations are modified from the original function
get.allele.cat <- function(genpop, n_ind){
  #--Set up categories for GLOBAL ALLELES
  allele_freqs<-colSums(genpop@tab, na.rm = TRUE)/(sum(n_ind)*2)
  
  total <- as.vector(which(allele_freqs > 0))
  v_com <- as.vector(which(allele_freqs > 0.10))
  com <- as.vector(which(allele_freqs > 0.05))
  low_freq <- as.vector(which(allele_freqs < 0.10 & allele_freqs > 0.01))
  rare <- as.vector(which((allele_freqs < 0.01) & (allele_freqs > 0)))
  
  # Arrange into a list
  list(total, v_com, com, low_freq, rare)
}	

# %%%% QUAC %%%% ----
# ---- MSATS: COMPLETE ----
# Read in genind file (GCC_QUAC_ZAIN repo; QUAC_wK_garden_wild_clean.gen; includes Kessler population)
QUAC.MSAT.genpop.filePath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/Garden_Wild/"
setwd(QUAC.MSAT.genpop.filePath)
QUAC.MSAT.genind <- read.genepop("QUAC_wK_garden_wild_clean.gen", ncode = 3)
# Correct popNames: pop1 is Garden, pop2 is Wild
pop(QUAC.MSAT.genind) <- gsub("pop1", "garden", pop(QUAC.MSAT.genind))
pop(QUAC.MSAT.genind) <- gsub("pop2", "wild", pop(QUAC.MSAT.genind))

# Declare vector of results
QUAC.MSAT.repRates <- vector(length = 5)
# Create vector to capture the number of each allele category (of which there are 5)
alleles_existing_by_sp <- vector(length = 5)
# First population is all garden samples ("garden"), population 2 is "wild"
garden_p <- 1 ; wild_p <- 2
n_ind_W <- table(QUAC.MSAT.genind@pop)[wild_p]; n_ind_G <- table(QUAC.MSAT.genind@pop)[garden_p]; 
QUAC.MSAT.genpop <- genind2genpop(QUAC.MSAT.genind)
alleles_cap <- colSums(QUAC.MSAT.genpop[garden_p]@tab,na.rm=T)
# The below calls returns each alleles for each category specified
allele_cat_tot <- get.allele.cat(QUAC.MSAT.genpop[wild_p], n_ind_W)
# For each allele category, find out how many alleles there are
for (i in 1:5) alleles_existing_by_sp[i] <- (sum((allele_cat_tot[[i]])>0,na.rm=T))
# For each allele category, calculate the captured allele percentages 
# by dividing by total number of alleles in that category
for (l in 1:length(allele_cat_tot)) QUAC.MSAT.repRates[l] <- round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
# Assign names to values, and print
names(QUAC.MSAT.repRates)<-c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%","Rare (<1%)")
print(QUAC.MSAT.repRates*100)

# ---- SNPS: COMPLETE ----
# Read in genind file: Optimized de novo assembly; R0, NOMAF, first SNP/locus, 2 populations
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])

# Declare vector of results
QUAC.SNP.repRates <- vector(length = 5)
# Create vector to capture the number of each allele category (of which there are 5)
alleles_existing_by_sp <- vector(length = 5)
# First population is all garden samples ("garden"), population 2 is "wild"
garden_p <- 1 ; wild_p <- 2
n_ind_W <- table(QUAC.SNP.genind@pop)[wild_p]; n_ind_G <- table(QUAC.SNP.genind@pop)[garden_p]; 
QUAC.SNP.genpop <- genind2genpop(QUAC.SNP.genind)
alleles_cap <- colSums(QUAC.SNP.genpop[garden_p]@tab,na.rm=T)
# The below calls returns each alleles for each category specified
allele_cat_tot <- get.allele.cat(QUAC.SNP.genpop[wild_p], n_ind_W)
# For each allele category, find out how many alleles there are
for (i in 1:5) alleles_existing_by_sp[i] <- (sum((allele_cat_tot[[i]])>0,na.rm=T))
# For each allele category, calculate the captured allele percentages 
# by dividing by total number of alleles in that category
for (l in 1:length(allele_cat_tot)) QUAC.SNP.repRates[l] <- round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
# Assign names to values, and print
names(QUAC.SNP.repRates)<-c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%","Rare (<1%)")
print(QUAC.SNP.repRates*100)

# ---- MSATS AND SNPS: SUBSET ----
# MSAT: read in Tissue database names from GCC_QUAC_ZAIN repository, and rename MSAT genind matrix
QUAC.MSAT.tissueNames_filepath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Data_Frames/QUAC_allpop_clean_df.csv"
QUAC.MSAT.tissueNames <- unlist(read.csv2(QUAC.MSAT.tissueNames_filepath, header = TRUE, sep=",")[1])
rownames(QUAC.MSAT.genind@tab) <- QUAC.MSAT.tissueNames

# SNP: read in Tissue database names, and rename SNP genind matrix
# This file was created (by Austin K.), and can be found on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
QUAC.SNP.tissueNames_filepath <- paste0(SSRvSNP.wd,"exSituRepresentation/Resampling/QUAC_SNP_TissueNames.csv")
QUAC.SNP.tissueNames <- unlist(read.csv2(QUAC.SNP.tissueNames_filepath, header = TRUE, sep = ",")[3])
rownames(QUAC.SNP.genind@tab) <- QUAC.SNP.tissueNames

# Subset SNP sample names by those that are also seen within the MSAT samples 
QUAC_sharedSamples <- sort(QUAC.SNP.tissueNames[which(QUAC.SNP.tissueNames %in% QUAC.MSAT.tissueNames)])
# Subset MSAT and SNP genind matrices to strictly shared samples, dropping now absent alleles
QUAC.MSAT_subset.genind <- QUAC.MSAT.genind[QUAC_sharedSamples,, drop=TRUE]
QUAC.SNP_subset.genind <- QUAC.SNP.genind[QUAC_sharedSamples,, drop=TRUE]

# Allelic capture of subset MSATs ----
# Declare vector of results
QUAC.MSAT_subset.repRates <- vector(length = 5)
# Create vector to capture the number of each allele category (of which there are 5)
alleles_existing_by_sp <- vector(length = 5)
# First population is all garden samples ("garden"), population 2 is "wild"
garden_p <- 1 ; wild_p <- 2
n_ind_W <- table(QUAC.MSAT_subset.genind@pop)[wild_p]; n_ind_G <- table(QUAC.MSAT_subset.genind@pop)[garden_p]; 
QUAC.MSAT_subset.genpop <- genind2genpop(QUAC.MSAT_subset.genind)
alleles_cap <- colSums(QUAC.MSAT_subset.genpop[garden_p]@tab,na.rm=T)
# The below calls returns each alleles for each category specified
allele_cat_tot <- get.allele.cat(QUAC.MSAT_subset.genpop[wild_p], n_ind_W)
# For each allele category, find out how many alleles there are
for (i in 1:5) alleles_existing_by_sp[i] <- (sum((allele_cat_tot[[i]])>0,na.rm=T))
# For each allele category, calculate the captured allele percentages 
# by dividing by total number of alleles in that category
for (l in 1:length(allele_cat_tot)) QUAC.MSAT_subset.repRates[l] <- round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
# Assign names to values, and print
names(QUAC.MSAT_subset.repRates)<-c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%","Rare (<1%)")
print(QUAC.MSAT_subset.repRates*100)

# Allelic capture of subset SNPs ----
# Declare vector of results
QUAC.SNP_subset.repRates <- vector(length = 5)
# Create vector to capture the number of each allele category (of which there are 5)
alleles_existing_by_sp <- vector(length = 5)
# First population is all garden samples ("garden"), population 2 is "wild"
garden_p <- 1 ; wild_p <- 2
n_ind_W <- table(QUAC.SNP_subset.genind@pop)[wild_p]; n_ind_G <- table(QUAC.SNP_subset.genind@pop)[garden_p]; 
QUAC.SNP_subset.genpop <- genind2genpop(QUAC.SNP_subset.genind)
alleles_cap <- colSums(QUAC.SNP_subset.genpop[garden_p]@tab,na.rm=T)
# The below calls returns each alleles for each category specified
allele_cat_tot <- get.allele.cat(QUAC.SNP_subset.genpop[wild_p], n_ind_W)
# For each allele category, find out how many alleles there are
for (i in 1:5) alleles_existing_by_sp[i] <- (sum((allele_cat_tot[[i]])>0,na.rm=T))
# For each allele category, calculate the captured allele percentages 
# by dividing by total number of alleles in that category
for (l in 1:length(allele_cat_tot)) QUAC.SNP_subset.repRates[l] <- round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
# Assign names to values, and print
names(QUAC.SNP_subset.repRates)<-c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%","Rare (<1%)")
print(QUAC.SNP_subset.repRates*100)

# %%%% QUBO %%%% ----
# ---- MSATS: COMPLETE ----
# Read in genind file (GCC_QUBO_ZAIN repo; QUBO_wK_garden_wild_clean.gen; includes Kessler population)
# Read in genind file (Southeast Oaks repo; genetic_data/Qb_total.gen file)
genpop.filePath <- 
  "~/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
setwd(genpop.filePath)
QUBO.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_total.gen"), ncode=3)
# Correct popNames: last population (IMLS4_MP1_IMLS336_C05) is garden; rest (9) are wild
levels(QUBO.MSAT.genind@pop) <- c(rep("wild",9), "garden") 

# Declare vector of results
QUBO.MSAT.repRates <- vector(length = 5)
# Create vector to capture the number of each allele category (of which there are 5)
alleles_existing_by_sp <- vector(length = 5)
# First population is all garden samples ("garden"), population 2 is "wild"
garden_p <- 1 ; wild_p <- 2
n_ind_W <- table(QUBO.MSAT.genind@pop)[wild_p]; n_ind_G <- table(QUBO.MSAT.genind@pop)[garden_p]; 
QUBO.MSAT.genpop <- genind2genpop(QUBO.MSAT.genind)
alleles_cap <- colSums(QUBO.MSAT.genpop[garden_p]@tab,na.rm=T)
# The below calls returns each alleles for each category specified
allele_cat_tot <- get.allele.cat(QUBO.MSAT.genpop[wild_p], n_ind_W)
# For each allele category, find out how many alleles there are
for (i in 1:5) alleles_existing_by_sp[i] <- (sum((allele_cat_tot[[i]])>0,na.rm=T))
# For each allele category, calculate the captured allele percentages 
# by dividing by total number of alleles in that category
for (l in 1:length(allele_cat_tot)) QUBO.MSAT.repRates[l] <- round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
# Assign names to values, and print
names(QUBO.MSAT.repRates)<-c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%","Rare (<1%)")
print(QUBO.MSAT.repRates*100)

# ---- SNPS: COMPLETE ----
# Read in genind file (QUBO GSNAP4 alignment; R0, min-maf=0; 1 SNP/locus; 2 populations, garden and wild)
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUBO.SNP.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])

# Declare vector of results
QUBO.SNP.repRates <- vector(length = 5)
# Create vector to capture the number of each allele category (of which there are 5)
alleles_existing_by_sp <- vector(length = 5)
# First population is all garden samples ("garden"), population 2 is "wild"
garden_p <- 1 ; wild_p <- 2
n_ind_W <- table(QUBO.SNP.genind@pop)[wild_p]; n_ind_G <- table(QUBO.SNP.genind@pop)[garden_p]; 
QUBO.SNP.genpop <- genind2genpop(QUBO.SNP.genind)
alleles_cap <- colSums(QUBO.SNP.genpop[garden_p]@tab,na.rm=T)
# The below calls returns each alleles for each category specified
allele_cat_tot <- get.allele.cat(QUBO.SNP.genpop[wild_p], n_ind_W)
# For each allele category, find out how many alleles there are
for (i in 1:5) alleles_existing_by_sp[i] <- (sum((allele_cat_tot[[i]])>0,na.rm=T))
# For each allele category, calculate the captured allele percentages 
# by dividing by total number of alleles in that category
for (l in 1:length(allele_cat_tot)) QUBO.SNP.repRates[l] <- round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
# Assign names to values, and print
names(QUBO.SNP.repRates)<-c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%","Rare (<1%)")
print(QUBO.SNP.repRates*100)

# ---- MSATS AND SNPS: SUBSET ----
# MSAT: Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.MSAT.genind@tab), function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.MSAT.genind@tab) <- QUBO.MSAT.sampleNames

# SNP: Remove QUBO_W_ headers from sample names
QUBO.SNP.sampleNames <- gsub("QUBO_G_",replacement = "", row.names(QUBO.SNP.genind@tab))
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
# Rename sample matrix
rownames(QUBO.SNP.genind@tab) <- QUBO.SNP.sampleNames

# Subset SNP sample names by those that are also seen within the MSAT samples
QUBO_sharedSamples <- sort(QUBO.SNP.sampleNames[which(QUBO.SNP.sampleNames %in% QUBO.MSAT.sampleNames)])
# Subset MSAT and SNP wild matrix objects to strictly shared samples
QUBO.MSAT_subset.genind <- QUBO.MSAT.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP_subset.genind <- QUBO.SNP.genind[QUBO_sharedSamples,, drop=TRUE]

# Allelic capture of subset MSATs ----
# Declare vector of results
QUBO.MSAT_subset.repRates <- vector(length = 5)
# Create vector to capture the number of each allele category (of which there are 5)
alleles_existing_by_sp <- vector(length = 5)
# First population is all garden samples ("garden"), population 2 is "wild"
garden_p <- 1 ; wild_p <- 2
n_ind_W <- table(QUBO.MSAT_subset.genind@pop)[wild_p]; n_ind_G <- table(QUBO.MSAT_subset.genind@pop)[garden_p]; 
QUBO.MSAT_subset.genpop <- genind2genpop(QUBO.MSAT_subset.genind)
alleles_cap <- colSums(QUBO.MSAT_subset.genpop[garden_p]@tab,na.rm=T)
# The below calls returns each alleles for each category specified
allele_cat_tot <- get.allele.cat(QUBO.MSAT_subset.genpop[wild_p], n_ind_W)
# For each allele category, find out how many alleles there are
for (i in 1:5) alleles_existing_by_sp[i] <- (sum((allele_cat_tot[[i]])>0,na.rm=T))
# For each allele category, calculate the captured allele percentages 
# by dividing by total number of alleles in that category
for (l in 1:length(allele_cat_tot)) QUBO.MSAT_subset.repRates[l] <- round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
# Assign names to values, and print
names(QUBO.MSAT_subset.repRates)<-c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%","Rare (<1%)")
print(QUBO.MSAT_subset.repRates*100)

# Allelic capture of subset SNPs ----
# Declare vector of results
QUBO.SNP_subset.repRates <- vector(length = 5)
# Create vector to capture the number of each allele category (of which there are 5)
alleles_existing_by_sp <- vector(length = 5)
# First population is all garden samples ("garden"), population 2 is "wild"
garden_p <- 1 ; wild_p <- 2
n_ind_W <- table(QUBO.SNP_subset.genind@pop)[wild_p]; n_ind_G <- table(QUBO.SNP_subset.genind@pop)[garden_p]; 
QUBO.SNP_subset.genpop <- genind2genpop(QUBO.SNP_subset.genind)
alleles_cap <- colSums(QUBO.SNP_subset.genpop[garden_p]@tab,na.rm=T)
# The below calls returns each alleles for each category specified
allele_cat_tot <- get.allele.cat(QUBO.SNP_subset.genpop[wild_p], n_ind_W)
# For each allele category, find out how many alleles there are
for (i in 1:5) alleles_existing_by_sp[i] <- (sum((allele_cat_tot[[i]])>0,na.rm=T))
# For each allele category, calculate the captured allele percentages 
# by dividing by total number of alleles in that category
for (l in 1:length(allele_cat_tot)) QUBO.SNP_subset.repRates[l] <- round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
# Assign names to values, and print
names(QUBO.SNP_subset.repRates)<-c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%","Rare (<1%)")
print(QUBO.SNP_subset.repRates*100)
