# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% CONVERT MSAT SUBSET GENIND FILES TO STRUCTURE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script takes genind files that have been subset to include samples shared between SNP and 
# MSAT analyses, and then uses the genind2structure function to convert the subset genind objects
# to STRUCTURE files (.str), which are then passed onto the STRUCTURE software.

# Because of computational limits, only R80 datasets are analyzed, and only wild samples are included
# in the genind object (because we're interested in the structure of wild populations, and preliminary
# analyses have shown that including garden samples can lead to ambiguous results).

# For SNP markers: SNP genind2structure steps are commented out, because actual STRUCTURE files
# (i.e. what's really used in analyses) are generated using Stacks (and specifying the required samples
# in the popmap file). Using Stacks to subset samples, rather than subsetting the genind objects here in R,
# retains more loci. The SNP genind2structure commands are kept in this script nevertheless, but are relic code.

# For QUAC: MSAT samples use "Tissue name" format, while SNP samples use "DNA name" format. To subset,
# SNP sample names are converted to Tissue name format, but after Subset samples have been determined,
# MSAT sample names are converted to DNA format (for documentation/ease of interpretation)

library(adegenet)

# Set the working directory, in order to properly read in tissue and population names
SSRvSNP.wd <- "/home/akoontz/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)
source("exSituRep/functions_exSituRepresentation.R")
# Specify folders to save STRUCTURE files to
structFilesFolder <- paste0(SSRvSNP.wd, "popAnalyses/STRUCTURE/")

# %%%% FUNCTIONS ----
# Function for converting genind to STRUCTURE 
# (website: https://github.com/lvclark/R_genetics_conv/blob/master/genind2structure.R)
genind2structure <- function(obj, file="", pops=FALSE){
  if(!"genind" %in% class(obj)){
    warning("Function was designed for genind objects.")
  }
  
  # get the max ploidy of the dataset
  pl <- max(obj@ploidy)
  # get the number of individuals
  S <- adegenet::nInd(obj)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind=rep(indNames(obj), each=pl))
  # column of pop ids to write
  if(pops){
    popnums <- 1:adegenet::nPop(obj)
    names(popnums) <- as.character(unique(adegenet::pop(obj)))
    popcol <- rep(popnums[as.character(adegenet::pop(obj))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- adegenet::locNames(obj) 
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=adegenet::nLoc(obj),
                           dimnames=list(NULL,loci)))
  
  # begin going through loci
  for(L in loci){
    thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), 
                              dimnames(obj@tab)[[2]]), 
                        drop = FALSE] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  
  # export table
  names(tab)[1] <- ""
  write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
} 

# %%%% QUAC %%%% ----
# ---- READ IN GENIND FILES ----
# MICROSATELLITE
QUAC.MSAT.genpop.filePath <- 
  "/home/akoontz/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/"
QUAC.MSAT.GW.genind <- read.genepop(paste0(QUAC.MSAT.genpop.filePath, "QUAC_woK_allpop_clean.gen"), ncode = 3)
# Correct popNames: samples with popName pattern QAc-G- are garden 
levels(QUAC.MSAT.GW.genind@pop)[grep(pattern = "QAc-G-", levels(QUAC.MSAT.GW.genind @pop))] <- 
  rep("garden", length(grep(pattern = "QAc-G-", levels(QUAC.MSAT.GW.genind @pop))))
# Correct popNames: samples with popName pattern QAc-W- are wild
levels(QUAC.MSAT.GW.genind@pop)[grep(pattern = "QAc-W-", levels(QUAC.MSAT.GW.genind @pop))] <- 
  rep("wild", length(grep(pattern = "QAc-W-", levels(QUAC.MSAT.GW.genind @pop))))
# Create a vector of MSAT tissue names, for subsetting later
QUAC.MSAT.tissueNames <- indNames(QUAC.MSAT.GW.genind)
# Create a genind that is only wild individuals
QUAC.MSAT.W.genind <- QUAC.MSAT.GW.genind[which(pop(QUAC.MSAT.GW.genind)=="wild"),, drop=TRUE]
# Create a vector of wild MSAT tissue names, for subsetting later
QUAC.MSAT.W.tissueNames <- indNames(QUAC.MSAT.W.genind)

# SNP: DE NOVO, R80, GARDEN AND WILD
genpop.filePath <-
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_gardenProv_Ordered/"
QUAC.SNP.DN.GW.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.GW.R80.genind) <- factor(read.table(paste0(genpop.filePath, "QUAC_popmap_gardenProv_Ordered"), header=FALSE)[,2])
# Capture sample names, to give the samples in the final Subset MSAT file names identical to SNP names (rather than "Tissue" names)
QUAC.SNP.GW.Subset.sampleNames <- factor(read.table(paste0(genpop.filePath, "QUAC_popmap_gardenProv_Ordered"), header=FALSE)[,1])
# Stacks population module generates a .genepop output that is order by population, rather than the
# order specified in the popmap file. To address this, we have to rearrange genetic matrix rows 
# according to the order in the popmap file.
QUAC.SNP.DN.GW.R80.genind@tab <- 
  QUAC.SNP.DN.GW.R80.genind@tab[as.character(QUAC.SNP.GW.Subset.sampleNames),,drop=FALSE]

# SNP: DE NOVO, R80, WILD
genpop.filePath <-
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R80_NOMAF_1SNP_NoK/"
QUAC.SNP.DN.W.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.W.R80.genind) <- factor(read.table(paste0(genpop.filePath, "QUAC_popmap_wild_NoK"), header=FALSE)[,2])
# Capture sample names, to give the samples in the final Subset MSAT file names identical to SNP names (rather than "Tissue" names)
QUAC.SNP.W.sampleNames <- factor(read.table(paste0(genpop.filePath, "QUAC_popmap_wild_NoK"), header=FALSE)[,1])

# SNP: REFERENCE, R80, WILD
genpop.filePath <-
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_wild_R80_NOMAF_1SNP_NoK/"
QUAC.SNP.REF.W.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.REF.W.R80.genind) <- factor(read.table(paste0(genpop.filePath, "QUAC_popmap_wild_NoK"), header=FALSE)[,2])

# ---- SUBSET GENIND FILES ----
# SNP: read in Tissue database names
# This file was created (by Austin K.), and can be found on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
QUAC.SNP.sampleTable_filepath <- paste0(SSRvSNP.wd,"exSituRepresentation/QUAC_SNP_TissueNames.csv")
QUAC.SNP.sampleTable <- read.table(QUAC.SNP.sampleTable_filepath, header = TRUE, sep = ",")
# Create vector of tissue names: return values in column 3, according to matching column 1 and
# rownames of genind object
QUAC.SNP.tissueNames <- 
  QUAC.SNP.sampleTable[match(rownames(QUAC.SNP.DN.GW.R80.genind@tab), QUAC.SNP.sampleTable[,1]),3]
# Rename genind rows to tissue names (instead of sample names)
rownames(QUAC.SNP.DN.GW.R80.genind@tab) <- QUAC.SNP.tissueNames

# For wild only analyses, remove garden names from tissue names vector, and sort
QUAC.SNP.W.tissueNames <- sort(QUAC.SNP.tissueNames[-grep(pattern = "QAc-G-", QUAC.SNP.tissueNames)])
# Rename wild SNP genind objects
rownames(QUAC.SNP.DN.W.R80.genind@tab) <- rownames(QUAC.SNP.REF.W.R80.genind@tab) <- QUAC.SNP.W.tissueNames

# Subset SNP sample names by those that are also seen within the MSAT samples
# In order to keep the sample order seen in SNP garden Prov genind files, we do not sort these names
QUAC_sharedSamples <- QUAC.MSAT.tissueNames[match(QUAC.SNP.tissueNames, QUAC.MSAT.tissueNames)]
# Subset Wild SNP sample names by those that are also seen within the MSAT samples 
QUAC_Wild_sharedSamples <- sort(QUAC.SNP.W.tissueNames[which(QUAC.SNP.W.tissueNames %in% QUAC.MSAT.W.tissueNames)])
# Subset SNP and MSAT genind matrices to strictly shared samples, dropping now absent alleles
# SNP
# Garden and Wild
QUAC.SNP.DN.GW.R80_subset.genind <- QUAC.SNP.DN.GW.R80.genind[QUAC_sharedSamples,, drop=TRUE]
# Wild
QUAC.SNP.DN.W.R80_subset.genind <- QUAC.SNP.DN.W.R80.genind[QUAC_Wild_sharedSamples,, drop=TRUE]
QUAC.SNP.REF.W.R80_subset.genind <- QUAC.SNP.REF.W.R80.genind[QUAC_Wild_sharedSamples,, drop=TRUE]

# MSAT
# Garden and Wild
QUAC.MSAT.GW_subset.genind <- QUAC.MSAT.GW.genind[QUAC_sharedSamples,, drop=TRUE]
# Wild
QUAC.MSAT.W_subset.genind <- QUAC.MSAT.W.genind[QUAC_Wild_sharedSamples,, drop=TRUE]
# Rename samples in MSAT genind files, to match sample names in SNP geninds
indNames(QUAC.MSAT.GW_subset.genind) <- QUAC.SNP.GW.Subset.sampleNames
indNames(QUAC.MSAT.W_subset.genind) <- QUAC.SNP.W.sampleNames

# ---- CONVERT GENIND TO STRUCTURE ----
# Convert complete MSAT genind objects
genind2structure(QUAC.MSAT.GW.genind, file = paste0(structFilesFolder, "QUAC.MSAT.GardenAndWild.str"))
genind2structure(QUAC.MSAT.W.genind, file = paste0(structFilesFolder, "QUAC.MSAT.Wild.str"))
# Convert subset genind objects
genind2structure(QUAC.MSAT.GW_subset.genind, file = paste0(structFilesFolder, "QUAC.MSAT.GardenAndWild_Subset.str"))
genind2structure(QUAC.MSAT.W_subset.genind, file = paste0(structFilesFolder, "QUAC.MSAT.Wild_Subset.str"))
# SNP: Commands commented out, as SNP .str files are generated using Stacks populations module
# genind2structure(QUAC.SNP.DN.GW.R80_subset.genind, file = paste0(structFilesFolder,"QUAC.SNP.DN.R80_Subset.str"))
# genind2structure(QUAC.SNP.REF.R80_subset.genind, file = paste0(structFilesFolder,"QUAC.SNP.REF.R80_Subset.str"))

# %%%% QUBO %%%% ----
setwd(SSRvSNP.wd)
# ---- READ IN GENIND FILES ----
# MICROSATELLITE
genpop.filePath <- 
  "/home/akoontz/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
QUBO.MSAT.GW.genind <- read.genepop(paste0(genpop.filePath,"Qb_total.gen"), ncode=3)
# Correct popNames: last population (IMLS4_MP1_IMLS336_C05) is garden; rest (9) are wild
levels(QUBO.MSAT.GW.genind@pop) <- c(rep("wild",9), "garden") 
# Subset to only wild individuals
QUBO.MSAT.W.genind <- QUBO.MSAT.GW.genind[which(pop(QUBO.MSAT.GW.genind)=="wild"),, drop=TRUE]

# SNP: DE NOVO, R80, WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80_NOMAF_1SNP_Subset_Ordered2/"
QUBO.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R80.genind) <- factor(read.table(paste0(genpop.filePath, "QUBO_popmap_wild_Subset_Ordered2"), header=FALSE)[,2])
# Capture sample names in wild population order specified in Stacks
# This is to ensure Subset MSAT and SNP structure files have identical sample orders (useful for plotting)
QUBO.SNP.sampleNames <- as.character(read.table(paste0(genpop.filePath, "QUBO_popmap_wild_Subset_Ordered2"), header=FALSE)[,1])
# Order the samples in the matrix of the genind file according to the sample order specified in the popmap
# During the process of reading in the genepop files, sample order changes (samples become sorted according to name)
QUBO.SNP.DN.R80.genind@tab <- QUBO.SNP.DN.R80.genind@tab[QUBO.SNP.sampleNames,]

# SNP: REFERENCE, R80, WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80_NOMAF_1SNP_Subset_Ordered2/"
QUBO.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R80.genind) <- factor(read.table(paste0(genpop.filePath,"QUBO_popmap_wild_Subset_Ordered2"), header=FALSE)[,2])
# Order the samples in the matrix of the genind file according to the sample order specified in the popmap
# During the process of reading in the genepop files, sample order changes (samples become sorted according to name)
QUBO.SNP.REF.R80.genind@tab <- QUBO.SNP.REF.R80.genind@tab[QUBO.SNP.sampleNames,]

# ---- SUBSET GENIND FILES ----
# MSAT: Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.MSAT.W.genind@tab), function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.MSAT.W.genind@tab) <- QUBO.MSAT.sampleNames

# SNP: Remove QUBO_W_ headers from sample names
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
QUBO_sharedSamples <- QUBO.SNP.sampleNames[which(QUBO.SNP.sampleNames %in% QUBO.MSAT.sampleNames)]
# Subset MSAT and SNP wild matrix objects to strictly shared samples
# MSAT
QUBO.MSAT.W_subset.genind <- QUBO.MSAT.W.genind[QUBO_sharedSamples,, drop=TRUE]
# SNP
QUBO.SNP.DN.R80_subset.genind <- QUBO.SNP.DN.R80.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.REF.R80_subset.genind <- QUBO.SNP.REF.R80.genind[QUBO_sharedSamples,, drop=TRUE]

# ---- CONVERT GENIND TO STRUCTURE ----
# Convert complete MSAT genind objects
genind2structure(QUBO.MSAT.GW.genind, file = paste0(structFilesFolder, "QUBO.MSAT.GardenAndWild.str"))
genind2structure(QUBO.MSAT.W.genind, file = paste0(structFilesFolder, "QUBO.MSAT.Wild.str"))
# Convert subset genind objects
genind2structure(QUBO.MSAT.W_subset.genind, file = paste0(structFilesFolder, "QUBO.MSAT.Wild_Subset.str"))
# SNP: Commands commented out, as SNP .str files are generated using Stacks populations module
# genind2structure(QUBO.SNP.DN.R80_subset.genind, file = paste0(structFilesFolder, "QUBO.SNP.DN.R80_Subset.str"))
# genind2structure(QUBO.SNP.REF.R80_subset.genind, file = paste0(structFilesFolder, "QUBO.SNP.REF.R80_Subset.str"))
