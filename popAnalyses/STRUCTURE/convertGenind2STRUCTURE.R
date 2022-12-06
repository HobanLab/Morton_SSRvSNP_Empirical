# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% CONVERT SUBSET GENIND FILES TO STRUCTURE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script takes genind files that have been subset to include samples shared between SNP and 
# MSAT analyses, and then uses the genind2structure function to convert the subset genind objects
# to STRUCTURE files (.str), which are then passed onto the STRUCTURE software.

# Because of computational limits, only R80 datasets are analyzed, and only wild samples are included
# in the genind object (because we're interested in the structure of wild populations, and preliminary
# analyses have shown that including garden samples can lead to ambiguous results).

library(adegenet)

# Set the working directory, in order to properly read in tissue and population names
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)
source("exSituRepresentation/functions_exSituRepresentation.R")

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
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/Garden_Wild/"
setwd(QUAC.MSAT.genpop.filePath)
QUAC.MSAT.genind <- read.genepop("QUAC_wK_garden_wild_clean.gen", ncode = 3)
# Correct popNames: pop1 is Garden, pop2 is Wild
pop(QUAC.MSAT.genind) <- gsub("pop1", "garden", pop(QUAC.MSAT.genind))
pop(QUAC.MSAT.genind) <- gsub("pop2", "wild", pop(QUAC.MSAT.genind))
# Subset to only wild individuals
QUAC.MSAT.genind <- QUAC.MSAT.genind[which(pop(QUAC.MSAT.genind)=="wild"),, drop=TRUE]

# SNP: DE NOVO, R80
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops/"
setwd(genpop.filePath)
QUAC.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.DN.R80.genind) <- factor(read.table("QUAC_popmap_GardenWild", header=FALSE)[,2])
# Subset to only wild individuals
QUAC.SNP.DN.R80.genind <- QUAC.SNP.DN.R80.genind[which(pop(QUAC.SNP.DN.R80.genind)=="wild"),, drop=TRUE]

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
rownames(QUAC.SNP.DN.R80.genind@tab) <- rownames(QUAC.SNP.REF.R80.genind@tab) <- QUAC.SNP.tissueNames

# Subset SNP sample names by those that are also seen within the MSAT samples 
QUAC_sharedSamples <- sort(QUAC.SNP.tissueNames[which(QUAC.SNP.tissueNames %in% QUAC.MSAT.tissueNames)])
# Subset MSAT and SNP genind matrices to strictly shared samples, dropping now absent alleles
# MSAT
QUAC.MSAT_subset.genind <- QUAC.MSAT.genind[QUAC_sharedSamples,, drop=TRUE]
# SNP
QUAC.SNP.DN.R80_subset.genind <- QUAC.SNP.DN.R80.genind[QUAC_sharedSamples,, drop=TRUE]
QUAC.SNP.REF.R80_subset.genind <- QUAC.SNP.REF.R80.genind[QUAC_sharedSamples,, drop=TRUE]

# ---- CONVERT GENIND TO STRUCTURE ----
# Save STRUCTURE files to folder on GitHub repo
setwd(paste0(SSRvSNP.wd, "popAnalyses/STRUCTURE/"))
# Convert subset genind objects
genind2structure(QUAC.MSAT_subset.genind, file = "QUAC.MSAT_Subset.str")
genind2structure(QUAC.SNP.DN.R80_subset.genind, file = "QUAC.SNP.DN.R80_Subset.str")
genind2structure(QUAC.SNP.REF.R80_subset.genind, file = "QUAC.SNP.REF.R80_Subset.str")

# %%%% QUBO %%%% ----
setwd(SSRvSNP.wd)
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

# SNP: DE NOVO, R80, WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80_NOMAF_1SNP_Subset_Ordered/"
setwd(genpop.filePath)
QUBO.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R80.genind) <- factor(read.table("QUBO_popmap_wild_Subset_Ordered", header=FALSE)[,2])

# SNP: REFERENCE, R80, WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80_NOMAF_1SNP_Subset_Ordered//"
setwd(genpop.filePath)
QUBO.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R80.genind) <- factor(read.table("QUBO_popmap_wild_Subset_Ordered", header=FALSE)[,2])

# ---- SUBSET GENIND FILES ----
# MSAT: Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.MSAT.genind@tab), function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.MSAT.genind@tab) <- QUBO.MSAT.sampleNames

# SNP: Remove QUBO_W_ headers from sample names
QUBO.SNP.sampleNames <- gsub("QUBO_W_",replacement = "", row.names(QUBO.SNP.REF.R80.genind@tab))
QUBO.SNP.sampleNames <- gsub("QUBO_W_",replacement = "", row.names(QUBO.SNP.DN.R80.genind@tab))
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
QUBO.MSAT_subset.genind <- QUBO.MSAT.genind[QUBO_sharedSamples,, drop=TRUE]
# SNP
QUBO.SNP.DN.R80_subset.genind <- QUBO.SNP.DN.R80.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.REF.R80_subset.genind <- QUBO.SNP.REF.R80.genind[QUBO_sharedSamples,, drop=TRUE]

# ---- CONVERT GENIND TO STRUCTURE ----
# Save STRUCTURE files to folder on GitHub repo
setwd(paste0(SSRvSNP.wd, "popAnalyses/STRUCTURE/"))
# Convert subset genind objects
genind2structure(QUBO.MSAT_subset.genind, file = "QUBO.MSAT_Subset.str")
genind2structure(QUBO.SNP.DN.R80_subset.genind, file = "QUBO.SNP.DN.R80_Subset.str")
genind2structure(QUBO.SNP.REF.R80_subset.genind, file = "QUBO.SNP.REF.R80_Subset.str")
