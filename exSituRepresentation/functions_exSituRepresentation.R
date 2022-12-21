# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FUNCTIONS FOR EX SITU REPRESENTATION %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script declares the functions used for ex situ analyses: garden representation of
# wild alleles, and iteratively resampling allele matrices to determine minimum sample sizes.

# Load adegenet library, since some functions use the pop() accessor
library(adegenet)
# Load parallel library, for parallelized resampling functions
library(parallel)

# ---- REPRESENTATION FUNCTIONS ----
# Function for reporting representation rates, using a vector of allele frequencies and a sample matrix.
# This function assumes that the freqVector represents the absolute allele frequencies
# for the population of interest (typically, the entire wild population). Allele names 
# between the frequency vector and the sample matrix must correspond in order for values to be comparable. 
# First, the length of matches between garden alleles and wild alleles of a given category 
# is calculated (numerator). Then, the number of wild alleles of that category (denominator) 
# is calculated. From these 2 values, a percentage is calculated. This function returns (for each
# allele category) the number of alleles seen in garden samples (numerators), the number of 
# alleles seen in wild samples (denominators),and the proportion of the two (representation rates) in a matrix.
getAlleleCategories <- function(freqVector, sampleMat){
  # Determine how many total alleles in the sample matrix are found in the frequency vector 
  garden.total_Alleles <- length(which(names(freqVector) %in% colnames(sampleMat)))
  wild.total_Alleles <- length(freqVector)
  total_Percentage <- (garden.total_Alleles/wild.total_Alleles)*100
  # Very common alleles (greater than 10%)
  garden.vCom_Alleles <- length(which(names(which(freqVector > 10)) %in% colnames(sampleMat)))
  wild.vCom_Alleles <- length(which(freqVector > 10))
  vCom_Percentage <- (garden.vCom_Alleles/wild.vCom_Alleles)*100
  # Common alleles (greater than 5%)
  garden.com_Alleles <- length(which(names(which(freqVector > 5)) %in% colnames(sampleMat)))
  wild.com_Alleles <- length(which(freqVector > 5))
  com_Percentage <- (garden.com_Alleles/wild.com_Alleles)*100
  # Low frequency alleles (between 1% and 10%)
  garden.lowFr_Alleles <- length(which(names(which(freqVector < 10 & freqVector > 1)) %in% colnames(sampleMat)))
  wild.lowFr_Alleles <- length(which(freqVector < 10 & freqVector > 1))
  lowFr_Percentage <- (garden.lowFr_Alleles/wild.lowFr_Alleles)*100
  # Rare alleles (less than 1%)
  garden.rare_Alleles <- length(which(names(which(freqVector < 1)) %in% colnames(sampleMat)))
  wild.rare_Alleles <- length(which(freqVector < 1))
  rare_Percentage <- (garden.rare_Alleles/wild.rare_Alleles)*100
  # Concatenate values to vectors
  gardenAlleles <- c(garden.total_Alleles, garden.vCom_Alleles, garden.com_Alleles, garden.lowFr_Alleles, garden.rare_Alleles)
  wildAlleles <- c(wild.total_Alleles, wild.vCom_Alleles, wild.com_Alleles, wild.lowFr_Alleles, wild.rare_Alleles)
  repRates <- c(total_Percentage,vCom_Percentage,com_Percentage,lowFr_Percentage,rare_Percentage) 
  # Bind vectors to a matrix, name dimensions, and return
  exSituValues <- cbind(gardenAlleles, wildAlleles, repRates)
  rownames(exSituValues) <- c("Total","Very common (>10%)","Common (>5%)",
                              "Low frequency (1% -- 10%)","Rare (<1%)")
  colnames(exSituValues) <- c("Garden", "Wild", "Rate (%)")
  return(exSituValues)
}

# This function is a wrapper of getAlleleCategories, and takes as an argument a single genind object
# (containing both garden and wild samples to analyze). It processes that genind object to extract the
# objects it needs to calculate ex situ representation: a vector of wild allele frequencies, and a 
# sample matrix of garden samples. Missing alleles (those with allele frequencies or colSums of 0) 
# are removed from the wild allele frequency vector and the garden sample matrix, prior to 
# ex situ representation rates being calculated.
exSitu_Rep <- function(gen.obj){
  # Generate numerical vectors corresponding to garden and wild rows
  garden.Rows <- which(gen.obj@pop == "garden")
  wild.Rows <- which(gen.obj@pop == "wild")
  # Build the wild allele frequency vector: sum the allele counts, and divide by the number of wild samples
  # (which is equal to the number of wild rows in the genind matrix) times 2 (assuming diploid individuals)
  wildFreqs <- colSums(gen.obj@tab[wild.Rows,], na.rm = TRUE)/(length(wild.Rows)*2)*100
  # Remove any missing alleles (those with frequencies of 0) from the wild allele frequency vector
  wildFreqs <- wildFreqs[which(wildFreqs != 0)]
  # Generate the sample matrix to analyze (i.e. matrix of garden samples).
  gardenMat <- gen.obj@tab[garden.Rows,]
  # Remove any missing alleles (those with colSums of 0) from the matrix of garden samples
  gardenMat <- gardenMat[,which(colSums(gardenMat, na.rm = TRUE) != 0)]
  # Calculate how many alleles (of each category) the garden samples represent, and return
  repValues <- getAlleleCategories(freqVector=wildFreqs, sampleMat = gardenMat)
  return(repValues)
}

# RESAMPLING FUNCTIONS ----
# Ex situ sample function, which finds the level of ex situ representation of a sample of individuals.
# It takes a matrix (rows are samples, columns are loci) and an integer for the number of samples to 
# extract from that matrix. Using the sample and getAlleleCategories functions, it calculates 
# the allelic representation of randomly sampled rows from the  specified matrix. Prior to passing
# the required frequency vector and sample matrix to getAlleleCategories, missing loci are removed.
# Just the 3rd column of the matrix generated by getAlleleCategories (i.e. the representation rates) 
# is retained for resampling, since we don't use the absolute allele count values.
exSitu_Sample <- function(wildMat, numSamples){
  # Calculate a vector of allele frequencies, based on the total sample matrix
  freqVector <- colSums(wildMat, na.rm = TRUE)/(nrow(wildMat)*2)*100
  # Remove any missing alleles (those with frequencies of 0) from the frequency vector
  freqVector <- freqVector[which(freqVector != 0)]
  # From a matrix of individuals, select a set of random individuals (rows)
  samp <- wildMat[sample(nrow(wildMat), size=numSamples, replace = FALSE),]
  # Remove any missing alleles (those with colSums of 0) from the sample matrix
  samp <- samp[,which(colSums(samp, na.rm = TRUE) != 0)]
  # Calculate how many alleles (of each category) that sample represents
  repRates <- getAlleleCategories(freqVector, samp)
  # Subset matrix returned by getAlleleCategories to just 3rd column (representation rates), and return
  repRates <- repRates[,3]
  return(repRates)
}

# Wrapper for the exSitu_Sample function, iterating that function over the entire sample matrix
exSitu_Resample <- function(gen.obj){
  # Create a matrix of wild individuals (those with population "wild") from genind object
  wildMat <- gen.obj@tab[which(pop(gen.obj) == "wild"),]
  # Apply the exSitu_Sample function to all rows of the wild matrix
  # (except row 1, because we need at least 2 individuals to sample)
  # The resulting matrix needs to be transposed, in order to keep columns as different allele categories
  representationMatrix <- t(sapply(2:nrow(wildMat), function(x) exSitu_Sample(wildMat, x)))
  # Return the matrix of representation values
  return(representationMatrix)
}

# Wrapper for exSitu_Resample, which runs resampling in parallel over a specified cluster. Results
# are saved to a specified path
exSitu_Resample_Parallel <- function(gen.obj, cluster, reps, arrayFilepath="~/resamplingArray.Rdata"){
  # Run resampling in parallel, capturing results to an array
  resamplingArray <- 
    parSapply(cluster, 1:reps, function(a) exSitu_Resample(gen.obj = gen.obj), simplify = "array")
  # Save the resampling array object to disk, for later usage
  saveRDS(resamplingArray, file = arrayFilepath)
  # Return resampling array to global environment
  return(resamplingArray)
}

# From resampling array, calculate the mean minimum sample size to represent 95% of the total wild diversity
resample_min95_mean <- function(resamplingArray){
  # resampling array[,1,]: returns the Total column values for each replicate (3rd array dimension)
  # apply(resamplingArray[,1,],1,mean): calculates the average across replicates for each row
  # which(apply(resamplingArray[,1,],1,mean) > 95): returns the rows with averages greater than 95
  # min(which(apply(resamplingArray[,1,],1,mean) > 95)): the lowest row with an average greater than 95
  meanValue <- min(which(apply(resamplingArray[,1,],1, mean, na.rm=TRUE) > 95))
  return(meanValue)
}

# From resampling array, calculate the standard deviation, at the mean 95% value
resample_min95_sd <- function(resamplingArray){
  # Determine the mean value for representing 95% of allelic diversity
  meanValue <- resample_min95_mean(resamplingArray)
  # Calculate the standard deviation, at that mean value, and return
  sdValue <- apply(resamplingArray[,1,],1,sd)[meanValue]
  return(sdValue)
}

# From resampling array, calculate the mean values (across replicates) for each allele frequency category
resample_meanValues <- function(resamplingArray){
  # Declare a matrix to receive average values
  meanValue_mat <- matrix(nrow=nrow(resamplingArray), ncol=ncol(resamplingArray))
  # For each column in the array, average results across replicates (3rd array dimension)
  meanValue_mat[,1] <- apply(resamplingArray[,1,], 1, mean, na.rm=TRUE)
  meanValue_mat[,2] <- apply(resamplingArray[,2,], 1, mean, na.rm=TRUE)
  meanValue_mat[,3] <- apply(resamplingArray[,3,], 1, mean, na.rm=TRUE)
  meanValue_mat[,4] <- apply(resamplingArray[,4,], 1, mean, na.rm=TRUE)
  meanValue_mat[,5] <- apply(resamplingArray[,5,], 1, mean, na.rm=TRUE)
  # Give names to meanValue_mat columns, and return
  colnames(meanValue_mat) <- c("Total","Very common","Common","Low frequency","Rare")
  return(meanValue_mat)
}

# EXPLORATORY FUNCTIONS ----
# Function for generating a vector of wild allele frequencies from a genind object
getWildFreqs <- function(gen.obj){
  # Build a vector of rows corresponding to wild individuals (those that do not have a population of "garden")
  wildRows <- which(pop(gen.obj)!="garden")
  # Build the wild allele frequency vector: colSums of alleles (removing NAs), divided by number of haplotypes (Ne*2)
  wildFreqs <- colSums(gen.obj@tab[wildRows,], na.rm = TRUE)/(length(wildRows)*2)*100
  return(wildFreqs)
}

# Function for generating a vector of total allele frequencies from a genind object
getTotalFreqs <- function(gen.obj){
  # Build allele frequency vector: colSums of alleles (removing NAs), divided by number of haplotypes (Ne*2)
  totalFreqs <- colSums(gen.obj@tab, na.rm = TRUE)/(nInd(gen.obj)*2)*100
  return(totalFreqs)
}

# Exploratory function for reporting the proprtion of alleles of each category, from a (wild) frequency vector
getWildAlleleFreqProportions <- function(gen.obj){
  # Build the wild allele frequency vector, using the getWildFreqs function
  wildFreqs <- getWildFreqs(gen.obj)
  # Very common
  veryCommonAlleles <- wildFreqs[which(wildFreqs > 10)]
  veryCommon_prop <- (length(veryCommonAlleles)/length(wildFreqs))*100
  # Low frequency
  lowFrequencyAlleles <- wildFreqs[which(wildFreqs < 10 & wildFreqs > 1)]
  lowFrequency_prop <- (length(lowFrequencyAlleles)/length(wildFreqs))*100
  # Rare
  rareAlleles <- wildFreqs[which(wildFreqs < 1)]
  rare_prop <- (length(rareAlleles)/length(wildFreqs))*100
  # Build list of proportions, and return
  freqProportions <- c(veryCommon_prop, lowFrequency_prop, rare_prop)
  names(freqProportions) <- c("Very common (>10%)","Low frequency (1% -- 10%)","Rare (<1%)")
  return(freqProportions)
}

# Exploratory function for reporting the proprtion of alleles of each category, 
# from a frequency vector (of ALL alleles--garden AND wild)
getTotalAlleleFreqProportions <- function(gen.obj){
  # Build the wild allele frequency vector, using the getWildFreqs function
  totalFreqs <- getTotalFreqs(gen.obj)
  # Very common
  veryCommonAlleles <- totalFreqs[which(totalFreqs > 10)]
  veryCommon_prop <- (length(veryCommonAlleles)/length(totalFreqs))*100
  # Low frequency
  lowFrequencyAlleles <- totalFreqs[which(totalFreqs < 10 & totalFreqs > 1)]
  lowFrequency_prop <- (length(lowFrequencyAlleles)/length(totalFreqs))*100
  # Rare
  rareAlleles <- totalFreqs[which(totalFreqs < 1)]
  rare_prop <- (length(rareAlleles)/length(totalFreqs))*100
  # Build list of proportions, and return
  freqProportions <- c(veryCommon_prop, lowFrequency_prop, rare_prop)
  names(freqProportions) <- c("Very common (>10%)","Low frequency (1% -- 10%)","Rare (<1%)")
  return(freqProportions)
}

# %%%% ARCHIVED %%%% ----
# These functions were used in the data exploration phase of the SSRvSNP Empirical project, 
# but aren't utilized in final analyses.
# # ALTERNATIVE EX SITU REPRESENTATION FUNCTIONS ----
# 
# # This is an altered, outdated version of the getAlleleCategories function 
# # Characters following the underscore of allele names are stripped 
# # Therefore, only whole loci are compared across garden and wild individuals (i.e. not SNP positions/values)
# # This is why the unique() function needs to be added to each calculation
# getAlleleCategories_Partial <- function(freqVector, sampleMat){
#   # Determine how many total alleles in the sample (i.e. greater than 0) are found in the frequency vector 
#   total <- length(which(unique(names(which(freqVector > 0))) %in% unique(names(which(colSums(sampleMat, na.rm = TRUE) > 0)))))/length(unique(names(which(freqVector > 0))))*100
#   # Very common alleles (greater than 10%)
#   v_com <- length(which(unique(names(which(freqVector > 10))) %in% unique(names(which(colSums(sampleMat, na.rm = TRUE) > 0)))))/length(unique(names(which(freqVector > 10))))*100
#   # Common alleles (greater than 5%)
#   com <- length(which(unique(names(which(freqVector > 5))) %in% unique(names(which(colSums(sampleMat, na.rm = TRUE) > 0)))))/length(unique(names(which(freqVector > 5))))*100
#   # Low frequency alleles (between 1% and 10%)
#   low_freq <- length(which(unique(names(which(freqVector < 10 & freqVector > 1))) %in% unique(names(which(colSums(sampleMat, na.rm = TRUE) > 0)))))/length(unique(names(which(freqVector < 10 & freqVector > 1))))*100
#   # Rare alleles (less than 1%)
#   rare <- length(which(unique(names(which(freqVector < 1 & freqVector > 0))) %in% unique(names(which(colSums(sampleMat, na.rm = TRUE) > 0)))))/length(unique(names(which(freqVector < 1 & freqVector > 0))))*100
#   # Concatenate values to a vector, name that vector, and return
#   exSituRepRates <- c(total,v_com,com,low_freq,rare) 
#   names(exSituRepRates) <- c("Total","Very common (>10%)","Common (>5%)",
#                              "Low frequency (1% -- 10%","Rare (<1%)")
#   return(exSituRepRates)
# }
# 
# # This function is a wrapper of getAlleleCategories_Partial, using a single genind object
# # Alleles are renamed to drop information on SNP position
# reportAllelicRepresentation_Together_Partial <- function(gen.obj){
#   # Generate numerical vectors corresponding to garden and wild rows, for later calculations
#   garden.Rows <- which(gen.obj@pop == "garden")
#   wild.Rows <- which(gen.obj@pop == "wild")
#   # Rename alleles, dropping the portion of allele names following the underscore (SNP information)
#   colnames(gen.obj@tab) <- gsub(pattern = "_[0-9]{1,4}.[0-9]{1,2}", replacement = "", colnames(gen.obj@tab))
#   # Build the wild allele frequency vector
#   wildFreqs <- colSums(gen.obj@tab[wild.Rows,], na.rm = TRUE)/(length(wild.Rows)*2)*100
#   # Calculate how many alleles (of each category) the garden samples represent, and return
#   repRates <- getAlleleCategories_Partial(wildFreqs, sampleMat = gen.obj@tab[garden.Rows,])
#   return(repRates)
# }
# 
# # This function is a wrapper of getAlleleCategories, using two separate genind objects (one garden, one wild)
# # This function was motivated by the concern of grouping garden and wild samples together influencing
# # the genotyping (i.e. SNP calling) process
# reportAllelicRepresentation_Separate <- function(gen.obj.garden, gen.obj.wild){
#   # Build the wild allele frequency vector
#   wildFreqs <- colSums(gen.obj.wild@tab, na.rm = TRUE)/(nInd(gen.obj.wild)*2)*100
#   # Calculate how many alleles (of each category) the garden samples represent, and return
#   repRates <- getAlleleCategories(wildFreqs, sampleMat = gen.obj.garden@tab)
#   return(repRates)
# }
# 
# # This function is a wrapper of getAlleleCategories_Partial, using two genind objects (one garden, one wild)
# # Alleles are renamed to drop information on SNP position
# reportAllelicRepresentation_Separate_Partial <- function(gen.obj.garden, gen.obj.wild){
#   # Build the wild allele frequency vector
#   wildFreqs <- colSums(gen.obj.wild@tab, na.rm = TRUE)/(nInd(gen.obj.wild)*2)*100
#   # Rename wild frequencies, dropping the portion of allele names following the underscore
#   names(wildFreqs) <- gsub(pattern = "_[0-9]{1,4}.[0-9]{1,2}", replacement = "", names(wildFreqs))
#   # Rename garden colnames, dropping the portion of allele names following the underscore
#   colnames(gen.obj.garden@tab) <- gsub(pattern = "_[0-9]{1,4}.[0-9]{1,2}", replacement = "", colnames(gen.obj.garden@tab))
#   # Calculate how many alleles (of each category) the garden samples represent, and return
#   repRates <- getAlleleCategories_Partial(wildFreqs, sampleMat = gen.obj.garden@tab)
#   return(repRates)
# }
