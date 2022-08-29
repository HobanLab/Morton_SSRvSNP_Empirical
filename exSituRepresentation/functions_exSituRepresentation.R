# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FUNCTIONS FOR EX SITU REPRESENTATION %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script declares the functions used for ex situ analyses: garden representation of
# wild alleles, and iteratively resampling allele matrices to determine minimum sample sizes.

# ---- FUNCTIONS ----
# Function for reporting representation rates, using a sample matrix and a vector of allele frequencies
# This function assumes that the freqVector represents the absolute allele frequencies, for the entire population
getAlleleCategories <- function(freqVector, sampleMat){
  # Determine how many total alleles in the sample (i.e. greater than 0) are found in the frequency vector 
  total <- length(which(names(which(freqVector > 0)) %in% names(which(colSums(sampleMat, na.rm = TRUE) > 0))))/length(which(freqVector > 0))*100
  # Very common alleles (greater than 10%)
  v_com <- length(which(names(which(freqVector > 10)) %in% names(which(colSums(sampleMat, na.rm = TRUE) > 0))))/length(which(freqVector > 10))*100
  # Common alleles (greater than 5%)
  com <- length(which(names(which(freqVector > 5)) %in% names(which(colSums(sampleMat, na.rm = TRUE) > 0))))/length(which(freqVector > 5))*100
  # Low frequency alleles (between 1% and 10%)
  low_freq <- length(which(names(which(freqVector < 10 & freqVector > 1)) %in% names(which(colSums(sampleMat, na.rm = TRUE) > 0))))/length(which(freqVector < 10 & freqVector > 1))*100
  # Rare alleles (less than 1%)
  rare <- length(which(names(which(freqVector < 1 & freqVector > 0)) %in% names(which(colSums(sampleMat, na.rm = TRUE) > 0))))/length(which(freqVector < 1 & freqVector > 0))*100
  # Concatenate values to a vector, name that vector, and return
  exSituRepRates <- c(total,v_com,com,low_freq,rare) 
  names(exSituRepRates) <- c("total","v_com","com","low_freq","rare")
  return(exSituRepRates)
}

# Ex situ sample function, which finds the level of ex situ representation of a sample of individuals
# (using the getAlleleCategories function above)
exSitu_Sample <- function(wildMat, numSamples){
  # Calculate a vector of allele frequencies, based on the total sample matrix
  freqVector <- colSums(wildMat, na.rm = TRUE)/(nrow(wildMat)*2)*100
  # From a matrix of individuals, select a set of random individuals (rows)
  samp <- wildMat[sample(nrow(wildMat), size=numSamples, replace = FALSE),]
  # Calculate how many alleles (of each category) that sample captures, and return
  repRates <- getAlleleCategories(freqVector, samp)
  return(repRates)
}

# Wrapper for the exSituSample function, iterating that function over the entire sample matrix
exSitu_Resample <- function(wildMat){
  # Apply the exSituSample function to all rows of the sample matrix
  # (except row 1, because we need at least 2 individuals to sample)
  # Resulting matrix needs to be transposed in order to keep columns as different allele categories
  representationMatrix <- t(sapply(2:nrow(wildMat), function(x) exSitu_Sample(wildMat, x)))
  # Name columns according to categories of allelic representation, and return matrix
  colnames(representationMatrix) <- c("Total","Very common","Common","Low frequency","Rare")
  return(representationMatrix)
}

# Wrapper for exSitu_Resample: Function version that works from genind object, rather than wild matrix
ExSituResample <- function(gen.obj){
  # Checks: argument is a genind object, with only 2 populations, whose names are "garden" and "wild"
  
}
