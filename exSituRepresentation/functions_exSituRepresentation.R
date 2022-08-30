# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FUNCTIONS FOR EX SITU REPRESENTATION %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script declares the functions used for ex situ analyses: garden representation of
# wild alleles, and iteratively resampling allele matrices to determine minimum sample sizes.

# ---- FUNCTIONS ----
# Function for reporting representation rates, using a sample matrix and a vector of allele frequencies
# This function assumes that the freqVector represents the absolute allele frequencies, 
# for the entire population of interest
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
  names(exSituRepRates) <- c("Total","Very common (>10%)","Common (>5%)",
                             "Low frequency (1% -- 10%","Rare (<1%)")
  return(exSituRepRates)
}

# This is an altered, outdated version of the getAlleleCategories function 
# Characters following the underscore of allele names are stripped 
# Therefore, only whole RAD loci are compared across garden and wild individuals.
# This is why the unique() function needs to be added to each calculation
getAlleleCategories_Partial <- function(freqVector, sampleMat){
  # Determine how many total alleles in the sample (i.e. greater than 0) are found in the frequency vector 
  total <- length(which(unique(names(which(freqVector > 0))) %in% unique(names(which(colSums(sampleMat, na.rm = TRUE) > 0)))))/length(unique(names(which(freqVector > 0))))*100
  # Very common alleles (greater than 10%)
  v_com <- length(which(unique(names(which(freqVector > 10))) %in% unique(names(which(colSums(sampleMat, na.rm = TRUE) > 0)))))/length(unique(names(which(freqVector > 10))))*100
  # Common alleles (greater than 5%)
  com <- length(which(unique(names(which(freqVector > 5))) %in% unique(names(which(colSums(sampleMat, na.rm = TRUE) > 0)))))/length(unique(names(which(freqVector > 5))))*100
  # Low frequency alleles (between 1% and 10%)
  low_freq <- length(which(unique(names(which(freqVector < 10 & freqVector > 1))) %in% unique(names(which(colSums(sampleMat, na.rm = TRUE) > 0)))))/length(unique(names(which(freqVector < 10 & freqVector > 1))))*100
  # Rare alleles (less than 1%)
  rare <- length(which(unique(names(which(freqVector < 1 & freqVector > 0))) %in% unique(names(which(colSums(sampleMat, na.rm = TRUE) > 0)))))/length(unique(names(which(freqVector < 1 & freqVector > 0))))*100
  # Concatenate values to a vector, name that vector, and return
  exSituRepRates <- c(total,v_com,com,low_freq,rare) 
  names(exSituRepRates) <- c("Total","Very common (>10%)","Common (>5%)",
                             "Low frequency (1% -- 10%","Rare (<1%)")
  return(exSituRepRates)
}

# REPORT ALLELIC REPRESENTATION FUNCTIONS----
# This function is a wrapper of getAlleleCategories, using a single genind object
reportAllelicRepresentation_Together <- function(gen.obj){
  # Generate numerical vectors corresponding to garden and wild rows, for later calculations
  garden.Rows <- seq_len(length(which(pop(gen.obj)=="garden")))
  wild.Rows <- seq(from=length(which(pop(gen.obj)=="garden"))+1, to=nInd(gen.obj))
  # Build the wild allele frequency vector
  wildFreqs <- colSums(gen.obj@tab[wild.Rows,], na.rm = TRUE)/(length(wild.Rows)*2)*100
  # Calculate how many alleles (of each category) the garden samples represent, and return
  repRates <- getAlleleCategories(wildFreqs, sampleMat = gen.obj@tab[garden.Rows,])
  return(repRates)
}

# This function is a wrapper of getAlleleCategories_Partial, using a single genind object
# Alleles are renamed to drop information on SNP position
reportAllelicRepresentation_Together_Partial <- function(gen.obj){
  # Generate numerical vectors corresponding to garden and wild rows, for later calculations
  garden.Rows <- seq_len(length(which(pop(gen.obj)=="garden")))
  wild.Rows <- seq(from=length(which(pop(gen.obj)=="garden"))+1, to=nInd(gen.obj))
  # Rename alleles, dropping the portion of allele names following the underscore (SNP information)
  colnames(gen.obj@tab) <- gsub(pattern = "_[0-9]{1,4}.[0-9]{1,2}", replacement = "", colnames(gen.obj@tab))
  # Build the wild allele frequency vector
  wildFreqs <- colSums(gen.obj@tab[wild.Rows,], na.rm = TRUE)/(length(wild.Rows)*2)*100
  # Calculate how many alleles (of each category) the garden samples represent, and return
  repRates <- getAlleleCategories_Partial(wildFreqs, sampleMat = gen.obj@tab[garden.Rows,])
  return(repRates)
}

# This function is a wrapper of getAlleleCategories, using two genind objects (one garden, one wild)
reportAllelicRepresentation_Separate <- function(gen.obj.garden, gen.obj.wild){
  # Build the wild allele frequency vector
  wildFreqs <- colSums(gen.obj.wild@tab, na.rm = TRUE)/(nInd(gen.obj.wild)*2)*100
  # Calculate how many alleles (of each category) the garden samples represent, and return
  repRates <- getAlleleCategories(wildFreqs, sampleMat = gen.obj.garden@tab)
  return(repRates)
}

# This function is a wrapper of getAlleleCategories_Partial, using two genind objects (one garden, one wild)
# Alleles are renamed to drop information on SNP position
reportAllelicRepresentation_Separate_Partial <- function(gen.obj.garden, gen.obj.wild){
  # Build the wild allele frequency vector
  wildFreqs <- colSums(gen.obj.wild@tab, na.rm = TRUE)/(nInd(gen.obj.wild)*2)*100
  # Rename wild frequencies, dropping the portion of allele names following the underscore
  names(wildFreqs) <- gsub(pattern = "_[0-9]{1,4}.[0-9]{1,2}", replacement = "", names(wildFreqs))
  # Rename garden colnames, dropping the portion of allele names following the underscore
  colnames(gen.obj.garden@tab) <- gsub(pattern = "_[0-9]{1,4}.[0-9]{1,2}", replacement = "", colnames(gen.obj.garden@tab))
  # Calculate how many alleles (of each category) the garden samples represent, and return
  repRates <- getAlleleCategories_Partial(wildFreqs, sampleMat = gen.obj.garden@tab)
  return(repRates)
}

# RESAMPLING FUNCTIONS ----
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
