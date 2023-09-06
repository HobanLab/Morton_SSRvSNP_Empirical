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
## Function for reporting representation rates, using a vector of allele frequencies and a sample matrix.
# Assumes that freqVector represents the absolute allele frequencies for the population of interest 
# (typically, entire wild population). Allele names between frequency vector and sample matrix must match! 
# 1. The length of matches between garden and wild alleles is calculated (numerator). 
# 2. The complete number of wild alleles of that category (denominator) is calculated. 
# 3. From these 2 values, a percentage is calculated. 
# This function returns the numerators, denominators, and the proportion (representation rates) in a matrix.
getAlleleCategories <- function(freqVector, sampleMat){
  # Determine how many Total alleles in the sample matrix are found in the frequency vector 
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

# Wrapper of getAlleleCategories: processes a genind object (containing wild and garden samples) 
# to extract  a vector of wild allele frequencies and a matrix of garden samples, and then calculate
# ex situ representation. Missing alleles (those with allele frequencies/colSums of 0) 
# are removed (from frequency vector and the sample matrix), prior to representation rates being calculated.
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
# Wrapper of getAlleleCategories: provided a matrix of wild samples and an integer
# specifying the number of samples to randomly draw from that matrix, will calculate a frequency
# vector (based on entire sample matrix) and a smaller matrix of samples (representing a simulated
# ex situ dataset). These objects are passed to the getAlleleCategories function, to calculate 
# representation rates (just the 3rd column is retained, no numerators or denominators).
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

# Wrapper of exSitu_Sample: iterates that function over the entire sample matrix
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

# Wrapper of exSitu_Resample: runs resampling in parallel over a specified cluster. Results
# are saved to a specified file path.
exSitu_Resample_Parallel <- function(gen.obj, cluster, reps, arrayFilepath="~/resamplingArray.Rdata"){
  # Run resampling in parallel, capturing results to an array
  resamplingArray <- 
    parSapply(cluster, 1:reps, function(a) exSitu_Resample(gen.obj = gen.obj), simplify = "array")
  # Save the resampling array object to disk, for later usage
  saveRDS(resamplingArray, file = arrayFilepath)
  # Return resampling array to global environment
  return(resamplingArray)
}

# From resampling array, calculate the mean minimum sample size to represent 95% of the Total wild diversity
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
  for(i in 1:ncol(resamplingArray)){
    meanValue_mat[,i] <- apply(resamplingArray[,i,], 1, mean, na.rm=TRUE)
  }
  # Give names to meanValue_mat columns, and return
  colnames(meanValue_mat) <- c("Total","Very common","Common","Low frequency","Rare")
  return(meanValue_mat)
}

# From resampling array, calculate the 95% quantile (across replicates) for just the Total allele frequency category
resample_quantiles <- function(resamplingArray, CI){
  # Calculate the quantiles (across number of samples, or rows) from the array, and return
  # CI argument specifies the interval used to calculate quantiles
  # Since we're only interested in Total allelic representation, subset columns to just the 1st
  quantiles_95 <- apply(resamplingArray[,1,], 1, quantile, CI)
  return(quantiles_95)
}

# From resampling array, plot the results of the resampling analysis and save to a PDF file
resample_singlePlot_PDF <- function(arrayPath, imagePath="~/", colors, xLeg, yLeg, minSampleLineDist, mainText){
  # Create two vectors for colors. This is to show points on the graph and in the legend clearly
  fullColors <- colors
  fadedColors <- c(colors[1], alpha(colors[2:5], 0.4))
  # Read in the resampling array, based on the array path argument
  resamplingArray <- readRDS(file=arrayPath)
  # Generate the average values (across replicates) for each allele frequency category 
  averageValueMat <- resample_meanValues(resamplingArray)
  # Generate the minimum sample size to represent 95% of allelic diversity (across replicates)
  min95_Value <- resample_min95_mean(resamplingArray)
  # Call pdf command, to save resampling plot to disk. 
  pdf(file = imagePath, width = 9, height = 7.5)
  # Use the matplot function to plot the matrix of average values, with specified settings
  matplot(averageValueMat, ylim=c(0,110), col=fadedColors, pch=16, ylab="Allelic Representation (%)")
  # Add title and x-axis labels to the graph
  title(main=mainText, line=0.5)
  mtext(text="Number of individuals", side=1, line=2.4)
  # Mark the 95% threshold line, as well as the 95% minimum sampling size
  abline(h=95, col="black", lty=3); abline(v=min95_Value, col="black")
  # Add text for the minimum sampling size line. Location based on min 95 value and function argument
  mtext(text=paste0("Minimum sampling size (95%) = ", min95_Value),
        side=1, line=-1.5, at=min95_Value-minSampleLineDist, cex=1)
  # Add legend
  legend(x=xLeg, y=yLeg, inset = 0.05,
         legend = c("Total","Very common","Common","Low frequency", "Rare"),
         col=fullColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty="n", y.intersp = 1)
  # Turn off plotting device
  dev.off()
}

# Plots the results of two different resampling analyses (usually R0 and R80), and saves to a PDF file
resample_doublePlot_PDF <- function(arrayPath1, arrayPath2, imagePath="~/",
                                    colors, xLeg, yLeg, minSampleLineDist, mainText1, mainText2){
  # Create two vectors for colors. This is to show points on the graph and in the legend clearly
  fullColors <- colors
  fadedColors <- c(colors[1], alpha(colors[2:5], 0.4))
  # Read in and process resampling arrays
  # %%% FIRST ARRAY
  resamplingArray1 <- readRDS(file=arrayPath1)
  averageValueMat1 <- resample_meanValues(resamplingArray1)
  min95_Value1 <- resample_min95_mean(resamplingArray1)
  # %%% SECOND ARRAY
  resamplingArray2 <- readRDS(file=arrayPath2)
  averageValueMat2 <- resample_meanValues(resamplingArray2)
  min95_Value2 <- resample_min95_mean(resamplingArray2)
  # Call pdf command, to save resampling plot to disk. 
  pdf(file = imagePath, width = 9, height = 7.5)
  # Set plotting window to stack 2 graphs vertically
  par(mfcol=c(2,1), oma=rep(0.1,4), mar=c(3,4,2,1))
  
  # %%% FIRST ARRAY
  # Use the matplot function to plot the matrix of average values, with specified settings
  matplot(averageValueMat1, ylim=c(0,110), col=fadedColors, pch=16, ylab="Allelic Representation (%)")
  # Add title and x-axis labels to the graph
  title(main=mainText1, line=0.5)
  mtext(text="Number of individuals", side=1, line=1.8)
  # Mark the 95% threshold line, as well as the 95% minimum sampling size
  abline(h=95, col="black", lty=3); abline(v=min95_Value1, col="black")
  # Add text for the minimum sampling size line. Location based on min 95 value and function argument
  mtext(text=paste0("Minimum sampling size (95%) = ", min95_Value1),
        side=1, line=-1.5, at=min95_Value1-minSampleLineDist, cex=1)
  # Add legend
  legend(x=xLeg, y=yLeg, inset = 0.05,
         legend = c("Total","Very common","Common","Low frequency", "Rare"),
         col=fullColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty="n", y.intersp = 1)
  
  # %%% SECOND ARRAY
  # Use the matplot function to plot the matrix of average values, with specified settings
  matplot(averageValueMat2, ylim=c(0,110), col=fadedColors, pch=16, ylab="Allelic Representation (%)")
  # Add title and x-axis labels to the graph
  title(main=mainText2, line=0.5)
  mtext(text="Number of individuals", side=1, line=1.8)
  # Mark the 95% threshold line, as well as the 95% minimum sampling size
  abline(h=95, col="black", lty=3); abline(v=min95_Value2, col="black")
  # Add text for the minimum sampling size line. Location based on min 95 value and function argument
  mtext(text=paste0("Minimum sampling size (95%) = ", min95_Value2),
        side=1, line=-1.5, at=min95_Value2-minSampleLineDist, cex=1)
  # Add legend
  legend(x=xLeg, y=yLeg, inset = 0.05,
         legend = c("Total","Very common","Common","Low frequency", "Rare"),
         col=fullColors, pch = c(20,20,20), cex=0.9, pt.cex = 2, bty="n", y.intersp = 1)
  
  # Turn off plotting device
  dev.off()
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
