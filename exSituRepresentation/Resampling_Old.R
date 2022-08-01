# %%%%%%%%%%%%%%%%%%%%%%
# %%% RESAMPLING_OLD %%%
# %%%%%%%%%%%%%%%%%%%%%%

# This script uses the approach taken in the QUAC in situ/ex situ
# analyses to generate Resampling arrays, and plots their results.

library(adegenet)
library(RColorBrewer)
setwd("/home/akoontz/Documents/SSRvSNP/Code/exSituCapture/")

# %%%% FUNCTIONS %%%% ----
colMax <- function(data) sapply(data, max, na.rm = TRUE)
# get.allele.cat function, adapted from QUAC_analyses/RScripts/Fa_sample_funcs.R
get.allele.cat<-function(UK_genpop, n_ind_p_pop, n_drop=0){
  
  n_pops<-length(n_ind_p_pop)
  #--Set up categories for GLOBAL ALLELES
  allele_freqs<-colSums(UK_genpop@tab, na.rm = TRUE)/(sum(n_ind_p_pop)*2)	
  min_freq<-(n_drop/(2*sum(n_ind_p_pop)))
  glob<-as.vector(which(allele_freqs>min_freq))
  glob_v_com<-as.vector(which(allele_freqs>0.10))
  glob_com<-as.vector(which(allele_freqs>0.05))
  glob_lowfr<-as.vector(which(allele_freqs<0.10&allele_freqs>0.01))
  glob_rare<-as.vector(which((allele_freqs<0.01)&(allele_freqs>min_freq)))
  
  list(glob, glob_v_com, glob_com, glob_lowfr, glob_rare)
}

# %%%% QUAC %%%% ----
# ---- SNPS ----
# READ IN GENIND FILE (QUAC DNFA; R0, min-maf=0; 1 SNP/locus) ----

# ---- MSATS ----
# %%%% QUBO %%%% ----
# ---- SNPS ----
# READ IN GENIND FILE (QUBO GSNAP4 alignment; R0, min-maf=0; 1 SNP/locus; only 2 pops) ----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_TwoPops//"
setwd(genpop.filePath)
QUBO.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF.genind) <- factor(read.table("QUBO_popmap_GardenWild", header=FALSE)[,2])
# Create genind object of only wild individuals
QUBO.genind.sep <- seppop(QUBO.R0_NOMAF.genind)
QUBO.R0_NOMAF.Wild <- QUBO.genind.sep$wild

# GENERATE RESAMPLING RESULTS ARRAY ----
# Set number of resampling reps 
num_reps<-5
# List out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare")

# Create summary arrays for n to drop = 0 and n to drop = 2
summ_results_tree_ndrop2 <-array(dim=c(nrow(QUBO.R0_NOMAF.Wild@tab)-1, 5, num_reps))
summ_results_tree_ndrop0 <- array(dim=c(nrow(QUBO.R0_NOMAF.Wild@tab)-1, 5, num_reps)) 
# List of total number of individuals
n_total_indivs <- length(QUBO.R0_NOMAF.Wild@tab[,1])
#list of all populations
n_ind_p_pop<-table(QUBO.R0_NOMAF.Wild@pop)
#calculate total allele frequencies
allele_freqs <- colSums(QUBO.R0_NOMAF.Wild@tab, na.rm = TRUE)/(n_total_indivs*2) 
#create list for inclusion of very rare alleles or not
n_drop <- c(0,2)

# loop to resample wild trees and determine at what sample size alleles are calculated 
# Outer loop runs it for inclusion of very rare alleles or not (n to drop = 0 includes all rare alleles)
# n to drop = 2 removes very rare alleles 
# Inner loop calculates sample size for allelic resampling 
for(ndrop in n_drop){
  # Repeat the resampling based on specified number of replicates
  for (nrep in 1:num_reps){ 
    # Calculate the allelic category 
    allele_cat <- get.allele.cat(QUBO.R0_NOMAF.Wild, n_ind_p_pop, n_drop=ndrop)
    
    # Create empty matrix
    alleles_samp <- matrix(nrow=nrow(QUBO.R0_NOMAF.Wild@tab)-1,ncol=length(allele_freqs))
    
    # This loop will sample trees from t = 2 to the total number of trees
    for (t in 2:(nrow(QUBO.R0_NOMAF.Wild@tab)-1)){ ##minus one because our number of trees is 172, have to compare between 2 trees 
      
      # Create a sample of trees of length t, by using 'sample()' which randomly samples rows
      alleles_samp <- colSums(QUBO.R0_NOMAF.Wild@tab[sample(1:nrow(QUBO.R0_NOMAF.Wild@tab), t),],na.rm=T)
      
      # Add a loop to store results 
      if(ndrop == 2){
        
        # Store allele frequency summaries in an array - for dropping very rare alleles
        for (l in 1:length(allele_cat)) summ_results_tree_ndrop2[(t),(l),nrep]<-sum(alleles_samp[allele_cat[[l]]]>0, na.rm=T)
        
      }else{
        # Store allele frequency summaries in an array - including very rare alleles 
        for (l in 1:length(allele_cat)) summ_results_tree_ndrop0[(t),(l),nrep]<-sum(alleles_samp[allele_cat[[l]]]>0, na.rm=T)
      }
    }
  }  
}    

# Divide by numbers of alleles in each category to calculate % of each frequency of alleles captured
for (n in 1:num_reps) summ_results_tree_ndrop2[,,n]<-t(t(summ_results_tree_ndrop2[,,n])/summ_results_tree_ndrop2[length(summ_results_tree_ndrop2[,1,1]),,n])
for (n in 1:num_reps) summ_results_tree_ndrop0[,,n]<-t(t(summ_results_tree_ndrop0[,,n])/summ_results_tree_ndrop0[length(summ_results_tree_ndrop0[,1,1]),,n])
# Convert fractions to percentages, to make comparable to new approach
summ_results_tree_ndrop0 <- summ_results_tree_ndrop0*100
summ_results_tree_ndrop2 <- summ_results_tree_ndrop2*100

# CALCULATE MEANS AND PLOT ----

# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min_95 <- min(which(apply(summ_results_tree_ndrop0[,1,],1,mean) > 95)); min_95

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(summ_results_tree_ndrop0[,1,], 1, mean)
total_sd <- apply(summ_results_tree_ndrop0[,1,], 1, sd)

v.com_means <- apply(summ_results_tree_ndrop0[,2,], 1, mean)
v.com_sd <- apply(summ_results_tree_ndrop0[,2,], 1, sd)

com_means <- apply(summ_results_tree_ndrop0[,3,], 1, mean)
com_sd <- apply(summ_results_tree_ndrop0[,3,], 1, sd)

lowfr_means <- apply(summ_results_tree_ndrop0[,4,], 1, mean)
lowfr_sd <- apply(summ_results_tree_ndrop0[,4,], 1, sd)

rare_means <- apply(summ_results_tree_ndrop0[,5,], 1, mean)
rare_sd <- apply(summ_results_tree_ndrop0[,5,], 1, sd)
# Plots all sets of points onto single graph, as well as 95% threshold line
plotColors <- c("red","firebrick","darkorange3","coral","deeppink4")
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture",
     main="Q. boyntonii Resampling (Old): SNPs")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend("bottomright", inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.2)
# Lines for 95% threshold
abline(h=95, col="black", lty=3)
abline(v=min_95, col="black")

# ---- MSATS ----
# READ IN GENIND FILE (Southeast Oaks repo; Qb_total.gen file) ----
genpop.filePath <- 
  "~/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
setwd(genpop.filePath)
QUBO.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_total.gen"), quiet = TRUE, ncode=3)
# Create a genind object of only wild populations (which are all populations, except the last one)
QUBO.MSAT.genind.sep <- seppop(QUBO.MSAT.genind)
QUBO.MSAT.Wild.genind <- repool(QUBO.MSAT.genind.sep$IMLS1.2_MP1_IMLS002_A02,QUBO.MSAT.genind.sep$IMLS1.2_MP1_IMLS032_C10,
                                QUBO.MSAT.genind.sep$IMLS1.2_MP1_IMLS048_E02, QUBO.MSAT.genind.sep$IMLS1.2_MP1_IMLS068_F10,
                                QUBO.MSAT.genind.sep$IMLS1.2_MP1_IMLS138_G08, QUBO.MSAT.genind.sep$IMLS3_MP1_IMLS280_E04,
                                QUBO.MSAT.genind.sep$IMLS3_MP1_IMLS244_B03, QUBO.MSAT.genind.sep$IMLS3_MP1_IMLS307_G07,
                                QUBO.MSAT.genind.sep$IMLS3_MP1_IMLS312_G12)
pop(QUBO.MSAT.Wild.genind) <- rep("wild", nInd(QUBO.MSAT.Wild.genind))



# GENERATE RESAMPLING RESULTS ARRAY ----
# Set number of resampling reps 
num_reps<-50
# List out allele categories
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare")

# Create summary arrays for n to drop = 0 and n to drop = 2
summ_results_tree_ndrop2 <-array(dim=c(nrow(QUBO.MSAT.Wild.genind@tab)-1, 5, num_reps))
summ_results_tree_ndrop0 <- array(dim=c(nrow(QUBO.MSAT.Wild.genind@tab)-1, 5, num_reps)) 
# List of total number of individuals
n_total_indivs <- length(QUBO.MSAT.Wild.genind@tab[,1])
#list of all populations
n_ind_p_pop<-table(QUBO.MSAT.Wild.genind@pop)
#calculate total allele frequencies
allele_freqs <- colSums(QUBO.MSAT.Wild.genind@tab, na.rm = TRUE)/(n_total_indivs*2) 
#create list for inclusion of very rare alleles or not
n_drop <- c(0,2)

# loop to resample wild trees and determine at what sample size alleles are calculated 
# Outer loop runs it for inclusion of very rare alleles or not (n to drop = 0 includes all rare alleles)
# n to drop = 2 removes very rare alleles 
# Inner loop calculates sample size for allelic resampling 
for(ndrop in n_drop){
  # Repeat the resampling based on specified number of replicates
  for (nrep in 1:num_reps){ 
    # Calculate the allelic category 
    allele_cat <- get.allele.cat(QUBO.MSAT.Wild.genind, n_ind_p_pop, n_drop=ndrop)
    
    # Create empty matrix
    alleles_samp <- matrix(nrow=nrow(QUBO.MSAT.Wild.genind@tab)-1,ncol=length(allele_freqs))
    
    # This loop will sample trees from t = 2 to the total number of trees
    for (t in 2:(nrow(QUBO.MSAT.Wild.genind@tab)-1)){ ##minus one because our number of trees is 172, have to compare between 2 trees 
      
      # Create a sample of trees of length t, by using 'sample()' which randomly samples rows
      alleles_samp <- colSums(QUBO.MSAT.Wild.genind@tab[sample(1:nrow(QUBO.MSAT.Wild.genind@tab), t),],na.rm=T)
      
      # Add a loop to store results 
      if(ndrop == 2){
        
        # Store allele frequency summaries in an array - for dropping very rare alleles
        for (l in 1:length(allele_cat)) summ_results_tree_ndrop2[(t),(l),nrep]<-sum(alleles_samp[allele_cat[[l]]]>0, na.rm=T)
        
      }else{
        # Store allele frequency summaries in an array - including very rare alleles 
        for (l in 1:length(allele_cat)) summ_results_tree_ndrop0[(t),(l),nrep]<-sum(alleles_samp[allele_cat[[l]]]>0, na.rm=T)
      }
    }
  }  
}    

# Divide by numbers of alleles in each category to calculate % of each frequency of alleles captured
for (n in 1:num_reps) summ_results_tree_ndrop2[,,n]<-t(t(summ_results_tree_ndrop2[,,n])/summ_results_tree_ndrop2[length(summ_results_tree_ndrop2[,1,1]),,n])
for (n in 1:num_reps) summ_results_tree_ndrop0[,,n]<-t(t(summ_results_tree_ndrop0[,,n])/summ_results_tree_ndrop0[length(summ_results_tree_ndrop0[,1,1]),,n])
# Convert fractions to percentages, to make comparable to new approach
summ_results_tree_ndrop0 <- summ_results_tree_ndrop0*100
summ_results_tree_ndrop2 <- summ_results_tree_ndrop2*100

# CALCULATE MEANS AND PLOT ----

# Average results across replicates (slices) of the sampling array, to determine
# the minimum number of samples required to capture 95% wild genetic diversity
# (We average samplingResults[,1,], since this column contains the total genetic diversity)
min_95 <- min(which(apply(summ_results_tree_ndrop0[,1,],1,mean) > 95)); min_95

# Calculate means and standard deviations, for each capture rate category
total_means <- apply(summ_results_tree_ndrop0[,1,], 1, mean)
total_sd <- apply(summ_results_tree_ndrop0[,1,], 1, sd)

v.com_means <- apply(summ_results_tree_ndrop0[,2,], 1, mean)
v.com_sd <- apply(summ_results_tree_ndrop0[,2,], 1, sd)

com_means <- apply(summ_results_tree_ndrop0[,3,], 1, mean)
com_sd <- apply(summ_results_tree_ndrop0[,3,], 1, sd)

lowfr_means <- apply(summ_results_tree_ndrop0[,4,], 1, mean)
lowfr_sd <- apply(summ_results_tree_ndrop0[,4,], 1, sd)

rare_means <- apply(summ_results_tree_ndrop0[,5,], 1, mean)
rare_sd <- apply(summ_results_tree_ndrop0[,5,], 1, sd)
# Plots all sets of points onto single graph, as well as 95% threshold line
plotColors <- c("red","firebrick","darkorange3","coral","deeppink4")
plot(total_means, ylim=c(0,110), col=plotColors[1], pch=16, 
     xlab="Number of Individuals", ylab="Percent Diversity Capture",
     main="Q. boyntonii Resampling (Old): Microsatellites")
points(v.com_means, col=plotColors[2], pch=16)
points(com_means, col=plotColors[3], pch=16)
points(lowfr_means, col=plotColors[4], pch=16)
points(rare_means, col=plotColors[5], pch=16)
legend("bottomright", inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 0.2)
# Lines for 95% threshold
abline(h=95, col="black", lty=3)
abline(v=min_95, col="black")
