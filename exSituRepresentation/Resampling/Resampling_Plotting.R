# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% RESAMPLING: PLOTTING %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in resampling arrays generated from the Resampling.R script, for 
# Quercus acerifolia (QUAC) and Quercus boyntonii (QUBO). It then generates plots for 
# those resampling analyses, which are used in the manuscript.

# There are 20 datasets that are read in, at the beginning of the script: 
# 10 for QUAC (5 Complete, 5 Subset), and 10 for QUBO (5 Complete, 5 Subset)
# Resampling plots are generated for each dataset.

library(adegenet)
library(RColorBrewer)
library(scales)
library(parallel)

# %%%% FUNCTIONS AND VARIABLES %%%% ----
# Read in relevant functions required for resampling analyses
SSRvSNP.wd <- "/home/akoontz/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)
source("exSituRepresentation/functions_exSituRepresentation.R")
# Specify path to the directory (on the lab server), where resampling arrays are held
resamplingDataDir <- paste0(SSRvSNP.wd, "exSituRepresentation/Resampling/resamplingData/")
# Specify path to the directory (on the lab server), where resampling plots (PDFs) will be saved
imageOutDir <- "/home/akoontz/Documents/SSRvSNP/Documentation/Images/MolEcol_202305_Images/"
# Plotting colors (for all plots!)
plotColors <- c("red","red4","darkorange3","coral","purple")

# %%%% MAIN TEXT PLOTS %%%% ----

# %%%% SUPPLEMENTAL PLOTS %%%% ----
# For each dataset, input (resampling arrays) and output (PDF of graph) filepaths are declared.
# Then, using this input/output, a plot is generated using the resample_plot_PDF function.
# Plots are generated according to their order in the supplement
# %%%% QUAC ----
# MSAT ----
# Complete
QUAC.MSAT.Complete.arrayPath <- paste0(resamplingDataDir,"QUAC.MSAT.Complete_resampArr.Rdata")
QUAC.MSAT.Complete.imagePath <- paste0(imageOutDir,"supplement/Fig_S8.pdf")
# Generate resampling plot
resample_singlePlot_PDF(QUAC.MSAT.Complete.arrayPath, QUAC.MSAT.Complete.imagePath,  
                  colors = plotColors, xLeg=120, minSampleLineDist=40, 
                  mainText = "QUAC, Microsatellites (Complete: 164 samples; 5,000 Replicates)")
# Subset
QUAC.MSAT.Subset.arrayPath <- paste0(resamplingDataDir,"QUAC.MSAT.Subset_resampArr.Rdata")
QUAC.MSAT.Subset.imagePath <- paste0(imageOutDir,"supplement/Fig_S9.pdf")
# Generate resampling plot
resample_singlePlot_PDF(QUAC.MSAT.Subset.arrayPath, QUAC.MSAT.Subset.imagePath,  
                  colors = plotColors, xLeg=65, minSampleLineDist=25, 
                  mainText = "QUAC, Microsatellites (Subset: 91 samples; 5,000 Replicates)")
# SNP ---
# De novo
# Complete
QUAC.SNP.DN.R0.Complete.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.DN.R0.Complete_resampArr.Rdata")
QUAC.SNP.DN.R80.Complete.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.DN.R80.Complete_resampArr.Rdata")
QUAC.SNP.DN.Complete.imagePath <- paste0(imageOutDir,"supplement/Fig_S10.pdf")
# Generate resampling plot (R0 and R80)
resample_doublePlot_PDF(QUAC.SNP.DN.R0.Complete.arrayPath, QUAC.SNP.DN.R80.Complete.arrayPath, 
                        QUAC.SNP.DN.Complete.imagePath, colors = plotColors, xLeg=120, minSampleLineDist=40,
                        mainText1 = "QUAC, SNPs: De novo, R0 (Complete: 91 samples; 5,000 Replicates)",
                        mainText2 = "QUAC, SNPs: De novo, R80 (Complete: 91 samples; 5,000 Replicates)")



# ---- SNPS (COMPLETE) ----
# %%% DE NOVO ----
# R0 ----
# *** PLOTTING ----
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUAC.SNP.DN.R0-R80.Complete.png"), width = plotWidth, height = plotHeight)
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4))
# R0
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.DN.R0.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUAC, SNPs: De novo, R0 (Complete: 91 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.DN.R0.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.DN.R0.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.DN.R0.min95_mean-10)
# R80
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.DN.R80.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUAC, SNPs: De novo, R80 (Complete: 91 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.DN.R80.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.DN.R80.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.DN.R80.min95_mean-10)
# Turn off plotting device
dev.off()

# %%% REFERENCE ----
# *** PLOTTING ----
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUAC.SNP.REF.R0-R80.Complete.png"), width = plotWidth, height = plotHeight)
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4))
# R0
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.REF.R0.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUAC, SNPs: Reference, R0 (Complete: 91 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.REF.R0.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.REF.R0.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.REF.R0.min95_mean-10)
# R80
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.REF.R80.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUAC, SNPs: Reference, R80 (Complete: 91 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.REF.R80.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.REF.R80.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.REF.R80.min95_mean-10)
# Turn off plotting device
dev.off()

# *** PLOTTING ----
# MSAT
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUAC.MSAT.Subset.png"), width = plotWidth, height = plotHeight)
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.MSAT.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUAC, Microsatellites (Subset: 91 samples; ", num_reps, " Replicates)"))
legend(x=80, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.MSAT.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.MSAT.Subset.min95_mean), 
      side=1, line=-1.5, at=QUAC.MSAT.Subset.min95_mean-15)
# Turn off plotting device
dev.off()

# SNP SUBSET
# DE NOVO
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUAC.SNP.DN.R0-R80.Subset.png"), width = plotWidth, height = plotHeight)
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4))
# R0
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.DN.R0.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUAC, SNPs: De novo, R0 (Subset: 91 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.DN.R0.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.DN.R0.Subset.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.DN.R0.Subset.min95_mean-10)
# R80
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.DN.R80.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUAC, SNPs: De novo, R80 (Subset: 91 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.DN.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.DN.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.DN.R80.Subset.min95_mean-10)
# Turn off plotting device
dev.off()

# REFERENCE
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUAC.SNP.REF.R0-R80.Subset.png"), width = plotWidth, height = plotHeight)
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4))
# R0
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.REF.R0.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUAC, SNPs: Reference, R0 (Subset: 91 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.REF.R0.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.REF.R0.Subset.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.REF.R0.Subset.min95_mean-10)
# R80
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.REF.R80.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUAC, SNPs: Reference, R80 (Subset: 91 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.REF.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.REF.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.REF.R80.Subset.min95_mean-10)
# Turn off plotting device
dev.off()

# ALL 3 SUBSET DATASETS (MSAT, SNP DE NOVO R80, SNP REFERENCE R80)
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUAC.MSAT.SNP_R80.Subset.png"), width = 1100, height = 586)
# Set plotting window to stack 3 graphs vertically
# Margin values are changed prior to each plot, in order to allow for sufficient margins at the top/bottom
par(mfcol=c(3,1), oma=rep(0.2,4), mar=c(2,4.5,4,1.5)+0.1)
# MSAT
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.MSAT.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=19,
        xlab="", ylab="")
title("Microsatellites", line=0.15, cex=2)
legend(x=75, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1.5, pt.cex = 2, bty="n", y.intersp = 0.7)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.MSAT.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.MSAT.Subset.min95_mean), 
      side=1, line=-1.5, at=QUAC.MSAT.Subset.min95_mean-10)
# Graph title 
mtext(text="QUAC: Subset Datasets (91 Samples; 5,000 Replicates; R80 SNP Loci)", 
      side=3, line=2, cex=2)
# SNP DN R80
par(mar=c(3,4.5,3,1.5)+0.1)
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.DN.R80.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=19,
        xlab="", ylab="")
title("SNPs: De novo", line=0.15, cex=2)
legend(x=75, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1.5, pt.cex = 2, bty="n", y.intersp = 0.7)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.DN.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.DN.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.DN.R80.Subset.min95_mean-10)
# Y-axis
mtext(text="Wild allelic diversity (%)", side=2, line=3, cex=1.2, srt=90)
# SNP REF R80
par(mar=c(4,4.5,2,1.5)+0.1)
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.REF.R80.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=19,
        xlab="", ylab="")
title("SNPs: Reference", line=0.15, cex=2)
legend(x=75, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1.5, pt.cex = 2, bty="n", y.intersp = 0.7)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.REF.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.REF.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.REF.R80.Subset.min95_mean-10)
# X-axis
mtext(text="Number of individuals", side=1, line=2.5, cex=1.2)
# Turn off plotting device
dev.off()

# %%%% QUBO %%%% ----
# ---- MSATS (COMPLETE) ----
# Read in genind file (Southeast Oaks repo; genetic_data/Qb_total.gen file)
genpop.filePath <- 
  "~/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
QUBO.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_total.gen"), ncode=3)
# Correct popNames: last population (IMLS4_MP1_IMLS336_C05) is garden; rest (9) are wild
levels(QUBO.MSAT.genind@pop) <- c(rep("wild",9), "garden") 
# Rename samples: Split sample names on underscore, and return 3rd element. This is the "IMLS" portion
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.MSAT.genind@tab),function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.MSAT.genind@tab) <- QUBO.MSAT.sampleNames

# CREATE RESAMPLING ARRAY AND CALCULATE SUMMARY STATISTICS
# Export genind object
clusterExport(cl, varlist = "QUBO.MSAT.genind")
# Run resampling in parallel, to generate an array that's saved to disc
QUBO.MSAT.resamplingResults <- 
  exSitu_Resample_Parallel(QUBO.MSAT.genind, cluster = cl, reps = num_reps,
                           arrayFilepath=paste0(resamplingDataDir,"QUBO.MSAT.Complete_resampArr.Rdata"))
# Close cores
# stopCluster(cl)

# Average results across replicates (slices) of the resampling array
QUBO.MSAT.meanValuesMat <- resample_meanValues(QUBO.MSAT.resamplingResults)
# Calculate and report mean minimum 95% sample size (and standard deviation)
QUBO.MSAT.min95_mean <- resample_min95_mean(QUBO.MSAT.resamplingResults)
QUBO.MSAT.min95_sd <- resample_min95_sd(QUBO.MSAT.resamplingResults)
print(c(QUBO.MSAT.min95_mean, QUBO.MSAT.min95_sd))

# *** PLOTTING ----
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUBO.MSAT.Complete.png"), width = plotWidth, height = plotHeight)

# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.MSAT.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUBO, Microsatellites (Complete: 245 samples; ", num_reps, " Replicates)"))
legend(x=200, y=86.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.MSAT.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.MSAT.min95_mean), 
      side=1, line=-1.5, at=QUBO.MSAT.min95_mean-23)
# Turn off plotting device
dev.off()

# ---- SNPS (COMPLETE) ----
# %%% DE NOVO ----
# R0 ----
# READ IN GENIND FILE (Optimized de novo assembly; R0, min-maf=0, first SNP/locus, 2 populations (garden and wild))
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R0_NOMAF_1SNP_2Pops/"
QUBO.SNP.DN.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R0.genind) <- 
  factor(read.table(paste0(genpop.filePath, "QUBO_popmap_GardenWild"), header=FALSE)[,2])
# Rename samples: get rid of "QUBO_G/W_" prefix. This will just leave the "IMLS" portion
QUBO.SNP.DN.R0.sampleNames <- gsub("QUBO_G_",replacement = "", row.names(QUBO.SNP.DN.R0.genind@tab))
QUBO.SNP.DN.R0.sampleNames <- gsub("QUBO_W_",replacement = "", QUBO.SNP.DN.R0.sampleNames)
# Replace SH-Q names in SNP list with IMLS names
# These were determined by Austin K., and are outlined on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
# Only 1 of the 11 SH_Q garden samples has an IMLS sample name (SHQ2177); others are unshared 
QUBO.SNP.DN.R0.sampleNames <- gsub("SH_Q2177",replacement = "IMLS338", QUBO.SNP.DN.R0.sampleNames)
QUBO.SNP.DN.R0.sampleNames <- gsub("SH_Q2178",replacement = "IMLS312", QUBO.SNP.DN.R0.sampleNames)
QUBO.SNP.DN.R0.sampleNames <- gsub("SH_Q2179",replacement = "IMLS062", QUBO.SNP.DN.R0.sampleNames)
QUBO.SNP.DN.R0.sampleNames <- gsub("SH_Q2180",replacement = "IMLS051", QUBO.SNP.DN.R0.sampleNames)
QUBO.SNP.DN.R0.sampleNames <- gsub("SH_Q2181",replacement = "IMLS011", QUBO.SNP.DN.R0.sampleNames)
QUBO.SNP.DN.R0.sampleNames <- gsub("SH_Q2182",replacement = "IMLS144", QUBO.SNP.DN.R0.sampleNames)
QUBO.SNP.DN.R0.sampleNames <- gsub("SH_Q2183",replacement = "IMLS170", QUBO.SNP.DN.R0.sampleNames)
QUBO.SNP.DN.R0.sampleNames <- gsub("SH_Q2184",replacement = "IMLS005", QUBO.SNP.DN.R0.sampleNames)
QUBO.SNP.DN.R0.sampleNames <- gsub("SH_Q2186",replacement = "IMLS017", QUBO.SNP.DN.R0.sampleNames)
# Rename sample matrix
rownames(QUBO.SNP.DN.R0.genind@tab) <- QUBO.SNP.DN.R0.sampleNames

# CREATE RESAMPLING ARRAY AND CALCULATE SUMMARY STATISTICS
# Export genind object
clusterExport(cl, varlist = "QUBO.SNP.DN.R0.genind")
# Run resampling in parallel, to generate an array that's saved to disc
QUBO.SNP.DN.R0.resamplingResults <- 
  exSitu_Resample_Parallel(QUBO.SNP.DN.R0.genind, cluster = cl, reps = num_reps,
                           arrayFilepath=paste0(resamplingDataDir,"QUBO.SNP.DN.R0.Complete_resampArr.Rdata"))
# Close cores
# stopCluster(cl)

# Average results across replicates (slices) of the resampling array
QUBO.SNP.DN.R0.meanValuesMat <- resample_meanValues(QUBO.SNP.DN.R0.resamplingResults)
# Calculate and report mean minimum 95% sample size (and standard deviation)
QUBO.SNP.DN.R0.min95_mean <- resample_min95_mean(QUBO.SNP.DN.R0.resamplingResults)
QUBO.SNP.DN.R0.min95_sd <- resample_min95_sd(QUBO.SNP.DN.R0.resamplingResults)
print(c(QUBO.SNP.DN.R0.min95_mean, QUBO.SNP.DN.R0.min95_sd))

# R80 ----
# READ IN GENIND FILE (Optimized de novo assembly; R80, min-maf=0, first SNP/locus, 2 populations (garden and wild))
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_R80_NOMAF_1SNP_2Pops/"
QUBO.SNP.DN.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.DN.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath, "QUBO_popmap_GardenWild"), header=FALSE)[,2])
# Rename sample matrix, based on QUBO.SNP.DN.R0.sampleNames
rownames(QUBO.SNP.DN.R80.genind@tab) <- QUBO.SNP.DN.R0.sampleNames

# CREATE RESAMPLING ARRAY AND CALCULATE SUMMARY STATISTICS
# Export genind object
clusterExport(cl, varlist = "QUBO.SNP.DN.R80.genind")
# Run resampling in parallel, to generate an array that's saved to disc
QUBO.SNP.DN.R80.resamplingResults <- 
  exSitu_Resample_Parallel(QUBO.SNP.DN.R80.genind, cluster = cl, reps = num_reps,
                           arrayFilepath=paste0(resamplingDataDir,"QUBO.SNP.DN.R80.Complete_resampArr.Rdata"))
# Close cores
# stopCluster(cl)

# Average results across replicates (slices) of the resampling array
QUBO.SNP.DN.R80.meanValuesMat <- resample_meanValues(QUBO.SNP.DN.R80.resamplingResults)
# Calculate and report mean minimum 95% sample size (and standard deviation)
QUBO.SNP.DN.R80.min95_mean <- resample_min95_mean(QUBO.SNP.DN.R80.resamplingResults)
QUBO.SNP.DN.R80.min95_sd <- resample_min95_sd(QUBO.SNP.DN.R80.resamplingResults)
print(c(QUBO.SNP.DN.R80.min95_mean, QUBO.SNP.DN.R80.min95_sd))

# *** PLOTTING ----
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUBO.SNP.DN.R0-R80.Complete.png"), width = plotWidth, height = plotHeight)
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4))
# R0
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.DN.R0.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUBO, SNPs: De novo, R0 (Complete: 95 samples; ", num_reps, " Replicates)"))
legend(x=83, y=86.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.DN.R0.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.DN.R0.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.DN.R0.min95_mean-10)
# R80
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.DN.R80.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUBO, SNPs: De novo, R80 (Complete: 95 samples; ", num_reps, " Replicates)"))
legend(x=83, y=86.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.DN.R80.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.DN.R80.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.DN.R80.min95_mean-10)
# Turn off plotting device
dev.off()

# %%% REFERENCE ----
# R0 ----
# READ IN GENIND FILE (QUBO GSNAP4 alignment; R0, min-maf=0, first SNP/locus, 2 populations (garden and wild))
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_1SNP_2Pops/"
QUBO.SNP.REF.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R0.genind) <- 
  factor(read.table(paste0(genpop.filePath, "QUBO_popmap_GardenWild"), header=FALSE)[,2])
# Rename sample matrix, based on QUBO.SNP.DN.R0.sampleNames
rownames(QUBO.SNP.REF.R0.genind@tab) <- QUBO.SNP.DN.R0.sampleNames

# CREATE RESAMPLING ARRAY AND CALCULATE SUMMARY STATISTICS
# Export genind object
clusterExport(cl, varlist = "QUBO.SNP.REF.R0.genind")
# Run resampling in parallel, to generate an array that's saved to disc
QUBO.SNP.REF.R0.resamplingResults <- 
  exSitu_Resample_Parallel(QUBO.SNP.REF.R0.genind, cluster = cl, reps = num_reps,
                           arrayFilepath=paste0(resamplingDataDir,"QUBO.SNP.REF.R0.Complete_resampArr.Rdata"))
# Close cores
# stopCluster(cl)

# Average results across replicates (slices) of the resampling array
QUBO.SNP.REF.R0.meanValuesMat <- resample_meanValues(QUBO.SNP.REF.R0.resamplingResults)
# Calculate and report mean minimum 95% sample size (and standard deviation)
QUBO.SNP.REF.R0.min95_mean <- resample_min95_mean(QUBO.SNP.REF.R0.resamplingResults)
QUBO.SNP.REF.R0.min95_sd <- resample_min95_sd(QUBO.SNP.REF.R0.resamplingResults)
print(c(QUBO.SNP.REF.R0.min95_mean, QUBO.SNP.REF.R0.min95_sd))

# R80 ----
# READ IN GENIND FILE (QUBO GSNAP4 alignment; R80, min-maf=0, first SNP/locus, 2 populations (garden and wild))
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_1SNP_2Pops/"
QUBO.SNP.REF.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.REF.R80.genind) <- 
  factor(read.table(paste0(genpop.filePath, "QUBO_popmap_GardenWild"), header=FALSE)[,2])
# Rename sample matrix, based on QUBO.SNP.DN.R0.sampleNames
rownames(QUBO.SNP.REF.R80.genind@tab) <- QUBO.SNP.DN.R0.sampleNames

# CREATE RESAMPLING ARRAY AND CALCULATE SUMMARY STATISTICS
# Export genind object
clusterExport(cl, varlist = "QUBO.SNP.REF.R80.genind")
# Run resampling in parallel, to generate an array that's saved to disc
QUBO.SNP.REF.R80.resamplingResults <- 
  exSitu_Resample_Parallel(QUBO.SNP.REF.R80.genind, cluster = cl, reps = num_reps,
                           arrayFilepath=paste0(resamplingDataDir,"QUBO.SNP.REF.R80.Complete_resampArr.Rdata"))
# Close cores
# stopCluster(cl)

# Average results across replicates (slices) of the resampling array
QUBO.SNP.REF.R80.meanValuesMat <- resample_meanValues(QUBO.SNP.REF.R80.resamplingResults)
# Calculate and report mean minimum 95% sample size (and standard deviation)
QUBO.SNP.REF.R80.min95_mean <- resample_min95_mean(QUBO.SNP.REF.R80.resamplingResults)
QUBO.SNP.REF.R80.min95_sd <- resample_min95_sd(QUBO.SNP.REF.R80.resamplingResults)
print(c(QUBO.SNP.REF.R80.min95_mean, QUBO.SNP.REF.R80.min95_sd))

# *** PLOTTING ----
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUBO.SNP.REF.R0-R80.Complete.png"), width = plotWidth, height = plotHeight)
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4))
# R0
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.REF.R0.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUBO, SNPs: Reference, R0 (Complete: 95 samples; ", num_reps, " Replicates)"))
legend(x=83, y=86.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.REF.R0.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.REF.R0.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.REF.R0.min95_mean-10)
# R80
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.REF.R80.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUBO, SNPs: Reference, R80 (Complete: 95 samples; ", num_reps, " Replicates)"))
legend(x=83, y=86.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.REF.R80.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.REF.R80.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.REF.R80.min95_mean-10)
# Turn off plotting device
dev.off()

# ---- MSATS AND SNPS: SUBSET ----
# Subset SNP sample names by those that are also seen within the MSAT samples
QUBO_sharedSamples <- sort(QUBO.SNP.DN.R0.sampleNames[which(QUBO.SNP.DN.R0.sampleNames %in% QUBO.MSAT.sampleNames)])
# Subset MSAT and SNP wild matrix objects to strictly shared samples
QUBO.MSAT_subset.genind <- QUBO.MSAT.genind[QUBO_sharedSamples,, drop=TRUE]
# De novo
QUBO.SNP.DN.R0_subset.genind <- QUBO.SNP.DN.R0.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.DN.R80_subset.genind <- QUBO.SNP.DN.R80.genind[QUBO_sharedSamples,, drop=TRUE]
# Reference
QUBO.SNP.REF.R0_subset.genind <- QUBO.SNP.REF.R0.genind[QUBO_sharedSamples,, drop=TRUE]
QUBO.SNP.REF.R80_subset.genind <- QUBO.SNP.REF.R80.genind[QUBO_sharedSamples,, drop=TRUE]

# CREATE RESAMPLING ARRAY AND CALCULATE SUMMARY STATISTICS
# Export genind objects
clusterExport(cl, varlist = c("QUBO.MSAT_subset.genind", "QUBO.SNP.DN.R0_subset.genind",
                              "QUBO.SNP.DN.R80_subset.genind", "QUBO.SNP.REF.R0_subset.genind",
                              "QUBO.SNP.REF.R80_subset.genind"))
# Run resampling in parallel, to generate arrays
# MSAT
QUBO.MSAT.Subset.resamplingResults <- 
  exSitu_Resample_Parallel(QUBO.MSAT_subset.genind, cluster = cl, reps = num_reps,
                           arrayFilepath=paste0(resamplingDataDir,"QUBO.MSAT.Subset_resampArr.Rdata"))
# SNP: De novo
# R0
QUBO.SNP.DN.R0.Subset.resamplingResults <- 
  exSitu_Resample_Parallel(QUBO.SNP.DN.R0_subset.genind, cluster = cl, reps = num_reps,
                           arrayFilepath=paste0(resamplingDataDir,"QUBO.SNP.DN.R0.Subset_resampArr.Rdata"))
# R80
QUBO.SNP.DN.R80.Subset.resamplingResults <- 
  exSitu_Resample_Parallel(QUBO.SNP.DN.R80_subset.genind, cluster = cl, reps = num_reps,
                           arrayFilepath=paste0(resamplingDataDir,"QUBO.SNP.DN.R80.Subset_resampArr.Rdata"))
# SNP: Reference
# R0
QUBO.SNP.REF.R0.Subset.resamplingResults <- 
  exSitu_Resample_Parallel(QUBO.SNP.REF.R0_subset.genind, cluster = cl, reps = num_reps,
                           arrayFilepath=paste0(resamplingDataDir,"QUBO.SNP.REF.R0.Subset_resampArr.Rdata"))
# R80
QUBO.SNP.REF.R80.Subset.resamplingResults <- 
  exSitu_Resample_Parallel(QUBO.SNP.REF.R80_subset.genind, cluster = cl, reps = num_reps,
                           arrayFilepath=paste0(resamplingDataDir,"QUBO.SNP.REF.R80.Subset_resampArr.Rdata"))
# Close cores
stopCluster(cl)

# MSAT Subset
# Average results across replicates (slices) of the resampling array
QUBO.MSAT.Subset.meanValuesMat <- resample_meanValues(QUBO.MSAT.Subset.resamplingResults)
# Calculate and report mean minimum 95% sample size (and standard deviation)
QUBO.MSAT.Subset.min95_mean <- resample_min95_mean(QUBO.MSAT.Subset.resamplingResults)
QUBO.MSAT.Subset.min95_sd <- resample_min95_sd(QUBO.MSAT.Subset.resamplingResults)
print(c(QUBO.MSAT.Subset.min95_mean, QUBO.MSAT.Subset.min95_sd))

# SNP Subset: De novo
# R0
# Average results across replicates (slices) of the resampling array
QUBO.SNP.DN.R0.Subset.meanValuesMat <- resample_meanValues(QUBO.SNP.DN.R0.Subset.resamplingResults)
# Calculate and report mean minimum 95% sample size (and standard deviation)
QUBO.SNP.DN.R0.Subset.min95_mean <- resample_min95_mean(QUBO.SNP.DN.R0.Subset.resamplingResults)
QUBO.SNP.DN.R0.Subset.min95_sd <- resample_min95_sd(QUBO.SNP.DN.R0.Subset.resamplingResults)
print(c(QUBO.SNP.DN.R0.Subset.min95_mean, QUBO.SNP.DN.R0.Subset.min95_sd))
# R80
# Average results across replicates (slices) of the resampling array
QUBO.SNP.DN.R80.Subset.meanValuesMat <- resample_meanValues(QUBO.SNP.DN.R80.Subset.resamplingResults)
# Calculate and report mean minimum 95% sample size (and standard deviation)
QUBO.SNP.DN.R80.Subset.min95_mean <- resample_min95_mean(QUBO.SNP.DN.R80.Subset.resamplingResults)
QUBO.SNP.DN.R80.Subset.min95_sd <- resample_min95_sd(QUBO.SNP.DN.R80.Subset.resamplingResults)
print(c(QUBO.SNP.DN.R80.Subset.min95_mean, QUBO.SNP.DN.R80.Subset.min95_sd))

# SNP Subset: Reference
# R0
# Average results across replicates (slices) of the resampling array
QUBO.SNP.REF.R0.Subset.meanValuesMat <- resample_meanValues(QUBO.SNP.REF.R0.Subset.resamplingResults)
# Calculate and report mean minimum 95% sample size (and standard deviation)
QUBO.SNP.REF.R0.Subset.min95_mean <- resample_min95_mean(QUBO.SNP.REF.R0.Subset.resamplingResults)
QUBO.SNP.REF.R0.Subset.min95_sd <- resample_min95_sd(QUBO.SNP.REF.R0.Subset.resamplingResults)
print(c(QUBO.SNP.REF.R0.Subset.min95_mean, QUBO.SNP.REF.R0.Subset.min95_sd))
# R80
# Average results across replicates (slices) of the resampling array
QUBO.SNP.REF.R80.Subset.meanValuesMat <- resample_meanValues(QUBO.SNP.REF.R80.Subset.resamplingResults)
# Calculate and report mean minimum 95% sample size (and standard deviation)
QUBO.SNP.REF.R80.Subset.min95_mean <- resample_min95_mean(QUBO.SNP.REF.R80.Subset.resamplingResults)
QUBO.SNP.REF.R80.Subset.min95_sd <- resample_min95_sd(QUBO.SNP.REF.R80.Subset.resamplingResults)
print(c(QUBO.SNP.REF.R80.Subset.min95_mean, QUBO.SNP.REF.R80.Subset.min95_sd))

# *** PLOTTING ----
# MSAT
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUBO.MSAT.Subset.png"), width = plotWidth, height = plotHeight)
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.MSAT.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUBO, Microsatellites (Subset: 94 samples; ", num_reps, " Replicates)"))
legend(x=80, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.MSAT.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.MSAT.Subset.min95_mean), 
      side=1, line=-1.5, at=QUBO.MSAT.Subset.min95_mean-10)
# Turn off plotting device
dev.off()

# SNP SUBSET
# DE NOVO
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUBO.SNP.DN.R0-R80.Subset.png"), width = plotWidth, height = plotHeight)
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4))
# R0
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.DN.R0.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUBO, SNPs: De novo, R0 (Subset: 94 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.DN.R0.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.DN.R0.Subset.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.DN.R0.Subset.min95_mean-10)
# R80
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.DN.R80.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUBO, SNPs: De novo, R80 (Subset: 94 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.DN.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.DN.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.DN.R80.Subset.min95_mean-10)
# Turn off plotting device
dev.off()

# REFERENCE
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUBO.SNP.REF.R0-R80.Subset.png"), width = plotWidth, height = plotHeight)
# Set plotting window to stack 2 graphs vertically
par(mfcol=c(2,1), oma=rep(0.2,4))
# R0
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.REF.R0.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUBO, SNPs: Reference, R0 (Subset: 94 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.REF.R0.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.REF.R0.Subset.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.REF.R0.Subset.min95_mean-10)
# R80
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.REF.R80.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=16,
        xlab="Number of Individuals", ylab="Percent Diversity Representation",
        main=paste0("QUBO, SNPs: Reference, R80 (Subset: 94 samples; ", num_reps, " Replicates)"))
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1, pt.cex = 2, bty="n", y.intersp = 1)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.REF.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.REF.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.REF.R80.Subset.min95_mean-10)
# Turn off plotting device
dev.off()

# ALL 3 SUBSET DATASETS (MSAT, SNP DE NOVO R80, SNP REFERENCE R80)
# Call png command, to save resampling plots to disk
png(file = paste0(resamplingDataDir, "QUBO.MSAT.SNP_R80.Subset.png"), width = 1100, height = 586)
# Set plotting window to stack 3 graphs vertically
# Margin values are changed prior to each plot, in order to allow for sufficient margins at the top/bottom
par(mfcol=c(3,1), oma=rep(0.2,4), mar=c(2,4.5,4,1.5)+0.1)
# MSAT
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.MSAT.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=19,
        xlab="", ylab="")
title("Microsatellites", line=0.15, cex=2)
legend(x=80, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1.5, pt.cex = 2, bty="n", y.intersp = 0.7)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.MSAT.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.MSAT.Subset.min95_mean), 
      side=1, line=-1.5, at=QUBO.MSAT.Subset.min95_mean-10)
# Graph title 
mtext(text="QUBO: Subset Datasets (94 Samples; 5,000 Replicates; R80 SNP Loci)", 
      side=3, line=2, cex=2)
# SNP DN R80
par(mar=c(3,4.5,3,1.5)+0.1)
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.DN.R80.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=19,
        xlab="", ylab="")
title("SNPs: De novo", line=0.15, cex=1.5)
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1.5, pt.cex = 2, bty="n", y.intersp = 0.7)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.DN.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.DN.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.DN.R80.Subset.min95_mean-10)
# Y-axis
mtext(text="Wild allelic diversity (%)", side=2, line=3, cex=1.2, srt=90)
# SNP REF R80
par(mar=c(4,4.5,2,1.5)+0.1)
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.REF.R80.Subset.meanValuesMat, ylim=c(0,110), col=plotColors, pch=19,
        xlab="", ylab="")
title("SNPs: Reference", line=0.15, cex=2)
legend(x=83, y=70.13276, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=plotColors, pch = c(20,20,20), cex=1.5, pt.cex = 2, bty="n", y.intersp = 0.7)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.REF.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.REF.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.REF.R80.Subset.min95_mean-10)
# X-axis
mtext(text="Number of individuals", side=1, line=2.5, cex=1.2)
# Turn off plotting device
dev.off()
