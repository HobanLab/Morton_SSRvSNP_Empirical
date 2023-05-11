# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% RESAMPLING: PLOTTING %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in resampling arrays generated from the Resampling.R script, for 
# Quercus acerifolia (QUAC) and Quercus boyntonii (QUBO). It then generates plots for 
# those resampling analyses, which are used in the manuscript.

# There are 20 reasmpling datasets that are utilized in this script: 
# 10 for QUAC (5 Complete, 5 Subset), and 10 for QUBO (5 Complete, 5 Subset)
# Resampling plots are generated for each dataset. Some of these are used in the main
# text of the manuscript; all of them are presented in the supplemental material.

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
# Two plots are generated (one for QUAC, one for QUBO), in which there are 3 sets of resampling curves:
# one for microsatellites, one for SNP De novo (R80), and one for SNP Reference (R80). All are Subset.
# %%%% QUAC ----
# READ IN AND PROCESS DATASETS ----
# MSAT
QUAC.MSAT.Subset.arrayPath <- paste0(resamplingDataDir,"QUAC.MSAT.Subset_resampArr.Rdata")
QUAC.MSAT.Subset.resamplingArr <- readRDS(QUAC.MSAT.Subset.arrayPath)
# Average results across replicates (slices) of the resampling array, and mean minimum 95% sampling size
QUAC.MSAT.Subset.meanValuesMat <- resample_meanValues(QUAC.MSAT.Subset.resamplingArr)
QUAC.MSAT.Subset.min95_mean <- resample_min95_mean(QUAC.MSAT.Subset.resamplingArr)
# SNP, De novo (R80)
QUAC.SNP.DN.R80.Subset.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.DN.R80.Subset_resampArr.Rdata")
QUAC.SNP.DN.R80.Subset.resamplingArr <- readRDS(QUAC.SNP.DN.R80.Subset.arrayPath)
# Average results across replicates (slices) of the resampling array, and mean minimum 95% sampling size
QUAC.SNP.DN.R80.Subset.meanValuesMat <- resample_meanValues(QUAC.SNP.DN.R80.Subset.resamplingArr)
QUAC.SNP.DN.R80.Subset.min95_mean <- resample_min95_mean(QUAC.SNP.DN.R80.Subset.resamplingArr)
# SNP, Reference (R80)
QUAC.SNP.REF.R80.Subset.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.REF.R80.Subset_resampArr.Rdata")
QUAC.SNP.REF.R80.Subset.resamplingArr <- readRDS(QUAC.SNP.REF.R80.Subset.arrayPath)
# Average results across replicates (slices) of the resampling array, and mean minimum 95% sampling size
QUAC.SNP.REF.R80.Subset.meanValuesMat <- resample_meanValues(QUAC.SNP.REF.R80.Subset.resamplingArr)
QUAC.SNP.REF.R80.Subset.min95_mean <- resample_min95_mean(QUAC.SNP.REF.R80.Subset.resamplingArr)

# PLOTTING ----
# Call pdf command, to save resampling plots to disk
pdf(file = paste0(imageOutDir, "mainText/Fig_2.pdf"), width = 9, height = 7.5)
# Set plotting window to stack 3 graphs vertically
# Margin values are changed prior to each plot, in order to allow for sufficient margins at the top/bottom
par(mfcol=c(3,1), oma=rep(0.2,4), mar=c(2,4.5,4,1.5)+0.1)
# Create two vectors for colors. This is to show points on the graph and in the legend clearly
fullColors <- plotColors
fadedColors <- c(plotColors[1], alpha(plotColors[2:5], 0.5))
# MSAT
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.MSAT.Subset.meanValuesMat, ylim=c(0,110), col=fadedColors, pch=16, ylab="")
title("Microsatellites", line=0.15, cex=2)
legend(x=75, y=76, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=fullColors, pch = c(20,20,20), cex=1.2, pt.cex = 2, bty="n", y.intersp = 1.3)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.MSAT.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.MSAT.Subset.min95_mean), 
      side=1, line=-1.5, at=QUAC.MSAT.Subset.min95_mean-18, cex = 1)
# Graph title 
mtext(text="QUAC: Subset Datasets (91 Samples)", 
      side=3, line=2, cex=2)
# SNP DN R80
par(mar=c(3,4.5,3,1.5)+0.1)
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.DN.R80.Subset.meanValuesMat, ylim=c(0,110), col=fadedColors, pch=16, ylab="")
title("SNPs: De novo (R80)", line=0.15, cex=2)
legend(x=75, y=76, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=fullColors, pch = c(20,20,20), cex=1.2, pt.cex = 2, bty="n", y.intersp = 1.3)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.DN.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.DN.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.DN.R80.Subset.min95_mean-18, cex = 1)
# Y-axis
mtext(text="Wild allelic diversity (%)", side=2, line=3, cex=1.2, srt=90)
# SNP REF R80
par(mar=c(4,4.5,2,1.5)+0.1)
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUAC.SNP.REF.R80.Subset.meanValuesMat, ylim=c(0,110), col=fadedColors, pch=16, ylab="")
title("SNPs: Reference (R80)", line=0.15, cex=2)
legend(x=75, y=76, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=fullColors, pch = c(20,20,20), cex=1.2, pt.cex = 2, bty="n", y.intersp = 1.3)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUAC.SNP.REF.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUAC.SNP.REF.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUAC.SNP.REF.R80.Subset.min95_mean-18, cex = 1)
# X-axis
mtext(text="Number of individuals", side=1, line=2.5, cex=1.2)
# Turn off plotting device
dev.off()
# %%%% QUBO ----
# READ IN AND PROCESS DATASETS ----
# MSAT
QUBO.MSAT.Subset.arrayPath <- paste0(resamplingDataDir,"QUBO.MSAT.Subset_resampArr.Rdata")
QUBO.MSAT.Subset.resamplingArr <- readRDS(QUBO.MSAT.Subset.arrayPath)
# Average results across replicates (slices) of the resampling array, and mean minimum 95% sampling size
QUBO.MSAT.Subset.meanValuesMat <- resample_meanValues(QUBO.MSAT.Subset.resamplingArr)
QUBO.MSAT.Subset.min95_mean <- resample_min95_mean(QUBO.MSAT.Subset.resamplingArr)
# SNP, De novo (R80)
QUBO.SNP.DN.R80.Subset.arrayPath <- paste0(resamplingDataDir,"QUBO.SNP.DN.R80.Subset_resampArr.Rdata")
QUBO.SNP.DN.R80.Subset.resamplingArr <- readRDS(QUBO.SNP.DN.R80.Subset.arrayPath)
# Average results across replicates (slices) of the resampling array, and mean minimum 95% sampling size
QUBO.SNP.DN.R80.Subset.meanValuesMat <- resample_meanValues(QUBO.SNP.DN.R80.Subset.resamplingArr)
QUBO.SNP.DN.R80.Subset.min95_mean <- resample_min95_mean(QUBO.SNP.DN.R80.Subset.resamplingArr)
# SNP, Reference (R80)
QUBO.SNP.REF.R80.Subset.arrayPath <- paste0(resamplingDataDir,"QUBO.SNP.REF.R80.Subset_resampArr.Rdata")
QUBO.SNP.REF.R80.Subset.resamplingArr <- readRDS(QUBO.SNP.REF.R80.Subset.arrayPath)
# Average results across replicates (slices) of the resampling array, and mean minimum 95% sampling size
QUBO.SNP.REF.R80.Subset.meanValuesMat <- resample_meanValues(QUBO.SNP.REF.R80.Subset.resamplingArr)
QUBO.SNP.REF.R80.Subset.min95_mean <- resample_min95_mean(QUBO.SNP.REF.R80.Subset.resamplingArr)

# PLOTTING ----
# Call pdf command, to save resampling plots to disk
pdf(file = paste0(imageOutDir, "mainText/Fig_3.pdf"), width = 9, height = 7.5)
# Set plotting window to stack 3 graphs vertically
# Margin values are changed prior to each plot, in order to allow for sufficient margins at the top/bottom
par(mfcol=c(3,1), oma=rep(0.2,4), mar=c(2,4.5,4,1.5)+0.1)
# Create two vectors for colors. This is to show points on the graph and in the legend clearly
fullColors <- plotColors
fadedColors <- c(plotColors[1], alpha(plotColors[2:5], 0.5))
# MSAT
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.MSAT.Subset.meanValuesMat, ylim=c(0,110), col=fadedColors, pch=16, ylab="")
title("Microsatellites", line=0.15, cex=2)
legend(x=80, y=76, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=fullColors, pch = c(20,20,20), cex=1.2, pt.cex = 2, bty="n", y.intersp = 1.3)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.MSAT.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.MSAT.Subset.min95_mean), 
      side=1, line=-1.5, at=QUBO.MSAT.Subset.min95_mean-18, cex = 1)
# Graph title 
mtext(text="QUBO: Subset Datasets (94 Samples)", 
      side=3, line=2, cex=2)
# SNP DN R80
par(mar=c(3,4.5,3,1.5)+0.1)
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.DN.R80.Subset.meanValuesMat, ylim=c(0,110), col=fadedColors, pch=16, ylab="")
title("SNPs: De novo (R80)", line=0.15, cex=2)
legend(x=80, y=76, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=fullColors, pch = c(20,20,20), cex=1.2, pt.cex = 2, bty="n", y.intersp = 1.3)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.DN.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.DN.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.DN.R80.Subset.min95_mean-18, cex = 1)
# Y-axis
mtext(text="Wild allelic diversity (%)", side=2, line=3, cex=1.2, srt=90)
# SNP REF R80
par(mar=c(4,4.5,2,1.5)+0.1)
# Plots all sets of points onto single graph, as well as 95% threshold line
matplot(QUBO.SNP.REF.R80.Subset.meanValuesMat, ylim=c(0,110), col=fadedColors, pch=16, ylab="")
title("SNPs: Reference (R80)", line=0.15, cex=2)
legend(x=80, y=76, inset = 0.05, legend = c("Total","Very common","Common","Low frequency", "Rare"),
       col=fullColors, pch = c(20,20,20), cex=1.2, pt.cex = 2, bty="n", y.intersp = 1.3)
# Lines for 95% threshold
abline(h=95, col="black", lty=3); abline(v=QUBO.SNP.REF.R80.Subset.min95_mean, col="black")
# Text for number of individuals to capture 95% threshold, with position based on variable
mtext(text=paste0("Minimum sampling size (95%) = ", QUBO.SNP.REF.R80.Subset.min95_mean), 
      side=1, line=-1.5, at=QUBO.SNP.REF.R80.Subset.min95_mean-18, cex = 1)
# X-axis
mtext(text="Number of individuals", side=1, line=2.5, cex=1.2)
# Turn off plotting device
dev.off()

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
                  colors = plotColors, xLeg=120, yLeg=70, minSampleLineDist=40, 
                  mainText = "QUAC, Microsatellites (Complete: 164 samples)")
# Subset
QUAC.MSAT.Subset.arrayPath <- paste0(resamplingDataDir,"QUAC.MSAT.Subset_resampArr.Rdata")
QUAC.MSAT.Subset.imagePath <- paste0(imageOutDir,"supplement/Fig_S9.pdf")
# Generate resampling plot
resample_singlePlot_PDF(QUAC.MSAT.Subset.arrayPath, QUAC.MSAT.Subset.imagePath,  
                  colors = plotColors, xLeg=65, yLeg=70, minSampleLineDist=25, 
                  mainText = "QUAC, Microsatellites (Subset: 91 samples)")
# SNP ---
# De novo ----
# Complete
QUAC.SNP.DN.R0.Complete.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.DN.R0.Complete_resampArr.Rdata")
QUAC.SNP.DN.R80.Complete.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.DN.R80.Complete_resampArr.Rdata")
QUAC.SNP.DN.Complete.imagePath <- paste0(imageOutDir,"supplement/Fig_S10.pdf")
# Generate resampling plot (R0 and R80)
resample_doublePlot_PDF(QUAC.SNP.DN.R0.Complete.arrayPath, QUAC.SNP.DN.R80.Complete.arrayPath, 
                        QUAC.SNP.DN.Complete.imagePath, colors = plotColors, xLeg=76, yLeg=70, minSampleLineDist=20,
                        mainText1 = "QUAC, SNPs: De novo, R0 (Complete: 91 samples)",
                        mainText2 = "QUAC, SNPs: De novo, R80 (Complete: 91 samples)")
# Subset
QUAC.SNP.DN.R0.Subset.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.DN.R0.Subset_resampArr.Rdata")
QUAC.SNP.DN.R80.Subset.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.DN.R80.Subset_resampArr.Rdata")
QUAC.SNP.DN.Subset.imagePath <- paste0(imageOutDir,"supplement/Fig_S11.pdf")
# Generate resampling plot (R0 and R80)
resample_doublePlot_PDF(QUAC.SNP.DN.R0.Subset.arrayPath, QUAC.SNP.DN.R80.Subset.arrayPath, 
                        QUAC.SNP.DN.Subset.imagePath, colors = plotColors, xLeg=76, yLeg=70, minSampleLineDist=20,
                        mainText1 = "QUAC, SNPs: De novo, R0 (Subset: 91 samples)",
                        mainText2 = "QUAC, SNPs: De novo, R80 (Subset: 91 samples)")
# Reference ----
# Complete
QUAC.SNP.REF.R0.Complete.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.REF.R0.Complete_resampArr.Rdata")
QUAC.SNP.REF.R80.Complete.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.REF.R80.Complete_resampArr.Rdata")
QUAC.SNP.REF.Complete.imagePath <- paste0(imageOutDir,"supplement/Fig_S12.pdf")
# Generate resampling plot (R0 and R80)
resample_doublePlot_PDF(QUAC.SNP.REF.R0.Complete.arrayPath, QUAC.SNP.REF.R80.Complete.arrayPath, 
                        QUAC.SNP.REF.Complete.imagePath, colors = plotColors, xLeg=76, yLeg=70, minSampleLineDist=20,
                        mainText1 = "QUAC, SNPs: Reference, R0 (Complete: 91 samples)",
                        mainText2 = "QUAC, SNPs: Reference, R80 (Complete: 91 samples)")
# Subset
QUAC.SNP.REF.R0.Subset.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.REF.R0.Subset_resampArr.Rdata")
QUAC.SNP.REF.R80.Subset.arrayPath <- paste0(resamplingDataDir,"QUAC.SNP.REF.R80.Subset_resampArr.Rdata")
QUAC.SNP.REF.Subset.imagePath <- paste0(imageOutDir,"supplement/Fig_S13.pdf")
# Generate resampling plot (R0 and R80)
resample_doublePlot_PDF(QUAC.SNP.REF.R0.Subset.arrayPath, QUAC.SNP.REF.R80.Subset.arrayPath, 
                        QUAC.SNP.REF.Subset.imagePath, colors = plotColors, xLeg=76, yLeg=70, minSampleLineDist=20,
                        mainText1 = "QUAC, SNPs: Reference, R0 (Subset: 91 samples)",
                        mainText2 = "QUAC, SNPs: Reference, R80 (Subset: 91 samples)")
# %%%% QUBO ----
# MSAT ----
# Complete
QUBO.MSAT.Complete.arrayPath <- paste0(resamplingDataDir,"QUBO.MSAT.Complete_resampArr.Rdata")
QUBO.MSAT.Complete.imagePath <- paste0(imageOutDir,"supplement/Fig_S14.pdf")
# Generate resampling plot
resample_singlePlot_PDF(QUBO.MSAT.Complete.arrayPath, QUBO.MSAT.Complete.imagePath,  
                        colors = plotColors, xLeg=190, yLeg=70, minSampleLineDist=60, 
                        mainText = "QUBO, Microsatellites (Complete: 245 samples)")
# Subset
QUBO.MSAT.Subset.arrayPath <- paste0(resamplingDataDir,"QUBO.MSAT.Subset_resampArr.Rdata")
QUBO.MSAT.Subset.imagePath <- paste0(imageOutDir,"supplement/Fig_S15.pdf")
# Generate resampling plot
resample_singlePlot_PDF(QUBO.MSAT.Subset.arrayPath, QUBO.MSAT.Subset.imagePath,  
                        colors = plotColors, xLeg=78, yLeg=40, minSampleLineDist=25, 
                        mainText = "QUBO, Microsatellites (Subset: 94 samples)")
# SNP ---
# De novo ----
# Complete
QUBO.SNP.DN.R0.Complete.arrayPath <- paste0(resamplingDataDir,"QUBO.SNP.DN.R0.Complete_resampArr.Rdata")
QUBO.SNP.DN.R80.Complete.arrayPath <- paste0(resamplingDataDir,"QUBO.SNP.DN.R80.Complete_resampArr.Rdata")
QUBO.SNP.DN.Complete.imagePath <- paste0(imageOutDir,"supplement/Fig_S16.pdf")
# Generate resampling plot (R0 and R80)
resample_doublePlot_PDF(QUBO.SNP.DN.R0.Complete.arrayPath, QUBO.SNP.DN.R80.Complete.arrayPath, 
                        QUBO.SNP.DN.Complete.imagePath, colors = plotColors, xLeg=80, yLeg=65, minSampleLineDist=25,
                        mainText1 = "QUBO, SNPs: De novo, R0 (Complete: 95 samples)",
                        mainText2 = "QUBO, SNPs: De novo, R80 (Complete: 95 samples)")
# Subset
QUBO.SNP.DN.R0.Subset.arrayPath <- paste0(resamplingDataDir,"QUBO.SNP.DN.R0.Subset_resampArr.Rdata")
QUBO.SNP.DN.R80.Subset.arrayPath <- paste0(resamplingDataDir,"QUBO.SNP.DN.R80.Subset_resampArr.Rdata")
QUBO.SNP.DN.Subset.imagePath <- paste0(imageOutDir,"supplement/Fig_S17.pdf")
# Generate resampling plot (R0 and R80)
resample_doublePlot_PDF(QUBO.SNP.DN.R0.Subset.arrayPath, QUBO.SNP.DN.R80.Subset.arrayPath, 
                        QUBO.SNP.DN.Subset.imagePath, colors = plotColors, xLeg=80, yLeg=65, minSampleLineDist=25,
                        mainText1 = "QUBO, SNPs: De novo, R0 (Subset: 94 samples)",
                        mainText2 = "QUBO, SNPs: De novo, R80 (Subset: 94 samples)")
# Reference ----
# Complete
QUBO.SNP.REF.R0.Complete.arrayPath <- paste0(resamplingDataDir,"QUBO.SNP.REF.R0.Complete_resampArr.Rdata")
QUBO.SNP.REF.R80.Complete.arrayPath <- paste0(resamplingDataDir,"QUBO.SNP.REF.R80.Complete_resampArr.Rdata")
QUBO.SNP.REF.Complete.imagePath <- paste0(imageOutDir,"supplement/Fig_S18.pdf")
# Generate resampling plot (R0 and R80)
resample_doublePlot_PDF(QUBO.SNP.REF.R0.Complete.arrayPath, QUBO.SNP.REF.R80.Complete.arrayPath, 
                        QUBO.SNP.REF.Complete.imagePath, colors = plotColors, xLeg=80, yLeg=65, minSampleLineDist=25,
                        mainText1 = "QUBO, SNPs: Reference, R0 (Complete: 95 samples)",
                        mainText2 = "QUBO, SNPs: Reference, R80 (Complete: 95 samples)")
# Subset
QUBO.SNP.REF.R0.Subset.arrayPath <- paste0(resamplingDataDir,"QUBO.SNP.REF.R0.Subset_resampArr.Rdata")
QUBO.SNP.REF.R80.Subset.arrayPath <- paste0(resamplingDataDir,"QUBO.SNP.REF.R80.Subset_resampArr.Rdata")
QUBO.SNP.REF.Subset.imagePath <- paste0(imageOutDir,"supplement/Fig_S19.pdf")
# Generate resampling plot (R0 and R80)
resample_doublePlot_PDF(QUBO.SNP.REF.R0.Subset.arrayPath, QUBO.SNP.REF.R80.Subset.arrayPath, 
                        QUBO.SNP.REF.Subset.imagePath, colors = plotColors, xLeg=80, yLeg=65, minSampleLineDist=25,
                        mainText1 = "QUBO, SNPs: Reference, R0 (Subset: 94 samples)",
                        mainText2 = "QUBO, SNPs: Reference, R80 (Subset: 94 samples)")
