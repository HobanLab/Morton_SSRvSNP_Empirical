# %%%%%%%%%%%%%%%%%%%%%%%%%
# %%% REPORT FST VALUES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%

# This script generates Fst matrices from genind files and other Stacks outputs, 
# and plots these as heatmaps using the image command. Values are derived for wild populations only.

# Different methods are used to generate Fst values, as outlined below
# 1. Fst_Nei_report: uses a genind file, from which a pairwise Fst is calculated using Nei's 1987 method
# 2. Fst_WC_report: uses a genind file, from which a pairwise Fst is calculated using Weir and Cockham's 1984 method
# 3. Fst_Stacks_report: reads in the populations.fst_summary.tsv file generated from the Stacks populations module (--fstats flag)

library(adegenet)
library(hierfstat)
# Set working directory
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)

# %%%% FUNCTIONS %%%% ----
# Function for plotting Fst values using a heatmap, given a matrix of values (and a title)
# The lower half of the fst_mat argument needs to have NA values, as does the diagonal 
Fst_plot <- function(fst_mat, title="Fst Values"){
  # Use image to plot a heatmap. First two arguments specify the boundaries of the heatmap; z provides actual values
  # z is transposed in order to plot numeral values later on
  image(x=1:ncol(fst_mat), y=1:nrow(fst_mat), z=t(fst_mat), axes=FALSE, xlab="", ylab="", 
        main=title)
  # Add boundary lines
  grid(nx=ncol(fst_mat), ny=nrow(fst_mat), col="black", lty=1)
  # Include sample names, to understand the context of genetic distances
  axis(1, 1:ncol(fst_mat), colnames(fst_mat), cex.axis=1.2, tick=FALSE)
  text(1, c(1:nrow(fst_mat)), labels=rownames(fst_mat), cex=1.2)
  # Go through matrix, and plot values on each cell
  for(x in 1:ncol(fst_mat)){
    for(y in 1:nrow(fst_mat)){
      text(x, y, fst_mat[y,x], cex=1.5)
    }
  }
}

# Function for calculating and plotting Fst values (Nei, 1987), from genind object (using hierfstat package)
Fst_Nei_report <- function(gen.obj, title){
  # Convert genind object to hierfstat format
  hierfstat.obj <- genind2hierfstat(gen.obj)
  # Calculate pairwise Fst (Nei, 1987)
  fst_mat <- pairwise.neifst(hierfstat.obj)
  # pairwise.neifst returns a full matrix, whereas we need only the upper half
  # The line below gets the matrix into the format needed for plotting (lower triangle values removed)
  fst_mat[lower.tri(fst_mat, diag = TRUE)] <- NA
  
  # Plot using the Fst_plot command
  Fst_plot(fst_mat, title=title)
  # Return matrix
  return(fst_mat)
}

# Function for calculating and plotting Fst values (Weir & Cockham, 1984), from genind object (using hierfstat package)
Fst_WC_report <- function(gen.obj, title){
  # Convert genind object to hierfstat format
  hierfstat.obj <- genind2hierfstat(gen.obj)
  # Calculate pairwise Fst (Weir & Cockham, 1984)
  fst_mat <- pairwise.WCfst(hierfstat.obj)
  # pairwise.neifst returns a full matrix, whereas we need only the upper half
  # The line below gets the matrix into the format needed for plotting (lower triangle values removed)
  fst_mat[lower.tri(fst_mat, diag = TRUE)] <- NA
  # Round the values in the Fst matrix to the 4th digit, to make it comparable to Nei values
  fst_mat <- round(fst_mat, digits = 4)
  
  # Plot using the Fst_plot command
  Fst_plot(fst_mat, title=title)
  # Return matrix
  return(fst_mat)
}

# Function for calculating and plotting Fst values, using hierfstat package
Fst_Stacks_report <- function(filepath.fst_tab, title){
  # Read the populations.fst_summary.tsv file in the specified directory
  fst_mat <- as.matrix(read.table(paste0(filepath.fst_tab,"populations.fst_summary.tsv"), 
                                  header=TRUE, row.names=1, sep = "\t"))
  # Add a row at the bottom to make matrix symmetrical (nrow=ncol)
  fst_mat <- rbind(fst_mat, rep(NA, ncol(fst_mat)))
  # Update row names based on column names
  rownames(fst_mat) <- colnames(fst_mat)
  # Round the values in the Fst matrix to the 4th digit, to make it comparable to Nei values
  fst_mat <- round(fst_mat, digits = 4)
  
  # Plot using the Fst_plot command
  Fst_plot(fst_mat, title=title)
  # Return matrix
  return(fst_mat)
}

# %%%% QUAC %%%% ----
# R0 ----
# Read in genind file: Optimized de novo assembly; R0, min-maf=0, first SNP/locus, wild samples only
QUAC.SNP.R0.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF_1SNP/"
setwd(QUAC.SNP.R0.filePath)
QUAC.SNP.R0.genind <- read.genepop(paste0(QUAC.SNP.R0.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.R0.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# Generate matrix and plot hierfstat (Nei) Fst values--COMMENTED OUT BECAUSE OF VERY LONG PROCESSING TIME
# Fst_Nei_report(QUAC.SNP.R0.genind, title = "QUAC Fst (Nei) Values: De novo, R0, NOMAF (108,409 loci)")
# Generate matrix and plot hierfstat (Weir and Cockham) Fst values--COMMENTED OUT BECAUSE OF VERY LONG PROCESSING TIME
# Fst_WC_report(QUAC.SNP.R0.genind, title = "QUAC Fst (Weir and Cockham) Values: De novo, R0, NOMAF (108,409 loci)")
# Generate matrix and plot Stacks Fst values
Fst_Stacks_report(QUAC.SNP.R0.filePath, title = "QUAC Fst (Stacks) Values: De novo, R0, NOMAF (108,409 loci)")

# R80 ----
# Read in genind file: Optimized de novo assembly; R80, min-maf=0, first SNP/locus, wild samples only
QUAC.SNP.R80.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R80_NOMAF_1SNP/"
setwd(QUAC.SNP.R80.filePath)
QUAC.SNP.R80.genind <- read.genepop(paste0(QUAC.SNP.R80.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.SNP.R80.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# Generate matrix and plot hierfstat (Nei) Fst values
Fst_Nei_report(QUAC.SNP.R80.genind, title = "QUAC Fst (Nei) Values: De novo, R80, NOMAF (6,548 loci)")
# Generate matrix and plot hierfstat (Weir and Cockham) Fst values
Fst_WC_report(QUAC.SNP.R80.genind, title = "QUAC Fst (Weir and Cockham) Values: De novo, R80, NOMAF (6,548 loci)")
# Generate matrix and plot Stacks Fst values
Fst_Stacks_report(QUAC.SNP.R80.filePath, title = "QUAC Fst (Stacks) Values: De novo, R80, NOMAF (6,548 loci)")

# %%%% QUBO %%%% ----
# R0 ----
# Read in genind file: QUBO GSNAP4 alignment; R0, min-maf=0, first SNP/locus, wild samples only
QUBO.SNP.R0.filePath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0_NOMAF_1SNP/"
setwd(QUBO.SNP.R0.filePath)
QUBO.SNP.R0.genind <- read.genepop(paste0(QUBO.SNP.R0.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.R0.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# Generate matrix and plot hierfstat (Nei) Fst values--COMMENTED OUT BECAUSE OF VERY LONG PROCESSING TIME
# Fst_Nei_report(QUBO.SNP.R0.genind, title = "QUBO Fst (Nei) Values: De novo, R0, NOMAF (61,800 loci)")
# Generate matrix and plot hierfstat (Weir and Cockham) Fst values--COMMENTED OUT BECAUSE OF VERY LONG PROCESSING TIME
# Fst_WC_report(QUBO.SNP.R0.genind, title = "QUBO Fst (Weir and Cockham) Values: De novo, R0, NOMAF (61,800 loci)")
# Generate matrix and plot Fst values
Fst_Stacks_report(QUBO.SNP.R0.filePath, title = "QUBO Fst (Stacks) Values: Q. robur reference alignment, R0, NOMAF (61,800 loci)")

# R80 ----
# Read in genind file: QUBO GSNAP4 alignment; R80, min-maf=0, first SNP/locus, wild samples only
QUBO.SNP.R80.filePath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80_NOMAF_1SNP/"
setwd(QUBO.SNP.R80.filePath)
QUBO.SNP.R80.genind <- read.genepop(paste0(QUBO.SNP.R80.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.SNP.R80.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# Generate matrix and plot hierfstat (Nei) Fst values
Fst_Nei_report(QUBO.SNP.R80.genind, title = "QUBO Fst (Nei) Values: De novo, R80, NOMAF (6,174 loci)")
# Generate matrix and plot hierfstat (Weir and Cockham) Fst values
Fst_WC_report(QUBO.SNP.R80.genind, title = "QUBO Fst (Weir and Cockham) Values: De novo, R80, NOMAF (6,174 loci)")
# Generate matrix and plot Fst values
Fst_Stacks_report(QUBO.SNP.R80.filePath, title = "QUBO Fst (Stacks) Values: Q. robur reference alignment, R80, NOMAF (6,174 loci)")
