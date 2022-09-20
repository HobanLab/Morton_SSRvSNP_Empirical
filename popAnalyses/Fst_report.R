# %%%%%%%%%%%%%%%%%%%%%%%%%
# %%% PLOT FST HEATMAPS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%

# This plots Fst charts generated for de novo, reference, and hybrid Stacks datasets
# These Fst values are derived for wild populations only, 
# and are pulled from the populations.fst_summary.tsv Stacks files output from the populations module
# (generated when the --fstats flag is included). Fst values are calculated by Stacks, not in R

# Set working directory
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)

# %%%% FUNCTIONS %%%% ----
# Function for calculating and plotting Fst values, from the Stacks fst_summary.tsv file
Fst_report <- function(filepath.fst_tab, title){
  # CREATE FST MATRIX ----
  # Read the populations.fst_summary.tsv file in the specified directory
  fst_mat <- as.matrix(read.table(paste0(filepath.fst_tab,"populations.fst_summary.tsv"), header=TRUE, row.names=1, sep = "\t"))
  # Add a row at the bottom to make matrix symmetrical (nrow=ncol)
  fst_mat <- rbind(fst_mat, rep(NA, ncol(fst_mat)))
  # Update row names based on column names
  rownames(fst_mat) <- colnames(fst_mat) 
  
  # PLOT FST MATRIX ----
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
  # Return matrix
  return(fst_mat)
}

# %%%% QUAC %%%% ----
# Read in genind file: Optimized de novo assembly; R0, min-maf=0, first SNP/locus, wild samples only
QUAC.SNP.R0.filepath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF_1SNP/"
# Generate matrix and plot Fst values
Fst_report(QUAC.SNP.R0.filepath, title = "QUAC Fst Values: De novo, R0, NOMAF (108,409 loci)")

# %%%% QUBO %%%% ----
# Read in genind file: QUBO GSNAP4 alignment; R0, min-maf=0, first SNP/locus, wild samples only
QUBO.SNP.R0.filepath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0_NOMAF_1SNP/"
# Generate matrix and plot Fst values
Fst_report(QUBO.SNP.R0.filepath, title = "QUBO Fst Values: Q. robur reference alignment, R0, NOMAF (61,800 loci)")
