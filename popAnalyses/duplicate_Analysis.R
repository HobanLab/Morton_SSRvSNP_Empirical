# Duplicate analysis. Do the presence of duplicates effect the number of loci?
library(adegenet)

# QUAC
# With duplicates included
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_sum/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUAC.genind.dup <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUAC.genind.dup) <- factor(read.table("../../../QUAC_popmap2", header=FALSE)[,2])
# How many loci are present?
nLoc(QUAC.genind.dup) ; ncol(QUAC.genind.dup@tab)

# Without duplicates
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_noDuplicates/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUAC.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUAC.genind) <- factor(read.table("./QUAC_popmap_noDuplicates", header=FALSE)[,2])
# How many loci are present?
nLoc(QUAC.genind) ; ncol(QUAC.genind@tab)

# How many loci are different between the dataset with and without duplicates?
nLoc(QUAC.genind.dup) - nLoc(QUAC.genind)
# 98 more loci in the dataset with duplicates.

# QUBO
# With duplicates included
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP/output/populations_sum/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUBO.genind.dup <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUBO.genind.dup) <- factor(read.table("../../../QUBO_popmap2", header=FALSE)[,2])
# How many loci are present?
nLoc(QUBO.genind.dup) ; ncol(QUBO.genind.dup@tab)

# Without duplicates
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP/output/populations_noDuplicates/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUBO.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUBO.genind) <- factor(read.table("./QUBO_popmap_noDuplicates", header=FALSE)[,2])
# How many loci are present?
nLoc(QUBO.genind) ; ncol(QUBO.genind@tab)

# How many loci are different between the dataset with and without duplicates?
nLoc(QUBO.genind.dup) - nLoc(QUBO.genind)
# 98 more loci in the dataset with duplicates.