# Duplicate analysis. Do the presence of duplicates effect the number of loci?
library(adegenet)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% REFERENCE -- FILTERED, GSNAP %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# QUAC----
# With duplicates included
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_sum/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUAC.R.genind.dup <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUAC.R.genind.dup) <- factor(read.table("../../../QUAC_popmap2", header=FALSE)[,2])
# How many loci are present?
nLoc(QUAC.R.genind.dup) ; ncol(QUAC.R.genind.dup@tab)
# How many samples?
nInd(QUAC.R.genind.dup)

# Without duplicates
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_noDuplicates/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUAC.R.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUAC.R.genind) <- factor(read.table("./QUAC_popmap_noDuplicates", header=FALSE)[,2])
# How many loci are present?
nLoc(QUAC.R.genind) ; ncol(QUAC.R.genind@tab)
# How many samples?
nInd(QUAC.R.genind)

# How many loci are different between the dataset with and without duplicates?
nLoc(QUAC.R.genind.dup) - nLoc(QUAC.R.genind)
# 98 more loci in the dataset with duplicates

# QUBO----
# With duplicates included
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP/output/populations_sum/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUBO.R.genind.dup <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUBO.R.genind.dup) <- factor(read.table("../../../QUBO_popmap2", header=FALSE)[,2])
# How many loci are present?
nLoc(QUBO.R.genind.dup) ; ncol(QUBO.R.genind.dup@tab)
# How many samples?
nInd(QUBO.R.genind.dup)

# Without duplicates
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP/output/populations_noDuplicates/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUBO.R.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUBO.R.genind) <- factor(read.table("./QUBO_popmap_noDuplicates", header=FALSE)[,2])
# How many loci are present?
nLoc(QUBO.R.genind) ; ncol(QUBO.R.genind@tab)
# How many samples?
nInd(QUBO.R.genind)

# How many loci are different between the dataset with and without duplicates?
nLoc(QUBO.R.genind.dup) - nLoc(QUBO.R.genind)
# 97 fewer loci in the dataset with duplicates

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% REFERENCE -- FILTERED, GSNAP2 %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# QUAC----
# With duplicates included
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP2/output/populations_sum/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUAC.RG2.genind.dup <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUAC.RG2.genind.dup) <- factor(read.table("../../../QUAC_popmap2", header=FALSE)[,2])
# How many loci are present?
nLoc(QUAC.RG2.genind.dup) ; ncol(QUAC.RG2.genind.dup@tab)
# How many samples?
nInd(QUAC.RG2.genind.dup)

# Without duplicates
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP2/output/populations_noDuplicates/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUAC.RG2.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUAC.RG2.genind) <- factor(read.table("./QUAC_popmap_noDuplicates", header=FALSE)[,2])
# How many loci are present?
nLoc(QUAC.RG2.genind) ; ncol(QUAC.RG2.genind@tab)
# How many samples?
nInd(QUAC.RG2.genind)

# How many loci are different between the dataset with and without duplicates?
nLoc(QUAC.RG2.genind.dup) - nLoc(QUAC.RG2.genind)
# 28 more loci in the dataset with duplicates

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% REFERENCE -- FILTERED, GSNAP4 %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# QUAC----
# With duplicates included
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP4/output/populations_sum/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUAC.RG4.genind.dup <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUAC.RG4.genind.dup) <- factor(read.table("../../../QUAC_popmap2", header=FALSE)[,2])
# How many loci are present?
nLoc(QUAC.RG4.genind.dup) ; ncol(QUAC.RG4.genind.dup@tab)
# How many samples?
nInd(QUAC.RG4.genind.dup)

# Without duplicates
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP4/output/populations_noDuplicates/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUAC.RG4.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUAC.RG4.genind) <- factor(read.table("./QUAC_popmap_noDuplicates", header=FALSE)[,2])
# How many loci are present?
nLoc(QUAC.RG4.genind) ; ncol(QUAC.RG4.genind@tab)
# How many samples?
nInd(QUAC.RG4.genind)

# How many loci are different between the dataset with and without duplicates?
nLoc(QUAC.RG4.genind.dup) - nLoc(QUAC.RG4.genind)
# 98 more loci in the dataset with duplicates

# QUBO----
# With duplicates included
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_sum/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUBO.RG4.genind.dup <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUBO.RG4.genind.dup) <- factor(read.table("../../../QUBO_popmap2", header=FALSE)[,2])
# How many loci are present?
nLoc(QUBO.RG4.genind.dup) ; ncol(QUBO.RG4.genind.dup@tab)
# How many samples?
nInd(QUBO.RG4.genind.dup)

# Without duplicates
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_noDuplicates/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUBO.RG4.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUBO.RG4.genind) <- factor(read.table("./QUBO_popmap_noDuplicates", header=FALSE)[,2])
# How many loci are present?
nLoc(QUBO.RG4.genind) ; ncol(QUBO.RG4.genind@tab)
# How many samples?
nInd(QUBO.RG4.genind)

# How many loci are different between the dataset with and without duplicates?
nLoc(QUBO.RG4.genind.dup) - nLoc(QUBO.RG4.genind)
# 97 fewer loci in the dataset with duplicates

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% DE NOVO -- FINAL ASSEMBLIES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# QUAC----
# With duplicates included
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_sum/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUAC.DN.genind.dup <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUAC.DN.genind.dup) <- factor(read.table("../../QUAC_popmap2", header=FALSE)[,2])
# How many loci are present?
nLoc(QUAC.DN.genind.dup) ; ncol(QUAC.DN.genind.dup@tab)
# How many samples?
nInd(QUAC.DN.genind.dup)

# Without duplicates
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_noDuplicates/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUAC.DN.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUAC.DN.genind) <- factor(read.table("./QUAC_popmap2_noDuplicates", header=FALSE)[,2])
# How many loci are present?
nLoc(QUAC.DN.genind) ; ncol(QUAC.DN.genind@tab)
# How many samples?
nInd(QUAC.DN.genind)

# How many loci are different between the dataset with and without duplicates?
nLoc(QUAC.DN.genind.dup) - nLoc(QUAC.DN.genind)
# 20 more loci in the dataset with duplicates

# QUBO----
# With duplicates included
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_sum/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUBO.DN.genind.dup <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUBO.DN.genind.dup) <- factor(read.table("../../QUBO_popmap2", header=FALSE)[,2])
# How many loci are present?
nLoc(QUBO.DN.genind.dup) ; ncol(QUBO.DN.genind.dup@tab)
# How many samples?
nInd(QUBO.DN.genind.dup)

# Without duplicates
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_noDuplicates/"
setwd(genpop.filePath)
# Reading in and processing .gen file
QUBO.DN.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
pop(QUBO.DN.genind) <- factor(read.table("./QUBO_popmap2_noDuplicates", header=FALSE)[,2])
# How many loci are present?
nLoc(QUBO.DN.genind) ; ncol(QUBO.DN.genind@tab)
# How many samples?
nInd(QUBO.DN.genind)

# How many loci are different between the dataset with and without duplicates?
nLoc(QUBO.DN.genind.dup) - nLoc(QUBO.DN.genind)
# 27 fewer loci in the dataset with duplicates