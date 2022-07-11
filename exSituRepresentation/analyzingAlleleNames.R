# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SCRIPT ANALYZING LOCI/ALLELE NAMES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script is exploratory, looking at the names given to loci/alleles built using Stacks

library(adegenet)

# QUBO R0 Loci (NOMAF, All SNPs; garden and wild samples Together)----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF_SNPs/"
setwd(genpop.filePath)
QUBO.R80.NOMAF.SNPs.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80.NOMAF.SNPs.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# Number of loci and alleles
nLoc(QUBO.R80.NOMAF.SNPs.genind)
ncol(QUBO.R80.NOMAF.SNPs.genind@tab)
# Names of loci
head(locNames(QUBO.R80.NOMAF.SNPs.genind))
# Names of loci with alleles
head(locNames(QUBO.R80.NOMAF.SNPs.genind, withAlleles=TRUE))
tail(locNames(QUBO.R80.NOMAF.SNPs.genind, withAlleles=TRUE))

# QUBO R80 Loci (garden and wild samples Together)----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80/"
setwd(genpop.filePath)
QUBO.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# Number of loci and alleles
nLoc(QUBO.R80.genind)
ncol(QUBO.R80.genind@tab)
# Names of loci
head(locNames(QUBO.R80.genind))
# Names of loci with alleles
head(locNames(QUBO.R80.genind, withAlleles=TRUE))
tail(locNames(QUBO.R80.genind, withAlleles=TRUE))

# Printing out polymorphic sites
# 78_9.04 and 78_9.02
QUBO.R80.genind@tab[,1:2] # Individual QUBO_W_IMLS280 is heterozygous

# 108018_4.02 and 108018_4.04
QUBO.R80.genind@tab[,763:764] # Individual QUBO_G_IMLS313 is heterozygous

# 23291_17.03 and 23291_17.01
QUBO.R80.genind@tab[,191:192] # Individual QUBO_G_IMLS316 is heterozygous

# QUAC R80 Loci (garden and wild samples Together)----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80/"
setwd(genpop.filePath)
QUAC.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# Number of loci and alleles
nLoc(QUAC.R80.genind)
ncol(QUAC.R80.genind@tab)
# Names of loci
head(locNames(QUAC.R80.genind))
# Names of loci with alleles
head(locNames(QUAC.R80.genind, withAlleles=TRUE))
tail(locNames(QUAC.R80.genind, withAlleles=TRUE))

# Printing out polymorphic sites
# 5_33.03 and 5_33.01
QUAC.R80.genind@tab[,1:2] # Individual QUAC_W_SH_Q2051 is heterozygous

# 10820_22.02 and 10820_22.04
QUAC.R80.genind@tab[,1839:1840] # Individual QUAC_G_SH_Q1485 is heterozygous

# 5424_3.04 and 5424_3.02
QUAC.R80.genind@tab[,973:974] # Individual QUAC_G_SH_Q1504 is heterozygous

# QUBO R80 Loci (garden and wild samples Separately)----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/"
setwd(genpop.filePath)
QUBO.R80.garden.genind <- read.genepop(paste0(genpop.filePath,"populations_garden_R80/populations.snps.gen"), quiet = TRUE)
QUBO.R80.wild.genind <- read.genepop(paste0(genpop.filePath,"populations_wild_R80/populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80.garden.genind) <- factor(read.table("populations_garden_R80/QUBO_popmap_garden", header=FALSE)[,2])
pop(QUBO.R80.wild.genind) <- factor(read.table("populations_wild_R80/QUBO_popmap_wild", header=FALSE)[,2])
# Number of loci and alleles
nLoc(QUBO.R80.garden.genind)
ncol(QUBO.R80.garden.genind@tab)
nLoc(QUBO.R80.wild.genind)
ncol(QUBO.R80.wild.genind@tab)
# Names of loci
head(locNames(QUBO.R80.garden.genind))
head(locNames(QUBO.R80.wild.genind))
# Names of loci with alleles
head(locNames(QUBO.R80.garden.genind, withAlleles=TRUE))
head(locNames(QUBO.R80.wild.genind, withAlleles=TRUE))

# Printing out polymorphic sites
# 1421_8.04 and 1421_8.03
QUBO.R80.garden.genind@tab[,7:8] # Individual QUBO_G_IMLS314 is heterozygous
QUBO.R80.wild.genind@tab[,9:10] # Individual QUBO_W_IMLS002 is heterozygous

# QUAC R80 Loci (garden and wild samples Separately)----
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/"
setwd(genpop.filePath)
QUAC.R80.garden.genind <- read.genepop(paste0(genpop.filePath,"populations_garden_R80/populations.snps.gen"), quiet = TRUE)
QUAC.R80.wild.genind <- read.genepop(paste0(genpop.filePath,"populations_wild_R80/populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80.garden.genind) <- factor(read.table("populations_garden_R80/QUAC_popmap_garden", header=FALSE)[,2])
pop(QUAC.R80.wild.genind) <- factor(read.table("populations_wild_R80/QUAC_popmap_wild", header=FALSE)[,2])
# Number of loci and alleles
nLoc(QUAC.R80.garden.genind)
ncol(QUAC.R80.garden.genind@tab)
nLoc(QUAC.R80.wild.genind)
ncol(QUAC.R80.wild.genind@tab)
# Names of loci
head(locNames(QUAC.R80.garden.genind), n=30L)
head(locNames(QUAC.R80.wild.genind), n=30L)

# Names of loci with alleles
head(locNames(QUAC.R80.garden.genind, withAlleles=TRUE))
head(locNames(QUAC.R80.wild.genind, withAlleles=TRUE))

# Printing out polymorphic sites
# 346_2.03 and 346_2.01
QUAC.R80.garden.genind@tab[,105:106] # Individual QUAC_G_SH_Q1500 is heterozygous
QUAC.R80.wild.genind@tab[,107:108] # Individual QUAC_W_SH_Q1885 is heterozygous
