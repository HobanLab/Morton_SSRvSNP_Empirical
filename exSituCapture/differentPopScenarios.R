# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% ANALYZING CAPTURE RATES IN DIFFERENT SCENARIOS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(adegenet)

# ---- FUNCTIONS ----
# Function for reporting allelic capture rates, using a single genind object
reportAllelicCapture_Together <- function(gen.obj){
  # Generate numberical vectors corresponding to garden and wild rows, for later calculations
  garden.Rows <- seq_len(length(which(pop(gen.obj)=="garden")))
  wild.Rows <- seq(from=length(which(pop(gen.obj)=="garden"))+1, to=nInd(gen.obj))
  garden.N <- length(garden.Rows)
  wild.N <- length(wild.Rows)
  # Build the wild allele frequency vector
  wildFreqs <- colSums(gen.obj@tab[wild.Rows,], na.rm = TRUE)/(wild.N*2)*100
  # Calculate capture rates
  # Total
  total <- length(which(names(which(wildFreqs > 0)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 0))*100
  # Very common
  veryCommon <- length(which(names(which(wildFreqs > 10)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 10))*100
  # Common
  common <- length(which(names(which(wildFreqs > 5)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 5))*100
  # Low frequency
  lowFrequency <- length(which(names(which(wildFreqs < 10 & wildFreqs > 1)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs < 10 & wildFreqs > 1))*100
  # Rare
  rare <- length(which(names(which(wildFreqs < 1 & wildFreqs > 0)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs < 1 & wildFreqs > 0))*100
  # Build list of rates
  captureRates <- c(total,veryCommon,common,lowFrequency,rare)
  names(captureRates) <- c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%)","Rare (<1%)")
  # Print capture rates and return
  return(captureRates)
}

# Function for reporting allelic capture rates, using two genind objects (one for garden, one for wild)
reportAllelicCapture_Separate <- function(gen.obj.garden, gen.obj.wild){
  # Build the wild allele frequency vector
  wildFreqs <- colSums(gen.obj.wild@tab, na.rm = TRUE)/(nInd(gen.obj.wild)*2)*100
  # Calculate capture rates
  # Total
  total <- length(which(names(which(wildFreqs > 0)) %in% colnames(gen.obj.garden@tab)))/length(which(wildFreqs > 0))*100
  # Very common
  veryCommon <- length(which(names(which(wildFreqs > 10)) %in% colnames(gen.obj.garden@tab)))/length(which(wildFreqs > 10))*100
  # Common
  common <- length(which(names(which(wildFreqs > 5)) %in% colnames(gen.obj.garden@tab)))/length(which(wildFreqs > 5))*100
  # Low frequency
  lowFrequency <- length(which(names(which(wildFreqs < 10 & wildFreqs > 1)) %in% colnames(gen.obj.garden@tab)))/length(which(wildFreqs < 10 & wildFreqs > 1))*100
  # Rare
  rare <- length(which(names(which(wildFreqs < 1 & wildFreqs > 0)) %in% colnames(gen.obj.garden@tab)))/length(which(wildFreqs < 1 & wildFreqs > 0))*100
  # Build list of rates
  captureRates <- c(total,veryCommon,common,lowFrequency,rare)
  names(captureRates) <- c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%)","Rare (<1%)")
  # Print capture rates and return
  return(captureRates)
}

# %%%% QUAC ----
# ---- SINGLE GENIND ----
# R0 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0/"
setwd(genpop.filePath)
QUAC.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R0 Capture rates
reportAllelicCapture_Together(QUAC.R0.genind)
nLoc(QUAC.R0.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80/"
setwd(genpop.filePath)
QUAC.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R80 Capture rates
reportAllelicCapture_Together(QUAC.R80.genind)
nLoc(QUAC.R80.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R100/"
setwd(genpop.filePath)
QUAC.R100.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R100 Capture rates
reportAllelicCapture_Together(QUAC.R100.genind)
nLoc(QUAC.R100.genind)

# ----TWO GENINDS ----
# R0 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R0/"
setwd(genpop.filePath)
QUAC.R0.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0/"
setwd(genpop.filePath)
QUAC.R0.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R0 Capture rates
reportAllelicCapture_Separate(QUAC.R0.garden.genind, QUAC.R0.wild.genind)
nLoc(QUAC.R0.garden.genind)
nLoc(QUAC.R0.wild.genind)

# R80 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R80/"
setwd(genpop.filePath)
QUAC.R80.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R80/"
setwd(genpop.filePath)
QUAC.R80.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R80 Capture rates
reportAllelicCapture_Separate(QUAC.R80.garden.genind, QUAC.R80.wild.genind)
nLoc(QUAC.R80.garden.genind)
nLoc(QUAC.R80.wild.genind)

# R100 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R100/"
setwd(genpop.filePath)
QUAC.R100.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R100/"
setwd(genpop.filePath)
QUAC.R100.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R100 Capture rates
reportAllelicCapture_Separate(QUAC.R100.garden.genind, QUAC.R100.wild.genind)
nLoc(QUAC.R100.garden.genind)
nLoc(QUAC.R100.wild.genind)

# ---- NO MINOR ALLELE FREQUENCY (NOMAF) ----
# ---- SINGLE GENIND ----
# R0 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R0_NOMAF/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R0_NOMAF Capture rates
reportAllelicCapture_Together(QUAC.R0_NOMAF.genind)
nLoc(QUAC.R0_NOMAF.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF/"
setwd(genpop.filePath)
QUAC.R80_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80_NOMAF.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R80_NOMAF Capture rates
reportAllelicCapture_Together(QUAC.R80_NOMAF.genind)
nLoc(QUAC.R80_NOMAF.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R100_NOMAF/"
setwd(genpop.filePath)
QUAC.R100_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100_NOMAF.genind) <- factor(read.table("QUAC_popmap2", header=FALSE)[,2])
# R100_NOMAF Capture rates
reportAllelicCapture_Together(QUAC.R100_NOMAF.genind)
nLoc(QUAC.R100_NOMAF.genind)

# ----TWO GENINDS ----
# R0 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R0_NOMAF/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R0_NOMAF/"
setwd(genpop.filePath)
QUAC.R0_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R0_NOMAF.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R0 Capture rates
reportAllelicCapture_Separate(QUAC.R0_NOMAF.garden.genind, QUAC.R0_NOMAF.wild.genind)
nLoc(QUAC.R0_NOMAF.garden.genind)
nLoc(QUAC.R0_NOMAF.wild.genind)

# R80 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R80_NOMAF/"
setwd(genpop.filePath)
QUAC.R80_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80_NOMAF.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R80_NOMAF/"
setwd(genpop.filePath)
QUAC.R80_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80_NOMAF.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R80 Capture rates
reportAllelicCapture_Separate(QUAC.R80_NOMAF.garden.genind, QUAC.R80_NOMAF.wild.genind)
nLoc(QUAC.R80_NOMAF.garden.genind)
nLoc(QUAC.R80_NOMAF.wild.genind)

# R100 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R100_NOMAF/"
setwd(genpop.filePath)
QUAC.R100_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100_NOMAF.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R100_NOMAF/"
setwd(genpop.filePath)
QUAC.R100_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100_NOMAF.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R100 Capture rates
reportAllelicCapture_Separate(QUAC.R100_NOMAF.garden.genind, QUAC.R100_NOMAF.wild.genind)
nLoc(QUAC.R100_NOMAF.garden.genind)
nLoc(QUAC.R100_NOMAF.wild.genind)

# TROUBLESHOOTING ----
# TWO GENINDS 
# R80 
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R80/"
setwd(genpop.filePath)
QUAC.R80.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R80/"
setwd(genpop.filePath)
QUAC.R80.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R80.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R80 Capture rates
reportAllelicCapture_Separate(QUAC.R80.garden.genind, QUAC.R80.wild.genind)

# Wild allele frequency vector
QUAC.R80_wildFreqs <- colSums(QUAC.R80.wild.genind@tab, na.rm = TRUE)/(nInd(QUAC.R80.wild.genind)*2)*100
QUAC.R80_wildFreqs[1:10]

min(QUAC.R80_wildFreqs)
which.min(QUAC.R80_wildFreqs)
max(QUAC.R80_wildFreqs)
which.max(QUAC.R80_wildFreqs)
hist(QUAC.R80_wildFreqs)

# Lowest frequency allele across wild samples
QUAC.R80.wild.genind@tab[,which.min(QUAC.R80_wildFreqs)]
# Which samples contain lowest frequency allele
which(QUAC.R80.wild.genind@tab[,which.min(QUAC.R80_wildFreqs)] != 0)
# Prevalence of lowest frequency allele in samples that contain (i.e. homozygote or heterozygote?)
QUAC.R80.wild.genind@tab[which(QUAC.R80.wild.genind@tab[,which.min(QUAC.R80_wildFreqs)] != 0), which.min(QUAC.R80_wildFreqs)]

# Total
length(which(names(which(QUAC.R80_wildFreqs > 0)) %in% colnames(QUAC.R80.garden.genind@tab)))/length(which(QUAC.R80_wildFreqs > 0))*100
# Very common
length(which(names(which(QUAC.R80_wildFreqs > 10)) %in% colnames(QUAC.R80.garden.genind@tab)))/length(which(QUAC.R80_wildFreqs > 10))*100
# Common
length(which(names(which(QUAC.R80_wildFreqs > 5)) %in% colnames(QUAC.R80.garden.genind@tab)))/length(which(QUAC.R80_wildFreqs > 5))*100
# Low frequency
length(which(names(which(QUAC.R80_wildFreqs < 10 & QUAC.R80_wildFreqs > 1)) %in% colnames(QUAC.R80.garden.genind@tab)))/length(which(QUAC.R80_wildFreqs < 10 & QUAC.R80_wildFreqs > 1))*100
# Rare
length(which(names(which(QUAC.R80_wildFreqs < 1)) %in% colnames(QUAC.R80.garden.genind@tab)))/length(which(QUAC.R80_wildFreqs < 1))*100

# R100 
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_garden_R100/"
setwd(genpop.filePath)
QUAC.R100.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100.garden.genind) <- factor(read.table("QUAC_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R100/"
setwd(genpop.filePath)
QUAC.R100.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUAC.R100.wild.genind) <- factor(read.table("QUAC_popmap_wild", header=FALSE)[,2])

# R100 Capture rates
reportAllelicCapture_Separate(QUAC.R100.garden.genind, QUAC.R100.wild.genind)

# Wild allele frequency vector
QUAC.R100_wildFreqs <- colSums(QUAC.R100.wild.genind@tab, na.rm = TRUE)/(nInd(QUAC.R100.wild.genind)*2)*100
QUAC.R100_wildFreqs[1:10]

min(QUAC.R100_wildFreqs)
which.min(QUAC.R100_wildFreqs)
max(QUAC.R100_wildFreqs)
which.max(QUAC.R100_wildFreqs)
hist(QUAC.R100_wildFreqs)

# Lowest frequency allele across wild samples
QUAC.R100.wild.genind@tab[,which.min(QUAC.R100_wildFreqs)]
# Which samples contain lowest frequency allele
which(QUAC.R100.wild.genind@tab[,which.min(QUAC.R100_wildFreqs)] != 0)
# Prevalence of lowest frequency allele in samples that contain (i.e. homozygote or heterozygote?)
QUAC.R100.wild.genind@tab[which(QUAC.R100.wild.genind@tab[,which.min(QUAC.R100_wildFreqs)] != 0), which.min(QUAC.R100_wildFreqs)]

# Total
length(which(names(which(QUAC.R100_wildFreqs > 0)) %in% colnames(QUAC.R100.garden.genind@tab)))/length(which(QUAC.R100_wildFreqs > 0))*100
# Very common
length(which(names(which(QUAC.R100_wildFreqs > 10)) %in% colnames(QUAC.R100.garden.genind@tab)))/length(which(QUAC.R100_wildFreqs > 10))*100
# Common
length(which(names(which(QUAC.R100_wildFreqs > 5)) %in% colnames(QUAC.R100.garden.genind@tab)))/length(which(QUAC.R100_wildFreqs > 5))*100
# Low frequency
length(which(names(which(QUAC.R100_wildFreqs < 10 & QUAC.R100_wildFreqs > 1)) %in% colnames(QUAC.R100.garden.genind@tab)))/length(which(QUAC.R100_wildFreqs < 10 & QUAC.R100_wildFreqs > 1))*100
# Rare
length(which(names(which(QUAC.R100_wildFreqs < 1)) %in% colnames(QUAC.R100.garden.genind@tab)))/length(which(QUAC.R100_wildFreqs < 1))*100

# Different R values for different populations
reportAllelicCapture_Separate(QUAC.R0.garden.genind, QUAC.R100.wild.genind)
reportAllelicCapture_Separate(QUAC.R0.garden.genind, QUAC.R80.wild.genind)

# %%%% QUBO ----
# ---- SINGLE GENIND ----
# R0 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0/"
setwd(genpop.filePath)
QUBO.R0.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R0 Capture rates
reportAllelicCapture_Together(QUBO.R0.genind)
nLoc(QUBO.R0.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80/"
setwd(genpop.filePath)
QUBO.R80.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R80 Capture rates
reportAllelicCapture_Together(QUBO.R80.genind)
nLoc(QUBO.R80.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R100/"
setwd(genpop.filePath)
QUBO.R100.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R100 Capture rates
reportAllelicCapture_Together(QUBO.R100.genind)
nLoc(QUBO.R100.genind)

# ----TWO GENINDS ----
# R0 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R0/"
setwd(genpop.filePath)
QUBO.R0.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0/"
setwd(genpop.filePath)
QUBO.R0.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R0 Capture rates
reportAllelicCapture_Separate(QUBO.R0.garden.genind, QUBO.R0.wild.genind)
nLoc(QUBO.R0.garden.genind)
nLoc(QUBO.R0.wild.genind)

# R80 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R80/"
setwd(genpop.filePath)
QUBO.R80.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80/"
setwd(genpop.filePath)
QUBO.R80.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R80 Capture rates
reportAllelicCapture_Separate(QUBO.R80.garden.genind, QUBO.R80.wild.genind)
nLoc(QUBO.R80.garden.genind)
nLoc(QUBO.R80.wild.genind)

# R100 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R100/"
setwd(genpop.filePath)
QUBO.R100.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R100/"
setwd(genpop.filePath)
QUBO.R100.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R100 Capture rates
reportAllelicCapture_Separate(QUBO.R100.garden.genind, QUBO.R100.wild.genind)
nLoc(QUBO.R100.garden.genind)
nLoc(QUBO.R100.wild.genind)

# ---- NO MINOR ALLELE FREQUENCY (NOMAF) ----
# ---- SINGLE GENIND ----
# R0 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R0_NOMAF/"
setwd(genpop.filePath)
QUBO.R0_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R0_NOMAF Capture rates
reportAllelicCapture_Together(QUBO.R0_NOMAF.genind)
nLoc(QUBO.R0_NOMAF.genind)

# R80 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF/"
setwd(genpop.filePath)
QUBO.R80_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80_NOMAF.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R80_NOMAF Capture rates
reportAllelicCapture_Together(QUBO.R80_NOMAF.genind)
nLoc(QUBO.R80_NOMAF.genind)

# R100 ----
# Read in genind file
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R100_NOMAF/"
setwd(genpop.filePath)
QUBO.R100_NOMAF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100_NOMAF.genind) <- factor(read.table("QUBO_popmap2", header=FALSE)[,2])
# R100_NOMAF Capture rates
reportAllelicCapture_Together(QUBO.R100_NOMAF.genind)
nLoc(QUBO.R100_NOMAF.genind)

# ----TWO GENINDS ----
# R0 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R0_NOMAF/"
setwd(genpop.filePath)
QUBO.R0_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R0_NOMAF/"
setwd(genpop.filePath)
QUBO.R0_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R0_NOMAF.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R0_NOMAF Capture rates
reportAllelicCapture_Separate(QUBO.R0_NOMAF.garden.genind, QUBO.R0_NOMAF.wild.genind)
nLoc(QUBO.R0_NOMAF.garden.genind)
nLoc(QUBO.R0_NOMAF.wild.genind)

# R80 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R80_NOMAF/"
setwd(genpop.filePath)
QUBO.R80_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80_NOMAF.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80_NOMAF/"
setwd(genpop.filePath)
QUBO.R80_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R80_NOMAF.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R80_NOMAF Capture rates
reportAllelicCapture_Separate(QUBO.R80_NOMAF.garden.genind, QUBO.R80_NOMAF.wild.genind)
nLoc(QUBO.R80_NOMAF.garden.genind)
nLoc(QUBO.R80_NOMAF.wild.genind)

# R100 ----
# GARDEN
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_garden_R100_NOMAF/"
setwd(genpop.filePath)
QUBO.R100_NOMAF.garden.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100_NOMAF.garden.genind) <- factor(read.table("QUBO_popmap_garden", header=FALSE)[,2])

# WILD
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R100_NOMAF/"
setwd(genpop.filePath)
QUBO.R100_NOMAF.wild.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"), quiet = TRUE)
# Correct popNames
pop(QUBO.R100_NOMAF.wild.genind) <- factor(read.table("QUBO_popmap_wild", header=FALSE)[,2])

# R100_NOMAF Capture rates
reportAllelicCapture_Separate(QUBO.R100_NOMAF.garden.genind, QUBO.R100_NOMAF.wild.genind)
nLoc(QUBO.R100_NOMAF.garden.genind)
nLoc(QUBO.R100_NOMAF.wild.genind)
