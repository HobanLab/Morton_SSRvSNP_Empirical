# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS (DAPC) %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses the adegenet package to run DAPC on two datasets
# For QUAC, the de novo final assembly (DNFA) dataset is used
# For QUBO, the dataset aligned to the Quercus robur reference genome (GSNAP4) is used
# In both analyses, garden and wild samples are included 

library(adegenet)
library(scales)
# Set the working directory
SSRvSNP.wd <- "~/Documents/SSRvSNP/Code/"
setwd(SSRvSNP.wd)

# TO DO: Emily wants you to run a normal PCA on the dataset (i.e. not through adegenet; instructions
# in adegenet Basics tutorial)

# %%%% QUAC %%%% ----
# %%%% COMPLETE: GARDEN & WILD ----
# To be populated
# %%%% COMPLETE: WILD ONLY ----
# ---- READ IN AND PROCESS GENIND FILES ----
# MSAT
QUAC.MSAT.genpop.filePath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/"
QUAC.MSAT.genind <- read.genepop(paste0(QUAC.MSAT.genpop.filePath, "QUAC_woK_allpop_clean.gen"), ncode = 3)
# Specify filepath to GCC_QUAC_ZAIN dataframe, containing sample names and population names
QUAC.MSAT.dataframe_filepath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Data_Frames/QUAC_woK_allpop_clean_df.csv"
# Assign sample names: read in Tissue database names from GCC_QUAC_ZAIN repository
QUAC.MSAT.tissueNames <- unlist(read.csv2(QUAC.MSAT.dataframe_filepath, header = TRUE, sep=",")[1])
rownames(QUAC.MSAT.genind@tab) <- QUAC.MSAT.tissueNames
# Correct population names: read in a dataframe containing population values
pop(QUAC.MSAT.genind) <- unlist(read.csv2(QUAC.MSAT.dataframe_filepath, header = TRUE, sep=",")[2])
# Remove garden samples from tissueNames vector
QUAC.MSAT.tissueNames <- QUAC.MSAT.tissueNames[-grep(pattern = "QAc-G-", QUAC.MSAT.tissueNames)]
# Subset Complete MSAT genind object to just wild populations (last 5 populations)
QUAC.W.MSAT.genind <- QUAC.MSAT.genind[QUAC.MSAT.tissueNames,,drop=TRUE]

# ---- RUN AND PLOT DAPC ----
# MSAT ----
# Use k-means clustering to find a number of groups which maximizes the variation between groups
# (after transforming the data using PCA, in order to increase computation times)
QUAC.W.MSAT.grp <- find.clusters(QUAC.W.MSAT.genind)
# QUAC.W.MSAT.grp <- find.clusters(QUAC.W.MSAT.genind, max.n.clust=20, n.pca = 120, n.clust = 4)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 120
# Clusters (seeking to minimize the Bayesian Information Criterion value): 4
# Conduct DAPC: transform the data using PCA, then run discriminant analysis on retained principal components,
# using the inferred groupings from the k-means clustering step above
QUAC.W.MSAT.dapc1 <- dapc(x=QUAC.W.MSAT.genind, pop=QUAC.W.MSAT.grp$grp)
# QUAC.W.MSAT.dapc1 <- dapc(x=QUAC.W.MSAT.genind, pop=QUAC.W.MSAT.grp$grp, n.pca = 50, n.da = 3)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 50
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 3
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUAC.W.MSAT.dapc1$grp, as.character(pop(QUAC.W.MSAT.genind)), rownames(QUAC.W.MSAT.genind@tab))

# Show DAPC as a scatterplot
scatter(QUAC.W.MSAT.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"),
        txt.leg = c("Pryor", "Porter", "Magazine", "Sugarloaf"))
mtext("QUAC MSAT: Complete Wild", adj=0.07)

# %%%% SUBSET: GARDEN & WILD ----
# To be populated

# %%%% SUBSET: WILD ONLY ----
# ---- READ IN AND PROCESS GENIND FILES ----
# MSAT
QUAC.MSAT.genpop.filePath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Adegenet_Files/"
QUAC.MSAT.genind <- read.genepop(paste0(QUAC.MSAT.genpop.filePath, "QUAC_woK_allpop_clean.gen"), ncode = 3)
# Specify filepath to GCC_QUAC_ZAIN dataframe, containing sample names and population names
QUAC.MSAT.dataframe_filepath <- 
  "~/Documents/peripheralProjects/GCC_QUAC_ZAIN/Data_Files/Data_Frames/QUAC_woK_allpop_clean_df.csv"
# Assign sample names: read in Tissue database names from GCC_QUAC_ZAIN repository
QUAC.MSAT.tissueNames <- unlist(read.csv2(QUAC.MSAT.dataframe_filepath, header = TRUE, sep=",")[1])
rownames(QUAC.MSAT.genind@tab) <- QUAC.MSAT.tissueNames
# Correct population names: read in a dataframe containing population values
pop(QUAC.MSAT.genind) <- unlist(read.csv2(QUAC.MSAT.dataframe_filepath, header = TRUE, sep=",")[2])

# SNP
# De novo
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_wild_R80_NOMAF_1SNP_NoK/"
QUAC.Subset.W.SNP.DN.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.Subset.W.SNP.DN.genind) <- factor(read.table(paste0(genpop.filePath, "QUAC_popmap_wild_NoK"), header=FALSE)[,2])
# Reference
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/Q_rubra/output/populations_wild_R80_NOMAF_1SNP_NoK/"
QUAC.Subset.W.SNP.REF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUAC.Subset.W.SNP.REF.genind) <- factor(read.table(paste0(genpop.filePath, "QUAC_popmap_wild_NoK"), header=FALSE)[,2])
# Read in Tissue database names, for subsetting genetic matrices later
# This file was created (by Austin K.), and can be found on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
QUAC.SNP.tissueNames_filepath <- paste0(SSRvSNP.wd,"exSituRepresentation/QUAC_SNP_TissueNames.csv")
QUAC.SNP.tissueNames <- unlist(read.csv2(QUAC.SNP.tissueNames_filepath, header = TRUE, sep = ",")[3])
# Remove garden samples from SNP tissue database names vector
QUAC.SNP.tissueNames <- QUAC.SNP.tissueNames[-grep(pattern = "QAc-G-", QUAC.SNP.tissueNames)]
# Rename SNP De novo and Reference genind matrices
rownames(QUAC.Subset.W.SNP.DN.genind@tab) <- QUAC.SNP.tissueNames
rownames(QUAC.Subset.W.SNP.REF.genind@tab) <- QUAC.SNP.tissueNames

# SUBSET MSAT GENIND
QUAC_sharedSamples <- sort(QUAC.SNP.tissueNames[which(QUAC.SNP.tissueNames %in% QUAC.MSAT.tissueNames)])
# Subset MSAT matrix to strictly shared samples, dropping now absent alleles
QUAC.Subset.W.MSAT.genind <- QUAC.MSAT.genind[QUAC_sharedSamples,, drop=TRUE]

# ---- RUN AND PLOT DAPC ----
# MSAT ----
# Use k-means clustering to find a number of groups which maximizes the variation between groups
# (after transforming the data using PCA, in order to increase computation times)
# QUAC.Subset.W.MSAT.grp <- find.clusters(QUAC.Subset.W.MSAT.genind, max.n.clust=20, n.pca = 90, n.clust = 4)
QUAC.Subset.W.MSAT.grp <- find.clusters(QUAC.Subset.W.MSAT.genind)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 90
# Clusters (seeking to minimize the Bayesian Information Criterion value): 4
# Conduct DAPC: transform the data using PCA, then run discriminant analysis on retained principal components,
# using the inferred groupings from the k-means clustering step above
# QUAC.Subset.W.MSAT.dapc1 <- dapc(x=QUAC.Subset.W.MSAT.genind, pop=QUAC.Subset.W.MSAT.grp$grp, n.pca = 40, n.da = 3)
QUAC.Subset.W.MSAT.dapc1 <- dapc(x=QUAC.Subset.W.MSAT.genind, pop=QUAC.Subset.W.MSAT.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 40
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 3
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUAC.Subset.W.MSAT.dapc1$grp, as.character(pop(QUAC.Subset.W.MSAT.genind)), rownames(QUAC.Subset.W.MSAT.genind@tab))

# Show DAPC as a scatterplot
scatter(QUAC.Subset.W.MSAT.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"),
        txt.leg = c("Magazine", "Sugarloaf", "Pryor", "Porter"))
mtext("QUAC MSAT: Subset Wild", adj=0.07)

# SNP, DE NOVO (R80) ----
# SUPPORTED CLUSTERS ----
# K-means clustering step
# QUAC.Subset.W.SNP.DN.grp <- find.clusters(QUAC.Subset.W.SNP.DN.genind, max.n.clust=20, n.pca = 150, n.clust = 2)
QUAC.Subset.W.SNP.DN.grp <- find.clusters(QUAC.Subset.W.SNP.DN.genind)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# Clusters (seeking to minimize the Bayesian Information Criterion value): 2
# PCA step
# QUAC.Subset.W.SNP.DN.dapc1 <- dapc(x=QUAC.Subset.W.SNP.DN.genind, pop=QUAC.Subset.W.SNP.DN.grp$grp, n.pca = 70, n.da = 1)
QUAC.Subset.W.SNP.DN.dapc1 <- dapc(x=QUAC.Subset.W.SNP.DN.genind, pop=QUAC.Subset.W.SNP.DN.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 70
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 1
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUAC.Subset.W.SNP.DN.dapc1$grp, as.character(pop(QUAC.Subset.W.SNP.DN.genind)), rownames(QUAC.Subset.W.SNP.DN.genind@tab))

# Show DAPC as a scatterplot
scatter(QUAC.Subset.W.SNP.DN.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"),
        txt.leg = c("Pryor", "Sugarloaf/Magazine/Porter"))
mtext("QUAC SNP De novo (R80): Subset Wild (Supported Clusters)", adj=0.15)

# GEOGRAPHIC CLUSTERS ----
# K-means clustering step
# QUAC.Subset.W.SNP.DN.grp <- find.clusters(QUAC.Subset.W.SNP.DN.genind, max.n.clust=20, n.pca = 150, n.clust = 4)
QUAC.Subset.W.SNP.DN.grp <- find.clusters(QUAC.Subset.W.SNP.DN.genind)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# Clusters (seeking to minimize the Bayesian Information Criterion value): 4
# PCA step
# QUAC.Subset.W.SNP.DN.dapc1 <- dapc(x=QUAC.Subset.W.SNP.DN.genind, pop=QUAC.Subset.W.SNP.DN.grp$grp, n.pca = 70, n.da = 3)
QUAC.Subset.W.SNP.DN.dapc1 <- dapc(x=QUAC.Subset.W.SNP.DN.genind, pop=QUAC.Subset.W.SNP.DN.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 70
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 3
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUAC.Subset.W.SNP.DN.dapc1$grp, as.character(pop(QUAC.Subset.W.SNP.DN.genind)), rownames(QUAC.Subset.W.SNP.DN.genind@tab))

# Show DAPC as a scatterplot
scatter(QUAC.Subset.W.SNP.DN.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"),
        txt.leg = c("Sugarloaf", "Porter", "Magazine", "Pryor"))
mtext("QUAC SNP De novo (R80): Subset Wild (Geographic Clusters)", adj=0.07, line=2.5)

# SNP, REFERENCE (R80) ----
# SUPPORTED CLUSTERS ----
# K-means clustering step
# QUAC.Subset.W.SNP.REF.grp <- find.clusters(QUAC.Subset.W.SNP.REF.genind, max.n.clust=20, n.pca = 150, n.clust = 2)
QUAC.Subset.W.SNP.REF.grp <- find.clusters(QUAC.Subset.W.SNP.REF.genind)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# Clusters (seeking to minimize the Bayesian Information Criterion value): 2
# PCA step
# QUAC.Subset.W.SNP.REF.dapc1 <- dapc(x=QUAC.Subset.W.SNP.REF.genind, pop=QUAC.Subset.W.SNP.REF.grp$grp, n.pca = 70, n.da = 1)
QUAC.Subset.W.SNP.REF.dapc1 <- dapc(x=QUAC.Subset.W.SNP.REF.genind, pop=QUAC.Subset.W.SNP.REF.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 70
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 1
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUAC.Subset.W.SNP.REF.dapc1$grp, as.character(pop(QUAC.Subset.W.SNP.REF.genind)), rownames(QUAC.Subset.W.SNP.REF.genind@tab))

# Show DAPC as a scatterplot
scatter(QUAC.Subset.W.SNP.REF.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"),
        txt.leg = c("Pryor", "Sugarloaf/Magazine/Porter"))
mtext("QUAC SNP Reference (R80): Subset Wild (Supported Clusters)", adj=0.15)

# GEOGRAPHIC CLUSTERS ----
# K-means clustering step
# QUAC.Subset.W.SNP.REF.grp <- find.clusters(QUAC.Subset.W.SNP.REF.genind, max.n.clust=20, n.pca = 150, n.clust = 4)
QUAC.Subset.W.SNP.REF.grp <- find.clusters(QUAC.Subset.W.SNP.REF.genind)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# Clusters (seeking to minimize the Bayesian Information Criterion value): 4
# PCA step
# QUAC.Subset.W.SNP.REF.dapc1 <- dapc(x=QUAC.Subset.W.SNP.REF.genind, pop=QUAC.Subset.W.SNP.REF.grp$grp, n.pca = 70, n.da = 3)
QUAC.Subset.W.SNP.REF.dapc1 <- dapc(x=QUAC.Subset.W.SNP.REF.genind, pop=QUAC.Subset.W.SNP.REF.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 70
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 3
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUAC.Subset.W.SNP.REF.dapc1$grp, as.character(pop(QUAC.Subset.W.SNP.REF.genind)), rownames(QUAC.Subset.W.SNP.REF.genind@tab))

# Show DAPC as a scatterplot
scatter(QUAC.Subset.W.SNP.REF.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"),
        txt.leg = c("Magazine 2", "Pryor", "Sugarloaf", "Porter/Magazine 1"))
mtext("QUAC SNP Reference (R80): Subset Wild (Geographic Clusters)", adj=0.15, line=2.5)

# Show DAPC as a STRUCTURE-like plot
compoplot(QUAC.Subset.W.SNP.REF.dapc1, posi="bottomright", 
          txt.leg = c("Porter/Magazine 1", "Magazine 2", "Sugarloaf", "Pryor"), lab="",
          ncol=1, xlab="Individuals")

# %%%% QUBO %%%% ----
# %%%% COMPLETE: GARDEN & WILD ----
# %%%% COMPLETE: WILD ONLY ----
# ---- READ IN AND PROCESS GENIND FILES ----
# MSAT
genpop.filePath <- 
  "~/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
QUBO.W.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_wild_w_ALL.gen"), ncode=3)
# Correct popNames: read in a dataframe which contains population names, in order of 
# samples in genind file. This document was created manually by Austin Koontz, and can be 
# found on the Hoban Lab Drive
QUBO.MSAT.dataframe_filepath <- 
  "~/Documents/SSRvSNP/Code/popAnalyses/QUBO_MSAT_Complete_WildPops.csv"
# Correct population names: read in a dataframe containing population values
pop(QUBO.W.MSAT.genind) <- unlist(read.csv2(QUBO.MSAT.dataframe_filepath, header = TRUE, sep=",", skip = 1)[2])
# Rename MSAT samples 
# Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.W.MSAT.genind@tab), function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.W.MSAT.genind@tab) <- QUBO.MSAT.sampleNames

# ---- RUN AND PLOT DAPC ----
# MSAT ----
# Use k-means clustering to find a number of groups which maximizes the variation between groups
# (after transforming the data using PCA, in order to increase computation times)
# QUBO.W.MSAT.grp <- find.clusters(QUBO.W.MSAT.genind, max.n.clust=20, n.pca = 150, n.clust = 6)
QUBO.W.MSAT.grp <- find.clusters(QUBO.W.MSAT.genind)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# Clusters (seeking to minimize the Bayesian Information Criterion value): 6
# Conduct DAPC: transform the data using PCA, then run discriminant analysis on retained principal components,
# using the inferred groupings from the k-means clustering step above
# QUBO.W.MSAT.dapc1 <- dapc(x=QUBO.W.MSAT.genind, pop=QUBO.W.MSAT.grp$grp, n.pca = 40, n.da = 5)
QUBO.W.MSAT.dapc1 <- dapc(x=QUBO.W.MSAT.genind, pop=QUBO.W.MSAT.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 40
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 5
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUBO.W.MSAT.dapc1$grp, as.character(pop(QUBO.W.MSAT.genind)), rownames(QUBO.W.MSAT.genind@tab))

# Show DAPC as a scatterplot
scatter(QUBO.W.MSAT.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=1.3, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5", "#D95F02","#8C510A"),
        txt.leg = c("Worldsong/MossRock/EBSCO/BTKC", "Peavine/EBSCO/Hinds", "Worldsong/Irondale/EBSCO/Hinds", 
                    "EBSCO/Hinds/Pop11", "EBSCO/Hinds/Peavine", "EBSCO/BTKC/Hinds/Peavine"))
mtext("QUBO MSAT: Complete Wild", adj=0.15)

# %%%% SUBSET: WILD ONLY ----
# ---- READ IN AND PROCESS GENIND FILES ----
# MSAT
genpop.filePath <- 
  "~/Documents/peripheralProjects/SE_oaks_genetics/genetic_data/"
QUBO.W.MSAT.genind <- read.genepop(paste0(genpop.filePath,"Qb_wild_w_ALL.gen"), ncode=3)
# Correct popNames: read in a dataframe which contains population names, in order of 
# samples in genind file. This document was created manually by Austin Koontz, and can be 
# found on the Hoban Lab Drive
QUBO.MSAT.dataframe_filepath <- 
  "~/Documents/SSRvSNP/Code/popAnalyses/QUBO_MSAT_Complete_WildPops.csv"
# Correct population names: read in a dataframe containing population values
pop(QUBO.W.MSAT.genind) <- unlist(read.csv2(QUBO.MSAT.dataframe_filepath, header = TRUE, sep=",", skip = 1)[2])
# Rename MSAT samples 
# Split sample names on underscore, and return 3rd element. Rename the sample matrix 
QUBO.MSAT.sampleNames <- unlist(lapply(rownames(QUBO.W.MSAT.genind@tab), function(x) strsplit(x, "_")[[1]][3]))
rownames(QUBO.W.MSAT.genind@tab) <- QUBO.MSAT.sampleNames

# SNP
# De novo
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUBO/output/populations_wild_R80_NOMAF_1SNP_Subset_Ordered2/"
QUBO.Subset.W.SNP.DN.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.Subset.W.SNP.DN.genind) <- factor(read.table(paste0(genpop.filePath, "QUBO_popmap_wild_Subset_Ordered2"), header=FALSE)[,2])
# Reference
genpop.filePath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_wild_R80_NOMAF_1SNP_Subset_Ordered2/"
QUBO.Subset.W.SNP.REF.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Correct popNames
pop(QUBO.Subset.W.SNP.REF.genind) <- factor(read.table(paste0(genpop.filePath, "QUBO_popmap_wild_Subset_Ordered2"), header=FALSE)[,2])
# Remove QUBO_W_ headers from sample names
QUBO.SNP.sampleNames <- gsub("QUBO_W_",replacement = "", row.names(QUBO.Subset.W.SNP.DN.genind@tab))
# Replace SH-Q names in SNP list with IMLS names
# These were determined by Austin K., and are outlined on the Hoban Lab Drive ("MSATcomparisons_TissueNames")
# Only 1 of the 11 SH_Q garden samples has an IMLS sample name (SHQ2177); others are unshared 
QUBO.SNP.sampleNames <- gsub("SH_Q2177",replacement = "IMLS338", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2178",replacement = "IMLS312", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2179",replacement = "IMLS062", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2180",replacement = "IMLS051", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2181",replacement = "IMLS011", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2182",replacement = "IMLS144", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2183",replacement = "IMLS170", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2184",replacement = "IMLS005", QUBO.SNP.sampleNames)
QUBO.SNP.sampleNames <- gsub("SH_Q2186",replacement = "IMLS017", QUBO.SNP.sampleNames)
# Rename sample matrices
rownames(QUBO.Subset.W.SNP.DN.genind@tab) <- QUBO.SNP.sampleNames
rownames(QUBO.Subset.W.SNP.REF.genind@tab) <- QUBO.SNP.sampleNames

# SUBSET MSAT GENIND
QUBO_sharedSamples <- sort(QUBO.SNP.sampleNames[which(QUBO.SNP.sampleNames %in% QUBO.MSAT.sampleNames)])
# Subset MSAT matrix to strictly shared samples, dropping now absent alleles
QUBO.Subset.W.MSAT.genind <- QUBO.W.MSAT.genind[QUBO_sharedSamples,, drop=TRUE]

# ---- RUN AND PLOT DAPC ----
# MSAT ----
# Use k-means clustering to find a number of groups which maximizes the variation between groups
# (after transforming the data using PCA, in order to increase computation times)
QUBO.Subset.W.MSAT.grp <- find.clusters(QUBO.Subset.W.MSAT.genind, max.n.clust=20, n.pca = 150, n.clust = 3)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# Clusters (seeking to minimize the Bayesian Information Criterion value): 3
# Conduct DAPC: transform the data using PCA, then run discriminant analysis on retained principal components,
# using the inferred groupings from the k-means clustering step above
QUBO.Subset.W.MSAT.dapc1 <- dapc(x=QUBO.Subset.W.MSAT.genind, pop=QUBO.Subset.W.MSAT.grp$grp, n.pca = 40, n.da = 2)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 40
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 2
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUBO.Subset.W.MSAT.dapc1$grp, as.character(pop(QUBO.Subset.W.MSAT.genind)), rownames(QUBO.Subset.W.MSAT.genind@tab))

# Show DAPC as a scatterplot
scatter(QUBO.Subset.W.MSAT.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=1.3, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"),
        txt.leg = c("Worldsong/Irondale/EBSCO/Hinds", "Irondale/Pop11/BTKC/Hinds", 
                    "EBSCOWattsville/Peavine/Worldsong"))
mtext("QUBO MSAT: Subset Wild", adj=0.65)

# SNP: DE NOVO (R80) ----
# SUPPORTED CLUSTERS ----
# K-means clustering step
# QUBO.Subset.W.SNP.DN.grp <- find.clusters(QUBO.Subset.W.SNP.DN.genind, max.n.clust=20, n.pca = 150, n.clust = 2)
QUBO.Subset.W.SNP.DN.grp <- find.clusters(QUBO.Subset.W.SNP.DN.genind)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# Clusters (seeking to minimize the Bayesian Information Criterion value): 2
# PCA step
# QUBO.Subset.W.SNP.DN.dapc1 <- dapc(x=QUBO.Subset.W.SNP.DN.genind, pop=QUBO.Subset.W.SNP.DN.grp$grp, n.pca = 70, n.da = 1)
QUBO.Subset.W.SNP.DN.dapc1 <- dapc(x=QUBO.Subset.W.SNP.DN.genind, pop=QUBO.Subset.W.SNP.DN.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 70
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 1
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUBO.Subset.W.SNP.DN.dapc1$grp, as.character(pop(QUBO.Subset.W.SNP.DN.genind)), rownames(QUBO.Subset.W.SNP.DN.genind@tab))

# Show DAPC as a scatterplot
scatter(QUBO.Subset.W.SNP.DN.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=1.3, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"),
        txt.leg = c("Oakbrook/Irondale/Hinds/Pop11", "Worldsong/EBSCO/BTKC/Wattsville/Peavine"))
mtext("QUBO SNP De novo (R80): Subset Wild (Supported Clusters)", adj=0.15)

# GEOGRAPHIC CLUSTERS (K=6) ----
# K-means clustering step
# QUBO.Subset.W.SNP.DN.grp <- find.clusters(QUBO.Subset.W.SNP.DN.genind, max.n.clust=20, n.pca = 150, n.clust = 6)
QUBO.Subset.W.SNP.DN.grp <- find.clusters(QUBO.Subset.W.SNP.DN.genind)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# Clusters (seeking to minimize the Bayesian Information Criterion value): 6
# PCA step
# QUBO.Subset.W.SNP.DN.dapc1 <- dapc(x=QUBO.Subset.W.SNP.DN.genind, pop=QUBO.Subset.W.SNP.DN.grp$grp, n.pca = 70, n.da = 5)
QUBO.Subset.W.SNP.DN.dapc1 <- dapc(x=QUBO.Subset.W.SNP.DN.genind, pop=QUBO.Subset.W.SNP.DN.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 70
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 5
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUBO.Subset.W.SNP.DN.dapc1$grp, as.character(pop(QUBO.Subset.W.SNP.DN.genind)), rownames(QUBO.Subset.W.SNP.DN.genind@tab))

# Show DAPC as a scatterplot
scatter(QUBO.Subset.W.SNP.DN.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=1.3, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5", "#D95F02","#8C510A"),
        txt.leg = c("Hinds/MossRock/Oakbrook", "Worldsong/Wattsville", "Peavine", "EBSCO/BTKC", "Irondale", "Pop11"))
mtext("QUBO SNP De novo (R80): Subset Wild (Geographic Clusters)", adj=0.15)

# SNP: REFERENCE (R80) ----
# SUPPORTED CLUSTERS ----
# K-means clustering step
# QUBO.Subset.W.SNP.REF.grp <- find.clusters(QUBO.Subset.W.SNP.REF.genind, max.n.clust=20, n.pca = 150, n.clust = 2)
QUBO.Subset.W.SNP.REF.grp <- find.clusters(QUBO.Subset.W.SNP.REF.genind)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# Clusters (seeking to minimize the Bayesian Information Criterion value): 2
# PCA step
# QUBO.Subset.W.SNP.REF.dapc1 <- dapc(x=QUBO.Subset.W.SNP.REF.genind, pop=QUBO.Subset.W.SNP.REF.grp$grp, n.pca = 70, n.da = 1)
QUBO.Subset.W.SNP.REF.dapc1 <- dapc(x=QUBO.Subset.W.SNP.REF.genind, pop=QUBO.Subset.W.SNP.REF.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 70
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 1
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUBO.Subset.W.SNP.REF.dapc1$grp, as.character(pop(QUBO.Subset.W.SNP.REF.genind)), rownames(QUBO.Subset.W.SNP.REF.genind@tab))

# Show DAPC as a scatterplot
scatter(QUBO.Subset.W.SNP.REF.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=1.3, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"),
        txt.leg = c("Worldsong/EBSCO/BTKC/Peavine/Wattsville", "Oakbrook/Irondale/Hinds/Pop11/Moss Rock"))
mtext("QUBO SNP Reference (R80): Subset Wild (Geographic Clusters)", adj=0.15)

# GEOGRAPHIC CLUSTERS (K=6) ----
# K-means clustering step
# QUBO.Subset.W.SNP.REF.grp <- find.clusters(QUBO.Subset.W.SNP.REF.genind, max.n.clust=20, n.pca = 150, n.clust = 6)
QUBO.Subset.W.SNP.REF.grp <- find.clusters(QUBO.Subset.W.SNP.REF.genind)
# Retained PCs (here, there is no cost to retaining a lot of PCs, even greater than xlim): 150
# Clusters (seeking to minimize the Bayesian Information Criterion value): 6
# PCA step
# QUBO.Subset.W.SNP.REF.dapc1 <- dapc(x=QUBO.Subset.W.SNP.REF.genind, pop=QUBO.Subset.W.SNP.REF.grp$grp, n.pca = 70, n.da = 5)
QUBO.Subset.W.SNP.REF.dapc1 <- dapc(x=QUBO.Subset.W.SNP.REF.genind, pop=QUBO.Subset.W.SNP.REF.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 70
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 5
# List grouping of each individual (helpful for determining correct coloration)
cbind(QUBO.Subset.W.SNP.REF.dapc1$grp, as.character(pop(QUBO.Subset.W.SNP.REF.genind)), rownames(QUBO.Subset.W.SNP.REF.genind@tab))

# Show DAPC as a scatterplot
scatter(QUBO.Subset.W.SNP.REF.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=1.3, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5", "#D95F02","#8C510A"),
        txt.leg = c("Pop11", "EBSCO/BTKC", "Worldsong/Peavine", "Wattsville", "Irondale", "EBSCO/BTKC"))
mtext("QUBO SNP Reference (R80): Subset Wild (Geographic Clusters)", adj=0.5, line=2.5)

# %%%% SUBSET: GARDEN AND WILD ----
# To be populated
