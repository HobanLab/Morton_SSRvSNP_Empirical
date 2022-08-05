# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS (DAPC) %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses the adegenet package to run DAPC on two datasets
# For QUAC, the de novo final assembly (DNFA) dataset is used
# For QUBO, the dataset aligned to the Quercus robur reference genome (GSNAP4) is used
# In both analyses, garden and wild samples are included 

library(adegenet)
library(scales)

# TO DO: Emily wants you to run a normal PCA on the dataset (i.e. not through adegenet; instructions
# in adegenet Basics tutorial)

# %%%% QUAC %%%% ----
# Specify filepath for finalized QUAC dataset, with both garden and wild samples, and the following filters:
# R80 (no filter for missing data), no minor allele frequency filter, first SNP per locus, 2 populations specified (garden and wild)
QUAC.genpop.filepath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/output/populations_R80_NOMAF_1SNP_2Pops/"
# Read in QUAC genind file
QUAC.genind <- read.genepop(paste0(QUAC.genpop.filepath,"populations.snps.gen"))

# Use k-means clustering to find a number of groups which maximizes the variation between groups
# (after transforming the data using PCA, in order to increase computation times)
QUAC.grp <- find.clusters(QUAC.genind, max.n.clust=20)
# Retained PCs (here, there is no cost to retaining a lot of PCs): 200
# Clusters (seeking to minimize the Bayesian Information Criterion value): 4
QUAC.grp

# Conduct DAPC: transform the data using PCA, then run discriminant analysis on retained principal components,
# using the inferred groupings from the k-means clustering step above
QUAC.dapc1 <- dapc(x=QUAC.genind, pop=QUAC.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting/unstable membership probability): 50
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 2

# List grouping of each individual (helpful for determining correct coloration)
QUAC.dapc1$grp
QUAC.dapc1$grp[97:193]
read.table("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/QUAC_popmap_wild")[,2]
cbind(QUAC.dapc1$grp[97:193], read.table("/RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/QUAC/QUAC_popmap_wild")[,2])

# Plot DAPC as a scatterplot
scatter(QUAC.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"),
        txt.leg = c("Pryor Mt.","Magazine Mt. 2", "Porter Mt.", "Kessler/Sugarloaf/Magazine Mt. 1"))
mtext("QUAC R80: Garden and Wild", adj=0.02)

# %%%% QUBO %%%% ----
# Specify filepath for finalized QUBO dataset, with both garden and wild samples, and the following filters:
# R0 (no filter for missing data), no minor allele frequency filter, first SNP per locus, 2 populations specified (garden and wild)
QUBO.genpop.filepath <- 
  "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP4/output/populations_R80_NOMAF_1SNP_2Pops/"
# Read in QUAC genind file
QUBO.genind <- read.genepop(paste0(QUBO.genpop.filepath,"populations.snps.gen"))

# Use k-means clustering to find a number of groups which maximizes the variation between groups
# (after transforming the data using PCA, in order to increase computation times)
QUBO.grp <- find.clusters(QUBO.genind, max.n.clust=20)
# Retained PCs (here, there is no cost to retaining a lot of PCs): 200
# Clusters (seeking to minimize the Bayesian Information Criterion value): 3
QUBO.grp

# Conduct DAPC: transform the data using PCA, then run discriminant analysis on retained principal components,
# using the inferred groupings from the k-means clustering step above
QUBO.dapc1 <- dapc(x=QUBO.genind, pop=QUBO.grp$grp)
# Retained PCs (specifying a large value can lead to overfitting, and unstable membership probability): 50
# Discriminant functions to retain (for less than 10 clusters, all eigenvalues can be retained): 4

# List grouping of each individual (helpful for determining correct coloration)
QUBO.dapc1$grp
QUBO.dapc1$grp[86:180]
read.table("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/QUBO_popmap_wild")[,2]
cbind(QUBO.dapc1$grp[86:180], read.table("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/QUBO_popmap_wild")[,2])

QUBO.dapc1$grp
read.table("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/QUBO_popmap_GardenWild")[1:85,2]
c(read.table("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/QUBO_popmap_GardenWild")[1:85,2],
  read.table("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/QUBO_popmap_wild")[,2])
cbind(QUBO.dapc1$grp, c(read.table("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/QUBO_popmap_GardenWild")[1:85,2],
                        read.table("/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/QUBO_popmap_wild")[,2]))

# Plot DAPC as a scatterplot
scatter(QUBO.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"), 
        txt.leg = c("Moss Rock/Oakbrook","EBSCO/Blue Trail Kings Chair/Peavine/Wattsville",
                    "Hinds Road/Pop 11/Wattsville"))
mtext("QUBO R80: Garden and Wild", adj=0.9)
