# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% POPULATION DIFFERENTIATION %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script quantifies the differences between wild populations (using DAPC and STRUCTURE)
# for the two study species of the SSRvSNP study: Quercus acerifolia (QUAC) and Q. boyntonii (QUBO)

library(adegenet)
library(pegas)
library(hierfstat)
library(RColorBrewer)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ---- QUERCUS ACERIFOLIA ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%

# READ IN/PROCESS GEN FILE----
genpop.filePath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUAC/GSNAP/output/populations_sum/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
QUAC.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Read in the Stacks popmap values, and use these to replace @pop values 
# (using the pops accessor; this is necessary because original pop names are incorrect)
pop(QUAC.genind) <- factor(read.table("../../../QUAC_popmap2", header=FALSE)[,2])

# Using Sean's modified repool function, from the SE_oaks_genetics repo
# https://github.com/smhoban/SE_oaks_genetics/blob/main/3_oak_pop_gen.R
repool_new<- function(genind_obj,vect_pops){
  # Setting reppop to use drop=TRUE leads to the error "Error in table(temp) : all arguments must have the same length"
  genind_obj_sep<-seppop(genind_obj)
  genind_obj_merge<-genind_obj_sep[[vect_pops[1]]]
  for (i in 1:(length(vect_pops)-1)) genind_obj_merge<-repool(genind_obj_merge,genind_obj_sep[[vect_pops[i+1]]])
  genind_obj_merge
}
# Repool samples, to generate a genind object of solely wild individuals
QUAC.genind.wild <- repool_new(QUAC.genind,c(2:6))

# Note the differences between the values below. The wild genind is missing loci specific to gardens
length(which(is.na(colSums(QUAC.genind@tab))))
length(which(is.na(colSums(QUAC.genind.wild@tab))))

# Calling standard seppop and repool functions fails to work
# QUAC.genind.wild <- seppop(QUAC.genind, drop=TRUE)[2:6]
# QUAC.genind.wild <- repool(QUAC.genind.wild$porterMt, QUAC.genind.wild$magazineMt, QUAC.genind.wild$pryorMt,
#                            QUAC.genind.wild$sugarloaf_midlandPeak, QUAC.genind.wild$kessler_shaleBarrenRidge)
# Error in table(temp) : attempt to make a table with >= 2^31 elements

# Find number of clusters
QUAC_DAPC.grp <- find.clusters(QUAC.genind.wild, max.n.clust=90)
# Number of retained PCs (here, there is no cost to retaining a lot of PCs): 90
# Number of clusters: 4 --- but note that THE MOST SUPPORTED NUMBER OF CLUSTERS IS 2 (AND, PROBABLY, 1!)
QUAC_DAPC.grp

# Conduct DAPC
QUAC_DAPC <- dapc(x=QUAC.genind.wild, pop=QUAC_DAPC.grp$grp)

# Original DAPC image
scatter(QUAC_DAPC, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = brewer.pal(n=4, name="Dark2"), 
        txt.leg = c("Magazine Mt.","Pryor Mt.","Porter/Kessler","Sugar Loaf"))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# ---- QUERCUS BOYNTONII ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# READ IN/PROCESS GEN FILE----
genpop.filePath <- "/RAID1/IMLS_GCCO/Analysis/Stacks/reference_filteredReads/QUBO/GSNAP/output/populations_sum/"
setwd(genpop.filePath)
# Pull in genepop object. For RADseq data, ncode argument = 2 (default) 
QUBO.genind <- read.genepop(paste0(genpop.filePath,"populations.snps.gen"))
# Read in the Stacks popmap values, and use these to replace @pop values 
# (using the pops accessor; this is necessary because original pop names are incorrect)
pop(QUBO.genind) <- factor(read.table("../../../QUBO_popmap2", header=FALSE)[,2])
# Repool samples, to generate a genind object of solely wild individuals
QUBO.genind.wild <- repool_new(QUBO.genind,c(2:12))

# Note the differences between the values below. The wild genind is missing loci specific to gardens
length(which(is.na(colSums(QUBO.genind@tab))))
length(which(is.na(colSums(QUBO.genind.wild@tab))))

# Find number of clusters
QUBO_DAPC.grp <- find.clusters(QUBO.genind.wild, max.n.clust=90)
# Number of retained PCs (here, there is no cost to retaining a lot of PCs): 90
# Number of clusters: 3 --- but note that THE MOST SUPPORTED NUMBER OF CLUSTERS IS 2 (AND, PROBABLY, 1!)
QUBO_DAPC.grp

# Conduct DAPC
QUBO_DAPC <- dapc(x=QUBO.genind.wild, pop=QUBO_DAPC.grp$grp)

# Original DAPC image
scatter(QUBO_DAPC, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = brewer.pal(n=3, name="Dark2"), 
        txt.leg = c("Worldsong/BTKC","Pop 11", "Oakbrook/Irondale/Hinds Road"))
