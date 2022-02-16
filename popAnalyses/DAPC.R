library(adegenet)
library(scales)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% DE NOVO FINAL ASSEMBLIES (DNFA) %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Working directory set to /RAID1/IMLS_GCCO/Analysis/Stacks/denovo_finalAssemblies/

# QUAC de novo final assembly (DNFA)----
# Working with wild samples only (populations_wild dataset)
# Pull in genepop object (with file suffix updated--Stacks writes as ".genepop", but needs to be ".gen")
QUAC_DNFA.genind <- read.genepop("QUAC/output/populations_wild/populations.snps.gen")
# Find number of clusters
QUAC_DNFA.grp <- find.clusters(QUAC_DNFA.genind, max.n.clust=20)
# Number of retained PCs (here, there is no cost to retaining a lot of PCs): 100
# Number of clusters: 6
QUAC_DNFA.grp

# Conduct DAPC
QUAC_DNFA.dapc1 <- dapc(x=QUAC_DNFA.genind, pop=QUAC_DNFA.grp$grp)
# Number of retained PCs: 100
# Number of discriminant functions to retain: 5
str(QUAC_DNFA.dapc1)
# List grouping of each individual (helpful for determining correct coloration)
QUAC_DNFA.dapc1$grp

# # Make a scatterplot of initial DAPC
# scatter(dapc1)
# scatter(dapc1, label = NULL)
# scatter(dapc1, clab = 0.5, label = letters[1:15])
# scatter(dapc1, clab = 0.5, scree.da = F)

# Original DAPC image
scatter(QUBO_DNFA.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"), 
        txt.leg = c("Porter/Kessler/Sugarloaf","Magazine Mt. 1","Pryor Mt.","Magazine Mt. 2"))

# QUBO de novo final assembly (DNFA)----
# Working with wild samples only (populations_wild dataset)
# Pull in genepop object (with file suffix updated--Stacks writes as ".genepop", but needs to be ".gen")
QUBO_DNFA.genind <- read.genepop("QUBO/output/populations_Sum/populations.snps.gen")
# Find number of clusters
QUBO_DNFA.grp <- find.clusters(QUBO_DNFA.genind, max.n.clust=20)
# Number of retained PCs (here, there is no cost to retaining a lot of PCs): 175
# Number of clusters: 4
QUBO_DNFA.grp

# Conduct DAPC
QUBO_DNFA.dapc1 <- dapc(x=QUBO_DNFA.genind, pop=QUBO_DNFA.grp$grp)
# Number of retained PCs: 150
# Number of discriminant functions to retain: 3
str(QUBO_DNFA.dapc1)
# List grouping of each individual (helpful for determining correct coloration)
QUBO_DNFA.dapc1$grp

# Original DAPC image
scatter(QUBO_DNFA.dapc1, scree.da=F, bg="white", pch=20, cell=0, cstar=0, solid=0.6, clab=0, legend=T,
        posi.leg=locator(n=1), cleg=1.0, cex=2, inset.solid=1,
        col = c("#66A61E","#E7298A","#7570B3","#2171B5"), 
        txt.leg = c("Moss Rock Preserve","Irondale/Wattsville/Pop11","Oakbrook","Worldsong/Kings Chair/EBSCO/Peavine"))
