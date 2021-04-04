library(adegenet)
library(pegas)
library(hierfstat)

#setwd("D:/Iowa State University/Debinski Lab/Parnassius genetics/ParnassiusGenetics")

### Read in file - uncomment the dataset you want to work with

source("analysis/adegenet_analysis/ReadUnfilImp.R")  #unfiltered-imputed, small file
#source("analysis/adegenet_analysis/ReadFiltNotimp.r") #filtered-not imputed, medium-sized file
#source("analysis/adegenet_analysis/ReadUnfiltNotimp.r")  #unfiltered, not imputed, large file


### Explore population structure with adegenet

# Find genetic clusters/groups
#grp <- find.clusters(PCdata1, max.n.clusters = 40)  #interactive
grp <- find.clusters(PCdata1, n.pca = 140, n.clust = 3, max.n.clust=40)
names(grp)
grp$Kstat
grp$stat
grp$grp #group membership for each individual
grp$size

# DAPC of data with existing populations
#dapc.pop <- dapc(PCdata1, grp$grp) #interactive
dapc.pop <- dapc(PCdata1, grp$grp,  n.pca = 45, n.da = 4)
dapc.pop

# Plots
#scatter(dapc.pop,1,1,bg="white", scree.da=FALSE, legend=TRUE, solid=.4)

scatter(dapc.pop, xax = 1, yax = 2, bg="white", pch=20, cell=0, cstar=0, solid=.4, cex=1, clab=0, scree.pca=TRUE, posi.pca="topright", ratio.da = .15, ratio.pca = .15, leg=TRUE, txt.leg=paste("Pop",1:17), posi.leg="bottomleft")

scatter(dapc.pop, xax = 1, yax = 2, bg="white", pch=20, cell=0, cstar=0, solid=.4, cex=1, clab=0,scree.da = FALSE, scree.pca=FALSE, posi.pca="topright", ratio.da = .15, ratio.pca = .15, leg=FALSE, txt.leg=paste("Pop",1:17), posi.leg="bottomleft")

# Loadings (i.e., variable contributions)
set.seed(4)
contrib.pop <- loadingplot(dapc.pop$var.contr, axis=2, thres=.07, lab.jitter=1)

# Probabilities for assignment of individuals into groups
round(head(dapc.pop$posterior),3)
#round(dapc.pop$posterior,3)
summary(dapc.pop)  #$assign.per.pop is probabilities of correct assignment into each group

# plot probabilities
assignplot(dapc.pop)
assignplot(dapc.pop, only.grp = 1)
assignplot(dapc.pop, only.grp = 2)
assignplot(dapc.pop, only.grp = 3)
#assignplot(dapc.pop, only.grp = 4)
#assignplot(dapc.pop, only.grp = 5)

compoplot(dapc.pop, posi="bottomright", txt.leg=paste("Population", 1:5), lab="", ncol=1, xlab="individuals", col=funky(7))
compoplot(dapc.pop, posi="bottomright", txt.leg=paste("Population", 1:5), lab="", ncol=1, xlab="individuals", col=funky(7), only.grp = 1)
compoplot(dapc.pop, posi="bottomright", txt.leg=paste("Population", 1:5), lab="", ncol=1, xlab="individuals", col=funky(7), only.grp = 2)
compoplot(dapc.pop, posi="bottomright", txt.leg=paste("Population", 1:5), lab="", ncol=1, xlab="individuals", col=funky(7), only.grp = 3)
#compoplot(dapc.pop, posi="bottomright", txt.leg=paste("Population", 1:5), lab="", ncol=1, xlab="individuals", col=funky(7), only.grp = 4)
#compoplot(dapc.pop, posi="bottomright", txt.leg=paste("Population", 1:5), lab="", ncol=1, xlab="individuals", col=funky(7), only.grp = 5)

# Prior vs post memberships
table(pop(PCdata1), grp$grp)

