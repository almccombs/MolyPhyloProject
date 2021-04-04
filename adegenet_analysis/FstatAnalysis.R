library(adegenet)
library(pegas)
library(hierfstat)

#setwd("D:/Iowa State University/Debinski Lab/Parnassius genetics/ParnassiusGenetics")

### Read in file - uncomment the dataset you want to work with

source("analysis/adegenet_analysis/ReadUnfilImp.R")  #unfiltered-imputed, small file

### Basic pop gen stats

my.Fst <- Fst(as.loci(PCdata1))
my.Fst

my.fstat <- fstat(PCdata1)
my.fstat

Fst.pairs <- pairwise.fst(PCdata1)
Fst.pairs

Gtest <- gstat.randtest(PCdata1, method = "global", nsim = 9999) #G-statistic test
 Gtest
plot(Gtest)

### Basic statistics with hierfstat

PCdata2 <- genind2hierfstat(PCdata1) # Create hierfstat object

basicstat <- basic.stats(PCdata2) #Not working - non-conformable arrays
names(basicstat)
