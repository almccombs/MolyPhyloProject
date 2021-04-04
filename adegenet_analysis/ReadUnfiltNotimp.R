library(adegenet)
library(pegas)
library(hierfstat)

source("analysis/adegenet_analysis/PopKey.R")
popkey <- popkey[,-c(3:4)]

pcdata <- read.vcf("../ParnassiusGeneticsData/parnassius_clodius_unfiltered_notimputed.vcf", which.loci = 1:1e6)
is.na(pcdata) <- pcdata == "./."

pcdata$SampleID <- rownames(pcdata)
PCdata <- merge(pcdata, popkey, by = "SampleID")
rm(popkey)
rm(pcdata)

ind <- PCdata$SampleID
population <- PCdata$SiteID
rownames(PCdata) <- PCdata[,1]
PCdata <- PCdata[, -c(1,42217)]

PCdata1 <- df2genind(PCdata, ploidy = 2, ind.names = ind, pop = population, sep = "/")
#PCdata1 <- df2genind(PCdata, ploidy = 2, sep = "/")
rm(PCdata)
rm(ind)
rm(population)
PCdata1
