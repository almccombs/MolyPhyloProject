library(adegenet)
library(pegas)
library(hierfstat)

source("analysis/adegenet_analysis/PopKey.R")
popkey <- popkey[,-c(3:4)]

pcdata <- read.vcf("../ParnassiusGeneticsData/parnassius_clodius_unfiltered_imputed.vcf", which.loci = 1:1e6)
is.na(pcdata) <- pcdata == "./."

pcdata$SampleID <- rownames(pcdata)
PCdata <- merge(pcdata, popkey, by = "SampleID")
rm(popkey)
rm(pcdata)

mydf <- as.data.frame(table(PCdata$SiteID))
bad.sites <- as.character(mydf$Var1[which(mydf$Freq <3)])
PCdata <- PCdata[which(PCdata$SiteID != "Bear paw lake intersection" & PCdata$SiteID != "Buffalo fork" & PCdata$SiteID != "Death Canyon Trail" & PCdata$SiteID != "Dump Road" & PCdata$SiteID != "Grand view 1" & PCdata$SiteID != "Hidden Falls Trail" & PCdata$SiteID != "Lozier Road" & PCdata$SiteID != "Lupine Meadow" & PCdata$SiteID != "Surprise fall down" & PCdata$SiteID != "Two Ocean Road 1"),]
PCdata <- droplevels(PCdata)
rm(bad.sites)
rm(mydf)

ind <- PCdata$SampleID
population <- PCdata$SiteID
rownames(PCdata) <- PCdata[,1]
PCdata <- PCdata[, -c(1,1003)]

PCdata1 <- df2genind(PCdata, ploidy = 2, ind.names = ind, pop = population, sep = "/")
#PCdata1 <- df2genind(PCdata, ploidy = 2, sep = "/")
rm(PCdata)
rm(ind)
rm(population)
PCdata1
