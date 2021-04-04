library(adegenet)
library(pegas)
library(hierfstat)

#setwd("D:/Iowa State University/Debinski Lab/Parnassius genetics/ParnassiusGenetics")

Mydata <- read.table("analysis/Master_Pinus_data_genotype.txt", header = TRUE)
dim(Mydata)

locus <- Mydata[, -c(1, 2, 3, 4, 17:3086)]    
colnames(locus) <- gsub("\\.", "_", colnames(locus)) # locus names can't have "."
ind <- as.character(Mydata$tree_id) # labels of the individuals
population <- as.character(Mydata$state) # labels of the populations
Mydata1 <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep = "")
Mydata1

nAll(Mydata1) # Number of alleles per locus

Mydata2 <- genind2hierfstat(Mydata1) # Create hierfstat object

div <- summary(Mydata1)
div
names(div)

plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")

plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")





#Parnassius data
source("analysis/LEA_analysis/import_only.R") #Sean's script for reading the .vcf file

source("analysis/adegenet_analysis/PopKey.r")

mydata <- as.data.frame(t(mark))
rm(geno)
rm(mark)
rm(dat3)
rm(colzeros)

mydata$SpecimenID <- rownames(mydata)
Mydata <- merge(mydata, popkey, by = "SpecimenID")
rm(popkey)
rm(mydata)

#colnames(Mydata) <- gsub("NA.", "", colnames(Mydata))
colnames(Mydata) <- gsub("_", "-", colnames(Mydata))
colnames(Mydata) <- gsub("\\.", "_", colnames(Mydata))

ind <- Mydata$SpecimenID
population <- Mydata$PopulationCode
rownames(Mydata) <- Mydata[,1]
Mydata <- Mydata[, -c(1,1003)]

mydata1 <- df2genind(Mydata, ploidy = 2, ind.names = ind, pop = population, sep = "")
rm(Mydata)
mydata1

nAll(mydata1)

mydata2 <- genind2hierfstat(mydata1) # Create hierfstat object

div <- summary(mydata1)
div
names(div)

plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")

plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")

bartlett.test(list(div$Hexp, div$Hobs)) # a test : H0: Hexp = Hobs

# These two functions are throwing an error: Error in sHo/2/n : non-conformable arrays
basicstat <- basic.stats(mydata2, diploid = TRUE, digits = 2) 
names(basicstat) 
boot.ppfis(mydata2) 

x <- indpca(mydata2) 
plot(x, cex = 0.7)

hw.test(mydata1, B = 1000)


