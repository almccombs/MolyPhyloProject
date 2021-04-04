library(adegenet)
library(pegas)

mydata <- read.vcf("SNPdata/parnassius_clodius_unfiltered_imputed.vcf")

mydata <- mydata[,2:(ncol(mydata))]

mydata.gen <- df2genind(mydata, ploidy = 2, sep = "/")
summary(mydata.gen)
