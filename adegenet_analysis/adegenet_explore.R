library(adegenet)
library(pegas)
library(hierfstat)

#setwd("D:/Iowa State University/Debinski Lab/Parnassius genetics/ParnassiusGenetics")

### Read in file - uncomment the dataset you want to work with

source("analysis/adegenet_analysis/ReadUnfilImp.R")  #unfiltered-imputed, small file
#source("analysis/adegenet_analysis/ReadFiltNotimp.r") #filtered-not imputed, medium-sized file
#source("analysis/adegenet_analysis/ReadUnfiltNotimp.r")  #unfiltered, not imputed, large file

### Basic pop gen stats

PCdata <- PCdata1

mydf <- as.data.frame(table(PCdata@pop))
good.sites <- as.numeric(mydf$Var1[which(mydf$Freq > 2)])

PCdata1 <- PCdata[good.sites,]
PCsumry <- summary(PCdata1)
PCsumry

plot(PCsumry$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", main="Observed heterozygosity per locus")

plot(PCsumry$Hobs,PCsumry$Hexp, xlab="Hobs", ylab="Hexp", main="Hexp as a function of Hobs per locus")

plot(PCsumry$n.by.pop, PCsumry$pop.n.all, xlab="Colonies sample size", ylab="Number of alleles",main="Alleles numbers and sample sizes", type="n")
text(PCsumry$n.by.pop,PCsumry$pop.n.all,lab=names(PCsumry$n.by.pop))

barplot(PCsumry$loc.n.all, ylab="Number of alleles", main="Number of alleles per locus")

barplot(PCsumry$Hexp-PCsumry$Hobs, main="Heterozygosity: expected-observed", ylab="Hexp - Hobs")

par(mar=c(10,4,4,2))
barplot(PCsumry$n.by.pop, main="Sample sizes per population", ylab="Number of genotypes",las=2)


#  Is mean observed H significantly lower than mean expected H?

bartlett.test(list(PCsumry$Hexp,PCsumry$Hobs))
#variances are very not equal (Ho = same variance)

t.test(PCsumry$Hexp,PCsumry$Hobs,pair=T,var.equal=F,alter="greater")
#mean observed H is not significantly lower than expected

#Hardy-Weinberg quilibrium
PCsumry.hwt <- hw.test(PCdata1, B=1000)
PCsumry.hwt
hist(PCsumry.hwt[,"Pr.exact"], xlab = "p-values", main = "p-values of by-locus HW test\nn=1000", breaks = 25)
abline(v = 0.05, col = "red")

#F statistics, 0 = panmixia; 1 = no shared genetic diversity (i.e., complete isolation)
my.fstat <- fstat(PCdata1)
my.fstat
#Panmixia


my.Fst <- Fst(as.loci(PCdata1))
my.Fst  #this is throwing "NaNs" - not sure what's up
#Monte Carlo test for genetic structuring of individuals given pop.
Gtest <- gstat.randtest(PCdata1, method = "global", nsim = 999) #G-statistic test
Gtest
plot(Gtest)
#p-value = 0.01, Ha -> proportions are greater than expected under HW equilibrium


popFst <- pairwise.fst(PCdata1[1:10,])  #first 10 individuals
popFst
is.euclid(popFst)

# Inbreeding
obj <- seppop(PCdata1)
dc <- repool(obj$`Death canyon Rangers`, obj$`Death canyon phelps lake`, obj$`Death Canyon Trail`)
pc <- repool(obj$`Aimees Meadow`, obj$`AMK Ranch`, obj$`AMK Road`, obj$`Cygnet Pond`, obj$`Dump Road`, obj$`Pilgrim Creek`)

temp <- inbreeding(PCdata1, N=1000) #throwing warnings about lack of precision
inbreed <- sapply(temp, mean)
hist(inbreed, main = "All populations")

temp <- inbreeding(dc, N=1000) #throwing warnings about lack of precision
inbreed <- sapply(temp, mean)  
hist(inbreed, main = "Death Canyon sites")

temp <- inbreeding(pc, N=1000) #throwing warnings about lack of precision
inbreed <- sapply(temp, mean)
hist(inbreed, main = "Pilgrim Creek area")

### Basic statistics with hierfstat

PCdata2 <- genind2hierfstat(PCdata1) # Create hierfstat object

#basicstat <- basic.stats(PCdata2) #Not working - non-conformable arrays
#names(basicstat)

mypca <- indpca(PCdata1)
plot(mypca, cex = 0.7)
