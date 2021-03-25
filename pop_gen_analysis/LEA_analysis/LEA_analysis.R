#Use read.vcf.markerstats3 for STACKS v1.44

setwd("C:/Users/Audrey McCombs/Desktop/ParnassiusGenetics")

##### For reading attributes from a VCF where only MAF is recorded (no entry for major allele frequency)
read.vcf.markerstats3 <- function(filename,max.marker) {
  con <- file(filename,"r") #open file for reading
  temp <- readLines(con,1)  #read one line
  comment.line <- 0
  while(substr(temp,1,2)=="##") {  #skip comment lines
    temp <- readLines(con,1)
    comment.line <- comment.line+1
  }
  header <- strsplit(temp,split="\t",fixed=TRUE)[[1]]
  
  #FORMAT is position 9
  n.sample <- length(header) - 9
  sample.names <- apply(array(header[-c(1:9)]),1,function(x){y=strsplit(x,split="\t",fixed=TRUE)[[1]][1];return(y)})
  close(con)
  
  marker.stats <- data.frame(Chrom=rep("",max.marker),Pos=rep(0,max.marker),Name=rep("",max.marker),NSamp=rep(0,max.marker),MAF=rep(0.0,max.marker),stringsAsFactors=F,Allele1=rep("",max.marker),Allele2=rep("",max.marker))
  GT <- matrix(0,max.marker,n.sample)
  colnames(GT) <- sample.names
  DP <- matrix(0,max.marker,n.sample)
  colnames(DP) <- sample.names
  
  con <- file(filename,"r")
  temp <- readLines(con,comment.line+1)
  m <- 0
  
  while ((m < max.marker) & (length(temp)>0)) {
    temp <- readLines(con,1)
    if (length(temp) > 0) {
      temp2 <- strsplit(temp,split="\t",fixed=TRUE)[[1]]
      if ((length(grep("/",temp2[5],fixed=TRUE))==0)&(temp2[5]!="-")&(temp2[4]!="-")&(length(grep(",",temp2[4],fixed=TRUE))==0)) {
        #only process bi-allelic SNPs (remove tri-allelic and indels)
        #for MAF
        temp3 <- strsplit(temp2[8],split=";",fixed=T)[[1]]
        #for NSamp
        temp4 <- strsplit(temp2[8],split=";",fixed=TRUE)[[1]]
        m <- m + 1
        marker.stats[m,"Chrom"] <- temp2[1]
        marker.stats[m,"Pos"] <- as.integer(temp2[2])
        marker.stats[m,"Name"] <- paste(temp2[2],temp2[3],sep=".")
        marker.stats[m,"NSamp"] <- as.integer(strsplit(temp4[1],split="=",fixed=T)[[1]][2])
        marker.stats[m,"MAF"] <- as.numeric(strsplit(temp3[2],split="=",fixed=T)[[1]][2])
        marker.stats[m,"Allele1"] <- temp2[4]
        marker.stats[m,"Allele2"] <- temp2[5]
      }
    }
    print(paste("Marker",m,"done"))
  }
  close(con) 
  
  marker.stats <- marker.stats[1:m,]
  rownames(marker.stats) <- marker.stats$Name
  
  return(marker.stats)
}  #end read.vcf


##command for Parnassius clodius data file, number is arbitrary but large
dat3 <- read.vcf.markerstats3(filename='SNPdata/parnassius_clodius_unfiltered_imputed.vcf',10000) 


read.vcf.genotypes <- function(filename,max.marker) {
  
  get_counts <- function(x,GT.pos) {
    if (length(grep(":",x,fixed=TRUE))>0) {
      y <- strsplit(x,split=":",fixed=TRUE)[[1]][GT.pos]
      z <- strsplit(y,split="/",fixed=TRUE)[[1]]
      return(as.numeric(z))
    } else {
      return(c(0,0))
    }
  }
  
  con <- file(filename,"r") #open file for reading
  temp <- readLines(con,1)  #read one line
  comment.line <- 0
  while(substr(temp,1,2)=="##") {  #skip comment lines
    temp <- readLines(con,1)
    comment.line <- comment.line+1
  }
  header <- strsplit(temp,split="\t",fixed=TRUE)[[1]]
  
  #format is position 9
  n.sample <- length(header) - 9
  
  sample.names <- apply(array(header[-c(1:9)]),1,function(x){y=strsplit(x,split=":",fixed=TRUE)[[1]][1];return(y)})
  
  temp <- readLines(con,1)
  temp2 <- strsplit(temp,split="\t",fixed=TRUE)[[1]]
  temp3 <- strsplit(temp2[9],split=":",fixed=TRUE)[[1]]
  GT.pos <- match("GT",temp3)
  close(con)
  
  print("Got genotype positions")
  
  genotypes1 <- matrix(NA,max.marker,n.sample)
  genotypes2 <- matrix(NA,max.marker,n.sample)
  genotypes <- matrix(NA,max.marker,n.sample)
  colnames(genotypes1) <- sample.names
  colnames(genotypes2) <- sample.names
  colnames(genotypes) <- sample.names
  marker.names <- array(rep("",max.marker))
  
  con <- file(filename,"r")
  temp <- readLines(con,comment.line+1)
  m <- 0
  
  while ((m < max.marker) & (length(temp)>0)) {
    temp <- readLines(con,1)
    if (length(temp) > 0) {
      temp2 <- strsplit(temp,split="\t",fixed=TRUE)[[1]]
      temp3 <- strsplit(temp2[1],split="d",fixed=TRUE)[[1]] #splitting by "d" because the scaffold number is preceded by "scaffold"
      if ((length(grep(",",temp2[5],fixed=TRUE))==0)&(temp2[5]!="-")&(temp2[4]!="-")&(length(grep(",",temp2[4],fixed=TRUE))==0)) {
        
        #only process bi-allelic SNPs (remove tri-allelic and indels)
        counts <- apply(array(temp2[-c(1:9)]),1,get_counts,GT.pos)
          m <- m + 1
          genotypes1[m,] <- counts[1,]
          genotypes2[m,] <- counts[2,]
          marker.names[m] <- paste(temp3[2],temp2[2],temp2[3],sep=".")
      }
    }
    print(paste("Marker",m,"done",sep=" "))
  }
  close(con) 
  
  i <- seq(1,n.sample,by=1)
  j <- seq(1,max.marker,by=1)
  for (i in 1:n.sample){
    for (j in 1:max.marker){
      genotypes[j,i] <- sum(genotypes1[j,i],genotypes2[j,i])
    }
  }
  genotypes[1:m,]
  rownames(genotypes) <- marker.names
  
  return(genotypes)
  
}  #end read.vcf



#command for Parnassius clodius, number refers to loci
geno <- read.vcf.genotypes(filename="SNPdata/parnassius_clodius_unfiltered_imputed.vcf",1001)
geno <- t(geno)

#tabulate frequencies
apply(geno[,1:5],2,table)
table(sapply(apply(geno,2,table), length))

#Conduct PCA and SNMF clustering
source("http://bioconductor.org/biocLite.R")
#biocLite("LEA")
library('LEA')
library('maps')
#install.packages('RColorBrewer')
library('RColorBrewer')
colors <- brewer.pal(5,"Accent")

#Remove invariant SNPs
colzeros <- apply(geno,2,sd)==0
mark <- geno[,colzeros==F]
rm(colzeros)

#Commenting this out because every time you run this it takes forever to commit. They should be good as is, unless we're analyzing a new dataset.
#write.geno(mark,"analysis/LEA_analysis/ParaFiles/para.geno")
#write.lfmm(mark,"analysis/LEA_analysis/ParaFiles/para.lfmm")

pc <- pca("analysis/LEA_analysis/ParaFiles/para.lfmm",scale=TRUE)

#Tracey Widom test.  From R doc: Perform tracy-widom tests on a set of eigenvalues to determine the number of significative eigenvalues and calculate the percentage of variance explained by each principal component.

tw <- tracy.widom(pc)
barplot(tw$pvalues[1:10],ylab="P value",xlab="PC",names.arg=1:10)
abline(h=0.05,lty=2)
plot(pc$sdev[1:10]^2/sum(pc$sdev^2),xlab="PC",ylab="Fraction Variation Explained")
plot(tw$percentage[1:10], xlim = c(0,10))  #3 major genetic clusters in the data
plot(pc$projections)

#Compute admixture
snmf2 <- snmf("analysis/LEA_analysis/ParaFiles/para.geno",K=1:10,ploidy=2,entropy=T,alpha=100,project="new")
plot(snmf2,col="blue4",cex=1.4,pch=19)  #minimum at K=3
K=3
snmf1 = snmf("analysis/LEA_analysis/ParaFiles/para.geno", K = K, alpha = 100, project = "new")
qmatrix = Q(snmf1, K = K)
barplot(t(qmatrix),col=brewer.pal(K,"Paired"),border=NA,space=0,xlab="Individuals",ylab="Admixture coefficients")

rm(K)

# create df for clodius
source("analysis/adegenet_analysis/PopKey.R")
coord <- popkey[,c(1,2,4,3)]
names(coord)[3:4] <- c("Latitude", "Longitude")
rm(popkey)

  #make df with only the good samples
ids <- noquote(rownames(mark))
samples_used <- (coord$SampleID %in% ids)
coords <- coord[samples_used,]
rm(ids)
rm(coord)
rm(samples_used)

  #get actual sample size for each pop
samp.size.all <- as.data.frame(table(coords$SiteID))
names(samp.size.all) <- c("SiteID","n")
samp.size <- samp.size.all[samp.size.all$n != 0,]
rm(samp.size.all)

  #full df for all 146 samples: site ID, lat-long, admixture
pops <- unique(coords$SiteID)
Npop = length(unique(coords$SiteID))
Npop
q3.samples <- cbind(coords,qmatrix)
rm(pops)
rm(Npop)
rm(coords)

  #reduced df, aggregated across sites
pop.means <- aggregate(.~SiteID, data=q3.samples, mean)
qpop <- data.matrix(pop.means[,5:ncol(pop.means)]) #just admixtures
coord.pop <- cbind(pop.means[,4],pop.means[,3]) #just lat-long
pop.means <- merge(pop.means, samp.size, by = "SiteID")
rm(samp.size)
pop.means <- pop.means[,c(1, 8, 3:7)]

#check that all pops' ancestry proportions sum to 1
unique(apply(qpop,1,function(x){round(sum(x),2)==1}))

#write file for later plotting
write.csv(x = pop.means, file = "analysis/LEA_analysis/PopMeans.csv", sep = ",", col.names = T, row.names = F)
