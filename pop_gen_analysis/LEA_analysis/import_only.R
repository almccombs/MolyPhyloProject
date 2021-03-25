#Use read.vcf.markerstats3 for STACKS v1.44

#setwd("D:/Iowa State University/Debinski Lab/Parnassius genetics/ParnassiusGenetics")

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
geno <- read.vcf.genotypes(filename='SNPdata/parnassius_clodius_unfiltered_imputed.vcf',1001)
geno <- t(geno)
