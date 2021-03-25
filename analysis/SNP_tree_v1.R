# https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#principal-component-analysis-pca

library(gdsfmt)
library(SNPRelate)


#biallelic by default
snpgdsVCF2GDS("data/parnassius_clodius_unfiltered_imputed.vcf", "data/parnassius_clodius_unfiltered_imputed.gds")
snpgdsSummary("data/parnassius_clodius_unfiltered_imputed.gds")
genofile <- snpgdsOpen("data/parnassius_clodius_unfiltered_imputed.gds")

#LD based SNP pruning
set.seed(1000)
snpset <- snpgdsLDpruning(genofile,ld.threshold = 1, autosome.only = F)
snp.id=unlist(snpset)

# Run PCA
pca_run <- snpgdsPCA(genofile, num.thread=2, autosome.only = F)
pc_var <- pca_run$varprop
head(round(pc_var, 2))
tab <- data.frame(sample.id = pca_run$sample.id,
                  EV1 = pca_run$eigenvect[,1],    # the first eigenvector
                  EV2 = pca_run$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

# distance matrix - use IBS
dissMatrix  =  snpgdsIBS(genofile, autosome.only = F)
snpgdsClose(genofile)

# Cluster and draw tree
snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)

cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL, 
                        col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, label.Z=TRUE, 
                        verbose=TRUE)

snpgdsDrawTree(cutTree, main = "Dataset 1",edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
               y.label.kinship=T,leaflab="perpendicular")


# Estimate IBD coefficients
set.seed(100)
ibd <- snpgdsIBDMLE(genofile, num.thread=2, autosome.only = F)
ibd_coeff <- snpgdsIBDSelection(ibd)
plot(ibd_coeff$k0, ibd_coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MLE)")
lines(c(0,1), c(1,0), col="red", lty=2)
