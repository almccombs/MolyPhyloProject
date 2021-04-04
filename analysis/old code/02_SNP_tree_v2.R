# https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#principal-component-analysis-pca

library(gdsfmt)
library(SNPRelate)
library(ggplot2)

load("data/pop_data_n146.Rdata")

# Read in SNP dataset

#snpgdsVCF2GDS("data/parnassius_clodius_unfiltered_imputed.vcf", "data/parnassius_clodius_unfiltered_imputed.gds")

snpgdsSummary("data/parnassius_clodius_unfiltered_imputed.gds")

#snpgdsClose(genofile)
#genofile <- snpgdsOpen("data/parnassius_clodius_unfiltered_imputed.gds", readonly = FALSE)
#genofile

# add population data
#id_order <- read.gdsn(index.gdsn(genofile, "sample.id"))
#ordered_pops <- mydat$SiteID[order(match(mydat$SampleID,id_order))]
  # spot-check
#tail(read.gdsn(index.gdsn(genofile, "sample.id")), 10)
#tail(ordered_pops, 10)
  # add in the correct order
#add.gdsn(genofile, "sample.annot", ordered_pops)
#rm(id_order, ordered_pops)
#genofile
#read.gdsn(index.gdsn(genofile, "sample.annot"), start = 1, count = 5)
#snpgdsClose(genofile)


# Open read-only copy
genofile <- snpgdsOpen("data/parnassius_clodius_unfiltered_imputed.gds", readonly = TRUE)
genofile
read.gdsn(index.gdsn(genofile, "sample.annot"), start = 1, count = 5)
read.gdsn(index.gdsn(genofile, "sample.id"), start = 1, count = 5)
read.gdsn(index.gdsn(genofile, "snp.allele"), start = 1, count = 5)
read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count=c(5,10))
# 0 = two B alleles, 1 = one A and one B allele, 2 = two A alleles, 3 = missing
sum(read.gdsn(index.gdsn(genofile, "genotype")) == 0)
sum(read.gdsn(index.gdsn(genofile, "genotype")) == 1)
sum(read.gdsn(index.gdsn(genofile, "genotype")) == 2)
sum(read.gdsn(index.gdsn(genofile, "genotype")) == 3)

pop_code <- read.gdsn(index.gdsn(genofile, path = "sample.annot"))
length(unique(pop_code))
table(pop_code)

sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# # LD-based SNP pruning
# set.seed(926)
# snpset <- snpgdsLDpruning(genofile,ld.threshold = 0.2, autosome.only = F)
# snp_id <- unlist(snpset)
# str(snpset)
# snpset_id <- unlist(unname(snpset))


# Run PCA
pca_run <- snpgdsPCA(genofile, snp.id = snpset_id, num.thread=2, autosome.only = F)
pca_run
pca_run$varprop
head(round(pc_var, 2))
pca_df <- data.frame(sample.id = pca_run$sample.id,
                     pop = pop_code,
                     EV1 = pca_run$eigenvect[,1],    # the first eigenvector
                     EV2 = pca_run$eigenvect[,2],    # the second eigenvector
                     stringsAsFactors = FALSE)
head(pca_df)

ggplot(pca_df, aes(x = EV1, y = EV2, color = pop)) + geom_point()



# IBD analysis for Pilgrim Creek
pc_id <- sample_id[pop_code == "Pilgrim Creek"]

  # MoM
ibd_mom <- snpgdsIBDMoM(genofile, sample.id = pc_id, autosome.only = FALSE)
ibd_mom_coef <- snpgdsIBDSelection(ibd_mom)
head(ibd_mom_coef)
  # The coefficient of kinship between two individuals B and C is the probability that homologous genese on gametes segregating from B and C are ibd.
  # k0 is the probability that B and C have exactly 0 alleles IBD
  # k1 is the probability that B and C have exactly 1 allele IBD
ggplot(ibd_mom_coef, aes(x = k0, y = k1)) + geom_point()

    # MLE
ibd_mle <- snpgdsIBDMLE(genofile, sample.id=pc_id, snp.id=snp_id,
                    num.thread=3, autosome.only = FALSE)
ibd_mle_coef <- snpgdsIBDSelection(ibd_mle)
head(ibd_mle_coef)

ggplot(ibd_mle_coef, aes(x = k0, y = k1)) + geom_point()


# IBS analysis
snpgdsIBSNum(genofile, autosome.only = FALSE)
# e.g., IBS0 is an nxn matrix that reports the number of sites (for a pair of samples) that both have 0.  E.g., IBS0[1,146] = 20, which means that sample 1 and sample 146 have twenty 0's at the same site (same SNP)

ibs <- snpgdsIBS(genofile, autosome.only = FALSE)

order_pop_code <- order(pop_code)
image(ibs$ibs[order_pop_code, order_pop_code], col = terrain.colors(16))


# multidimensional scaling
loc <- cmdscale(1-ibs$ibs, k = 2)
x <- loc[,1]; y <- loc[,2]
population <- as.factor(pop_code)

ggplot(data = NULL, aes(x = x, y = y, color = population)) + geom_point()

# Cluster and draw tree
set.seed(926)
ibs_hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread = 2, autosome.only = FALSE))

cut_tree <- snpgdsCutTree(ibs_hc, verbose=TRUE)
cut_tree_pop <- snpgdsCutTree(ibs_hc, samp.group = as.factor(pop_code), verbose=TRUE)

plot(cut_tree_pop$dendrogram, leaflab = "none")

snpgdsDrawTree(cut_tree_pop,
               y.label.kinship=T,
               leaflab="perpendicular")


# Estimate IBD coefficients
set.seed(100)
ibd <- snpgdsIBDMLE(genofile, num.thread=2, autosome.only = F)
ibd_coeff <- snpgdsIBDSelection(ibd)
plot(ibd_coeff$k0, ibd_coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MLE)")
lines(c(0,1), c(1,0), col="red", lty=2)
