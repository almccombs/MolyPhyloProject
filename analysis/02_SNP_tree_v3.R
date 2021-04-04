#https://markravinet.github.io/Chapter9.html
# file:///C:/Users/AUDREY~1/AppData/Local/Temp/practical-introphylo.1.0.pdf

library(adegenet)
library(pegas)
library(ggplot2)

load("data/pop_data_n146.Rdata")

snpdata <- read.vcf("data/parnassius_clodius_unfiltered_imputed.vcf")
snpdata_gen <- df2genind(snpdata, ploidy = 2, sep = "/")
#summary(snpdata_gen)
snpdata_df <- genind2df(snpdata_gen)
#snpdata_matrix <- as.matrix(snpdata_df)

# convert to genpop object including population information
id_order <- rownames(snpdata)
ordered_pops <- mydat$SiteID[order(match(mydat$SampleID,id_order))]
snpdata_genpop <- genind2genpop(x = snpdata_gen, pop = ordered_pops, quiet = FALSE)
snpdata_df <- as.data.frame(cbind(ordered_pops, snpdata_df))
rm(id_order, ordered_pops)

# calculate population distances
pop_dist <- dist.genpop(snpdata_genpop, method = 1, diag = T, upper = T)
  # method = 1 is nei's distance
pop_dist_df <- as.data.frame(as.matrix(pop_dist))
pop_dist_matrix <- as.matrix(pop_dist)

  # Plot a heat map
heatmap_df <- reshape2::melt(pop_dist_df)
names(heatmap_df) <- c("pop1", "dist")
heatmap_df$pop2 <- rep(unique(heatmap_df$pop1), 29)
heatmap_df <- heatmap_df[,c("pop1", "pop2", "dist")]
ggplot(heatmap_df, aes(x = pop1, y = pop2, fill = dist)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "#008080") +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90))
rm(heatmap_df, pop_dist_df)

# Make a NJ tree
nj_tree <- nj(pop_dist)
png(filename = "output/nj_tree_unrooted.png", width = 10, height = 7.5, units = "in", res = 300)
plot(nj_tree, type = "unrooted", use.edge.length = T, cex = .5)
dev.off()
  # write in Newick format
write.tree(nj_tree, file = "data/nj_tree.txt")

  # Use cophenetic to check
x <- as.vector(pop_dist)
y <- as.vector(as.dist(cophenetic(nj_tree)))
qplot(x, y)


# Make a UPGMA tree
upgma_tree <- as.phylo(hclust(pop_dist, method = "average"))
png(filename = "output/upgma_tree.png", width = 10, height = 7.5, units = "in", res = 300)
plot(upgma_tree)
dev.off()
  # write in Newick format
write.tree(upgma_tree, file = "data/nupgma_tree.txt")

  # Use cophenetic to check
y <- as.vector(as.dist(cophenetic(upgma_tree)))
qplot(x, y)
rm(x,y)


## Bootstrap
boot_nj <- boot.phylo(nj_tree, pop_dist_matrix, FUN = nj, B = 1000, rooted = F)
boot_nj
write.table(boot_nj, file = "output/nj_tree_boot_support.txt")

rm(upgma_tree, boot_nj)

