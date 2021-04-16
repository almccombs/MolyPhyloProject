#https://markravinet.github.io/Chapter9.html
# file:///C:/Users/AUDREY~1/AppData/Local/Temp/practical-introphylo.1.0.pdf

library(adegenet)
library(pegas)
library(ggplot2)
library(gridExtra)

load("data/pop_data_n142.Rdata")

snpdata <- read.vcf("data/parnassius_clodius_unfiltered_imputed.vcf")
snpdata_gen <- df2genind(snpdata, ploidy = 2, sep = "/")
#summary(snpdata_gen)
#snpdata_df <- genind2df(snpdata_gen)
#snpdata_matrix <- as.matrix(snpdata_df)

# remove 4 individuals from sites with only 1 or 2 samples
rm_ids <- rownames(snpdata)[!rownames(snpdata) %in% mydat$SampleID]
snpdata <- snpdata[which(!rownames(snpdata) %in% rm_ids),]
snpdata_gen <- snpdata_gen[!rownames(snpdata_gen@tab) %in% rm_ids]

# convert to genpop object including population information
id_order <- rownames(snpdata)
ordered_pops <- mydat$SiteID[order(match(mydat$SampleID,id_order))]
snpdata_genpop <- genind2genpop(x = snpdata_gen, pop = ordered_pops, quiet = FALSE)
rm(id_order, ordered_pops, rm_ids)


#### Nei's distance (drift + mutation)
pop_dist_nei <- dist.genpop(snpdata_genpop, method = 1, diag = T, upper = T)
  # method = 1 is nei's distance
pop_dist_nei_df <- as.data.frame(as.matrix(pop_dist_nei, labels = T))
pop_dist_nei_matrix <- as.matrix(pop_dist_nei, labels = T)

  # save distance matrix for later analysis
save(pop_dist_nei_matrix, file = "output/nei_popdist_matrix.Rdata")
  # save site names to match order later
nei_site_order <- colnames(pop_dist_nei_matrix)
save(nei_site_order, file = "output/nei_site_order.Rdata")
rm(nei_site_order)


  # Plot a heat map
heatmap_df <- reshape2::melt(pop_dist_nei_df)
names(heatmap_df) <- c("pop1", "dist")
heatmap_df$pop2 <- rep(unique(heatmap_df$pop1), 26)
heatmap_df <- heatmap_df[,c("pop1", "pop2", "dist")]

heatmap_nei <- ggplot(heatmap_df, aes(x = pop1, y = pop2, fill = dist)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "#008080") +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Population distances, Nei's") +
  labs(fill = "Distance")
heatmap_nei
save(heatmap_nei, file = "output/nei_heatmap_plot_object.Rdata")


png(filename = "output/nei_heatmap_plot.png", width = 7.5, height = 7.5, units = "in", res = 300)
heatmap_nei
dev.off()
rm(heatmap_df, pop_dist_nei_df)

  # Make a NJ tree
nj_tree <- nj(pop_dist_nei)
plot(nj_tree, type = "unrooted", use.edge.length = T, cex = .5)

  # Use cophenetic to check
x <- as.vector(pop_dist_nei)
y_nj <- as.vector(as.dist(cophenetic(nj_tree)))
plot_df <- as.data.frame(cbind(x, y_nj))
coph_nj_nei <- ggplot(data = plot_df, aes(x, y_nj)) + geom_point() +
  xlab("population distance") + ylab("tree distance") +
  labs(title = "Cophenetic Plot: Nei's Distance",
       subtitle = "Neighbor-Joining Tree")
  # This plots the population distance against the distance given by the dendrogram
coph_nj_nei
save(coph_nj_nei, file = "output/nei_nj_cophenetic_plot_object.Rdata")

  # Make a UPGMA tree
upgma_tree <- as.phylo(hclust(pop_dist_nei, method = "average"))
plot(upgma_tree)
  # write in Newick format
write.tree(upgma_tree, file = "output/nei_upgma_tree_newick.Rdata")
write.tree(upgma_tree, file = "output/nei_upgma_tree_newick.txt")

  # Use cophenetic to check
y_upgma <- as.vector(as.dist(cophenetic(upgma_tree)))
plot_df <- as.data.frame(cbind(x, y_upgma))
coph_upgma_nei <- ggplot(data = plot_df, aes(x, y_upgma)) + geom_point() +
  xlab("population distance") + ylab("tree distance") +
  labs(title = "Cophenetic Plot: Nei's Distance",
       subtitle = "UPGMA Tree")
  # This plots the population distance against the distance given by the dendrogram
coph_upgma_nei
save(coph_upgma_nei, file = "output/nei_upgma_cophenetic_plot_object.Rdata")

grid.arrange(coph_nj_nei, coph_upgma_nei, ncol = 2)

  # Bootstrap
boot_nj <- boot.phylo(nj_tree, pop_dist_nei_matrix, FUN = nj, B = 1000, rooted = F)
boot_nj
nj_tree$node.label <- boot_nj

# write in Newick format, with bootstrap support
write.tree(nj_tree, file = "output/nei_nj_tree_newick.txt")
write.tree(nj_tree, file = "output/nei_nj_tree_newick.Rdata")

rm(plot_df, coph_nj_nei, coph_upgma_nei, heatmap_nei, nj_tree, pop_dist_nei_matrix, upgma_tree, boot_nj, pop_dist_nei, x, y_nj, y_upgma)



### Reynolds' distance (drift only)
  # https://rdrr.io/cran/adegenet/man/dist.genpop.html

pop_dist_reynolds <- dist.genpop(snpdata_genpop, method = 3, diag = T, upper = T)
  # method = 3 is Reynolds's distance
pop_dist_reynolds_df <- as.data.frame(as.matrix(pop_dist_reynolds))
pop_dist_reynolds_matrix <- as.matrix(pop_dist_reynolds)

  # save distance matrix for later analysis
save(pop_dist_reynolds_matrix, file = "output/reynolds_popdist_matrix.Rdata")
  # save site names to match order later
reynolds_site_order <- colnames(pop_dist_reynolds_matrix)
save(reynolds_site_order, file = "output/reynolds_site_order.Rdata")
rm(reynolds_site_order) # the order of sites is the same for the Nei and Reynolds dist matrices

  # Plot a heat map
heatmap_df <- reshape2::melt(pop_dist_reynolds_df)
names(heatmap_df) <- c("pop1", "dist")
heatmap_df$pop2 <- rep(unique(heatmap_df$pop1), 26)
heatmap_df <- heatmap_df[,c("pop1", "pop2", "dist")]

heatmap_reynolds <- ggplot(heatmap_df, aes(x = pop1, y = pop2, fill = dist)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "#008080") +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Population distances, Reynolds'") +
  labs(fill = "Distance")
heatmap_reynolds
save(heatmap_reynolds, file = "output/reynolds_heatmap_plot_object.Rdata")


png(filename = "output/reynolds_heatmap_plot.png", width = 7.5, height = 7.5, units = "in", res = 300)
heatmap_reynolds
dev.off()
rm(heatmap_df, pop_dist_reynolds_df)

  # Make a NJ tree
nj_tree <- nj(pop_dist_reynolds)
plot(nj_tree, type = "unrooted", use.edge.length = T, cex = .5)

# Use cophenetic to check
x <- as.vector(pop_dist_reynolds)
y_nj <- as.vector(as.dist(cophenetic(nj_tree)))
plot_df <- as.data.frame(cbind(x, y_nj))
coph_nj_reynolds <- ggplot(data = plot_df,
                           aes(x, y_nj)) + geom_point() +
  xlab("population distance") + ylab("tree distance") +
  labs(title = "Cophenetic Plot: Reynolds' Distance",
       subtitle = "Neighbor-Joining Tree")
  # This plots the population distance against the distance given by the dendrogram
coph_nj_reynolds
save(coph_nj_reynolds, file = "output/reynolds_nj_cophenetic_plot_object.Rdata")

  # Make a UPGMA tree
upgma_tree <- as.phylo(hclust(pop_dist_reynolds, method = "average"))
plot(upgma_tree)
  # write in Newick format
write.tree(upgma_tree, file = "output/reynolds_upgma_tree_newick.Rdata")
write.tree(upgma_tree, file = "output/reynolds_upgma_tree_newick.txt")

# Use cophenetic to check
y_upgma <- as.vector(as.dist(cophenetic(upgma_tree)))
plot_df <- as.data.frame(cbind(x, y_upgma))
coph_upgma_reynolds <- ggplot(data = plot_df,
                              aes(x, y_upgma)) +
  geom_point() +
  xlab("population distance") + ylab("tree distance") +
  labs(title = "Cophenetic Plot: Reynolds' Distance",
       subtitle = "UPGMA Tree")
# This plots the population distance against the distance given by the dendrogram
coph_upgma_reynolds
save(coph_upgma_reynolds, file = "output/reynolds_upgma_cophenetic_plot_object.Rdata")

grid.arrange(coph_nj_reynolds, coph_upgma_reynolds, ncol = 2)

  # Bootstrap
boot_nj <- boot.phylo(nj_tree, pop_dist_reynolds_matrix, FUN = nj, B = 1000, rooted = F)
boot_nj
nj_tree$node.label <- boot_nj
# write in Newick format
write.tree(nj_tree, file = "output/reynolds_nj_tree_newick.Rdata")
write.tree(nj_tree, file = "output/reynolds_nj_tree_newick.txt")

rm(plot_df, coph_nj_reynolds, coph_upgma_reynolds, heatmap_reynolds, nj_tree, pop_dist_reynolds_matrix, upgma_tree, boot_nj, pop_dist_reynolds, x, y_nj, y_upgma)

