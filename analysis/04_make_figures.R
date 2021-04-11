library(ggplot2)
library(gridExtra)
library(ape)
library(phytools)


# Heatmap figure
load("output/nei_heatmap_plot_object.Rdata")
load("output/reynolds_heatmap_plot_object.Rdata")

heatmap_nei <- heatmap_nei + guides(fill = FALSE) +
  ggtitle("Nei's") + theme(aspect.ratio = 1)
heatmap_reynolds <- heatmap_reynolds +
  theme(axis.text.y = element_blank(), aspect.ratio = 1) +
  ggtitle("Reynolds'")

png(filename = "figures/heatmaps.png", width = 15, height = 7.5, units = "in", res = 300)
grid.arrange(heatmap_nei, heatmap_reynolds, nrow = 1)
dev.off()

rm(heatmap_nei, heatmap_reynolds)

# Cophenetic figure
load("output/nei_nj_cophenetic_plot_object.Rdata")
load("output/nei_upgma_cophenetic_plot_object.Rdata")
load("output/reynolds_nj_cophenetic_plot_object.Rdata")
load("output/reynolds_upgma_cophenetic_plot_object.Rdata")

coph_nj_reynolds <- coph_nj_reynolds + labs(title = "Reynolds' Distance",
                                  subtitle = "Neighbor-Joining Tree")
coph_upgma_nei <- coph_upgma_nei + labs(title = "", subtitle = "UPGMA Tree") +
  ylab("")
coph_upgma_reynolds <- coph_upgma_reynolds + labs(title = "", subtitle = "UPGMA Tree") +
  ylab("")

png(filename = "figures/cophenetics.png", width = 7.5, height = 7.5, units = "in", res = 300)
grid.arrange(coph_nj_nei, coph_upgma_nei, coph_nj_reynolds, coph_upgma_reynolds, nrow = 2)
dev.off()

rm(coph_nj_nei, coph_nj_reynolds, coph_upgma_nei, coph_upgma_reynolds)

# Trees
nei_nj <- read.tree("output/nei_nj_tree_newick.Rdata")
nei_upgma <- read.tree("output/nei_upgma_tree_newick.Rdata")
reynolds_nj <- read.tree("output/reynolds_nj_tree_newick.Rdata")
reynolds_upgma <- read.tree("output/reynolds_upgma_tree_newick.Rdata")

png(filename = "figures/Rtrees/nei_nj.png", width = 10, height = 7.5, units = "in", res = 300)
plot(nei_nj, type = "unrooted", cex = .75, use.edge.length = F, lab4ut = "axial", no.margin =T)
dev.off()

cairo_ps(filename = "figures/Rtrees/nei_nj.eps", width = 10, height = 7.5)
plot(nei_nj, type = "unrooted", cex = .75, use.edge.length = F, lab4ut = "axial", no.margin =T)
dev.off()

png(filename = "figures/Rtrees/nei_upgma.png", width = 10, height = 7.5, units = "in", res = 300)
plot(nei_upgma, type = "unrooted", cex = .75, use.edge.length = F, lab4ut = "axial", no.margin =T)
dev.off()

png(filename = "figures/Rtrees/reynolds_nj.png", width = 10, height = 7.5, units = "in", res = 300)
plot(reynolds_nj, type = "unrooted", cex = .75, use.edge.length = F, lab4ut = "axial", no.margin =T)
dev.off()

cairo_ps(filename = "figures/Rtrees/reynolds_nj.eps", width = 10, height = 7.5)
plot(reynolds_nj, type = "unrooted", cex = .75, use.edge.length = F, lab4ut = "axial", no.margin =T)
dev.off()


png(filename = "figures/Rtrees/reynolds_upgma.png", width = 10, height = 7.5, units = "in", res = 300)
plot(reynolds_upgma, type = "unrooted", cex = .75, use.edge.length = F, lab4ut = "axial", no.margin =T)
dev.off()


nei_both <- cophylo(nei_nj, nei_upgma)
plot(nei_both, link.type = "curved", link.lty = "solid", link.col = make.transparent("blue", 0.25), fsize = 0.8)

nj_both <- cophylo(nei_nj, reynolds_nj)
plot(nj_both, link.type = "curved", link.lty = "solid", link.col = make.transparent("blue", 0.25), fsize = 0.8)
