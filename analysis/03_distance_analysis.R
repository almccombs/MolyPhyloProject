library(geodist)
library(ape)
library(ggplot2)
library(gridExtra)

load("data/pop_data_n146.Rdata")
load("output/nei_site_order.Rdata")
sites <- mydat[match(unique(mydat$SiteID), mydat$SiteID),]
sites <- sites[match(nei_site_order, sites$SiteID),]
sites$SiteID == nei_site_order
site_dist <- geodist(x = sites, sequential = F, measure = "geodesic")
save(site_dist, file = "output/site_geodistance_matrix.Rdata")
rm(nei_site_order, sites)

load("output/nei_popdist_matrix.Rdata")
load("output/reynolds_popdist_matrix.Rdata")

nei_dist <- pop_dist_nei_matrix[lower.tri(pop_dist_nei_matrix)]
reynolds_dist <- pop_dist_reynolds_matrix[lower.tri(pop_dist_reynolds_matrix)]
geo_dist <- site_dist[lower.tri(site_dist)]

# Mantel Null: The distances among objects in a matrix of response variables are not linearly correlated with another matrix of explanatory variables.

  # Nei to geo
nei_to_geo <- ggplot(data = NULL, aes(x = nei_dist, y = geo_dist)) +
  geom_point() + xlab("Nei's distance") + ylab("geodesic distance")
nei_to_geo
mantel.test(m1 = pop_dist_nei_matrix, m2 = site_dist, nperm = 9999)

  # Reynolds to geo
reynolds_to_geo <- ggplot(data = NULL,
                          aes(x = reynolds_dist, y = geo_dist)) +
  geom_point() + xlab("Reynolds' distance") + ylab("geodesic distance")
reynolds_to_geo
mantel.test(m1 = pop_dist_reynolds_matrix, m2 = site_dist, nperm = 9999)
  
  # Nei to Reynolds
nei_to_reynolds <- ggplot(data = NULL,
                          aes(x = nei_dist, y = reynolds_dist)) +
  geom_point() + xlab("Nei's distance") + ylab("Reynolds' distance")
nei_to_reynolds
mantel.test(m1 = pop_dist_nei_matrix, m2 = pop_dist_reynolds_matrix, nperm = 9999)

# Save figure
png(filename = "figures/distance_comparisons.png", width = 7.5, height = 7.5, units = "in", res = 300)
grid.arrange(nei_to_geo, reynolds_to_geo, nei_to_reynolds, nrow = 2)
dev.off()

summary(lm(reynolds_dist ~ nei_dist))
  # Adj R^2 = 0.8832