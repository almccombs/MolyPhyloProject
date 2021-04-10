library(geodist)

load("data/pop_data_n146.Rdata")
sites <- mydat[match(unique(mydat$SiteID), mydat$SiteID),]

site_dist <- geodist(x = sites, sequential = F)

save(site_dist, file = "output/site_geodistance_matrix.Rdata")
