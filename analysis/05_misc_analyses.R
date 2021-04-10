load("output/nei_popdist_matrix.Rdata")
load("output/reynolds_popdist_matrix.Rdata")

nei_avg <- apply(pop_dist_nei_matrix, 1, mean)
reynolds_avg <- apply(pop_dist_reynolds_matrix, 1, mean)

avg_df <- as.data.frame(cbind(nei_avg, reynolds_avg))
vegan::kendall.global(avg_df, nperm = 9999)
