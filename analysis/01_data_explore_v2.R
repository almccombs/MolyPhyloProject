mydat <- read.csv("data/TetonSampleSites_working.csv", header = T)
load("data/sample_names_n146.Rdata")

#mydat <- mydat[which(mydat$SampleID %in% samples_small & mydat$PopulationCode != ""),]
mydat <- mydat[which(mydat$SampleID %in% samples_small),]
change_ids <- mydat$SampleID[which(mydat$SiteID == "Lozier hill meadow (lozier hill east)", arr.ind = T)]
mydat$SiteID[which(mydat$SampleID %in% change_ids)] <- "Lozier hill"
rm(change_ids)

unique(mydat$SiteID)
sample_counts <- as.data.frame(with(mydat, table(SiteID)))
sample_counts[order(sample_counts$Freq),]

# Remove sites with only 1 or 2 samples
mydat <- subset(mydat, SiteID != "Bear paw lake intersection" & SiteID != "Lozier Road" & SiteID != "Buffalo fork")
sample_counts <- as.data.frame(with(mydat, table(SiteID)))
sample_counts[order(sample_counts$Freq),]

save(mydat, file = "data/pop_data_n142.Rdata")
