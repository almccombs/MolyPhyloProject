mydat <- read.csv("data/TetonSampleSites.csv", header = T)
load("data/sample_names_n146.Rdata")

mydat <- mydat[which(mydat$SpecimenID %in% samples_small & mydat$PopulationCode != ""),]
change_ids <- mydat$SpecimenID[which(mydat$PopulationCode == "Lozier hill meadow (lozier hill east)", arr.ind = T)]
mydat$PopulationCode[which(mydat$SpecimenID %in% change_ids)] <- "Lozier hill"
rm(change_ids)

unique(mydat$PopulationCode)
with(mydat, table(PopulationCode))

