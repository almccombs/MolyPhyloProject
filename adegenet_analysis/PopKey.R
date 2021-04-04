## Read in population key

popkey <- read.csv("analysis/LocDataCleaned.csv", header = T)

popkey$SiteID <- gsub("road", "Road", popkey$SiteID, ignore.case = FALSE)
popkey$SiteID <- gsub("Death Canyon. UTM 514946, 4833685 - on death canyon trail, not in designated meadow", "Death Canyon Trail", popkey$SiteID, ignore.case = FALSE)
popkey$SiteID <- gsub("Death Canyon. UTM 0516038, 4833635 - on death canyon trail, not in a designated meadow", "Death Canyon Trail", popkey$SiteID, ignore.case = FALSE)
popkey$SiteID <- gsub("hill meadow \\(lozier hill east\\)", "Hill Meadow", popkey$SiteID, ignore.case = T)
popkey$SiteID <- gsub("Near Wilderness Rd 3 but not at a real site", "Wilderness Rd 3", popkey$SiteID, ignore.case = T)
popkey$SiteID <-  as.factor(popkey$SiteID)
levels(popkey$SiteID)
