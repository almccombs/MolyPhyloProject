library(sp)
library(RgoogleMaps)
library(mapplots)  #for draw.pie function

setwd("D:/Iowa State University/Debinski Lab/Parnassius genetics/ParnassiusGenetics")

op <- par()

# Set up data structures

pop.means <- read.csv("analysis/LEA_analysis/PopMeans.csv", header = T)

xmin <- min(pop.means$Longitude)
xmax <- max(pop.means$Longitude)
xmid <- (xmin + xmax)/2
ymin <- min(pop.means$Latitude)
ymax <- max(pop.means$Latitude)
ymid <- (ymin + ymax)/2

coords.sp <- pop.means
coordinates(coords.sp) <- c("Longitude", "Latitude")
proj4string(coords.sp) <- CRS("+proj=longlat +datum=WGS84")

#Plot the sample sites
gtnp <- GetMap(center = c(ymid, xmid), zoom = 10, maptype = "terrain")
PlotOnStaticMap(gtnp, coords.sp@coords[,2], coords.sp@coords[,1], pch = 19, col = 4)

#Add dicentra sites
load("analysis/gstudio/objects/search_sites.Rdata")
dicentra <- subset(pop.means, SiteID %in% search.sites, select = c(SiteID, Latitude, Longitude))
dicentra
dicentra.sp <- dicentra
coordinates(dicentra.sp) <- c("Longitude", "Latitude")
proj4string(dicentra.sp) <- CRS("+proj=longlat +datum=WGS84")

PlotOnStaticMap(gtnp, dicentra.sp@coords[,2], dicentra.sp@coords[,1], pch = 19, col = "red", add = TRUE)
# saved this plot as DicentraSites.png in the Maps folder