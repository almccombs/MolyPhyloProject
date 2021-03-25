library(sp)
library(RgoogleMaps)
library(mapplots)  #for draw.pie function

op <- par()

pop.means <- read.csv("analysis/LEA_analysis/PopMeans.csv", header = T)
pies <- as.matrix(pop.means[,c("V1","V2","V3")])
n <- as.vector(pop.means[,2])

coords.sp <- pop.means
coordinates(coords.sp) <- c("Longitude", "Latitude")
proj4string(coords.sp) <- CRS("+proj=longlat +datum=WGS84")


#Make smaller submaps interactively
G <- select.spatial(coords.sp, digitize = F)

group1 <- pop.means[G,]
group1

xminG1 <- min(group1$Longitude)
xmaxG1 <- max(group1$Longitude)
xmidG1 <- (xminG1 + xmaxG1)/2
yminG1 <- min(group1$Latitude)
ymaxG1 <- max(group1$Latitude)
ymidG1 <- (yminG1 + ymaxG1)/2

coords.spG1 <- group1
coordinates(coords.spG1) <- c("Longitude", "Latitude")
proj4string(coords.spG1) <- CRS("+proj=longlat +datum=WGS84")

group1.map <- GetMap(center = c(ymidG1, xmidG1), zoom = 12, maptype = "terrain")

par(op)
googG1 <- LatLon2XY.centered(group1.map, coords.spG1@coords[,2], coords.spG1@coords[,1], zoom = 12)
PlotOnStaticMap(group1.map)
draw.pie(z = pies[G,], x = googG1$newX, y = googG1$newY, radius = (sqrt(n/pi)*10), labels = "")

