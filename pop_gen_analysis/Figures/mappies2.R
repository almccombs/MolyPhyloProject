library(sp)
library(RgoogleMaps)
library(mapplots)  #for draw.pie function

setwd("H:/Iowa State University/Debinski Lab/Parnassius genetics/ParnassiusGenetics")

op <- par()

# Set up data structures

pop.means <- read.csv("analysis/LEA_analysis/PopMeans.csv", header = T)

xmin <- min(pop.means$Longitude)
xmax <- max(pop.means$Longitude)
xmid <- (xmin + xmax)/2
ymin <- min(pop.means$Latitude)
ymax <- max(pop.means$Latitude)
ymid <- (ymin + ymax)/2

pies <- as.matrix(pop.means[,c("V1","V2","V3")])
coords <- as.matrix(pop.means[,c("Longitude","Latitude")])
n <- as.vector(pop.means[,2])

coords.sp <- pop.means
coordinates(coords.sp) <- c("Longitude", "Latitude")
proj4string(coords.sp) <- CRS("+proj=longlat +datum=WGS84")

#Plot the sample sites
gtnp <- GetMap(center = c(ymid, xmid), zoom = 10, maptype = "terrain")
sites.map <- PlotOnStaticMap(gtnp, coords.sp@coords[,2], coords.sp@coords[,1], pch = 19, col = 4)
TextOnStaticMap(sites.map, lat = 44.08, lon = -110.39, labels = "Sample Populations", add = T)

#Plot the pies using the mapplots package
par(op)
goog2 <- LatLon2XY.centered(gtnp, coords.sp@coords[,2], coords.sp@coords[,1], zoom = 10)
gtnp.map <- PlotOnStaticMap(gtnp)
TextOnStaticMap(gtnp.map, lat = 44.06, lon = -110.39, labels = "Genetic admixture for each population\nArea of pie proportional to\nsample size at each site", add = T)
draw.pie(z = pies, x = goog2$newX, y = goog2$newY, radius = (sqrt(n/pi)*10), labels = "")

#Make smaller submaps
  #Northern sites
ns <- c(1,2,3,8,12,13,14,18,26,27,28,29) #indicies for northern sites
north.sites <- pop.means[ns,]
n.ns <- north.sites[,2]

xmin.ns <- min(north.sites$Longitude)
xmax.ns <- max(north.sites$Longitude)
xmid.ns <- (xmin.ns + xmax.ns)/2
ymin.ns <- min(north.sites$Latitude)
ymax.ns <- max(north.sites$Latitude)
ymid.ns <- (ymin.ns + ymax.ns)/2

coords.sp.ns <- north.sites
coordinates(coords.sp.ns) <- c("Longitude", "Latitude")
proj4string(coords.sp.ns) <- CRS("+proj=longlat +datum=WGS84")

north.sites.map <- GetMap(center = c(ymid.ns, xmid.ns), zoom = 12, maptype = "terrain")

par(op)
goog.ns <- LatLon2XY.centered(north.sites.map, coords.sp.ns@coords[,2], coords.sp.ns@coords[,1], zoom = 12)
ns.map <- PlotOnStaticMap(north.sites.map)
TextOnStaticMap(gtnp.map, lat = 44.06, lon = -110.39, labels = "Genetic admixture for northern sites", add = T)
draw.pie(z = pies[ns,], x = goog.ns$newX, y = goog.ns$newY, radius = (sqrt(n.ns/pi)*10), labels = "")

  #Western sites
ws <- c(4,5,7,16,17,21,22,23,24,25) #indicies for western sites
west.sites <- pop.means[ws,]
n.ws <- west.sites[,2]

xmin.ws <- min(west.sites$Longitude)
xmax.ws <- max(west.sites$Longitude)
xmid.ws <- (xmin.ws + xmax.ws)/2
ymin.ws <- min(west.sites$Latitude)
ymax.ws <- max(west.sites$Latitude)
ymid.ws <- (ymin.ws + ymax.ws)/2

coords.sp.ws <- west.sites
coordinates(coords.sp.ws) <- c("Longitude", "Latitude")
proj4string(coords.sp.ws) <- CRS("+proj=longlat +datum=WGS84")

west.sites.map <- GetMap(center = c(ymid.ws, xmid.ws), zoom = 12, maptype = "terrain")

par(op)
goog.ws <- LatLon2XY.centered(west.sites.map, coords.sp.ws@coords[,2], coords.sp.ws@coords[,1], zoom = 12)
PlotOnStaticMap(west.sites.map)
TextOnStaticMap(gtnp.map, lat = 44.06, lon = -110.39, labels = "Genetic admixture for western sites", add = T)
draw.pie(z = pies[ws,], x = goog.ws$newX, y = goog.ws$newY, radius = (sqrt(n.ws/pi)*10), labels = "")


#Bubble plots for the three genetic groups
radV1 <- sqrt(pop.means$V1/pi)
radV2 <- sqrt(pop.means$V2/pi)
radV3 <- sqrt(pop.means$V3/pi)

par(mfrow = c(2,2))

PlotOnStaticMap(gtnp)
goog2 <- LatLon2XY.centered(gtnp, coords.sp@coords[,2], coords.sp@coords[,1], zoom = 10)
draw.pie(z = pies, x = goog2$newX, y = goog2$newY, radius = (sqrt(n/pi)*10), labels = "")

symbols(x = coords[,1], y = coords[,2], circles = radV1, inches = 0.1, bg = "brown", main = "\nAdmixture 1", asp = 1)
symbols(x = coords[,1], y = coords[,2], circles = radV2, inches = 0.1, bg = "red", main = "\nAdmixture 2", asp = 1)
symbols(x = coords[,1], y = coords[,2], circles = radV3, inches = 0.1, bg = "green", main = "\nAdmixture 3", asp = 1)

coords.sp.bp <- coords.sp
coords.sp.bp@data$V1 <- coords.sp@data$V1 * 100
coords.sp.bp@data$V2 <- coords.sp@data$V2 * 100
coords.sp.bp@data$V3 <- coords.sp@data$V3 * 100

par(op)
bubbleMap(coords.sp.bp, coords = c(goog2$newX, goog2$newY), map = gtnp, zcol = "V1", col = "brown")
bubbleMap(coords.sp.bp, coords = c(goog2$newX, goog2$newY), map = gtnp, zcol = 'V2', col = "red")
bubbleMap(coords.sp.bp, coords = c(goog2$newX, goog2$newY), map = gtnp, zcol = "V3", col = "darkgreen", do.sqrt = T)

