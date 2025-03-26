## Load Libraries --------------
<<<<<<< HEAD
=======
library(rgeos)
library(raster)
library(ggmap)
library(rgdal)
library(RColorBrewer)
library(broom)
library(maps) # tool for maps
library(mapdata) # all your basemaps are here
library(marmap) # for bathymetry if needed
library(mapplots) # for add.pie
library(gplots) # for colour range
library(rworldmap)
library(lattice)
library(rnaturalearth)
library(sf)
library(sp)
library(dplyr)


##Snow Crab Locations Map#####
#source("PlotBathy2.R")

pops<-read.csv("data/Pop_Coordinates_alt.csv",header = T)

latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"

=======
library(maptools)
library(lattice)

setwd("R:/Science/CESD/HES_MPAGroup/Projects/Snow Crab Genetics/")
##Snow Crab Locations Map#####
source("PlotBathy2.R")

pops<-read.csv("data/Pop_Coordinates.csv",header = T)
>>>>>>> e9f7d1387d6708209891eee2af834d4ea8df08d0
#Get colour palettes
blues <- colorRampPalette(c("black","darkblue","darkslateblue","cyan2"))(1000)
greys <- colorRampPalette(c(grey(0.4),grey(0.99)))
blues<-colorRampPalette(colors=c("midnightblue","navy","darkblue","blue","cornflowerblue","lightslateblue","lightskyblue","cyan"),bias=0.5)(500)

#OR
ocean.pal <-   c("#000000", "#000413", "#000728", "#002650", "#005E8C",
                 "#0096C8", "#45BCBB", "#8AE2AE", "#BCF8B9", "#DBFBDC")


land.pal <-  c("#467832", "#887438", "#B19D48", "#DBC758", "#FAE769",
               "#FAEB7E", "#FCED93", "#FCF1A7", "#FCF6C1", "#FDFAE0")
#Add NAFO divisions
<<<<<<< HEAD
# NAFO<-shapefile("R:/Science ERD/HEDMG/MPA_Group/Projects/Fundian Channel/data/NAFO_Divisions_Shapefiles/Divisions.shp")
# NAFO= spTransform(NAFO, CRS("+proj=longlat +datum=WGS84"))
# 
# 
# MARCMAS<-shapefile("R:/Science ERD/HEDMG/MPA_Group/Projects/Snow Crab Genetics/data/CrabManagementAreas.shp")
# MARCMAS= spTransform(MARCMAS, CRS("+proj=longlat +datum=WGS84"))
# 
# NLCMAS<-shapefile("R:/Science ERD/HEDMG/MPA_Group/Projects/Snow Crab Genetics/data/NL_CMAs/CMA.shp")
# NLCMAS= spTransform(NLCMAS, CRS("+proj=longlat +datum=WGS84"))
# 
# 
# MARCFAS<-shapefile("R:/Science ERD/HEDMG/MPA_Group/Projects/Snow Crab Genetics/data/CrabFishingAreas_2016.shp")
# MARCFAS= spTransform(MARCFAS, CRS("+proj=longlat +datum=WGS84"))

#Basic map of Atlantic Canada

basemap_atlantic <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
                            dplyr::select(name_en,geometry)%>%
                            st_as_sf()%>%
                            st_union()%>%
                            st_transform(latlong)%>%
                            st_as_sf()%>%
                            mutate(country="Canada"),
                          ne_states(country = "United States of America",returnclass = "sf")%>%
                            dplyr::select(name_en,geometry)%>%
                            st_as_sf()%>%
                            st_union()%>%
                            st_transform(latlong)%>%
                            st_as_sf()%>%
                            mutate(country="USA"))
## Sampling range ####
Sample.Lat.lim=c(43,56)
Sample.Long.lim=c(-68.5,-46)
=======
NAFO<-shapefile("R:/Science ERD/HEDMG/MPA_Group/Projects/Fundian Channel/data/NAFO_Divisions_Shapefiles/Divisions.shp")
NAFO= spTransform(NAFO, CRS("+proj=longlat +datum=WGS84"))


MARCMAS<-shapefile("R:/Science ERD/HEDMG/MPA_Group/Projects/Snow Crab Genetics/data/CrabManagementAreas.shp")
MARCMAS= spTransform(MARCMAS, CRS("+proj=longlat +datum=WGS84"))

NLCMAS<-shapefile("R:/Science ERD/HEDMG/MPA_Group/Projects/Snow Crab Genetics/data/NL_CMAs/CMA.shp")
NLCMAS= spTransform(NLCMAS, CRS("+proj=longlat +datum=WGS84"))


MARCFAS<-shapefile("R:/Science ERD/HEDMG/MPA_Group/Projects/Snow Crab Genetics/data/CrabFishingAreas_2016.shp")
MARCFAS= spTransform(MARCFAS, CRS("+proj=longlat +datum=WGS84"))

#Basic map of Nova Scotia/Maritime region
## Sampling range ####
Sample.Lat.lim=c(41,50)
Sample.Long.lim=c(-67,-46)
>>>>>>> e9f7d1387d6708209891eee2af834d4ea8df08d0
map("worldHires", xlim=Sample.Long.lim, ylim=Sample.Lat.lim, 
    col="grey", fill=TRUE, resolution=0);map.axes()


<<<<<<< HEAD
bathydata <- getNOAA.bathy(-46,-68.5,57,43, res=1,keep=T)
=======
bathydata <- getNOAA.bathy(-45,-68,60,41, res=1,keep=T)
>>>>>>> e9f7d1387d6708209891eee2af834d4ea8df08d0
#bathydata<- read.csv("data/marmap_coord_-68;42;-46;54_res_1.csv")

png(filename = "UpdatedCrabMap.png",width = 4200,height = 3800)
plot(bathydata, image = T, land = TRUE, lwd = 0.01, col = "grey", 
     bpal = list(c(0, max(bathydata), "grey"), c(min(bathydata), 0, blues)))
plot(bathydata, lwd = 3.6, deep = 0, shallow = 0, step = 0, add = TRUE) # 
<<<<<<< HEAD
plot.bathy(bathydata, image=F,lwd = 1, deep = -3000, shallow = 0, step = 200, n=30, drawlabel=T, add = TRUE,col="black") 
#plot(NAFO,lwd=3, add=T)
#plot(MARCMAS,lwd=5,add=T)
#plot(NLCMAS,lwd=5,add=T)
#plot(MARCFAS,lwd=3,add=T)
points(pops$Long,pops$Lat,pch=19,cex=1.5,col="red")
dev.off()

autoplot.bathy(bathydata, geom=c("tile","contour")) +
  scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="cyan") +
  geom_point(data = pops, aes(x = Long, y = Lat),
             colour = 'red', size = 3, alpha = 1) +
  labs(y = "Latitude", x = "Longitude", fill = "Elevation") +
  coord_cartesian(expand = 0)
=======
plot.bathy(bathydata, image=F,lwd = 1, deep = -3000, shallow = 0, step = 200, n=30, drawlabel=T, add = TRUE,col="black") # 
#plot(NAFO,lwd=3, add=T)
plot(MARCMAS,lwd=5,add=T)
plot(NLCMAS,lwd=5,add=T)
#plot(MARCFAS,lwd=3,add=T)
points(pops$Long,pops$Lat,pch=19,cex=2,col="red")
dev.off()

>>>>>>> e9f7d1387d6708209891eee2af834d4ea8df08d0
