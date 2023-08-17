## Load Libraries --------------
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
library(maptools)
library(lattice)

setwd("R:/Science/CESD/HES_MPAGroup/Projects/Snow Crab Genetics/")
##Snow Crab Locations Map#####
source("PlotBathy2.R")

pops<-read.csv("data/Pop_Coordinates.csv",header = T)
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
map("worldHires", xlim=Sample.Long.lim, ylim=Sample.Lat.lim, 
    col="grey", fill=TRUE, resolution=0);map.axes()


bathydata <- getNOAA.bathy(-45,-68,60,41, res=1,keep=T)
#bathydata<- read.csv("data/marmap_coord_-68;42;-46;54_res_1.csv")

png(filename = "UpdatedCrabMap.png",width = 4200,height = 3800)
plot(bathydata, image = T, land = TRUE, lwd = 0.01, col = "grey", 
     bpal = list(c(0, max(bathydata), "grey"), c(min(bathydata), 0, blues)))
plot(bathydata, lwd = 3.6, deep = 0, shallow = 0, step = 0, add = TRUE) # 
plot.bathy(bathydata, image=F,lwd = 1, deep = -3000, shallow = 0, step = 200, n=30, drawlabel=T, add = TRUE,col="black") # 
#plot(NAFO,lwd=3, add=T)
plot(MARCMAS,lwd=5,add=T)
plot(NLCMAS,lwd=5,add=T)
#plot(MARCFAS,lwd=3,add=T)
points(pops$Long,pops$Lat,pch=19,cex=2,col="red")
dev.off()

