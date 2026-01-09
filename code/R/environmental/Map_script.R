#Map of Atlantic Canada with crab sampling locations included for Figure 1 of manuscript

# Load libraries ----------------------------------------------------------

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
library(ggplot2)
library(sf)
library(tidyr)
library(dplyr)
library(rnaturalearth)

#map projections ----
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"


#Crab fishing areas ----
cfas <- read_sf("~/Documents/GitHub/SnowCrabGenomics/data/shapefiles/dfo_east_snow_crab_areas_erase.shp")%>%
  st_transform(latlong)

#sample coordinates ----
crab_coords <- read.csv("~/Documents/GitHub/SnowCrabGenomics/data/Pop_Coords2025.csv")%>%
  st_as_sf(coords=c("Long","Lat"),crs=latlong)%>%
  st_transform(latlong)

# Read in bathymetry contour 

bathy <- read_sf("~/Documents/GitHub/SnowCrabGenomics/data/shapefiles/bathymetry/contour_250.shp") %>% 
  st_transform(latlong)


#Basemap ----
basemap <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
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

### set plot extent
plot_extent <- crab_coords%>%
  st_buffer(100*2000)%>%
  st_bbox()

# Add PC1 and PC2 axis values from pcadapt script to map

map.df2 <- left_join(crab_coords,map.df, by=c("SampleSite"="pca.with.pops.pops"))


# Plot the map with ggplot2 -----------------------------------------------

crab_map <- ggplot()+
  geom_sf(data=basemap,color="black")+
  #geom_sf(data=cfas,fill=NA)+
  geom_sf(data=bathy, fill=NA, color="lightgrey")+
  geom_sf(data=map.df2, aes(fill=mean_PC1),colour= "black",shape=21, size=4)+
  scale_fill_viridis_c()+
  #geom_sf_label(data=crab_coords, aes(label = SampleSite))+
  coord_sf(expand=0,xlim=plot_extent[c(1,3)],ylim=plot_extent[c(2,4)])+
  theme_bw();crab_map

ggsave(filename = "CrabMap_PC1.png",plot = crab_map, device = "png", 
       path = "~/Documents/GitHub/SnowCrabGenomics/figures/", 
       width = 10, height = 8, dpi = 320, units = "in")



#OR use this older way with marmap and some base plotting



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

