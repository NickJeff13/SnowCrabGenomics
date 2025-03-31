## code to make a 250 m contour for the snow crab stations

#load libraries --- 
library(sf)
library(tidyverse)
library(terra)
library(rnaturalearth)

#map projections
latlong <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load coordinates 
snow_crab <- read.csv("data/DTO extractions/Pop_Coords2024_with_GLORYS_info.csv")%>%
             st_as_sf(coords=c("Long","Lat"),crs=latlong)

snow_crab_extent <- snow_crab%>%
                    st_bbox()%>%
                    st_as_sfc()%>%
                    st_transform(CanProj)%>%
                    st_buffer(500*1000)%>% # ~ 500 km buffer
                    st_as_sf()%>%
                    st_transform(latlong)

#load canadawide GEBCO
can_gebco <- rast("r:/Science/CESD/HES_MPAGroup/Data/Bathymetry/GEBCO/gebco_2019_Canada.tif")

#crop to Atlantic canada
atlantic_gebco <- can_gebco%>%
                  terra::crop(.,snow_crab_extent%>%st_transform(st_crs(can_gebco)))
#create contour
contour_250 <- as.contour(atlantic_gebco,levels=-250)%>%
               st_as_sf()%>%
               st_transform(latlong)

st_write(contour_250,dsn = "data/bathymetry/contour_250.shp")

contour_250 <- st_read("data/bathymetry/contour_250.shp") #read it in 

#plot it up to make sure it looks right --------
#load basemap
basemap <- ne_states(country = "Canada",returnclass = "sf")%>%
  dplyr::select(name_en,geometry)%>%
  st_as_sf()%>%
  st_union()%>%
  st_transform(latlong)%>%
  st_as_sf()%>%
  mutate(country="Canada")%>%
  rbind(.,ne_states(country = "United States of America",returnclass = "sf")%>%
          dplyr::select(name_en,geometry)%>%
          st_as_sf()%>%
          st_union()%>%
          st_transform(latlong)%>%
          st_as_sf()%>%
          mutate(country="US"),
        ne_states(country = "Greenland",returnclass = "sf")%>%
          dplyr::select(name_en,geometry)%>%
          st_as_sf()%>%
          st_union()%>%
          st_transform(latlong)%>%
          st_as_sf()%>%
          mutate(country="Greenland"),
        ne_states(country = "Iceland",returnclass = "sf")%>%
          dplyr::select(name_en,geometry)%>%
          st_as_sf()%>%
          st_union()%>%
          st_transform(latlong)%>%
          st_as_sf()%>%
          mutate(country="Iceland"))%>%
  st_transform(CanProj)

plot_lims <- snow_crab%>%
             st_transform(CanProj)%>%
             st_buffer(100*1000)%>% # 100 km buffer
             st_bbox()

ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=snow_crab%>%st_transform(CanProj))+
  geom_sf(data=contour_250%>%st_transform(CanProj),linetype=2,lwd=0.3,col="grey40")+
  coord_sf(xlim=plot_lims[c(1,3)],ylim=plot_lims[c(2,4)],expand=0)+
  theme_bw()
