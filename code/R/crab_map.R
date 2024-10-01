## Mapping function 

#load libraries ----
library(ggplot2)
library(sf)
library(tidyr)
library(dplyr)
library(rnaturalearth)

#map projections ----
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"


#Crab fishing areas ----
cfas <- read_sf("data/shapefiles/dfo_east_snow_crab_areas_erase.shp")%>%
        st_transform(CanProj)

#sample coordinates ----
crab_coords <- read.csv("data/Pop_Coords2024.csv")%>%
               st_as_sf(coords=c("Long","Lat"),crs=latlong)%>%
               st_transform(CanProj)

#Basemap ----
basemap <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(CanProj)%>%
                   st_as_sf()%>%
                   mutate(country="Canada"),
                 ne_states(country = "United States of America",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(CanProj)%>%
                   st_as_sf()%>%
                   mutate(country="USA"),
                 ne_states(country = "Greenland",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(CanProj)%>%
                   st_as_sf()%>%
                   mutate(country="USA"))

### set plot extent
plot_extent <- crab_coords%>%
               st_buffer(100*1000)%>%
               st_bbox()

#make map

crab_map <- ggplot()+
            geom_sf(data=basemap)+
            geom_sf(data=cfas,fill=NA)+
            geom_sf(data=crab_coords, size=3)+
            coord_sf(expand=0,xlim=plot_extent[c(1,3)],ylim=plot_extent[c(2,4)])+
            theme_bw();crab_map

ggsave(filename = "CrapMap_2024_AllPops.png",plot = crab_map, device = "png", path = "figures/", width = 10, height=8, dpi = 320, units = "in")
