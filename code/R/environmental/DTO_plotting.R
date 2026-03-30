## DTO plotting

#load libraries ----------
library(tidyverse)
library(sf)
library(rnaturalearth)
library(stringr)
library(viridis)
library(marmap)
library(terra)
library(tidyterra)
library(ggspatial)
library(ggrepel)

#source functions
source("https://raw.githubusercontent.com/dfo-mar-mpas/MCRG_functions/refs/heads/main/code/marmap_to_isobath.R")
source("https://raw.githubusercontent.com/dfo-mar-mpas/MCRG_functions/refs/heads/main/code/trim_img_ws.R")

#Map projections ---------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#load the extracted data ---------
load("data/DTO extractions/process_extractions.RData")

#load Canadian EEZ
can_eez <- read_sf("c:/Users/stanleyr/Documents/GitHub/stannsbank_mpa/data/Shapefiles/can_eez.shp")%>%
           st_transform(CanProj)

#load the table matching ---------
table_match <- read.csv("data/pop_coords_table_match.csv")

#download a basemap --------
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

#set up plot limits
plot_lims <- mean_df%>%
            distinct(site,.keep_all=TRUE)%>%
            st_as_sf(coords=c("long","lat"),crs=latlong)%>%
            st_transform(CanProj)%>%
            st_bbox()%>% #get the bounding box
            st_as_sfc()%>%
            st_buffer(70000)%>% #create a buffer on that bounding box of 70 km - this is faster than doing a buffer on the polygon
            st_bbox()


plot_lims2 <- mean_df%>%
              distinct(site,.keep_all=TRUE)%>%
              st_as_sf(coords=c("long","lat"),crs=latlong)%>%
              st_transform(CanProj)%>%
              st_bbox()%>% #get the bounding box
              st_as_sfc()%>%
              st_buffer(300*1000)%>% #create a buffer on that bounding box of 70 km - this is faster than doing a buffer on the polygon
              st_bbox()

#download bathymetry data
bathy_lims <- mean_df%>%
              distinct(site,.keep_all=TRUE)%>%
              st_as_sf(coords=c("long","lat"),crs=latlong)%>%
              st_bbox()%>% #get the bounding box
              st_as_sfc()%>%
              st_buffer(325*1000)%>%#create a buffer on that bounding box of 70 km - this is faster than doing a buffer on the polygon
              st_bbox()
        

# noaabathy_plotregion <- getNOAA.bathy(bathy_lims[1]-1,bathy_lims[3]+1,bathy_lims[2]-1,bathy_lims[4]+1,resolution = 0.25,keep=T)
# save(noaabathy_plotregion,file = "data/noaabathy_plotregion.Rdata")
# bathy_xyz <- as.xyz(noaabathy_plotregion)#
# bathy_rast <- rast(x = bathy_xyz,type = "xyz",crs = "EPSG:4326")%>%project(CanProj)
# 
# writeRaster(bathy_rast,filename ="data/bathy_rast.tif")

bathy_rast <- rast("data/bathy_rast.tif")
deep_val <- -3500

# fill projection gaps
bathy_rast[is.na(bathy_rast)] <- deep_val

# ensure land values plot consistently
bathy_rast <- clamp(bathy_rast, upper = 0)
bathy_rast <- clamp(bathy_rast, lower = deep_val)

cont_250_plotregion <- marmap_to_isobath(noaabathy_plotregion,-250,latlong)%>%st_transform(CanProj)
cont_350_plotregion <- marmap_to_isobath(noaabathy_plotregion,-350,latlong)%>%st_transform(CanProj)



#process site mean data for the last two decades (2004-2024)
#temperature
# site_mean_df <- mean_df%>%
#   filter(var=="thetao",
#          grouping=="annual")%>%
#   mutate(year=as.numeric(time))%>%
#   filter(year>2003)%>%
#   group_by(site,season)%>%
#   summarise(lat=unique(lat),
#             long=unique(long),
#             prec_lower = quantile(mean,0.1),
#             prec_upper = quantile(mean,0.9),
#             mean=mean(mean),
#             min=min(min),
#             max=max(max),
#             sd=sd(mean))%>%
#   ungroup()%>%
#   data.frame()%>%
#   mutate(season=str_to_title(season))%>%
#   st_as_sf(coords=c("long","lat"),crs=latlong)%>%
#   st_transform(CanProj)
# 
# #salinity
# site_mean_df2 <- mean_df%>%
#   filter(var=="so",
#          grouping=="annual")%>%
#   mutate(year=as.numeric(time))%>%
#   filter(year>2003)%>%
#   group_by(site,season)%>%
#   summarise(lat=unique(lat),
#             long=unique(long),
#             prec_lower = quantile(mean,0.1),
#             prec_upper = quantile(mean,0.9),
#             mean=mean(mean),
#             min=min(min),
#             max=max(max),
#             sd=sd(mean))%>%
#   ungroup()%>%
#   data.frame()%>%
#   mutate(season=str_to_title(season))%>%
#   st_as_sf(coords=c("long","lat"),crs=latlong)%>%
#   st_transform(CanProj)

## annual differences plot

plot_mean_df2 <- mean_df%>%
  filter(var=="thetao",
         grouping=="annual")%>%
  mutate(year=as.numeric(time))%>%
  filter(year>2003)%>%
  group_by(site)%>%
  summarise(lat=unique(lat),
            long=unique(long),
            prec_lower = quantile(mean,0.1),
            prec_upper = quantile(mean,0.9),
            mean=mean(mean), #average climatology
            min=mean(min), #mean of the extremes among seasons and years
            max=mean(max),
            sd=sd(mean))%>%
  ungroup()%>%
  data.frame()%>%
  left_join(.,table_match%>%dplyr::select(site,table_name,table_code,table_region))%>%
  st_as_sf(coords=c("long","lat"),crs=latlong)%>%
  st_transform(CanProj)

#make some plots

##

p1 <- ggplot() +
  geom_spatraster(data = bathy_rast) +
  geom_sf(data=cont_350_plotregion,col="grey40",linewidth=0.25)+
  geom_sf(data = can_eez, fill = NA,col="grey90") +
  geom_sf(data = basemap) +
  geom_sf(data = basemap %>% filter(country == "Canada"), fill = "grey60") +
  geom_sf(data = plot_mean_df2, fill = "white", shape = 21, size = 3.5) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "inside",
    legend.position.inside = c(0.07,0.72),
    legend.background =  element_blank()
  ) +
  scale_fill_gradientn(
    colors = c("lightblue", "dodgerblue4"),
    values = scales::rescale(c(0, -50, -3500)),
    name = "Depth (m)",
    guide = guide_colorbar(reverse = FALSE,
                           frame.colour = "black",
                           ticks.colour = "black")
  )+
  coord_sf(
    expand = 0,
    xlim = plot_lims[c(1,3)],
    ylim = plot_lims[c(2,4)]
  )+
  annotation_scale(location="br")
  
plot_mean_df2_lab <- plot_mean_df2 %>%
  cbind(st_coordinates(.))

p1_labs <- ggplot() +
  geom_spatraster(data = bathy_rast) +
  geom_sf(data=cont_350_plotregion,col="grey40",linewidth=0.25)+
  geom_sf(data = can_eez, fill = NA) +
  geom_sf(data = basemap) +
  geom_sf(data = basemap %>% filter(country == "Canada"), fill = "grey60") +
  geom_sf(data = plot_mean_df2, fill = "white", shape = 21, size = 3.5) +
  geom_text_repel(
    data = plot_mean_df2_lab,
    aes(x = X, y = Y, label = table_code),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey40"
  )+
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    legend.position = "inside",
    legend.position.inside = c(0.07,0.72),
    legend.background =  element_blank()
  ) +
  scale_fill_gradientn(
    colors = c("lightblue", "dodgerblue4"),
    values = scales::rescale(c(0, -50, -3500)),
    name = "Depth (m)",
    guide = guide_colorbar(reverse = FALSE,
                           frame.colour = "black",
                           ticks.colour = "black")
  )+
  coord_sf(
    expand = 0,
    xlim = plot_lims[c(1,3)],
    ylim = plot_lims[c(2,4)]
  )+
  annotation_scale(location="br")

ggsave("figures/figure1_map.jpg",p1,height=6,width=6,units="in",dpi=300)
trim_img_ws("figures/figure1_map.jpg")
ggsave("figures/figure1_map_labs.jpg",p1_labs,height=6,width=6,units="in",dpi=300)

#make the globe inset

center_pt <- plot_lims%>%
  st_as_sfc()%>%
  st_transform(latlong)%>%
  st_centroid()

lon0 <- st_coordinates(center_pt)[1]
lat0 <- st_coordinates(center_pt)[2]

globe_crs <- sprintf("+proj=ortho +lat_0=%s +lon_0=%s",lat0, lon0)

#download the world globe basemap
world_globe <- ne_countries(
  scale = "medium",
  returnclass = "sf"
) %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>%
  st_transform(globe_crs)



#define the box denoting the study region you want to highlight
global_box <- plot_lims%>%
  st_as_sfc() %>%
  st_transform(globe_crs)

#make a circle to wrap the globe plot

globe_circle <- st_sfc(
  st_buffer(
    st_point(c(0, 0)),   # center of orthographic projection
    dist =  6378137  # meters
  ),
  crs = globe_crs
)

#crudgy way to make it so that the oceans are white in the plot
globe_disc <- st_sfc(
  st_point(c(0, 0)),  # center in projected coords
  crs = globe_crs
) %>%
  st_buffer(dist = 6378137) %>%   # Earth radius in meters
  st_as_sf()


global_inset2 <- ggplot() +
  
  #global background
  geom_sf(data = globe_disc, fill = "white", colour = "black", linewidth = 0.4)+
  
  # Land
  geom_sf(
    data = world_globe,
    colour = "grey20",
    linewidth = 0.2
  ) +
  
  geom_sf(
    data = world_globe%>%filter(formal_en == "Canada"),
    fill = "grey60",
    colour = "grey20",
    linewidth = 0.2
  ) +
  
  # Study region box
  geom_sf(
    data = global_box,
    fill = NA,
    colour = "black",
    linewidth = 0.9
  ) +
  
  geom_sf(data = globe_circle,
          fill = NA,
          colour = "grey30",
          linewidth = 0.4)+
  
  coord_sf(crs = globe_crs) +
  
  theme_void() +
  
  theme(
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background  = element_rect(fill = NA, colour = NA)
  )

ggsave(
  "figures/global_inset2.png",
  plot = global_inset2,
  width = 4,
  height = 4,
  dpi = 600,
  bg = "transparent"
)





#minimum temp
p1 <- ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=site_mean_df,aes(fill=min),shape=21,size=2.5)+
  facet_wrap(~season,ncol=2)+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        legend.position="bottom")+
  scale_fill_viridis(option="C")+
  coord_sf(expand=0,xlim=plot_lims[c(1,3)],ylim=plot_lims[c(2,4)])+
  labs(fill=expression("Temperature ("*degree*C*")"),
       title="Minimum bottom temperature")

#maximum temp
p2 <- ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=site_mean_df,aes(fill=max),shape=21,size=2.5)+
  facet_wrap(~season,ncol=2)+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        legend.position="bottom")+
  scale_fill_viridis(option="C")+
  coord_sf(expand=0,xlim=plot_lims[c(1,3)],ylim=plot_lims[c(2,4)])+
  labs(fill=expression("Temperature ("*degree*C*")"),
       title="Maximum bottom temperature")

#mean temp
p3 <- ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=site_mean_df,aes(fill=mean),shape=21,size=2.5)+
  facet_wrap(~season,ncol=2)+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        legend.position="bottom")+
  scale_fill_viridis(option="C")+
  coord_sf(expand=0,xlim=plot_lims[c(1,3)],ylim=plot_lims[c(2,4)])+
  labs(fill=expression("Temperature ("*degree*C*")"),
       title="Mean bottom temperature");p3

#save plots
ggsave("figures/minimum_bottom_temp.png",p1,height=6,width=6,units="in",dpi=300)
ggsave("figures/maximum_bottom_temp.png",p2,height=6,width=6,units="in",dpi=300)
ggsave("figures/mean_bottom_temp.png",p3,height=6,width=6,units="in",dpi=300)


#Drop the geometry values from the above data to make it easier to merge with coords file
site_temp <- site_mean_df %>% st_drop_geometry()
site_salin <- site_mean_df2 %>% st_drop_geometry()

write.csv(x = site_temp, file = "data/DTO extractions/Site_TemperatureData.csv", quote = F)
write.csv(x = site_salin, file = "data/DTO extractions/Site_SalinityData.csv", quote = F)
