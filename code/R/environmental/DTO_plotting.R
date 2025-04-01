## DTO plotting

#load libraries ----------
library(tidyverse)
library(sf)
library(rnaturalearth)
library(stringr)
library(viridis)

#Map projections ---------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#load the extracted data ---------
load("data/DTO extractions/process_extractions.RData")

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

#process site mean data for the last two decades (2004-2024)
#temperature
site_mean_df <- mean_df%>%
  filter(var=="thetao",
         grouping=="annual")%>%
  mutate(year=as.numeric(time))%>%
  filter(year>2003)%>%
  group_by(site,season)%>%
  summarise(lat=unique(lat),
            long=unique(long),
            prec_lower = quantile(mean,0.1),
            prec_upper = quantile(mean,0.9),
            mean=mean(mean),
            min=min(min),
            max=max(max),
            sd=sd(mean))%>%
  ungroup()%>%
  data.frame()%>%
  mutate(season=str_to_title(season))%>%
  st_as_sf(coords=c("long","lat"),crs=latlong)%>%
  st_transform(CanProj)

#salinity
site_mean_df2 <- mean_df%>%
  filter(var=="so",
         grouping=="annual")%>%
  mutate(year=as.numeric(time))%>%
  filter(year>2003)%>%
  group_by(site,season)%>%
  summarise(lat=unique(lat),
            long=unique(long),
            prec_lower = quantile(mean,0.1),
            prec_upper = quantile(mean,0.9),
            mean=mean(mean),
            min=min(min),
            max=max(max),
            sd=sd(mean))%>%
  ungroup()%>%
  data.frame()%>%
  mutate(season=str_to_title(season))%>%
  st_as_sf(coords=c("long","lat"),crs=latlong)%>%
  st_transform(CanProj)

#make some plots

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
