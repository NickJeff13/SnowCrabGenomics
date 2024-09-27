## This is a processing script for the Digital Twin Ocean (DTO) extractions

#load libraries ----- 
library(terra)
library(sf)
library(tidyverse)
library(ncdf4)
library(ncdf4.helpers)
library(stringr)
library(purrr)

#load functions
source("code/R/nc_open.R")# opening and (pre)processing scripts

#filepaths for all the nc datafiles
nc_paths <- dir("data/DTO extractions/locations_nc/",full.names = TRUE)

process_nc_file <- function(nc_path) {
        temp <- nc_process(nc_path)
        
        temp %>%
          filter(!is.na(value)) %>%
          group_by(year, var, season) %>%
          filter(depth == max(depth)) %>%  # Only for bottom temperature
          summarise(depth = unique(depth),
                    mean = mean(value), #mean of the value
                    min = min(value), #min of the value
                    max = max(value), #max of the value
                    sd = sd(value), #sd of the value
                    prec_lower = quantile(value,0.1), #10th percentile in case the 'min' is an outlier
                    perc_upper = quantile(value,0.9), #90th percentile in case the 'max' is an outlier
                    #count = n(), #number of observations -- note this is 3 for all except 2024, since it is only inclusive of part of the spring. 
                    .groups = "drop") %>%
          mutate(site = unique(temp$site))
}

#apply and group the dataset together. 
mean_df <- map_dfr(nc_paths, process_nc_file)

save(mean_df,file="data/DTO extractions/process_extractions.RData")

