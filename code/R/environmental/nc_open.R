## unload the data from the netcdf. Note that this code was translated from Mike Casey's Matlab code 
nc_unload <- function(filename, time_standard = TRUE, timename = "time", time_only = FALSE) {
  
  require(ncdf4)
  
  # Open the NetCDF file
  nc <- nc_open(filename)
  
  # Initialize the data list
  Data <- list()
  
  if (time_standard) {
    # Try to read the time variable
    varid_tcounter <- tryCatch({
      ncvar_get(nc, timename)
    }, error = function(e) {
      timename <- paste0(timename, "_counter")
      ncvar_get(nc, timename)
    })
    
    time_counter <- as.numeric(varid_tcounter)
    
    # Extract time units
    time_counter_units <- ncatt_get(nc, timename, "units")$value
    time_counter_units <- gsub("UTC|T|Z", "", time_counter_units)
    
    # Parse the time units for the reference date
    tempdat <- unlist(strsplit(time_counter_units, " "))
    
    if (length(tempdat) < 4) {
      time_base_dv <- as.POSIXlt(tempdat[3], format = "%Y-%m-%d")
    } else {
      time_base_dv <- as.POSIXlt(paste(tempdat[3], tempdat[4]), format = "%Y-%m-%d %H:%M:%S")
    }
    
    # Convert time based on units (days, hours, seconds)
    if (tempdat[1] == "days") {
      time_counter_dt <- time_base_dv + time_counter * 86400
    } else if (tempdat[1] == "hours") {
      time_counter_dt <- time_base_dv + time_counter * 3600
    } else if (tempdat[1] == "seconds") {
      time_counter_dt <- time_base_dv + time_counter
    } else {
      stop(paste("Unknown time origin:", tempdat[1]))
    }
    
    # Store time data in Data structure
    Data$time_counter_dt <- as.POSIXct(time_counter_dt, origin = "1970-01-01")
    Data$time_counter_unit <- time_counter_units
  }
  
  if (!time_only) {
    # Get information about the NetCDF variables
    tempinfo <- nc$var
    
    for (i in seq_along(tempinfo)) {
      # Get the variable name
      vname <- tempinfo[[i]]$name
      
      # Get the variable data
      tempdata <- ncvar_get(nc, vname)
      
      # Get attributes for the variable (units, long name, etc.)
      if (!is.character(tempdata)) {
        if ("units" %in% names(nc$var[[i]]$attributes)) {
          Data[[paste0(vname, "_unit")]] <- ncatt_get(nc, vname, "units")$value
        }
        if ("long_name" %in% names(nc$var[[i]]$attributes)) {
          Data[[paste0(vname, "_lname")]] <- ncatt_get(nc, vname, "long_name")$value
        }
        
        # Handle missing values
        if ("_FillValue" %in% names(nc$var[[i]]$attributes)) {
          missing_value <- ncatt_get(nc, vname, "_FillValue")$value
          tempdata[tempdata == missing_value] <- NA
        } else if ("missing_value" %in% names(nc$var[[i]]$attributes)) {
          missing_value <- ncatt_get(nc, vname, "missing_value")$value
          tempdata[tempdata == missing_value] <- NA
        }
        
        # Apply scale factor and offset
        if ("scale_factor" %in% names(nc$var[[i]]$attributes)) {
          scale_factor <- ncatt_get(nc, vname, "scale_factor")$value
          add_offset <- if ("add_offset" %in% names(nc$var[[i]]$attributes)) ncatt_get(nc, vname, "add_offset")$value else 0
          tempdata <- tempdata * scale_factor + add_offset
        }
      }
      
      # Store the variable data in Data structure
      Data[[vname]] <- tempdata
    }
  }
  
  # Close the NetCDF file
  nc_close(nc)
  
  return(Data)
}

#Process the data ------

##code to process monthly means variables from the bottom 

nc_process <- function(fname){
  
  #fname = file name with full or relative path
  
  #load libraries
  require(dplyr)
  require(lubridate)
  require(terra)
  require(ncdf4)
  require(purrr)
  
  source("code/R/environmental/nc_open.R") #this is needed to open the data. 
  
  #get metadata from netcdf
  nc <- nc_unload(fname)
  
  #extract time
  time <- nc[["time_counter_dt"]]%>%as.POSIXct()
  time_df <- data.frame(year=year(time),month = month(time),day = day(time))%>%
    mutate(id=1:n(),
           season=case_when(month %in% 1:3 ~ "winter",
                            month %in% 4:6 ~ "spring",
                            month %in% 7:9 ~ "summer",
                            month %in% 10:12 ~ "fall"))
  
  #extract depths
  nc_depth <- nc_open(fname)
  depths <- ncdf4::ncvar_get(nc_depth,"depth")
  
  #extract variables
  
  #just the variables we want averages of
  dto_vars <- names(nc)[!grepl("time",names(nc)) & !grepl("zos",names(nc))]
  
  #variables for the collapsing of matrices
  col_names <- time_df$id
  row_names <- depths
  
  #quick function to collapse the matrix using purr 
  flatten_matrix <- function(mat, mat_name) {
    as.data.frame(mat) %>%
      mutate(Row = row_names) %>%  # Add row names from vector
      pivot_longer(-Row, names_to = "Column", values_to = "Value") %>%
      mutate(Column = col_names[as.numeric(gsub("V", "", Column))]) %>%  # Add column names from vector
      mutate(Matrix = mat_name)  # Add matrix name
  }
  
  #site ID
  site <- basename(fname)%>%
    gsub("GLORYS_monthly_","",.)%>%
    gsub("_1993_2024.nc","",.)
  
  #flatten the dataframes and combine with the time dataframe. 
  flattened_df <- imap_dfr(nc[dto_vars], flatten_matrix)%>%
    data.frame()%>%
    rename(depth=Row,id=Column,value=Value,var=Matrix)%>%
    left_join(.,time_df)%>%
    dplyr::select(-id)%>%
    mutate(site=site)
  
  return(flattened_df)
  
} #end function
