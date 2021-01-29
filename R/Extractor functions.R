#' Extract Suspended Particulate Matter data from Globcolour
#'
#' This function takes a netcdf from Globcolour and returns the mean Suspended Particulate MAtter (SPM) from
#' within a lat/lon window.
#'
#' @param File a path to a GLobcolour netcdf file.
#' @param Year The year the file contains data for.
#' @param Month The month the file contains data for.
#' @param cropr The bounding box of an SF object to crop to, as returned by `st_bbox()`.
#' @return A dataframe containing point estimates within an area of SPM with time columns appended.
#' @export
get_SPM <- function(File, Year, Month, crop) {

  raw <- ncdf4::nc_open(File)                                               # Open a file
  row <- raw$dim$row$vals                                                   # Extract the data needed to calculate
  center_lat <- ncdf4::ncvar_get(raw, "center_lat")                         # Lat-lon coordinates of pixels
  center_lon <- ncdf4::ncvar_get(raw, "center_lon")
  lon_step <- ncdf4::ncvar_get(raw, "lon_step")
  col <- ncdf4::ncvar_get(raw, "col")                                       # Each row has a different number of columns

  data <- data.frame(Bin = raw$dim$bin$vals,                                # Pull pixel
                     SPM = ncdf4::ncvar_get(raw, "SPM-OC5_mean"),           # SPM
                     Year = Year,                                           # And attach date
                     Month = Month) %>%
    dplyr::mutate(index = (row[Bin] - row[1]+1),                            # Calculate Lat-lon positions
                  latitude = center_lat[index],
                  longitude = (center_lon[index] + col[Bin] * lon_step[index])) %>%
    dplyr::filter(dplyr::between(latitude, crop["ymin"], crop["ymax"]),     # Initial rough crop using the bounding box of a polygon
                  dplyr::between(longitude, crop["xmin"], crop["xmax"])) %>%
    dplyr::select(-index)                                                   # Ditch uneccessary column

  ncdf4::nc_close(raw)                                                      # Close file connection

  return(data) }

#' Get Rivers
#'
#' This function reads river NEMO-MEDUSA river runoff files and reshapes for StrathE2E.
#'
#' The river runoff file for a single year is imported. Points are checked to see if they fall in the model domain.
#' The points which do are summed by months.
#'
#' NOTE! The river run off files provided to MiMeMo have estimates in open ocean. Runoff is channeled to cells in
#' the vicinity of river mouths so that single grid cells don't get the full amount. Of course river run off really
#' comes from the coast, so you may need to modify your domain polygon.
#'
#' @param File The name (with path) of a netcdf file containing river runoff data.
#' @param Year The year the file contains data for.
#' @param domain The polygon used to capture point estimates relevant to the model domain.
#' @return A dataframe containing the total river runoff into the model domain by month with a year column.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_rivers <- function(File, Year, domain) {

  data <- raster::raster(File, varname = "nav_lat") %>%
    raster::as.data.frame(xy = TRUE) %>%                                                              # Get a dataframe of xy positions and latitudes
    dplyr::left_join(raster::as.data.frame(raster::raster(File, varname = "nav_lon"), xy = TRUE)) %>% # Bind to the same for longitudes
    dplyr::left_join(raster::as.data.frame(raster::brick(File, varname = "sorunoff"), xy = TRUE)) %>% # Bind to the river run off values for all months
    sf::st_as_sf(coords = c("nav_lon", "nav_lat"), crs = 4326) %>%                                    # Convert to sf
    dplyr::select(-c(x,y)) %>%                                                                        # Drop cell positions used for binding
    tidyr::drop_na() %>%                                                                              # Drop empty pixels
    sf::st_transform(sf::st_crs(domain)) %>%                                                          # Transform points to the same crs as domain polygon
    sf::st_join(domain) %>%                                                                           # Check which points are in the domain
    sf::st_drop_geometry() %>%                                                                        # Drop sf formatting
    tidyr::drop_na() %>%                                                                              # Drop points outside the domain
    colSums()                                                                                         # Total all river runoff in a month

  names(data) <- c(1:12)                                                                              # Label each column by month

  result <- data.frame(Year = Year, Month = 1:12, Runoff = data)                                      # Create a dataframe of runoff by month and add year column

  return(result)
}

#' Get Surface Irradiance & Air Temperature
#'
#' This function reads either title variable from a NEMO-MEDUSA model \strong{DRIVERS} file, and returns a monthly time series.
#'
#' The appropriate variable in the netcdf file is imported according to the `Type` parameter, only reading within an
#' x/y window specified in `Space`. The function then drops points outside the model domain before constructing the monthly
#' time series.
#'
#' Unlike other NEMO-MEDUSA related get_* functions, this function constructs the monthly time series directly. A wrapper function
#' to handle time isn't required, as each netcdf file for driving data contains 360 day steps for a single year. `stratify` also
#' isn't called, as these variables have no depth dimension.
#'
#' @param File The location of the netcdf file containing NEMO-MEDUSA driving data.
#' @param Type The variable contained in the netcdf file either "T150" (air temperature) or "SWF" (surface irradiance).
#' @param Year The year the necdf file contains data for.
#' @return A dataframe containing a monthly time series within a year of either average air temperature or surface
#' irradiance. Air temperature is also split by shore zone.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_air <- function(File, Type, Year) {

  #File <- Airtemp_files$File[1] ; Type <- Airtemp_files$Type[1] ; Year <- Airtemp_files$Year[1] # test
  if(Type == "SWF") months <- Light_months                                   # Get the file holding the months
  if(Type == "T150") months <- Airtemp_months                                # For the time steps of this data

  nc_raw <- ncdf4::nc_open(File)                                             # Open up a netcdf file to see it's raw contents (var names)
  nc_var <- ncdf4::ncvar_get(nc_raw, Type, c(Space$Limits$Lon_start, Space$Limits$Lat_start, 1),  # Extract the variable of interest
                             c(Space$Limits$Lon_count, Space$Limits$Lat_count, -1)) # cropped to window, with all time steps
  ncdf4::nc_close(nc_raw)                                                           # You must close an open netcdf file when finished to avoid data loss

  Data <- as.data.frame.table(nc_var, responseName = "Measured") %>%         # Reshape array as dataframe
    dplyr::rename(Longitude = Var1, Latitude = Var2, Time_step = Var3) %>%   # Name the columns
    dplyr::mutate(Longitude = rep(rep(Space$Lons,                            # Replace the factor levels with dimension values
                                      times = length(unique(Latitude))), times = length(unique(Time_step))),
                  Latitude = rep(rep(Space$Lats,
                                     each = length(unique(Longitude))), times = length(unique(Time_step))),
                  Time_step = rep(seq_along(unique(Time_step)),
                                  each = length(unique(Latitude)) * length(unique(Longitude)))) %>%
    dplyr::right_join(domains_mask) %>%                                      # Crop to domain
    dplyr::left_join(months) %>%                                             # Assign a month to each time step
    dplyr::mutate(Year = Year,                                               # Attach Year
                  Type = Type)                                                      # Attach variable name

  if(Type == "SWF") Data <- dplyr::group_by(Data, Month, Year, Type)         # We don't need to bother accounting for shore in light data
  if(Type == "T150") Data <- dplyr::group_by(Data, Month, Year, Type, Shore) # We care about shore for temperature, retain interesting columns

  Summary <- dplyr::summarise(Data, Measured = stats::weighted.mean(Measured, Cell_area)) # Average by time step.

  return(Summary)
}

#' Get Surface Irradiance & Air Temperature
#'
#' This function reads either title variable from a NEMO-MEDUSA model \strong{DRIVERS} file, and returns a monthly time series.
#'
#' The appropriate variable in the netcdf file is imported according to the `Type` parameter, only reading within an
#' x/y window specified in `Space`. The function then drops points outside the model domain before constructing the monthly
#' time series.
#'
#' Unlike other NEMO-MEDUSA related get_* functions, this function constructs the monthly time series directly. A wrapper function
#' to handle time isn't required, as each netcdf file for driving data contains 360 day steps for a single year. `stratify` also
#' isn't called, as these variables have no depth dimension.
#'
#' @param File The location of the netcdf file containing NEMO-MEDUSA driving data.
#' @param Type The variable contained in the netcdf file either "T150" (air temperature) or "SWF" (surface irradiance).
#' @param Year The year the necdf file contains data for.
#' @return A dataframe containing a monthly time series within a year of either average air temperature or surface
#' irradiance. Air temperature is also split by shore zone.
#' @importFrom data.table := as.data.table setnames
#' @family NEMO-MEDUSA variable extractors
#' @export
get_air_dt <- function(File, Type, Year) {

  #File <- all_files$File[1] ; Type <- all_files$Type[1] ; Year <- all_files$Year[1] # test
  if(Type == "SWF") months <- Light_months                                 # Get the file holding the months
  if(Type == "T150") months <- Airtemp_months                              # For the time steps of this data

  nc_raw <- ncdf4::nc_open(File)                                             # Open up a netcdf file to see it's raw contents (var names)
  nc_var <- ncdf4::ncvar_get(nc_raw, Type, c(Space$Limits$Lon_start, Space$Limits$Lat_start, 1),  # Extract the variable of interest
                             c(Space$Limits$Lon_count, Space$Limits$Lat_count, -1)) # cropped to window, with all time steps
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  DT <- as.data.table(nc_var, value.name = "Measured") %>%                   # Pull array
    setnames(old = c("V1", "V2", "V3"), new = c("Longitude", "Latitude", "Time_step")) %>% # Name the columns
    .[, ':='(Longitude = Space$Lons[Longitude],                              # read ':=' as mutate
             Latitude = Space$Lats[Latitude],                                # Replace the factor levels with dimension values
             Month = months[Time_step, "Month"],                             # Assign months to time steps
             Year = Year,                                                    # Add year
             Type = Type)]                                                   # Add variable name
  setkey(DT, Longitude, Latitude)                                            # Set key columns for speedy joins

  DT <- DT[domains_mask] #%>%
  #merge(domains_mask, all.y = TRUE)                                        # Crop to domain

  ## Variable specific summaries

  if(Type == "SWF") Data <- DT[,by = .(Month, Year, Type),                   # We don't need to bother accounting for shore in light data
                               .(Measured = weighted.mean(Measured, Cell_area))]  # Average by time ste, weighted by cell area
  if(Type == "T150") Data <- DT[,by= .(Month, Year, Type, Shore),            # We care about shore for temperature, retain interesting columns
                                .(Measured = weighted.mean(Measured, Cell_area))]  # Average by time ste, weighted by cell area
  return(Data)
}
