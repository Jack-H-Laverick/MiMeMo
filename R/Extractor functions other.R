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
                  Time_step = rep(1:length(unique(Time_step)),
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

#' Convert a U-V velocity field to speed and direction in degrees
#'
#' This function takes a vector of u and v velocities and calculates the direction and speed of the combined movement.
#'
#'This function was lifted from the `Rsenal` package, where it was originally used to calculate wind speeds. All I've done
#'is built a wrapper which accounts for different conventions when describing wind and flow directions.
#'
#' @param u A vector of Zonal currents (from West to East).
#' @param v A vector of Meridional currents (from South to North).
#' @return a dataframe of two columns is returned. Speed contains the composite speed of both velocities on the same scale.
#' Direction is the resolved direction of the flow in degrees, 0 heads north, 90 East, 180 South, 270 West.
#' @family NEMO-MEDUSA spatial tools
#' @export
vectors_2_direction <- function (u, v) {
  u <- -u                                        # This function was built to use wind direction
  v <- -v                                        # Winds  are "opposite", people care about where wind comes from, not where it goes

  # Lovingly lifted from the "Rsenal" package

  degrees <- function(radians) 180 * radians/pi
  mathdegs <- degrees(atan2(v, u))
  wdcalc <- ifelse(mathdegs > 0, mathdegs, mathdegs + 360)
  uvDirection <- ifelse(wdcalc < 270, 270 - wdcalc, 270 - wdcalc + 360)
  uvSpeed <- sqrt(u^2 + v^2)
  return(cbind(uvDirection, uvSpeed))
}

## used for Light and air temperature data which goes into NM, this data uses a different grid and has time stored differently

#' Get Indices to Use When Clipping netcdf Files at Import
#'
#' This function works out how much of a netcdf file to read, to capture the data between a given Lat-Lon window.
#'
#' The function reads in a vector for both latitudes and longitudes, and tests whether each entry is within the specified
#' window. The max and min position in these vectors where the condition == TRUE are taken to define the ends of the window
#' to import. The vectors of latitudes and longitudes between these limits are kept, so they can be added to the variables
#' of interest during extraction.
#'
#' @param file The full name of a netcdf file containing a longitude and latitude dimension.
#' @param w Degrees West to read from.
#' @param e Degrees East to read to.
#' @param s Degrees South to read from.
#' @param n Degrees North to read to.
#' @return A list of three elements:
#' \itemize{
#'  \item{\emph{Lats -}}{ A vector of latitudes from `s` to `n`.}
#'  \item{\emph{Lons -}}{ A vector of longitudes from `w` to `e`.}
#'  \item{\emph{Limits -}}{ A dataframe containing the index to start reading from (Lon_start, Lat_start)
#'  and the length of the vector to read (Lon_count, Lat_count.}
#'  }
#' @family NEMO-MEDUSA spatial tools
#' @export
Window <- function(file, w, e, s, n) {

  #file <- examples[1,]$File ; w = 0 ; e = 180 ; s = 0 ; n = 90

  raw <- ncdf4::nc_open(file)
  lon <- raw$dim$longitude$vals %>% between(w, e)
  W <- min(which(lon == TRUE))
  E <- max(which(lon == TRUE))

  lat <- raw$dim$latitude$vals %>% between(s, n)
  S <- min(which(lat == TRUE))
  N <- max(which(lat == TRUE))

  lons <- raw$dim$longitude$vals[W:E]
  lats <- raw$dim$latitude$vals[S:N]

  Limits <- data.frame("Lon_start" = W, "Lon_count" = E - W + 1, "Lat_start" = S, "Lat_count" = N - S + 1)

  Limits <- list(Lats = lats, Lons = lons, Limits = Limits)
  return(Limits)
}

#' Pull Coordinates from a Simple Feature Geometry Column
#'
#' This function takes an SF object, and adds two columns containing the coordinates in the geometry column.
#'
#' @param data An SF (Simple Feature) object.
#' @return The same object, now with two columns containing the coordinates in the geometry column.
#' @family NEMO-MEDUSA spatial tools
#' @export
sfc_as_cols <- function(x, names = c("x","y")) {
  stopifnot(inherits(x,"sf") && inherits(sf::st_geometry(x),"sfc_POINT"))
  ret <- sf::st_coordinates(x)
  ret <- as.data.frame(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}

#' Reproject from Latitude and Longitude to Project CRS
#'
#' This function takes a dataframe containing a latitude and longitude column, and replaces them with an X and Y column of coordinates
#' in a new CRS.
#'
#'The function converts a dataframe into an SF object and reprojects into a new CRS. Two coordinate columns are extracted from the
#'geometry column using `sfc_as_cols`, before the geometry column is dropped.
#'
#' @param data A dataframe containing Longitude and Latitude.
#' @param crs The new Coordinate Reference System  to project to.
#' @return A dataframe, now with an x and y column specifying the coordinates for points in the projects Coordiante Reference System.
#' @family NEMO-MEDUSA spatial tools
#' @export
reproj <- function(data, crs) {

  data %>%
    sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% # Specify original projection (crs)
    sf::st_transform(crs = crs) %>%                                   # Transform to crs specified in region file
    sfc_as_cols() %>%                                                 # Extract geometry column for geom_segment to work
    sf::st_set_geometry(NULL)                                         # Chuck geometry column
}
