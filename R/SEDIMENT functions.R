
#### Water movement tools    ####

#' @export
mync_get_pixel <- function(File, var, pixel_x, pixel_y) {

  ncvar_get(nc_open(File), var, c(pixel_x, pixel_y, 1, 1),  # Extract the variable of interest
            c(1, 1, -1, -1)) %>%                            # cropped to a target pixel, with all time steps and ensemble members
    colSums(na.rm = TRUE)                                   # Average across ensemble quickly
  # Setting na.rm = TRUE replaces iced NAs with 0s
}         # Import a variable clipped to Window

#' @export
get_waves_ts <- function(File, vars, Year, pixel_x, pixel_y) {

  #File <- all_files$File[1] ; Year <- all_files$Year[1]                           # test
  #pixel_lon <- ECMWF_mask$Longitude[1] ; pixel_lat <- ECMWF_mask$Latitude[1]  # test
  #pixel_x <- ECMWF_mask$x_index[1] ; pixel_y <- ECMWF_mask$y_index[1]         # test
  #vars <- c("swh", "mwd")                                                         # test

  Data <- map(vars, ~{mync_get_pixel(File, var = .x, pixel_x, pixel_y)})                   # Extract the variables of interest
  Data <- cbind(Data[[1]], Data[[2]], Data[[3]])
  colnames(Data) <- vars

  TS <- cbind(Data,
              Year = rep(Year, nrow(Data)),
              Month = rep(1:12, each = 8),
              Hour = rep(seq(0, 21, by = 3), times = 12))

  return(TS)
}    # Pull significant wave height monthly time series per zone in a year file

#' @export
sample_waves <- function(pixel) {

  months <- map2(wave_files$File, wave_files$Year, get_waves_ts, vars = c("swh", "mwd", "mwp"), # For each year file.
                 pixel_x = aligned$x_index[pixel], pixel_y = aligned$y_index[pixel]) %>%  # Extract the time series for a pixel
    do.call(rbind, .)                                                                  # Bind the years together quickly

  days <- matrix(rep(t(months), 28), ncol = ncol(months), byrow = TRUE) %>%            # Duplicate the months to get a daily tace (quickly)
    cbind(Day = rep(1:28, each = nrow(.)/28))

  colnames(days) <- c(colnames(months), "Day")

  Wave_step <- ISOdate(days[,"Year"], days[,"Month"], days[,"Day"], days[,"Hour"])     # Calculate time step

  days <- cbind(days[, !colnames(days) %in% c("Year", "Month", "Day", "Hour")], Wave_step)  # Replace time columns with time step
}

#' @export
get_tides_ts <- function(File, month, vars, pixel_x, pixel_y) {

  #File <- paste0(tide_files$path[1], tide_files$file[1]) ; Year <- tide_files$year[1] # test
  #pixel_x <- aligned$V1[1] ; pixel_y <- aligned$V2[1]                     # test
  #vars <- c("u_east", "v_north")                                                 # test

  Data <- map(vars, ~{mync_get_pixel(File, var = .x, pixel_x, pixel_y)})          # Extract the variables of interest
  Rescaled <- vectors_2_direction(Data[[1]], Data[[2]])    #*** warning, averaged across depth but not weighted right now by thickness.

  return(Rescaled)

}   # Pull the depth averaged currents from a file

#' @export
sample_tides_year <- function(packet, pixel) {

  #pixel <- 20000
  #months <- tide_files$month[1:11]
  #year <- tide_files$year[1]

  Rescaled <- map(paste0(packet$path, packet$file), get_tides_ts, # For each year file.
                  vars = c("u_east", "v_north"),                                    # Extract the time series for a pixel
                  pixel_x = aligned$V1[pixel], pixel_y = aligned$V2[pixel]) %>%
    do.call(rbind, .) %>%
    cbind(Tide_step = seq(ISOdate(packet$year[1], packet$month[1], 1, 0), by = "2 hours", length.out = nrow(.))) # Calculate time steps
  return(Rescaled)
}




#' Calculate bed shear stress at a location
#'
#' This function takes tidal data from SINMOD and wave data from ECMWF, aligns the two time series at a single location,
#' and then calculates the bed shear stress through time.
#'
#' @pixel a numeric value indicating the row number of a dataframe of target coordinates.
#' @tide_packets a list of dataframes containing instructions for processing all SINMOD files which share a year.
#' @return A dataframe containing a time series of bed shear stress metrics for a location.
#' @family water movements
#' @export
get_bed_shear_stress <- function(pixel, tide_packets, depth = 40) {

  # pixel <- 20000
  # tide_packets <- Instructions[1:2]

  ## PULL DATA ##

  wave_ts <- sample_waves(pixel)                                   # Extract wave time series
  wave_ts <- zoo::zoo(wave_ts[, colnames(wave_ts) != "Wave_step"],
                      wave_ts[, colnames(wave_ts) == "Wave_step"]) # Convert to zoo object to align time series
  tic()
  tide_ts <- purrr::map(tide_packets, sample_tides_year, pixel) %>%# Repeat sampling of tides for each year
    do.call(rbind, .)                                              # Bind the years together quickly
  tide_ts <- zoo::zoo(tide_ts[, colnames(tide_ts) != "Tide_step"],
                      tide_ts[, colnames(tide_ts) == "Tide_step"]) # Convert to zoo object to align time series
  toc()
  ## ALIGN TIME SERIES ##

  align <- merge(tide_ts, wave_ts)                                 # Match up time series (identify missing steps)

  align$mwd <- zoo::na.approx(align$mwd, rule = 2)                 # Interpolate wave direction onto tide time steps

  align$swh <- zoo::na.approx(align$mwd, rule = 2)                 # Interpolate wave height onto tide time steps

  align$mwp <- zoo::na.approx(align$mwd, rule = 2)                 # Interpolate wave period onto tide time steps

  align <- align[which(zoo::index(align) %in% zoo::index(tide_ts)),] # Keep only tide time steps

  ## CALCULATE BED SHEAR STRESS ##

  stress <- bedshear::shear_stress(
    bathymetry = depth,                                            # Depth to sea floor
    D50 = 0.000,                                                   # Nominal median grain size
    tidal_velocity = align[,"uvSpeed"],                            # Water movements
    tidal_direction = align[,"uvDirection"],
    wave_height = align[,"swh"],
    wave_period = align[,"mwp"],
    wave_direction = align[,"mwd"],
    switch = 0) %>%
    dplyr::mutate(Pixel = pixel,                                   # Add pixel ID
                  Time_Step = zoo::index(align))                   # Add time

  return(stress)

}

