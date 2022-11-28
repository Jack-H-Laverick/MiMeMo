
#### Water movement tools    ####

#' Subset a SINMOD file to depth-averaged water velocities in an area.
#'
#' This function is a wrapper for system calls to NCO functions. A single SINMOD netcdf file contains a
#' large grid covering the Arctic with extra variables we don't need. This function crops the spatial extent,
#' performs a weighted average across depth layers, and limits to variables `u_east` and `v_north`.
#'
#' Performing the reshaping directly on the netcdf files means I can host a small subset locally on my machine. This avoids
#' costly queries to the idrive later on. These files can also then be concatenated into a single netcdf file, allowing data
#' to be read into R pixel by pixel, unlike saving a clipped dataframe as a .rds object.
#'
#' @param path a character string directing to the source file location.
#' @param file the name of a netcdf file to process.
#' @param out_dir a character string directing to the folder the new files should be written to.
#' @param window a named list containing `xmin`, `xmax`, `ymin`, `ymax`, the indices of the window to crop to.
#' @return A netcdf file containing depth averaged UV water velocities cropped to the target window.
#' @family water movements
#' @export
reshape_SINMOD <- function(path, file, out_dir, window) {

  #path <- unique(all_files$path) ; file <- all_files$file[1]               # test

  ff_new <- paste0(out_dir, "SINMOD_", file)                             # Combine file name and new path

  if(!file.exists(ff_new)) {                                             # If the file doesn't exist, make a new one
    ff_dir <- dirname(ff_new)                                            # Get the directory of the new file
    if(!dir.exists(ff_dir))                                              # If that folder doesn't exist
      dir.create(ff_dir, recursive = TRUE)                               # Create the folder

    temp_file1 <- tempfile("dummy", tmpdir = ff_dir)
    temp_file2 <- tempfile("dummy", tmpdir = ff_dir)

    file.copy(paste0(path, file), ff_new)                                # Copy the i drive file to the new name/location

    ## Run command line functions
    # str_glue and {} allows us to programmatically buils character strings to run in Konsole.
    # Keep only these variables (-x -v means drop these variables),clip x and y dimensions
    # (names xc and yc specified in the file), from a file, saving result as new file.
    system(stringr::str_glue("ncea -v u_east,v_north,LayerDepths -d xc,{window$xmin},{window$xmax} -d yc,{window$ymin},{window$ymax} {ff_new} {temp_file1}"))

    # Weighted averages, over depth (zc), using layer depths, from a file, saving result
    system(stringr::str_glue("ncwa -a zc -w LayerDepths {temp_file1} {temp_file2}"))

    file.rename(temp_file2, ff_new)                                       # Name sensibly

    unlink(temp_file1)                                                    # Delete the intermediate step

  }
  usethis::ui_done("{usethis::ui_field(file)} cropped.")      # Announce finished file

}

#' Calculate bed shear stress at a location
#'
#' This function takes tidal data from SINMOD and wave data from ECMWF, aligns the two time series at a single location,
#' and then calculates the bed shear stress through time.
#'
#' @param pixel a numeric value indicating the row number of a dataframe of target coordinates.
#' @param tide_packets a list of dataframes containing instructions for processing all SINMOD files which share a year.
#' @return A dataframe containing a time series of bed shear stress metrics for a location.
#' @family water movements
#' @export
#' @keyword internal
mync_get_pixel <- function(File, var, pixel_x, pixel_y) {

  ncvar_get(nc_open(File), var, c(pixel_x, pixel_y, 1, 1),  # Extract the variable of interest
            c(1, 1, -1, -1)) %>%                            # cropped to a target pixel, with all time steps and ensemble members
    colSums(na.rm = TRUE)                                   # Average across ensemble quickly
  # Setting na.rm = TRUE replaces iced NAs with 0s
}         # Import a variable clipped to Window

#' Calculate bed shear stress at a location
#'
#' This function takes tidal data from SINMOD and wave data from ECMWF, aligns the two time series at a single location,
#' and then calculates the bed shear stress through time.
#'
#' @param pixel a numeric value indicating the row number of a dataframe of target coordinates.
#' @param tide_packets a list of dataframes containing instructions for processing all SINMOD files which share a year.
#' @return A dataframe containing a time series of bed shear stress metrics for a location.
#' @family water movements
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

#' Calculate bed shear stress at a location
#'
#' This function takes tidal data from SINMOD and wave data from ECMWF, aligns the two time series at a single location,
#' and then calculates the bed shear stress through time.
#'
#' @param pixel a numeric value indicating the row number of a dataframe of target coordinates.
#' @param tide_packets a list of dataframes containing instructions for processing all SINMOD files which share a year.
#' @return A dataframe containing a time series of bed shear stress metrics for a location.
#' @family water movements
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

#' Calculate bed shear stress at a location
#'
#' This function takes tidal data from SINMOD and wave data from ECMWF, aligns the two time series at a single location,
#' and then calculates the bed shear stress through time.
#'
#' @param pixel a numeric value indicating the row number of a dataframe of target coordinates.
#' @param tide_packets a list of dataframes containing instructions for processing all SINMOD files which share a year.
#' @return A dataframe containing a time series of bed shear stress metrics for a location.
#' @family water movements
#' @export
get_tides_ts <- function(File, month, vars, pixel_x, pixel_y) {

  #File <- paste0(tide_files$path[1], tide_files$file[1]) ; Year <- tide_files$year[1] # test
  #pixel_x <- aligned$V1[1] ; pixel_y <- aligned$V2[1]                     # test
  #vars <- c("u_east", "v_north")                                                 # test

  Data <- map(vars, ~{mync_get_pixel(File, var = .x, pixel_x, pixel_y)})          # Extract the variables of interest
  Rescaled <- vectors_2_direction(Data[[1]], Data[[2]])    #*** warning, averaged across depth but not weighted right now by thickness.

  return(Rescaled)

}   # Pull the depth averaged currents from a file

#' Calculate bed shear stress at a location
#'
#' This function takes tidal data from SINMOD and wave data from ECMWF, aligns the two time series at a single location,
#' and then calculates the bed shear stress through time.
#'
#' @param pixel a numeric value indicating the row number of a dataframe of target coordinates.
#' @param tide_packets a list of dataframes containing instructions for processing all SINMOD files which share a year.
#' @return A dataframe containing a time series of bed shear stress metrics for a location.
#' @family water movements
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
#' @param pixel a numeric value indicating the row number of a dataframe of target coordinates.
#' @param tide_packets a list of dataframes containing instructions for processing all SINMOD files which share a year.
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

#' Calculate sediment permeability
#'
#' This function takes a percent of mud/silt at a location and returns a measure of permeability. The default
#' parameterisation comes from Matt Pace's thesis for fine sediments.
#'
#' @param percent_mud The amount of all sediment at a location smaller than sand (mud/silt and finer) in percent.
#' @param scalar A parameter in the relationship
#' @param constant A parameter in the relationship.
#' @return A numeric vector containing estimates of sediment permeability in m^2 (can be length 1).
#' @family Sediment properties
#' @examples
#' data <- data.frame(Silt = seq(0.5, 100, length.out = 100)) %>%
#' mutate(Permeability = mud_to_permeability(Silt))
#'
#' ggplot(data) +
#' geom_line(aes(x = Silt, y = Permeability)) +
#' scale_y_continuous(trans = "log", breaks = c(0.000000000000001, 0.0000000000001, 0.00000000001))
#' @export
mud_to_permeability <- function(percent_mud, scalar = -2.171, constant = -10.232) {

  permeability <- 10^((scalar*log10(percent_mud)) + constant)
  return(permeability)
}

#' Calculate sediment porosity
#'
#' This function takes median grain size at a location and returns a measure of porosity. The default
#' parameterisation comes from Matt Pace's thesis for fine sediments. The numerical order of parameters
#' has been changed from Matt's to be consistent with StrathE2E.
#'
#' @param D50 The median grain size at a location.
#' @param p1 A parameter in the relationship.
#' @param p2 A parameter in the relationship.
#' @param p3 A parameter in the relationship.
#' @param p4 A parameter in the relationship.
#' @return A numeric vector containing estimates of sediment porosity in % (can be length 1).
#' @family Sediment properties
#' @examples
#' data <- data.frame(D50 = seq(0, 1, length.out = 1000)) %>%
#' mutate(Porosity = D50_to_porosity(D50))
#'
#' ggplot(data) +
#' geom_line(aes(x = D50, y = Porosity)) +
#' scale_x_continuous(trans = "log", breaks = c(0.001,0.01,0.1,1))
#' @export
D50_to_porosity <- function(D50, p1 = -1.035, p2 = -0.314, p3 = -0.435, p4 = 0.302) {

  complex <- 1+exp((-(log10(D50)-p1))/p2)

  porosity <- p3 + (p4*(1/(complex)))

  answer <- 10^porosity

  return(answer)
}
