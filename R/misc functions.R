#' Run an R Script in it's own Session Informatively
#'
#' This function is a wrapper for a few tweaks around sourcing an R script.
#'
#' The R script is run in it's own R session, to avoid crashing an R session handling the batch processing of scripts.
#' Activity is indicated with a spinner, and completed scripts are announced to help localise errors to a script. The time
#' taken for a script to complete is also tracked in a log file, which can be accessed with `tictoc::tic.log()`.
#'
#' @param filepath a character string indicating the R script you want to execute.
#' @return The R script is run for it's side effects and the run time added to a log by `tictoc`.
#' @export
execute <- function(filepath) {
  tictoc::tic(filepath)
  callr::r(function(filepath) source(filepath), args = list(filepath), spinner = TRUE)
  tictoc::toc(log = T, quiet = T)

  usethis::ui_done("{usethis::ui_field(filepath)} completed. {praise::praise('${Exclamation}!')}")}

#' Fill missing values in a vector
#'
#' This function takes a vector containing missing values ("" or NA), and overwrites empty positions with
#' the previous non-missing entry. This assumes the first value is not a missing value.
#'
#' The motivating example for this was to automatically fill in cells in BODC nutrient data where metadata
#' is not copied past the first row of a block of data representing a single CTD cast.
#'
#' @param data a vector containing "" or NA as missing values.
#' @return A filled in vector where "" or NA is replaced with the latest non-missing value.
#' @examples
#' # Create dummy data with missing entries
#' test <- c(1, 2, "", "", NA, 3, "")
#'
#' # Replace blanks
#' filled <- fill_in(test)
#' @export
fill_in <- function(data) {

  for (i in 2:length(data)) {               # From the second value onwards

    if(data[i] == "" | is.na(data[i])) data[i] <- data[i-1]  # If blank, overwrite with the previous value, otherwise do nothing
  }
return(data)}

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
#' @param shift True/False does the file use 0-360 longitudes? If yes, negative longitudes are corrected.
#' @param sf An sf polygon can be supplied instead of `w`,`e`,`s`, and `n`. The function takes these values from the bounding box.
#' @return A list of three elements:
#' \itemize{
#'  \item{\emph{Lats -}}{ A vector of latitudes from `s` to `n`.}
#'  \item{\emph{Lons -}}{ A vector of longitudes from `w` to `e`.}
#'  \item{\emph{Limits -}}{ A dataframe containing the index to start reading from (Lon_start, Lat_start)
#'  and the length of the vector to read (Lon_count, Lat_count.}
#'  }
#' @family NEMO-MEDUSA spatial tools
#' @export
Window <- function(file, w = NULL, e = NULL, s = NULL, n = NULL, shift = FALSE, sf = NULL){

  #file <- examples[1,]$File ; w = 0 ; e = 180 ; s = 0 ; n = 90

if (!is.null(sf)) {

  box <- st_transform(sf, 4326) %>%
    st_bbox()

  w <- box[1]; e <- box[3]; s <- box[2]; n <- box[4]

}

if (isTRUE(shift)) {        # If the data source uses longitudes 0-360 correct values

  if (w < 0) {w <- w+360}

  if (e < 0) {e <- e+360}

}

  if (w > e) {                              # Catch when w is larger than e which breaks between()
    if (isTRUE(shift)) { w <- 0 ; e <- 360 ; warning("W is larger than E. This may happen when crossing 0 longitude. Defaulting to all longitudes")} else{
      if (isFALSE(shift)) w <- -180 ; e <- 180 ; warning("W is larger than E. This may happen when crossing 0 longitude. Defaulting to all longitudes")}
  }

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

#' Zoom a ggplot to the bounding box of an sf object
#'
#' This function can be used with ggplot and uses an sf object to adjust the plotting window according to a bounding box.
#' If you want more breathing space use st_buffer when providing `sf`.
#'
#' @param sf an sf object that you want to crop the plot window to
#' @return Add to a ggplot with the standard `+` syntax to adjust the plotting window.
#' @export
box_zoom <- function(sf) {

  zoom <- sf::st_bbox(sf)
  ggplot2::coord_sf(xlim = zoom[c(1,3)], ylim = zoom[c(2,4)])

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
reproj <- function(data, crs) {

  data %>%
    sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% # Specify original projection (crs)
    sf::st_transform(crs = crs) %>%                                   # Transform to crs specified in region file
    sfc_as_cols() %>%                                                 # Extract geometry column for geom_segment to work
    sf::st_set_geometry(NULL)                                         # Chuck geometry column
}

#' Re-parameterise a Driver File for a Time Period in the MiMeMo Domain
#'
#' These function take a StrathE2E model and update the parameters in the driving data files to a specific time window.
#'
#' The functions are informative. Warnings are produced to highlight which parameters don't have data for the target time
#' window. In this case the old parameters are left in place. The function will refuse to update parameters based on data
#' which covers less than half of the target time interval. For parameters which are updated, the function returns the
#' years used.
#'
#' @param start The first year in the target time period.
#' @param last The last year in the target time period.
#' @param path A character string containing the path to a model variant, i.e. "./Models/Region/start-last".
#' @return Run this function for it's side effects. The appropriate driving file will be updated.
#' @name Update-drivers
NULL

#' @rdname Update-drivers
#' @export
update_boundary_period <- function(start, end, path){

  copied <- list.files(paste0(path,"/Driving"), full.names = T) %>%                            # List the driving files in the model
    .[stringr::str_detect(., "chemistry")]                                                     # Select the chemistry file

  Boundary_template <-  read.csv(copied)                                                       # Read in old boundary drivers

  My_DIN_fix <- readRDS("./Objects/Ammonia to DIN.rds")
  usethis::ui_warn("The correction from DIN to NH[4]and NO[3] comes from a fixed time period.")

  My_river_N <- readRDS("./Objects/River nitrate and ammonia.rds")

  if (nrow(My_river_N) != 12) {                                                                # If not already summarised
    My_river_N <- filter(My_river_N, between(Year, start, end)) %>%                            # Limit to reference period
      group_by(Month) %>%                                                                      # Average across years
      summarise(Ammonia = mean(Ammonia, na.rm = T),
                Nitrate = mean(Nitrate, na.rm = T)) %>%
      ungroup() %>%
      arrange(Month)                                                                           # Order months ascending
  }

  My_river_N <- dplyr::mutate(My_river_N, Ammonia = (Ammonia*(1/14.006720))*1e3,               # Convert mg/l to mmol/m^3
                  Nitrate = (Nitrate*(1/14.006720))*1e3)
  usethis::ui_warn("River nutrient concentrations come from a fixed time period.")

  #### With unchanging data sources ####

  Boundary_new <- dplyr::mutate(Boundary_template,
                                ## Rivers
                                RIV_nitrate = My_river_N$Nitrate,
                                RIV_ammonia = My_river_N$Ammonia,
                                RIV_detritus = 0)

  #### Update NEMO-MEDUSA data ####

  My_boundary_data <- readRDS("./Objects/Boundary measurements.rds") %>%                        # Import data
    dplyr::filter(between(Year, start, end))

  if(length(unique(My_boundary_data$Year)) > (end-start+1)*0.5){

    usethis::ui_info("Updating using NEMO-MEDUSA outputs from {usethis::ui_value(min(My_boundary_data$Year))} to {usethis::ui_value(max(My_boundary_data$Year))}.")

    My_boundary_data <- My_boundary_data %>%                                                        # Limit to reference period
      dplyr::group_by(Month, Compartment, Variable) %>%                                           # Average across years
      dplyr::summarise(Measured = mean(Measured, na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(Month) %>%                                                                   # Order months ascending
      dplyr::mutate(Compartment = factor(Compartment, levels = c("Inshore S", "Offshore S", "Offshore D"),
                                         labels = c("Inshore S" = "SI", "Offshore S" = "SO", "Offshore D" = "D"))) %>%
      tidyr::pivot_wider(names_from = c(Compartment, Variable), names_sep = "_", values_from = Measured) # Spread columns to match template

    Boundary_new <- dplyr::mutate(Boundary_new,
                                  SO_nitrate = My_boundary_data$SO_DIN * (1-dplyr::filter(My_DIN_fix, Depth_layer == "Shallow")$Proportion), # Multiply DIN by the proportion of total DIN as nitrate
                                  SO_ammonia = My_boundary_data$SO_DIN * dplyr::filter(My_DIN_fix, Depth_layer == "Shallow")$Proportion, # Multiply DIN by the proportion of total DIN as ammonium
                                  SO_phyt = My_boundary_data$SO_Phytoplankton,
                                  SO_detritus = My_boundary_data$SO_Detritus,
                                  D_nitrate = My_boundary_data$D_DIN * (1-dplyr::filter(My_DIN_fix, Depth_layer == "Deep")$Proportion), # Multiply DIN by the proportion of total DIN as nitrate
                                  D_ammonia = My_boundary_data$D_DIN * dplyr::filter(My_DIN_fix, Depth_layer == "Deep")$Proportion, # Multiply DIN by the proportion of total DIN as ammonium
                                  D_phyt = My_boundary_data$D_Phytoplankton,
                                  D_detritus = My_boundary_data$D_Detritus,
                                  SI_nitrate = My_boundary_data$SI_DIN * (1-dplyr::filter(My_DIN_fix, Depth_layer == "Shallow")$Proportion), # Multiply DIN by the proportion of total DIN as nitrate
                                  SI_ammonia = My_boundary_data$SI_DIN * dplyr::filter(My_DIN_fix, Depth_layer == "Shallow")$Proportion, # Multiply DIN by the proportion of total DIN as ammonium
                                  SI_phyt = My_boundary_data$SI_Phytoplankton,
                                  SI_detritus = My_boundary_data$SI_Detritus)

  } else {
    usethis::ui_warn("Did not update from NEMO-MEDUSA; fewer than half the target years are represented.")}

  My_atmosphere <- readRDS("./Objects/Atmospheric N deposition.rds") %>%
    filter(between(Year, start, end))                                   # Limit to reference period

  if(length(unique(My_atmosphere$Year)) > (end-start+1)*0.5){

    usethis::ui_info("Updating atmosphere using data from {usethis::ui_value(min(My_atmosphere$Year))} to {usethis::ui_value(max(My_atmosphere$Year))}.")

    My_atmosphere <-  My_atmosphere %>%
      dplyr::group_by(Month, Oxidation_state, Shore,  Year) %>%
      dplyr::summarise(Measured = sum(Measured, na.rm = T)) %>%                                  # Sum across deposition states
      dplyr::summarise(Measured = mean(Measured, na.rm = T)) %>%                                 # Average over years
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = c(Shore, Oxidation_state), values_from = Measured) %>%     # Spread to match template
      dplyr::arrange(Month)                                                                      # Order months ascending

    Boundary_new <- dplyr::mutate(Boundary_new,
                                  ## Atmosphere, daily deposition as monthly averages
                                  SO_ATM_nitrate_flux = My_atmosphere$Offshore_O,
                                  SO_ATM_ammonia_flux = My_atmosphere$Offshore_R,
                                  SI_ATM_nitrate_flux = My_atmosphere$Inshore_O,
                                  SI_ATM_ammonia_flux = My_atmosphere$Inshore_R,
                                  SI_other_nitrate_flux = 0,   # Can be used for scenarios
                                  SI_other_ammonia_flux = 0)} else {
                                    usethis::ui_warn("Did not update atmosphere; fewer than half the target years are represented.")}

  new <- stringr::str_split(path, "Models/")[[1]][2] %>%                 # Pull the text we need from the file path to name the new file
    stringr::str_replace("[[:punct:]]", "_") %>%                        # Remove the '/'
    stringr::str_replace(" ", "_") %>%                                  # Replace any spaces with '_'
    toupper() %>%                                                       # Capitalise
    paste0(path, "/Driving/chemistry_", ., ".csv")                      # Build new name

  write.csv(Boundary_new, file = new, row.names = F)                    # Save the new file
  unlink(copied)                                                        # Delete the old one

}

#' @rdname Update-drivers
#' @export
update_physics_period <- function(start, end, path){

  copied <- list.files(paste0(path,"/Driving"), full.names = T) %>%           # List the driving files in the model
    .[stringr::str_detect(., "physics")]                                      # Select the physics file

  Physics_template <-  read.csv(copied)                                       # Read in old physics drivers

  #### Static parameters ####

  My_scale <- readRDS("./Objects/Domains.rds") %>%                            # Calculate the volume of the three zones
    sf::st_drop_geometry() %>%
    dplyr::mutate(S = c(T, T),
                  D = c(F, T)) %>%
    tidyr::gather(key = "slab_layer", value = "Exists", S, D) %>%
    dplyr::filter(Exists == T) %>%
    dplyr::mutate(Elevation = c(Elevation[1], -60, Elevation[3] + 60)) %>%
    dplyr::mutate(Volume = area * abs(Elevation)) %>%
    dplyr::select(Shore, slab_layer, Volume)

  My_Waves <- readRDS("./Objects/Significant wave height.rds") %>%            #*2000 - 2010
    dplyr::arrange(Month)                                                     # Arrange to match template
  usethis::ui_warn("Significant wave height estimates come from a fixed time period.")

  My_Stress <- readRDS("./Objects/Habitat disturbance.rds") %>%
    dplyr::mutate(Month = factor(Month, levels = month.name)) %>%             # Set month as a factor for non-alphabetical ordering
    dplyr::arrange(Month)                                                     # Arrange to match template
  usethis::ui_warn("Habitat disturbance estimates come from a fixed time period.")

  Physics_new <- dplyr::mutate(Physics_template,
                               ## Daily proportion disturbed by natural bed shear stress
                               habS1_pdist = dplyr::filter(My_Stress, Shore == "Inshore", Habitat == "Silt")$Disturbance,
                               habS2_pdist = dplyr::filter(My_Stress, Shore == "Inshore", Habitat == "Sand")$Disturbance,
                               habS3_pdist = dplyr::filter(My_Stress, Shore == "Inshore", Habitat == "Gravel")$Disturbance,
                               habD1_pdist = dplyr::filter(My_Stress, Shore == "Offshore", Habitat == "Silt")$Disturbance,
                               habD2_pdist = dplyr::filter(My_Stress, Shore == "Offshore", Habitat == "Sand")$Disturbance,
                               habD3_pdist = dplyr::filter(My_Stress, Shore == "Offshore", Habitat == "Gravel")$Disturbance,
                               ## Monthly mean significant wave height inshore
                               Inshore_waveheight = My_Waves$Waves)

  #### Update light ####

  My_light <- readRDS("./Objects/Air temp and light.rds") %>%
    dplyr::filter(between(Year, start, end), grepl("Light", Type))              # Limit to reference period and variable

  if(length(unique(My_light$Year)) > (end-start+1)*0.5){

    usethis::ui_info("Updating surface irradiance using data from {usethis::ui_value(min(My_light$Year))} to {usethis::ui_value(max(My_light$Year))}.")

    My_light <-  My_light %>%
      dplyr::group_by(Month) %>%                                                  # Average across months
      dplyr::summarise(Measured = mean(Measured, na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(Month)                                                       # Order to match template

    Physics_new <- dplyr::mutate(Physics_new,
                                 SLight = My_light$Measured)} else {
                                   usethis::ui_warn("Did not update surface irradiance; fewer than half the target years are represented.")}

  #### Update air temperature ####

  My_AirTemp <- readRDS("./Objects/Air temp and light.rds") %>%
    dplyr::filter(between(Year, start, end), grepl("Air", Type))                # Limit to reference period and variable

  if(length(unique(My_AirTemp$Year)) > (end-start+1)*0.5){

    usethis::ui_info("Updating surface air temperature using data from {usethis::ui_value(min(My_AirTemp$Year))} to {usethis::ui_value(max(My_AirTemp$Year))}.")

    My_AirTemp <- My_AirTemp %>%
      dplyr::group_by(Month, Shore) %>%                                           # Average across months
      dplyr::summarise(Measured = mean(Measured, na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(Month)                                                       # Order to match template

    Physics_new <- dplyr::mutate(Physics_new,
                                 SO_AirTemp = dplyr::filter(My_AirTemp, Shore == "Offshore")$Measured,
                                 SI_AirTemp = dplyr::filter(My_AirTemp, Shore == "Inshore")$Measured)
  } else {
    usethis::ui_warn("Did not update surface air temperature; fewer than half the target years are represented.")}

  #### Update horizontal water exchange ####

  My_H_Flows <- readRDS("./Objects/H-Flows.rds") %>%
    dplyr::filter(between(Year, start, end))                                    # Limit to reference period

  if(length(unique(My_H_Flows$Year)) > (end-start+1)*0.5){

    usethis::ui_info("Updating horizontal flows using data from {usethis::ui_value(min(My_H_Flows$Year))} to {usethis::ui_value(max(My_H_Flows$Year))}.")

    My_H_Flows <- My_H_Flows %>%
      dplyr::group_by(dplyr::across(-c(Year, Flow))) %>%                          # Group over everything except year and variable of interest
      dplyr::summarise(Flow = mean(Flow, na.rm = T)) %>%                          # Average flows by month over years
      dplyr::ungroup() %>%
      dplyr::left_join(My_scale) %>%                                              # Attach compartment volumes
      dplyr::mutate(Flow = Flow/Volume) %>%                                       # Scale flows by compartment volume
      dplyr::mutate(Flow = abs(Flow * 86400)) %>%                                 # Multiply for total daily from per second, and correct sign for "out" flows
      dplyr::arrange(Month)                                                       # Order by month to match template

    Physics_new <- dplyr::mutate(Physics_new,
                                 ## Flows, should be proportions of volume per day
                                 SO_OceanIN = dplyr::filter(My_H_Flows, slab_layer == "S", Shore == "Offshore", Neighbour == "Ocean", Direction == "In")$Flow,
                                 D_OceanIN = dplyr::filter(My_H_Flows, slab_layer == "D", Shore == "Offshore", Neighbour == "Ocean", Direction == "In")$Flow,
                                 SI_OceanIN = dplyr::filter(My_H_Flows, slab_layer == "S", Shore == "Inshore", Neighbour == "Ocean", Direction == "In")$Flow,
                                 SI_OceanOUT = dplyr::filter(My_H_Flows, slab_layer == "S", Shore == "Inshore", Neighbour == "Ocean", Direction == "Out")$Flow,
                                 SO_SI_flow = dplyr::filter(My_H_Flows, slab_layer == "S", Shore == "Offshore", Neighbour == "Inshore", Direction == "Out")$Flow)
  } else {
    usethis::ui_warn("Did not update horizontal flows; fewer than half the target years are represented.")}

  #### Updating vertical water exchanges ####

  My_V_Flows <- readRDS("./Objects/vertical diffusivity.rds") %>%
    dplyr::filter(between(Year, start, end))                                    # Limit to reference period

  if(length(unique(My_V_Flows$Year)) > (end-start+1)*0.5){

    usethis::ui_info("Updating vertical diffusivity using data from {usethis::ui_value(min(My_V_Flows$Year))} to {usethis::ui_value(max(My_V_Flows$Year))}.")

    My_V_Flows <- My_V_Flows %>%
      dplyr::group_by(Month) %>%
      dplyr::summarise(V_diff = mean(Vertical_diffusivity, na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(Month)                                                       # Order by month to match template

    Physics_new <- dplyr::mutate(Physics_new,
                                 ## Vertical diffusivity
                                 log10Kvert = log10(My_V_Flows$V_diff))
  } else {
    usethis::ui_warn("Did not update vertical diffusivity; fewer than half the target years are represented.")}

  #### Update other volume based values ####

  My_volumes <- readRDS("./Objects/TS.rds") %>%
    dplyr::filter(between(Year, start, end))                                    # Limit to reference period

  if(length(unique(My_volumes$Year)) > (end-start+1)*0.5){

    usethis::ui_info("Updating water temperatures and ice (NEMO-MEDUSA) using data from {usethis::ui_value(min(My_volumes$Year))} to {usethis::ui_value(max(My_volumes$Year))}.")

    My_volumes <- My_volumes %>%
      dplyr::group_by(Compartment, Month) %>%                                     # By compartment and month
      dplyr::summarise(dplyr::across(Salinity_avg:Ice_conc_avg, mean, na.rm = T)) %>% # Average across years for multiple columns
      dplyr::ungroup() %>%
      dplyr::arrange(Month)                                                       # Order by month to match template

    Physics_new <- dplyr::mutate(Physics_new,
                                 ## Temperatures in volumes for each zone
                                 SO_temp = dplyr::filter(My_volumes, Compartment == "Offshore S")$Temperature_avg,
                                 D_temp =dplyr:: filter(My_volumes, Compartment == "Offshore D")$Temperature_avg,
                                 SI_temp = dplyr::filter(My_volumes, Compartment == "Inshore S")$Temperature_avg ,
                                 ## Cryo variables
                                 SO_IceFree = 1 - dplyr::filter(My_volumes, Compartment == "Offshore S")$Ice_pres,
                                 SI_IceFree = 1 - dplyr::filter(My_volumes, Compartment == "Inshore S")$Ice_pres,
                                 SO_IceCover = dplyr::filter(My_volumes, Compartment == "Offshore S")$Ice_conc_avg,
                                 SI_IceCover = dplyr::filter(My_volumes, Compartment == "Inshore S")$Ice_conc_avg,
                                 SO_IceThickness = dplyr::filter(My_volumes, Compartment == "Offshore S")$Ice_Thickness_avg,
                                 SI_IceThickness = dplyr::filter(My_volumes, Compartment == "Inshore S")$Ice_Thickness_avg,
                                 SO_SnowThickness = dplyr::filter(My_volumes, Compartment == "Offshore S")$Snow_Thickness_avg,
                                 SI_SnowThickness = dplyr::filter(My_volumes, Compartment == "Inshore S")$Snow_Thickness_avg)
  } else {
    usethis::ui_warn("Did not update water temperatures and ice (NEMO-MEDUSA); fewer than half the target years are represented.")}

  #### Update suspended particulate matter ####

  My_SPM <- readRDS("./Objects/Suspended particulate matter.rds") %>%
    dplyr::filter(between(Year, start, end))                                    # Limit to reference period

  if(length(unique(My_SPM$Year)) > (end-start+1)*0.5){

    usethis::ui_info("Updating suspended particulate matter using data from {usethis::ui_value(min(My_SPM$Year))} to {usethis::ui_value(max(My_SPM$Year))}.")

    My_SPM <- My_SPM %>%
      dplyr::group_by(Shore, Month) %>%
      dplyr::summarise(SPM = mean(SPM, na.rm = T)) %>%                            # Average by month across years
      dplyr::ungroup() %>%
      dplyr::arrange(Month)                                                       # Order by month to match template

    Physics_new <- dplyr::mutate(Physics_new,
                                 ## log e transformed suspended particulate matter concentration in zones
                                 SO_LogeSPM = log(dplyr::filter(My_SPM, Shore == "Offshore")$SPM),
                                 SI_LogeSPM = log(dplyr::filter(My_SPM, Shore == "Inshore")$SPM))
  } else {
    usethis::ui_warn("Did not update suspended particulate matter; fewer than half the target years are represented.")}

  #### Update river volumes ####

  My_Rivers <- readRDS("./Objects/River volume input.rds") %>%
    dplyr::filter(between(Year, start, end))                                    # Limit to reference period

  if(length(unique(My_Rivers$Year)) > (end-start+1)*0.5){

    usethis::ui_info("Updating river outflows using data from {usethis::ui_value(min(My_Rivers$Year))} to {usethis::ui_value(max(My_Rivers$Year))}.")

    My_Rivers <- My_Rivers %>%
      dplyr::group_by(Month) %>%
      dplyr::summarise(Runoff = mean(Runoff, na.rm = T)) %>%                      # Average by month across years
      dplyr::ungroup() %>%
      dplyr::arrange(as.numeric(Month))                                           # Order by month to match template

    Physics_new <- dplyr::mutate(Physics_new,
                                 ## River inflow,
                                 Rivervol_SI = My_Rivers$Runoff / dplyr::filter(My_scale, Shore == "Inshore")$Volume) # Scale as proportion of inshore volume
  } else {
    usethis::ui_warn("Did not update river outflows; fewer than half the target years are represented.")}

  #### Create new file ####

  new <- stringr::str_split(path, "Models/")[[1]][2] %>%                # Pull the text we need from the file path to name the new file
    stringr::str_replace("[[:punct:]]", "_") %>%                        # Remove the '/'
    stringr::str_replace(" ", "_") %>%                                  # Replace any spaces with '_'
    toupper() %>%                                                       # Capitalise
    paste0(path, "/Driving/physics_", ., ".csv")                        # Build new name

  write.csv(Physics_new, file = new, row.names = F)                     # Save the new file
  unlink(copied)                                                        # Delete the old one

}
