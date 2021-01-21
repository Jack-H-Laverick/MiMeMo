
#### Functions to average volumes once extracted ####

#' Prepare for Averaging by Decade
#'
#' This function cleans the saved NEMO-MEDUSA monthly summaries, for averaging into decades.
#'
#' @param saved A dataframe containing a summarised month from NEMO-MEDUSA model outputs.
#' @return A dataframe containing a summarised month of NEMO-MEDUSA output, gaining a decade column, and dropping columns
#' which aren't needed for spatial maps.
#' @family NEMO-MEDUSA averages
#' @export
decadal <- function(saved) {

  import <- readRDS(file = saved) %>%                                   # Read in wide format data file
    dplyr::select(-c(weights, Bathymetry)) %>%                          
    dplyr::rename(Decade = Year)

  stringr::str_sub(import$Decade, -1, -1) <- "0"                        # Overwite the 4th digit with a 0 to get the decade

  return(import)
}

#' Strip Snow and Ice Variables at Depth
#'
#' This function removes the snow and ice columns from a dataframe if depth = "D".
#'
#' Some variables are only relevant in the shallow zone of StrathE2E polar. There is no sea-ice 60 m under the sea.
#' This means, when dataframes containing both shallow and deep data are split by depth, empty columns can be introduced.
#' These empty columns can cause problems in downstream functions, such as plotting by column. This function removes the
#' empty columns.
#'
#' @param data A dataframe containing a summarised month from NEMO-MEDUSA model outputs, at a single depth.
#' @param dt Switch for using either data.table or dplyr methods (TRUE/FALSE respectively)
#' @return If the `data` contains shallow data, no action is taken. If `data` contains deep data, columns for variables only
#' relevant in the shallow zone are dropped.
#' @family NEMO-MEDUSA averages
#' @export
strip_ice <- function(data, dt) {

  if(dt == TRUE) {                                                                  # Run data.table method
    data <- data.table::setDT(data)

    if(data$Depth[1] == "D") data[, c("Ice_conc", "Ice_Thickness", "Snow_Thickness"):=NULL] else data
  } else{                                                                           # Run dplyr method

    if(data$Depth[1] == "D") {select(data, -c(starts_with("Ice"), Snow_Thickness))} else data}}

#' Average into Decadal Grids
#'
#' This function averages cleaned NEMO-MEDUSA monthly summaries into decadal grids.
#'
#' The function groups by all spatial variables (Longitude, Latitude, Depth, and Shore zone), and by decade and month.
#' The mean for every other variable is calculated within these groups.
#'
#' @param saved A dataframe containing a summarised month from NEMO-MEDUSA model outputs. It must contain the columns:
#' Longitude, Latitude, Decade, Month, Shore, and Depth.
#' @param dt Switch for using either data.table or dplyr methods (TRUE/FALSE respectively)
#' @return A dataframe containing a summarised decade of spatialy resolved NEMO-MEDUSA outputs.
#' @family NEMO-MEDUSA averages
#' @export
summarise_sp <- function(decade, dt) {

if(dt == TRUE){                                                               # Run data.table method
  # data.table::setDT(decade)                                                 # set as a data.table, not needed if decade is already a data.table
  Averaged <- decade[, lapply(.SD, mean, na.rm = TRUE),                       # Average data columns which aren't groups
                     by = c("Longitude", "Latitude", "Decade", "Month", "Shore", "Depth")] # Group by pixel and decade
} else{                                                                       # Run dplyr method
  Averaged <- decade %>%
    dplyr::group_by(Longitude, Latitude, Decade, Month, Shore, Depth) %>%     # Group by pixel and decade
    dplyr::summarise_all(mean, na.rm = TRUE) %>%                              # Average data columns
    dplyr::ungroup()                                                          # Ungroup
  }
  return(Averaged)
}

#' Average into Time Series
#'
#' This function averages NEMO-MEDUSA monthly summaries into time series for each model compartment.
#'
#' The function groups by model compartment (Depth and Shore zone) and time step (Month and Year).
#' The mean for every target variable is calculated within these groups.
#'
#' @param saved A dataframe containing a summarised month from NEMO-MEDUSA model outputs.
#' @return A dataframe containing a mean monthly time series of all target variables in NEMO-MEDUSA outputs.
#' @family NEMO-MEDUSA averages
#' @export
summarise_ts <- function(saved) {

  Groups <- readRDS(file = saved) %>%                                          # Read in wide format data file
    dplyr::filter(!weights < 0) %>%                                            # Drop points on land
    dplyr::mutate(weights = dplyr::na_if(weights, 0)) %>%                      # Replace 0 weights with NA so vector lengths match for weighted mean
    tidyr::drop_na(Year, Shore) %>%                                            # Drop points outside of the polygons
    dplyr::group_by(Shore, Year, Month, Depth)

  Ice <- dplyr::filter(Groups, Ice_pres > 0) %>%                               # Remove ice free pixels before averaging
    dplyr::summarise(Ice_Thickness_avg = mean(Ice_Thickness, na.rm = TRUE),    # Get monthly mean sea ice thickness
              Snow_Thickness_avg = mean(Snow_Thickness, na.rm = TRUE),         # Get monthly mean snow thickness
              Ice_conc_avg = mean(Ice_conc, na.rm = TRUE))                     # Get monthly mean sea ice concentration

  Averaged <- Groups %>%
    dplyr::summarise(Salinity_avg = stats::weighted.mean(Salinity, weights, na.rm = TRUE), # Get monthly mean salinity
              Temperature_avg = stats::weighted.mean(Temperature, weights, na.rm = TRUE),
              DIN_avg = stats::weighted.mean(DIN, weights, na.rm = TRUE),
              Detritus_avg = stats::weighted.mean(Detritus, weights, na.rm = TRUE),
              Phytoplankton_avg = stats::weighted.mean(Phytoplankton, weights, na.rm = TRUE),
              Ice_pres = mean(Ice_pres, na.rm = TRUE),                         # Proprtion of pixels covered by ice
              Vertical_diffusivity_avg = stats::weighted.mean(Vertical_diffusivity, weights, na.rm = TRUE),
              Vertical_velocity_avg = stats::weighted.mean(Vertical_velocity, weights, na.rm = TRUE),
              Meridional_avg = stats::weighted.mean(Meridional, weights, na.rm = TRUE),
              Zonal_avg = stats::weighted.mean(Zonal, weights, na.rm = TRUE)) %>%
      dplyr::left_join(Ice) %>%                                                # Add in ice and snow thicknesses
      dplyr::ungroup()

  return(Averaged) }
