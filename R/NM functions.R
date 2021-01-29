
#### Functions to summarise NEMO-MEDUSA model outputs ####

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

    if(data$slab_layer[1] == "D") data[, c("Ice_conc", "Ice_Thickness", "Snow_Thickness"):=NULL] else data
  } else{                                                                           # Run dplyr method

    if(data$slab_layer[1] == "D") {select(data, -c(starts_with("Ice"), Snow_Thickness))} else data}}

#' Summarise NEMO-MEDUSA Output into Decadal Grids
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
NM_decadal_summary <- function(decade, dt) {

if(dt == TRUE){                                                               # Run data.table method
  # data.table::setDT(decade)                                                 # set as a data.table, not needed if decade is already a data.table
  Averaged <- decade[, lapply(.SD, mean, na.rm = TRUE),                       # Average data columns which aren't groups
                     by = c("longitude", "latitude", "Decade", "Month", "Shore", "slab_layer")] # Group by pixel and decade
} else{                                                                       # Run dplyr method
  Averaged <- decade %>%
    dplyr::group_by(longitude, latitude, Decade, Month, Shore, slab_layer) %>%# Group by pixel and decade
    dplyr::summarise_all(mean, na.rm = TRUE) %>%                              # Average data columns
    dplyr::ungroup()                                                          # Ungroup
  }
  return(Averaged)
}

#' Summarise NEMO-MEDUSA Output into Time Series Within Model Compartments
#' This function averages NEMO-MEDUSA monthly summaries into time series for each model compartment.
#'
#' The function groups by model compartment (Depth and Shore zone) and time step (Month and Year).
#' The mean for every target variable is calculated within these groups.
#'
#' @param saved A dataframe containing a summarised month from NEMO-MEDUSA model outputs.
#' @return A dataframe containing a mean monthly time series of all target variables in NEMO-MEDUSA outputs.
#' @family NEMO-MEDUSA averages
#' @export
NM_volume_summary <- function(saved) {

  Groups <- readRDS(file = saved) %>%                                          # Read in wide format data file
#    dplyr::filter(!weights < 0) %>%                                            # Drop points on land
#    dplyr::mutate(weights = dplyr::na_if(weights, 0)) %>%                      # Replace 0 weights with NA so vector lengths match for weighted mean
    tidyr::drop_na(Year, Shore) %>%                                            # Drop points outside of the polygons
    dplyr::group_by(Shore, Year, Month, slab_layer)

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
              Ice_pres = mean(Ice_pres, na.rm = TRUE),                         # Proportion of pixels covered by ice
              Meridional_avg = stats::weighted.mean(Meridional, weights, na.rm = TRUE),
              Zonal_avg = stats::weighted.mean(Zonal, weights, na.rm = TRUE)) %>%
      dplyr::left_join(Ice) %>%                                                # Add in ice and snow thicknesses
      dplyr::ungroup()

  return(Averaged) }

#' Summarise NEMO-MEDUSA Output into Time Series Along Transects
#'
#' This function averages NEMO-MEDUSA monthly summaries into time series for the target boundaries of StrathE2E.
#'
#' The function subsets the NEMO-MEDUSA grid according to the transects object provided. Water exchanges between
#' model compartments are totaled. The boundary conditions of the model domain for variables needed by
#' StrathE2E are summarised as a flow-weighted mean, applying the flow rate at each transect.
#'
#' @param saved A dataframe containing a summarised month from NEMO-MEDUSA model outputs.
#' @param transects A dataframe containing the labelled transects along the model domain boundaries.
#' @param vars A character vector containing the column names to be summarised for boundary conditions. Defaults to the tragets for StrathE2E
#' @return A list containing two summaries.
#' \itemize{
#'  \item{Element 1}{A dataframe containing the total water exchanged between model compartments.}
#'  \item{Element 2}{A dataframe containing the flow-weighted boundary conditions around the model domain.}
#' @family NEMO-MEDUSA averages
#' @export
NM_boundary_summary <- function(saved, transects, vars = c("DIN", "Phytoplankton", "Temperature", "Salinity", "Detritus")) {

  Data <- readRDS(saved) %>%                                                  # Import a NM summary object
    dplyr::select(-c(Shore, weights))                                         # Drop duplicated columns which vonflict
  data.table::setDT(Data, key = c("x", "y", "slab_layer"))                    # Convert to a data.table keyed spatially for quick summaries.

  join <- Data[transects] %>%                                                 # 0.5% of the transects don't catch data (NA), it's because the deep offshore layer doesn't have a buffer on the shore side.
    dplyr::mutate(Flow = ifelse(current == "Zonal", Zonal, Meridional)) %>%   # Grab current perpendicular to the transect
    dplyr::mutate(Flow = ifelse(Flip == T, -1*Flow, Flow)) %>%                # Correct flow so that + always goes IN to a model box
    dplyr::mutate(Flow = Flow * weights,                                      # Weight by transect area to get the volume of water
           Direction = ifelse(Flow > 0, "In", "Out"))                         # Label direction based on flow rate

  ## Summarise water exchanges

  water <- join[, .(Flow = sum(Flow, na.rm = T)),                             # Tally up water movements
                by = c("Shore", "slab_layer", "Direction",                    # By exchanges we want to keep track of
                       "Neighbour", "Month", "Year")] %>%
    tidyr::drop_na()                                                          # The NA transects introduce a dead group, remove.

  ## Summarise boundary conditions

  #*# How do we weight by flow? some are negative
  boundary <- join[perimeter == T & Direction == "In",                        # For transects which bound the perimeter of the model domain
                   lapply(.SD, weighted.mean, w = Flow, na.rm = T),           # Flow-weighted mean of target variables
                   by = c("Shore", "slab_layer", "Neighbour",                 # By groups we want to keep track of
                          "Month", "Year"),
                   .SDcols = vars] %>%                                        # Specify target variables
    tidyr::drop_na() %>%                                                      # The NA transects introduce a dead group, remove.
    dplyr::mutate(Date = as.Date(paste(15, Month, Year, sep = "/"), format = "%d/%m/%Y"),
           Compartment = paste(Shore, slab_layer)) %>%
    tidyr::pivot_longer(eval(vars), names_to = "Variable", values_to = "Measured") # reshape

  result <- list(Flows = water, Boundarys = boundary)                         # Combine both summaries so they can be returned together

  return(result)
}
