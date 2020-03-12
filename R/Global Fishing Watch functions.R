
##**## A tidy place to keep functions and reduce clutter in programmes

#### Global Fishing Watch Data Extraction ####

#' Extract the values from a grid under transects along the external boundaries of the model domain
#'
#' This function reads in a datafile and attaches the values needed to transects.
#' It uses a precalculated set of indices of where transects intersect the grid for speed.
#'
#' @param Depth The depth layer to extract data from. Either "S" or "D"
#' @param Data The data object as provided by Sample_OOB
#' @param variables The variables to extract, provided by Sample_OOB
#' @return The function returns a dataframe of transects and their average DIN, chlorophyll, temperature, and salinity values by depth.
#' @family Global Fishing Watch functions
#' @export
Arctic_boats <- function (Data)                {
  Fish <- Data %>%                                                            # Take first csv file as a test
    mutate(Date = date) %>%                                            # We need Lon-lat to be the first two columns so move date column
    mutate(date = lon_bin/100, lat_bin = lat_bin/100)                 # Move Coordinate columns, and divide by 100 because raw data misses decimal point

  colnames(Fish)[2] <- "Latitude"                                             # Change column names for GIS functions
  colnames(Fish)[1] <- "Longitude"

  coordinates(Fish) <- ~ Longitude + Latitude                                 # Converts to a shapefile

  proj4string(Fish) <- proj4string(FAO_arctic)                                # Set the projection of the points to match the polygons

  # Filter the data to rows where polygons lie over the fishing data points
  inside.Arctic <- Fish[!is.na(over(Fish, as(FAO_arctic, "SpatialPolygons"))),] %>%
    fortify(inside.Arctic)                                    # Coerce for ggplot
  return(inside.Arctic) }     # Function to clip boat pings to the FAO regions of interest

#' Extract the values from a grid under transects along the external boundaries of the model domain
#'
#' This function reads in a datafile and attaches the values needed to transects.
#' It uses a precalculated set of indices of where transects intersect the grid for speed.
#'
#' @param Depth The depth layer to extract data from. Either "S" or "D"
#' @param Data The data object as provided by Sample_OOB
#' @param variables The variables to extract, provided by Sample_OOB
#' @return The function returns a dataframe of transects and their average DIN, chlorophyll, temperature, and salinity values by depth.
#' @family Global Fishing Watch functions
#' @export
FAO_boats <- function (Data, Clip)             {

  colnames(Data)[2] <- "Latitude"                                             # Change column names for GIS functions
  colnames(Data)[1] <- "Longitude"
  name <- deparse(substitute(Clip))                                           # Pull object name to attach as a column on the output

  coordinates(Data) <- ~ Longitude + Latitude                                 # Converts to a shapefile

  proj4string(Data) <- proj4string(Clip)                                      # Set the projection of the points to match the polygons

  # Filter the data to rows where polygons lie over the fishing data points
  inside.Arctic <- Data[!is.na(over(Data, as(Clip, "SpatialPolygons"))),] %>%
    fortify(inside.Arctic) %>%                                # Coerce for ggplot
    mutate(Region = name)
  return(inside.Arctic) }     # Function to clip boat pings to specific FAO regions
