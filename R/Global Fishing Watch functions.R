
##**## A tidy place to keep functions and reduce clutter in programmes

#### Global Fishing Watch Data Extraction ####

#' Import Fishing Effort Within Specified Areas
#'
#'
#'
#' @param file The full name to a .csv file containing Global Fishing Watch data.
#' @param area A Simple Feature object representing the areas to return fishing activity in.
#' @return A dataframe of cells from a 0.01 degree grid which contain fishing activity,
#' labelled by the FAO region they fall in.
#' @family Global Fishing Watch functions
#' @export
get_cropped_fishing <- function (file, area) {

  Fish <- utils::read.csv(file, header = TRUE) %>%                           # Take first csv file as a test
    dplyr::mutate(Longitude = lon_bin/100, Latitude = lat_bin/100) %>%       # Move Coordinate columns, and divide by 100 because raw data misses decimal point
    sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    sf::st_join(area) %>%
    tidyr::drop_na()
}

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
