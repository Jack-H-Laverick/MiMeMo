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

  usethis::ui_done("{usethis::ui_field(x)} completed. {praise::praise('${Exclamation}!')}")}

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
