
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

  raw <- ncdf4::nc_open(File)                                                 # Open a file
  row <- raw$dim$row$vals                                                   # Extract the data needed to calculate
  center_lat <- ncdf4::ncvar_get(raw, "center_lat")                         # Lat-lon coordinates of pixels
  center_lon <- ncdf4::ncvar_get(raw, "center_lon")
  lon_step <- ncdf4::ncvar_get(raw, "lon_step")
  col <- ncdf4::ncvar_get(raw, "col")                                       # Each row has a different number of columns

  data <- data.frame(Bin = raw$dim$bin$vals,                                # Pull pixel
                     SPM = ncvar_get(raw, "SPM-OC5_mean"),                  # SPM
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
