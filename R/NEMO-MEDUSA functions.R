
#### NEMO - MEDUSA data extraction    ####

#' Check Whether a Vector Contains Data
#'
#' This function provides a check for whether a vector of locations contains any data. This is useful
#' for checking whether a lat/lon pixel doesn't report any data for the deep zone. This is expected if the
#' sea floor at those coordinates is shallower than the maximum depth we set for the shallow zone.
#'
#' @param x A vector of values to check.
#' @return TRUE (if only NAs) or FALSE (if any entry is not NA).
#' @family NEMO-MEDUSA spatial tools
#' @export
empty <- function(x) all(is.na(x))

#' Calculate Water Layer Thicknesses Within an Array
#'
#' This function calculates the thickness of the water layer around each point in a NEMO-MEDUSA array. These thicknesses
#' are needed to calculate weighted averages across a depth window.
#'
#' The function starts by calculating the thickness of the depth window represented by the shallowest layer in the array.
#' This is the depth to halfway between the shallowest and the next layer, subtracting any depth to the top of the window, i.e.
#' If the top of the window is the sea surface, the first array layer is at 5m, and the second is at 10m, the water thickness about
#' the first point is 7.5m (mean(5, 10) - 0). The function then checks whether the depth of this layer of points is shallower than the depth
#' to the seafloor for any of the points. For these points which exist below the seafloor, thicknesses are replaced by the depth to the seafloor,
#' minus the depth to the top of the depth window.
#'
#' The function then populates the thicknesses of the mid-layers in a for-loop, working increasingly deeper. The thickness about a point is
#' calculated as the difference between the depth halfway to the point above, and the depth halfway to the point below. Checks and replacements
#' are performed, as explained above, to see whether any points are now beyond the limit of the bathymetry, or approaching the `top` or
#' `bottom` of the depth window. If they are, the claculations are performed again using the new max and min depths.
#'
#' The thicknesses for the deepest layer of points are calculated as the depth to the seafloor minus the mid-depth of the deepest two layers.
#'
#' If a weight is <= 0, it is replaced with NA as this indicates land.
#'
#' @param top The shallowest depth in the depth window.
#' @param bottom The deepest depth in the depth window.
#' @return An array of water layer thicknesses to match the dimensions of a NEMO-MEDUSA array.
#' @family NEMO-MEDUSA spatial tools
#' @export
get_weights <- function(top, bottom) {
  #top <- 200                                                                    # shallowest depth of the slice
  #bottom <- 1000                                                                # What's the bottom of the slice? for incorporating into a function

  weights <- array(NA, c(235,190,38))                                          # Initialise an array to hold weights
  first <- weights[,,1]
  first[] <- mean(Space$nc_depth[1:2]) - top %>% rep(times = 44650)            # Water above the first midpoint minus the depth at the top of the slice
  marks <- mask_bathy > mean(Space$nc_depth[1:2])                              # Which cells contain model outputs deeper than the seafloor?
  first[!marks] <- mask_bathy[!marks] - top                                    # Replace these with the depth to sea floor
  weights[,,1] <- first

  weights[,,38] <- mask_bathy - mean(Space$nc_depth[37:38])                    # The remaining water column thickness is the sea floor - the deepest midpoint.

  for (i in 2:37) {
    #i <- 23
    last_midpoint <- mean(Space$nc_depth[(i-1):i])                               # Find the mid depth to the layer above
    next_midpoint <- mean(Space$nc_depth[i:(i+1)])                               # Find the mid depth to the layer below

    if(top > last_midpoint) above <- top else above <- last_midpoint             # If the top of the slice is deeper than the previous midpoint, use the top of the slice
    if(bottom < next_midpoint) below <- bottom else below <- next_midpoint       # If the next midpoint is deeper than the bottom of the slice, use the bottom of the slice

    weights[,,i] <- below - above %>% rep(times = 44650)                         # Calculate layer thickness and repeat to fill the array

    marks <- mask_bathy > below                                                  # Is the seafloor deeper than the bottom of the layer?
    weights[,,i][!marks] <- mask_bathy[!marks] - above                           # If not, replace these with the depth to sea floor - the top of the water layer

  }                                                          # Roll through each matrix and calculate the water thickness using the next depth, bottom of the slice, or bathymetry, whichever is smaller
  no_weight <- weights[] <= 0; weights[no_weight] <- NA                        # Finally if a weight is <= 0 get NA

  return(weights)
}

#' Calculate Water Layer Thicknesses Within an Array (Velocities)
#'
#' This function is a variant of `get_weights` which applies to water velocities. Water velocities from NEMO-MEDUSA are on a different
#' vector of depths. This function calculates the thickness of the water layer around each point in a NEMO-MEDUSA array. These thicknesses
#' are needed to calculate weighted averages across a depth window.
#'
#' The function starts by calculating the thickness of the depth window represented by the shallowest layer in the array.
#' This is the depth to halfway between the shallowest and the next layer, subtracting any depth to the top of the window, i.e.
#' If the top of the window is the sea surface, the first array layer is at 5m, and the second is at 10m, the water thickness about
#' the first point is 7.5m (mean(5, 10) - 0). The function then checks whether the depth of this layer of points is shallower than the depth
#' to the seafloor for any of the points. For these points which exist below the seafloor, thicknesses are replaced by the depth to the seafloor,
#' minus the depth to the top of the depth window.
#'
#' The function then populates the thicknesses of the mid-layers in a for-loop, working increasingly deeper. The thickness about a point is
#' calculated as the difference between the depth halfway to the point above, and the depth halfway to the point below. Checks and replacements
#' are performed, as explained above, to see whether any points are now beyond the limit of the bathymetry, or approaching the `top` or
#' `bottom` of the depth window. If they are, the claculations are performed again using the new max and min depths.
#'
#' The thicknesses for the deepest layer of points are calculated as the depth to the seafloor minus the mid-depth of the deepest two layers.
#'
#' If a weight is <= 0, it is replaced with NA as this indicates land.
#'
#' @param top The shallowest depth in the depth window.
#' @param bottom The deepest depth in the depth window.
#' @return An array of water layer thicknesses to match the dimensions of a NEMO-MEDUSA array.
#' @family NEMO-MEDUSA spatial tools
#' @export
get_weights.W <- function(top, bottom) {

  weights <- array(NA, c(235,190,39))                                            # Initialise an array to hold weights
  first <- weights[,,1]
  first[] <- mean(DepthsW[1:2]) - top %>% rep(times = 44650)                     # Water above the first midpoint minus the depth at the top of the slice
  marks <- mask_bathy > mean(DepthsW[1:2])                                       # Which cells contain model outputs deeper than the seafloor?
  first[!marks] <- mask_bathy[!marks] - top                                      # Replace these with the depth to sea floor
  weights[,,1] <- first

  weights[,,39] <- mask_bathy - mean(DepthsW[38:39])                             # The remaining water column thickness is the sea floor - the deepest midpoint.

  for (i in 2:38) {
    #i <- 23
    last_midpoint <- mean(DepthsW[(i-1):i])                                      # Find the mid depth to the layer above
    next_midpoint <- mean(DepthsW[i:(i+1)])                                      # Find the mid depth to the layer below

    if(top > last_midpoint) above <- top else above <- last_midpoint             # If the top of the slice is deeper than the previous midpoint, use the top of the slice
    if(bottom < next_midpoint) below <- bottom else below <- next_midpoint       # If the next midpoint is deeper than the bottom of the slice, use the bottom of the slice

    weights[,,i] <- below - above %>% rep(times = 44650)                         # Calculate layer thickness and repeat to fill the array

    marks <- mask_bathy > below                                                  # Is the seafloor deeper than the bottom of the layer?
    weights[,,i][!marks] <- mask_bathy[!marks] - above                           # If not, replace these with the depth to sea floor - the top of the water layer

  }                                                          # Roll through each matrix and calculate the water thickness using the next depth, bottom of the slice, or bathymetry, whichever is smaller
  no_weight <- weights[] <= 0; weights[no_weight] <- NA                        # Finally if a weight is <= 0 get NA

  return(weights)
}

#' Calculate Water Layer Thicknesses Within an Array (Detritus)
#'
#' This function is a variant of `get_weights` which applies to the old NEMO-MEDUSA grid. Detrital nitrogen from NEMO-MEDUSA are on an old
#' spatial grid. This function calculates the thickness of the water layer around each point in a NEMO-MEDUSA array. These thicknesses
#' are needed to calculate weighted averages across a depth window.
#'
#' The function starts by calculating the thickness of the depth window represented by the shallowest layer in the array.
#' This is the depth to halfway between the shallowest and the next layer, subtracting any depth to the top of the window, i.e.
#' If the top of the window is the sea surface, the first array layer is at 5m, and the second is at 10m, the water thickness about
#' the first point is 7.5m (mean(5, 10) - 0). The function then checks whether the depth of this layer of points is shallower than the depth
#' to the seafloor for any of the points. For these points which exist below the seafloor, thicknesses are replaced by the depth to the seafloor,
#' minus the depth to the top of the depth window.
#'
#' The function then populates the thicknesses of the mid-layers in a for-loop, working increasingly deeper. The thickness about a point is
#' calculated as the difference between the depth halfway to the point above, and the depth halfway to the point below. Checks and replacements
#' are performed, as explained above, to see whether any points are now beyond the limit of the bathymetry, or approaching the `top` or
#' `bottom` of the depth window. If they are, the claculations are performed again using the new max and min depths.
#'
#' The thicknesses for the deepest layer of points are calculated as the depth to the seafloor minus the mid-depth of the deepest two layers.
#'
#' If a weight is <= 0, it is replaced with NA as this indicates land.
#'
#' @param top The shallowest depth in the depth window.
#' @param bottom The deepest depth in the depth window.
#' @return An array of water layer thicknesses to match the dimensions of a NEMO-MEDUSA array.
#' @family NEMO-MEDUSA spatial tools
#' @export
get_weights.old <- function(top, bottom) {
  #top <- 200                                                                    # shallowest depth of the slice
  #bottom <- 1000                                                                # What's the bottom of the slice? for incorporating into a function

  weights <- array(NA, c(nrow(Space$nc_lon),ncol(Space$nc_lon), 27))           # Initialise an array to hold weights
  first <- weights[,,1]
  first[] <- mean(Space$nc_depth[1:2]) - top %>% rep(times = nrow(Space$nc_lon) * ncol(Space$nc_lon)) # Water above the first midpoint minus the depth at the top of the slice
  marks <- mask_bathy > mean(Space$nc_depth[1:2])                              # Which cells contain model outputs deeper than the seafloor?
  first[!marks] <- mask_bathy[!marks] - top                                    # Replace these with the depth to sea floor
  weights[,,1] <- first

  weights[,,27] <- mask_bathy - mean(Space$nc_depth[26:27])                    # The remaining water column thickness is the sea floor - the deepest midpoint.

  for (i in 2:26) {
    #i <- 23
    last_midpoint <- mean(Space$nc_depth[(i-1):i])                               # Find the mid depth to the layer above
    next_midpoint <- mean(Space$nc_depth[i:(i+1)])                               # Find the mid depth to the layer below

    if(top > last_midpoint) above <- top else above <- last_midpoint             # If the top of the slice is deeper than the previous midpoint, use the top of the slice
    if(bottom < next_midpoint) below <- bottom else below <- next_midpoint       # If the next midpoint is deeper than the bottom of the slice, use the bottom of the slice

    weights[,,i] <- below - above %>% rep(times = nrow(Space$nc_lon) * ncol(Space$nc_lon)) # Calculate layer thickness and repeat to fill the array

    marks <- mask_bathy > below                                                  # Is the seafloor deeper than the bottom of the layer?
    weights[,,i][!marks] <- mask_bathy[!marks] - above                           # If not, replace these with the depth to sea floor - the top of the water layer

  }                                                          # Roll through each matrix and calculate the water thickness using the next depth, bottom of the slice, or bathymetry, whichever is smaller
  no_weight <- weights[] <= 0; weights[no_weight] <- NA                        # Finally if a weight is <= 0 get NA

  return(weights)
}

#' Get Latitudes, Longitudes, & Depths
#'
#' This function gets the latitudes, longitudes, and depths which define the spatial location of points in an array of NEMO-MEDUSA outputs.
#'
#' Each variable of interest in the netcdf file is imported, and then collected into a list.
#'
#' @param file The full name of a netcdf file.
#' @return A list of three elements:
#' \itemize{
#'  \item{\emph{nc_lat -}}{ A matrix of latitudes which maps onto the first and second dimension of a NEMO-MEDUSA array.}
#'  \item{\emph{nc_lon -}}{ A matrix of longitudes which maps onto the first and second dimension of a NEMO-MEDUSA array.}
#'  \item{\emph{nc_depth -}}{ A vector of depths which match the third dimension of a NEMO-MEDUSA array.}
#'  }
#' @family NEMO-MEDUSA variable extractors
#' @export
get_spatial <- function(file) {
  nc_raw <- nc_open(file)                              # Open up a netcdf file to see it's raw contents (var names)
  nc_lat <- ncvar_get(nc_raw, "nav_lat")               # Extract a matrix of all the latitudes
  nc_lon <- ncvar_get(nc_raw, "nav_lon")               # Extract a matrix of all the longitudes
  nc_depth <- ncvar_get(nc_raw, "deptht")              # Extract a matrix of depths
  nc_close(nc_raw)                                     # You must close an open netcdf file when finished to avoid data loss
  all <- list("nc_lat" = nc_lat, "nc_lon" = nc_lon, "nc_depth" = nc_depth)
  return(all)
}

#' Calculate the Domain Area per Grid Point
#'
#' This function takes an array of a variable, and an array of water thicknesses to perform a weighted average across depth. The depth
#' window to be averaged can be specified, so this function can be used to create both shallow and deep layers (or more for that matter).
#'
#' @param points A Simple Feature object og the grid points within the model domain.
#' @param area A Simple Feature object containing the model domain.
#' @return the `points` object is returned, but instead of points, the geometry column now contains polygons representing the area closest to each point. A column for the size of this area is also gained.
#' @family NEMO-MEDUSA spatial tools
#' @export
voronoi_grid <- function(points, area) {

  result <- purrr::map(1:nrow(area), ~{                            # For each polygon in area
    voronoi <- points %>%                                          # Take the grid points
      sf::st_geometry() %>%                                        # To get sfc from sf
      sf::st_union() %>%                                           # To get a sfc of MULTIPOINT type
      sf::st_voronoi(envelope = sf::st_geometry(area[.x,])) %>%    # Voronoi polygon for the area
      sf::st_collection_extract(type = "POLYGON") %>%              # A list of polygons
      sf::st_sf() %>%                                              # From list to sf object
      sf::st_join(points) %>%                                      # put names back
      sf::st_intersection(area[.x,]) %>%                           # Cut to shape of NC state
      dplyr::mutate(Cell_area = units::drop_units(sf::st_area(.))) # Area of each polygon
  }) %>%
    dplyr::bind_rows() %>%                                         # Combine the results from each area
    sf::st_sf(geomc = .$geometry, crs = 4326)                      # Reinstate attributes of the geometry column

}

#' Summarise Across Depths in a NEMO-MEDUSA Array
#'
#' This function takes an array of a variable, and an array of water thicknesses to perform a weighted average across depth. The depth
#' window to be averaged can be specified, so this function can be used to create both shallow and deep layers (or more for that matter).
#'
#' @param data An array containing estimates of a variable.
#' @param depth A TRUE FALSE vector indicating the depth layers to be extracted.
#' @param weights An array containing water thicknesses.
#' @return A matrix containing the weighted averages of the variable across the depth window of interest.
#' @family NEMO-MEDUSA spatial tools
#' @export
stratify  <- function(data, depth, weights) {

  # data <- nc_zonal ; depth <- Deep_mark ; weights <- dw                   # testing

  new <- data[,,depth] * weights[,,depth]                                      # Select slice of depths to average, multiply values by the weights
  empties <- apply(new, c(1,2), empty)                                         # Find pixels with all depths shown by NA (locations of fake 0s)

  new2 <- apply(new, c(1,2), sum, na.rm = TRUE)                                # Sum the weighted values at a pixel
  denominator <- apply(weights[,,depth], c(1,2), sum, na.rm = TRUE)            # Sum the weights
  weighted_mean <- new2/denominator                                            # Divide by the sum of the weights
  weighted_mean[empties] <- NA                                                 # Sum replaces an all NA dimension with 0, overwrite these by position
  return(weighted_mean)
}

#' Get Salinity, Temperature & Sea Ice Concentration
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Salinity, temperature & sea ice concentration can be found in files labelled "grid_T_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @param start3D
#' @param count3D
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_sea   <- function(path, file, grid, space) {

  print(stringr::str_glue("{file} Extracting Salinity, Temperature, and Sea Ice concentration"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_saline <- ncdf4::ncvar_get(nc_raw, "vosaline", space$start3D, space$count3D)          # Extract an array of salinities
  nc_temp <- ncdf4::ncvar_get(nc_raw, "votemper", space$start3D, space$count3D)            # Extract an array of temperatures
  nc_ice <- ncdf4::ncvar_get(nc_raw, "soicecov")                               # Extract a matrix of ice fractions
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  shallow <- grid %>%                                                          # Grab the tidied dataframe of lat-longs
    dplyr::mutate(Ice_conc = as.numeric(nc_ice),                               # Append new variable to coordinates (no depths for ice)
           Salinity = as.numeric(stratify(nc_saline, space$shallow, space$s.weights)),     # Collapse shallow salinity into 2D and convert to long format
           Temperature = as.numeric(stratify(nc_temp, space$shallow, space$s.weights)),     # Collapse shallow temperatures into 2D and convert to long format
           Depth = "S")                                                        # Introduce depth column

  deep <- grid %>%                                                             # Grab the tidied dataframe of lat-longs
    dplyr::mutate(Ice_conc = NA,                                               # Insert empty column to allow fast binding by position
           Salinity = as.numeric(stratify(nc_saline, space$deep, space$d.weights)),
           Temperature = as.numeric(stratify(nc_temp, space$deep, space$d.weights)),
           Depth = "D")                                                        # Collapse, reshape and append deepwater data

  all <- rbind(shallow, deep) %>%                                              # Bind both sets, this pipeline avoids computationally demanding reshaping
    dplyr::filter(Shore_dist > 0) %>%                                          # Remove points on land
    dplyr::mutate(Depth = as.factor(Depth))

  return(all)
}

#' Get Dissolved Inorganic Nitrogen & Chlorophyll a
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' DIN & chlorphyll a can be found in files labelled "ptrc_T_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_bio   <- function(path, file, grid, space) {

  print(stringr::str_glue("{file} Extracting Dissolved Inorganic Nitrogen and Chlorophyll"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_DIN <- ncdf4::ncvar_get(nc_raw, "DIN", space$start3D, space$count3D)                  # Extract an array for the variable
  nc_CHD <- ncdf4::ncvar_get(nc_raw, "CHD", space$start3D, space$count3D)
  nc_CHN <- ncdf4::ncvar_get(nc_raw, "CHN", space$start3D, space$count3D)
  nc_Chl <- nc_CHD + nc_CHN ; rm(nc_CHD, nc_CHN)
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  shallow <- grid %>%                                                          # Grab the tidied dataframe of lat-longs
    dplyr::mutate(DIN = as.numeric(stratify(nc_DIN, space$shallow, space$s.weights)),        # Collapse shallow DIN into 2D and convert to long format
           Chlorophyll = as.numeric(stratify(nc_Chl, space$shallow, space$s.weights)),       # Collapse shallow chlorophyll into 2D and convert to long format
           Depth = "S")                                                        # Introduce depth column

  deep <- grid %>%                                                             # Grab the tidied dataframe of lat-longs
    dplyr::mutate(DIN = as.numeric(stratify(nc_DIN, space$deep, space$d.weights)),
           Chlorophyll = as.numeric(stratify(nc_Chl, space$deep, space$d.weights)),
           Depth = "D")                                                        # Collapse, reshape and append deepwater data

  all <- rbind(shallow, deep) %>%                                              # Bind both sets, this pipeline avoids computationally demanding reshaping
    dplyr::filter(Shore_dist > 0) %>%                                          # Remove points on land
    dplyr::mutate(Depth = as.factor(Depth))

  return(all)
}

#' Get Ice Pesence & Thickness, & Snow Thickness
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Ice presence & thickness, & snow thickness can be found in files labelled "icemod_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_ice   <- function(path, file, grid, space) {

  print(stringr::str_glue("{file} Extracting Ice presence, and Ice and Snow thickness"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_Ice <- ncdf4::ncvar_get(nc_raw, "ice_pres")                               # Extract a matrix of ice presence
  nc_Ithick <- ncdf4::ncvar_get(nc_raw, "iicethic")                            # Extract ice thicknesses
  nc_Sthick <- ncdf4::ncvar_get(nc_raw, "isnowthi")                            # Extract snow thicknesses
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  shallow <- grid %>%                                                          # Grab the tidied dataframe of lat-longs
    dplyr::mutate(Ice_pres = as.numeric(nc_Ice),
           Ice_Thickness = as.numeric(nc_Ithick),
           Snow_Thickness = as.numeric(nc_Sthick),
           Depth = "S")                                                        # Introduce depth column

  deep <- grid %>%                                                             # Grab the tidied dataframe of lat-longs
    dplyr::mutate(Ice_pres = NA,
           Ice_Thickness = NA,
           Snow_Thickness = NA,
           Depth = "D")                                                        # Introduce depth column

  all <- rbind(shallow, deep) %>%                                              # Bind both sets, this pipeline avoids computationally demanding reshaping
    dplyr::filter(Shore_dist > 0) %>%                                          # Remove points on land
    dplyr::mutate(Depth = as.factor(Depth))

  return(all)
}

#' Get Vertical Velocity and Vertical Eddy Diffusivity
#'
#' This function reads in the title variables from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Vertical velocity and vertical eddy diffusivitiy can be found in files labelled "W".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for title variables
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_vertical   <- function(path, file, grid, space) {

  print(stringr::str_glue("{file} Extracting Vertical water movements"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                                 # Open up a netcdf file to see it's raw contents (var names)
  nc_vel <- ncvar_get(nc_raw, "vovecrtz", space$start3DW, space$count3DW)                  # Extract an array for the variable
  nc_dif <- ncvar_get(nc_raw, "votkeavt", space$start3DW, space$count3DW)
  ncdf4::nc_close(nc_raw)                                                      # You must close an open netcdf file when finished to avoid data loss

  shallow <- grid %>%                                                        # Grab the tidied dataframe of lat-longs
    dplyr::mutate(Vertical_velocity = as.numeric(stratify(nc_vel, space$shallow_W, space$s.wweights)), # Collapse shallow DIN into 2D and convert to long format
           Vertical_diffusivity = as.numeric(stratify(nc_dif, space$shallow_W, space$s.wweights)),     # Collapse shallow chlorophyll into 2D and convert to long format
           Depth = "S")                                                        # Introduce depth column

  deep <- grid %>%                                                           # Grab the tidied dataframe of lat-longs
    dplyr::mutate(Vertical_velocity = as.numeric(stratify(nc_vel, space$deep_W, space$d.wweights)),
           Vertical_diffusivity = as.numeric(stratify(nc_dif, space$deep_W, space$d.wweights)),
           Depth = "D")                                                        # Collapse, reshape and append deepwater data

  all <- rbind(shallow, deep) %>%                                              # Bind both sets, this pipeline avoids computationally demanding reshaping
    dplyr::filter(Shore_dist > 0) %>%                                          # Remove points on land
    dplyr::mutate(Depth = as.factor(Depth))

  return(all)
}

#' Get Meridional Currents
#'
#' This function reads in the title variable from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Meridional currents can be found in files labelled "grid_V_".
#'
#' Each variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variable.
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for the title variable
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_merid <- function(path, file, grid, space) {

  print(stringr::str_glue("{file} Extracting Meridional currents"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_merid <- ncdf4::ncvar_get(nc_raw, "vomecrty", space$start3D, space$count3D)         # Pull meridinal currents
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  shallow <- grid %>%                                                        # Grab the tidied dataframe of lat-longs
    mutate(Meridional = as.numeric(stratify(nc_merid, space$shallow, space$s.weights)),  # Collapse shallow meridional currents into 2D and convert to long format
           Depth = "S")                                                      # Introduce depth column

  deep <- grid %>%                                                           # Grab the tidied dataframe of lat-longs
    mutate(Meridional = as.numeric(stratify(nc_merid, space$deep, space$d.weights)),     # Collapse, reshape and append deepwater data
           Depth = "D")                                                      # Collapse, reshape and append deepwater data

  all <- rbind(shallow, deep) %>%                                            # Bind both sets, this pipeline avoids computationally demanding reshaping
    filter(Shore_dist > 0) %>%                                               # Remove points on land
    mutate(Depth = as.factor(Depth))

  return(all)
}

#' Get Zonal Currents
#'
#' This function reads in the title variable from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Zonal currents can be found in files labelled "grid_U_".
#'
#' The variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param path The path to the NEMO-MEDUSA model outputs.
#' @param file The name of a netcdf file containing the title variables.
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for the title variable
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_zonal <- function(path, file, grid, space) {

  print(stringr::str_glue("{file} Extracting Zonal currents"))
  nc_raw <- ncdf4::nc_open(paste0(path, file))                               # Open up a netcdf file to see it's raw contents (var names)
  nc_zonal <- ncdf4::ncvar_get(nc_raw, "vozocrtx", space$start3D, space$count3D)         # Pull zonal current
  nc_close(nc_raw)                                                           # You must close an open netcdf file when finished to avoid data loss

  shallow <- grid %>%                                                        # Grab the tidied dataframe of lat-longs
    dplyr::mutate(Zonal = as.numeric(stratify(nc_zonal, space$shallow, space$s.weights)),# Collapse shallow meridional currents into 2D and convert to long format
           Depth = "S")                                                      # Introduce depth column

  deep <- grid %>%                                                           # Grab the tidied dataframe of lat-longs
    dplyr::mutate(Zonal = as.numeric(stratify(nc_zonal, space$deep, space$d.weights)),   # Collapse, reshape and append deepwater data
           Depth = "D")                                                      # Collapse, reshape and append deepwater data

  all <- rbind(shallow, deep) %>%                                            # Bind both sets, this pipeline avoids computationally demanding reshaping
    dplyr::filter(Shore_dist > 0) %>%                                        # Remove points on land
    dplyr::mutate(Depth = as.factor(Depth))

  return(all)
}

#' Get Detritus
#'
#' This function reads in the title variable from NEMO-MEDUSA model outputs and reshapes for StrathE2E.
#' Detritus can be found in an older NEMO-MEDUSA data excerpt.
#'
#' The variable of interest in the netcdf file is imported, only reading within an x/y window specified with `start3D` and
#' `count3D`. The values are then passed to `stratify` to calculate two average matrices, one a weighted vertical average
#' of the shallow zone, and the same for the deep zone. Latitude and longitudes are also attached to each horizontal pixel.
#'
#' @param file The name of a netcdf file containing the title variables.
#' @param start The indices to start reading data from in 3 dimensions.
#' @param count The lengths of data reads in 3 dimensions.
#' @param grid A dataframe to bind the extracted values to.
#' @param shallow, A TRUE / FALSE vector indicating which depth layers should be combined into a shallow zone.
#' @param s.weights An array of water thicknesses to use as weights for collapsing layers into a shallow zone.
#' @param deep, A TRUE / FALSE vector indicating which depth layers should be combined into a deep zone.
#' @param d.weights An array of water thicknesses to use as weights for collapsing layers into a deep zone.
#' @return A dataframe containing points and lat/lon coordinates, with the average NEMO-MEDUSA model outputs for the title variable
#' in the shallow and deep zone. The dataframe contains the data for a single day.
#' @family NEMO-MEDUSA variable extractors
#' @export
get_detritus <- function(file, start, count, grid, shallow, s.weights, deep, d.weights) {

  print(stringr::str_glue("Extracting detrital nitrogen"))
  nc_raw <- ncdf4::nc_open(file)                                             # Open up a netcdf file to see it's raw contents (var names)
  nc_detritus <- ncdf4::ncvar_get(nc_raw, "DET", start, count)               # Pull detritus
  ncdf4::nc_close(nc_raw)                                                    # You must close an open netcdf file when finished to avoid data loss

  shallow <- grid %>%                                                        # Grab the tidied dataframe of lat-longs
    dplyr::mutate(Detritus = as.numeric(stratify(nc_detritus, shallow, s.weights)), # Collapse shallow meridional currents into 2D and convert to long format
           Depth = "S")                                                      # Introduce depth column

  deep <- grid %>%                                                           # Grab the tidied dataframe of lat-longs
    dplyr::mutate(Detritus = as.numeric(stratify(nc_detritus, deep, d.weights)), # Collapse, reshape and append deepwater data
           Depth = "D")                                                      # Collapse, reshape and append deepwater data

  all <- rbind(shallow, deep) %>%                                            # Bind both sets, this pipeline avoids computationally demanding reshaping
    dplyr::filter(Shore_dist > 0) %>%                                               # Remove points on land
    dplyr::mutate(Depth = as.factor(Depth))

  return(all)
}

#' Condense Daily netcdf Files into a Month by Type
#'
#' This function takes the metadata for multiple netcdf files and creates a spatial summary for the month.
#' \strong{For one file type only}.
#'
#' The function bridges the step between extracting data from a netcdf file, and creating an average dataframe
#' for a whole month of NEMO-MEDUSA model outputs.
#'
#' Different file types require a different get function. This function takes a collection of netcdf files of
#' the same time from the same month, and passes them to the correct `get_*` for data extraction. The results
#' are passed back to this function, before the estimates from different days at the same lat-lon-depth
#' combination are averaged to get a single number for the month.
#'
#' Creating intermediate monthly objects allows the terabytes of data to be reduced without encountering memory issues.
#' Also, working on independent monthly packets of data means we can parallelise any data processing for speed.
#'
#' @param data A dataframe containing the metadata of multiple netcdf files which share a type.
#' @param ... Additional arguments passed to the relevant `get_*` function.
#' @return The function returns a dataframe containing the monthly average shalllow and deep spatial grids for
#' variables of interest.
#' @family NEMO-MEDUSA variable extractors
#' @export
type_in_month <- function(data, ...) {

  Type <- data[1,3]                                                         # Pull type from the file

  if(Type == "grid_T_") get <- get_sea                                      # Change the extracting function based on file contents
  if(Type == "grid_U_") get <- get_zonal
  if(Type == "grid_V_") get <- get_merid
  if(Type == "grid_W_") get <- get_vertical
  if(Type == "icemod_") get <- get_ice
  if(Type == "ptrc_T_") get <- get_bio

  Month.type <- data %>%                                                    # Take the year
    dplyr::mutate(data = purrr::map2(Path, File, get, ...)) %>%             # Extract data from each file
    tidyr::unnest(data) %>%                                                 # Extract all encoded data
    dplyr::group_by(Longitude, Latitude, Depth, .drop=FALSE) %>%
    dplyr::summarise_if(is.numeric, mean)
  return(Month.type)
}

#' Condense Daily netcdf Files into a Monthly Summary
#'
#' This function takes the metadata for multiple netcdf files. It then creates a common spatial summary for all
#' the variables of interest for a month.
#'
#' The function takes the metadata for a collection of files which contain data from the same month. The files
#' are split into data packets which share the same file type, before being passed to `type_in_month` to be summarised.
#' `type_in_month` reduces the NEMO-MEDUSA model outputs from large arrays to effectively two matrices. The summaries for
#' each file type are returned to this function and get bound into a single dataframe. Points outside the project window are
#' removed before saving the dataframe for the month in "./Objects/Months/".
#'
#' Creating intermediate monthly objects allows the terabytes of data to be reduced without encountering memory issues.
#' Also, working on independent monthly packets of data means we can parallelise any data processing for speed.
#'
#' @param data A dataframe containing the metadata of multiple netcdf files from a common month.
#' @param ... Additional arguments to be passed to get_* functions.
#' @return The function returns a dataframe containing the monthly average shalllow and deep spatial grids for
#' \strong{all} the variables of interest in NEMO-MEDUSA outputs.
#' @family NEMO-MEDUSA variable extractors
#' @export
whole_month <- function(data, ...) {

  Month <- data[1,5] ; Year <- data[1,4]                                    # Pull date

  Month <- split(data, f = list(data$Type)) %>%                             # Split out the files for this month by type, so they can be averaged together
    purrr::map(type_in_month, ...) %>%                                      # Pull a whole month of data from a single file type
    purrr::reduce(dplyr::full_join) %>%                                     # Join together all the data packets
  # dplyr::right_join(spine) %>%     # a)                                   # Cut out rows outside of polygons and attach compartment labels
    dplyr::right_join(Window) %>%    # b)                                   # Cut out rows outside of plotting window
    saveRDS(., file = paste("./Objects/Months/NM", Month, Year, "rds", sep = "."))    # save out a data object for one whole month
}

#' Condense Daily netcdf Files into a Monthly Summary
#'
#' This function is a variant of `whole_month` which only considers detritus files. This is a convenience function.
#' Detritus data is available on it's own grid, so there's no reason to pass to `type_in_month` and specify the
#' required metadata.
#'
#' The function takes the metadata for a collection of files which contain data from the same month. The files
#' are passed to `get_detritus` to be summarised, reducing the NEMO-MEDUSA model outputs from large arrays to effectively
#' two matrices for the month. Points outside the project window are removed before saving the dataframe for the month in
#' "./Objects/Detritus/".
#'
#' Creating intermediate monthly objects allows the terabytes of data to be reduced without encountering memory issues.
#' Also, working on independent monthly packets of data means we can parallelise any data processing for speed.
#'
#' @param data A dataframe containing the metadata of multiple netcdf files from a common month.
#' @param targets A dataframe containing the points from the grid we want to retain data for.
#' @param ... Additional arguments passed on to `get_detritus`.
#' @return The function returns a dataframe containing the monthly average shalllow and deep spatial grids for
#' detrital nitrogen.
#' @family NEMO-MEDUSA variable extractors
#' @export
detritus_month <- function(data, targets, ...) {

  Month <- data[1,3] ; Year <- data[1,2]                                    # Pull date

  Month <- data %>%                                                         # Take the year
    dplyr::mutate(data = purrr::map(data$value, get_detritus, ...)) %>%     # Extract detritus data from each file
    tidyr::unnest(data) %>%                                                 # Extract all encoded data
    dplyr::group_by(Longitude, Latitude, Depth, .drop=FALSE) %>%
    dplyr::summarise_if(is.numeric, mean) %>%
    dplyr::right_join(targets) %>%                                                  # Cut out rows outside of plotting window
    saveRDS(., file = paste("./Objects/Detritus/Det", Month, Year, "rds", sep = "."))    # save out a data object for one whole month
}

#### NEMO - MEDUSA drivers extraction ####

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

  raw <- nc_open(file)
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

  nc_raw <- nc_open(File)                                                    # Open up a netcdf file to see it's raw contents (var names)
  nc_var <- ncvar_get(nc_raw, Type, c(Space$Limits$Lon_start, Space$Limits$Lat_start, 1),  # Extract the variable of interest
                      c(Space$Limits$Lon_count, Space$Limits$Lat_count, -1)) # cropped to window, with all time steps
  nc_close(nc_raw)                                                           # You must close an open netcdf file when finished to avoid data loss

  Data <- as.data.frame.table(nc_var, responseName = "Measured") %>%         # Reshape array as dataframe
    rename(Longitude = Var1, Latitude = Var2, Time_step = Var3) %>%          # Name the columns
    mutate(Longitude = rep(rep(Space$Lons,                                   # Replace the factor levels with dimension values
                               times = length(unique(Latitude))), times = length(unique(Time_step))),
           Latitude = rep(rep(Space$Lats,
                              each = length(unique(Longitude))), times = length(unique(Time_step))),
           Time_step = rep(1:length(unique(Time_step)),
                           each = length(unique(Latitude)) * length(unique(Longitude)))) %>%
    right_join(domains_mask) %>%                                             # Crop to domain
    left_join(months) %>%                                                    # Assign a month to each time step
    mutate(Year = Year,                                                      # Attach Year
           Type = Type)                                                      # Attach variable name

  if(Type == "SWF") Data <- group_by(Data, Month, Year, Type)                # We don't need to bother accounting for shore in light data
  if(Type == "T150") Data <- group_by(Data, Month, Year, Type, Shore)        # We care about shore for temperature, retain interesting columns

  Summary <- summarise(Data, Measured = weighted.mean(Measured, Cell_area))  # Average by time step.

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

  #File <- Airtemp_files$File[1] ; Type <- Airtemp_files$Type[1] ; Year <- Airtemp_files$Year[1] # test
  if(Type == "SWF") months <- Light_months                                 # Get the file holding the months
  if(Type == "T150") months <- Airtemp_months                              # For the time steps of this data

nc_raw <- nc_open(File)                                                    # Open up a netcdf file to see it's raw contents (var names)
nc_var <- ncvar_get(nc_raw, Type, c(Space$Limits$Lon_start, Space$Limits$Lat_start, 1),  # Extract the variable of interest
                    c(Space$Limits$Lon_count, Space$Limits$Lat_count, -1)) # cropped to window, with all time steps
nc_close(nc_raw)                                                           # You must close an open netcdf file when finished to avoid data loss

DT <- as.data.table(nc_var, value.name = "Measured") %>%       # Pull array
  setnames(old = c("V1", "V2", "V3"), new = c("Longitude", "Latitude", "Time_step")) %>% # Name the columns
  .[, ':='(Longitude = Space$Lons[Longitude],                              # read ':=' as mutate
           Latitude = Space$Lats[Latitude],                                # Replace the factor levels with dimension values
           Month = months[Time_step, "Month"],                             # Assign months to time steps
           Year = Year,                                                    # Add year
           Type = Type)] %>%                                               # Add variable name
  merge(as.data.table(domains_mask), all.y = TRUE)                         # Crop to domain

## Variable specific summaries

if(Type == "SWF") Data <- DT[,by = .(Month, Year, Type),                   # We don't need to bother accounting for shore in light data
                             .(Measured = weighted.mean(Measured, Cell_area))]  # Average by time ste, weighted by cell area
if(Type == "T150") Data <- DT[,by= .(Month, Year, Type, Shore),            # We care about shore for temperature, retain interesting columns
                             .(Measured = weighted.mean(Measured, Cell_area))]  # Average by time ste, weighted by cell area
return(Data)
}

#### Data averaging                   ####

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
 #   select(-c(geometry, weights, Day, Bathymetry, Shore_dist)) %>%
    dplyr::select(-c(weights, Day, Bathymetry, Shore_dist)) %>%
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
#' @return If the `data` contains shallow data, no action is taken. If `data` contains deep data, columns for variables only
#' relevant in the shallow zone are dropped.
#' @family NEMO-MEDUSA averages
#' @export
strip_ice <- function(data) {
  if(data$Depth[1] == "D") {select(data, -c(starts_with("Ice"), Snow_Thickness))} else data}

#' Average into Decadal Grids
#'
#' This function averages cleaned NEMO-MEDUSA monthly summaries into decadal grids.
#'
#' The function groups by all spatial variables (Longitude, Latitude, Depth, and Shore zone), and by decade and month.
#' The mean for every other variable is calculated within these groups.
#'
#' @param saved A dataframe containing a summarised month from NEMO-MEDUSA model outputs. It must contain the columns:
#' Longitude, Latitude, Decade, Month, Shore, and Depth.
#' @return A dataframe containing a summarised decade of spatialy resolved NEMO-MEDUSA outputs.
#' @family NEMO-MEDUSA averages
#' @export
summarise_sp <- function(decade) {

  Averaged <- decade %>%
    dplyr::group_by(Longitude, Latitude, Decade, Month, Shore, Depth) %>%     # Group by pixel and decade
    dplyr::summarise_all(mean, na.rm = TRUE) %>%                              # Average data columns
    dplyr::ungroup()                                                          # Ungroup
  return(Averaged)
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
  ret <- tibble::as_tibble(ret)
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

  # saved <- "./Objects/Months/NM.1.1981.rds"

  # test <- readRDS(file = saved) %>%                                    # Read in wide format data file
  #   filter(!weights < 0) %>%
  #   mutate(weights = na_if(weights, 0)) %>%                                        # Replace 0 weights with NA so vector lengths match for weighted mean
  #   drop_na(Year, Shore) %>%
  #   group_by(Shore, Year, Month, Depth)
  #
  # summarise(test, Salinity_avg = weighted.mean(Salinity, weights, na.rm = TRUE),  # Get monthly mean salinity
  #                 Salinity_sd = weighted.sd(Salinity, weights, na.rm = TRUE))  # Get monthly mean salinity

  ## Try and get the vertical SDs to work

  Groups <- readRDS(file = saved) %>%                                          # Read in wide format data file
    dplyr::filter(!weights < 0) %>%                                            # Drop points on land
    dplyr::mutate(weights = dplyr::na_if(weights, 0)) %>%                      # Replace 0 weights with NA so vector lengths match for weighted mean
    tidyr::drop_na(Year, Shore) %>%                                            # Drop points outside of the polygons
    dplyr::group_by(Shore, Year, Month, Depth)

  Ice <- dplyr::filter(Groups, Ice_pres > 0) %>%                               # Remove ice free pixels before averaging thicknesses
    dplyr::summarise(Ice_Thickness_avg = mean(Ice_Thickness, na.rm = TRUE),    # Get monthly mean sea ice thickness
              Snow_Thickness_avg = mean(Snow_Thickness, na.rm = TRUE),         # Get monthly mean snow thickness
              # SD
              Ice_Thickness_sd = stats::sd(Ice_Thickness, na.rm = TRUE),       # Get monthly mean sea ice thickness
              Snow_Thickness_sd = stats::sd(Snow_Thickness, na.rm = TRUE))     # Get monthly mean snow thickness

  Averaged <- Groups %>%
    dplyr::summarise(Salinity_avg = stats::weighted.mean(Salinity, weights, na.rm = TRUE), # Get monthly mean salinity
              Temperature_avg = stats::weighted.mean(Temperature, weights, na.rm = TRUE),
              DIN_avg = stats::weighted.mean(DIN, weights, na.rm = TRUE),
              Chlorophyll_avg = stats::weighted.mean(Chlorophyll, weights, na.rm = TRUE),
              Ice_pres = mean(Ice_pres, na.rm = TRUE),                         # Proprtion of pixels covered by ice
              Ice_conc_avg = mean(Ice_conc, na.rm = TRUE),                     # Get monthly mean sea ice concentration
              Vertical_diffusivity_avg = stats::weighted.mean(Vertical_diffusivity, weights, na.rm = TRUE),
              Vertical_velocity_avg = stats::weighted.mean(Vertical_velocity, weights, na.rm = TRUE),
              Meridional_avg = stats::weighted.mean(Meridional, weights, na.rm = TRUE),
              Zonal_avg = stats::weighted.mean(Zonal, weights, na.rm = TRUE),
              # SD
              Salinity_sd = radiant.data::weighted.sd(Salinity, weights, na.rm = TRUE), # Get monthly mean salinity
              Temperature_sd = radiant.data::weighted.sd(Temperature, weights, na.rm = TRUE),
              DIN_sd = radiant.data::weighted.sd(DIN, weights, na.rm = TRUE),
              Chlorophyll_sd = radiant.data::weighted.sd(Chlorophyll, weights, na.rm = TRUE),
              Ice_conc_sd = stats::sd(Ice_conc, na.rm = TRUE),                 # Get monthly mean sea ice concentration
              #Vertical_diffusivity_sd = weighted.sd(Vertical_diffusivity, weights, na.rm = TRUE),
              #Vertical_velocity_sd = weighted.sd(Vertical_velocity, weights, na.rm = TRUE),
              Meridional_sd = radiant.data::weighted.sd(Meridional, weights, na.rm = TRUE),
              Zonal_sd = radiant.data::weighted.sd(Zonal, weights, na.rm = TRUE)) %>%
      dplyr::left_join(Ice) %>%                                                # Add in ice and snow thicknesses
      dplyr::ungroup()

  return(Averaged) }

#' Average into Detrital Time Series
#'
#' This function is a variant of `summarise_ts` which only considers detritus files.
#'
#' The function groups by model compartment (Depth and Shore zone) and time step (Month and Year) before calculating
#' mean detrital nitrogen within groups.
#'
#' @param saved A dataframe containing a summarised month from NEMO-MEDUSA model outputs.
#' @return A dataframe containing a mean monthly time series of detrital nitrogen in NEMO-MEDUSA outputs.
#' @family NEMO-MEDUSA averages
#' @export
summarise_ts_detritus <- function(saved) {

  # saved <- "./Objects/Detritus/Det.1.1988.rds"

  Averaged <- readRDS(file = saved) %>%                                        # Read in wide format data file
    dplyr::filter(!weights < 0) %>%                                            # Drop points on land
    dplyr::mutate(weights = dplyr::na_if(weights, 0)) %>%                      # Replace 0 weights with NA so vector lengths match for weighted mean
    tidyr::drop_na() %>%                                                       # Drop points outside of the polygons and without weights
    dplyr::group_by(Shore, Year, Month, Depth) %>%
    dplyr::summarise(Detritus_avg = stats::weighted.mean(Detritus, weights, na.rm = TRUE),   # Get monthly mean detrital nitrogen
              Detritus_sd = radiant.data::weighted.sd(Detritus, weights, na.rm = TRUE)) %>%
    dplyr::ungroup()

  return(Averaged) }
