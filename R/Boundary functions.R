
##**## Functions needed to sample the exchanges between compartment bopundaries

#### Extracting vertical water exchanges ####

#' Interpolate vertical water movements to the boundary depth between shallow and deep compartments
#'
#' This function reads in a 3 dimensional array of water movements, and uses the `rcdo` package to interpolate
#' values to a plane at the boundary depth between shallow and deep compartments.
#'
#' The function uses `rcdo::nc_remap` to interpolate vertical eddy diffusivity and vertical water velocities at
#' the boundary depth between the shallow and deep model compartments. Points which fall outside the geographic
#' extent of model compartments are dropped. The retained points are passed back to `interp_month` for further averaging
#' as a dataframe.
#'
#' @param Path The path to a netcdf file containing target data, inherited from `interp_month`.
#' @param File The name of a netcdf file containing target data, inherited from `interp_month`.
#' @param Points An SF object containing the grid points which fall within the domain polygons, inherited from `interp_month`.
#' @param Boundary The boundary depth between shallow and deep compartments, defaults to 60 m for MiMeMo.
#' @return A data frame containing vertical eddy diffusivity and vertical velocity at the boundary depth
#' between shallow and deep compartments, for a single day as a spatial plane.
#' @family NEMO-MEDUSA variable extractors
# interp_vertical <- function(Path, File, Points, Boundary = 60) {
#
#   V_interface <- rcdo::nc_remap(paste0(Path, File), vars = c("vovecrtz", "votkeavt"), vert_depths = Boundary) %>% # Interpolate at 60 m
#     dplyr::filter(Latitude != 0) %>%                                                 # Remove incorrectlly labelled points (outside domain anyway)
#     dplyr::rename(Velocity = vovecrtz, `Eddy Diffusivity` = votkeavt) %>%            # Rename variables
#     dplyr::inner_join(Points, .)                                                     # Limit to points of interest
#
#   ## Check CDO has formed the grid correctly
#   # ggplot() + geom_sf(data = V_interface, aes(colour = Velocity))
#
#   Velocity <- dplyr::mutate(V_interface, Flow = ifelse(Velocity >= 0, "Upwelling", "Downwelling")) %>% # Label upwelling and downwelling
#     #  group_by(Flow) %>%                                                     # Group for averaging (splitting up and down welling)
#     dplyr::summarise(Value = mean(Velocity))                                         # Calculate mean vertical veolcities for each region
#
#   day <- dplyr::group_by(V_interface) %>%                                            # Repeat grouping for eddy diffusivitity
#     dplyr::summarise(Value = mean(`Eddy Diffusivity`)) %>%                           # Average by region
#     dplyr::mutate(Flow = "Eddy Diffusivity") %>%                                     # Attach a label (so we can group by variable when plotting timeseries)
#     dplyr::bind_rows(., Velocity)                                                    # Combine summaries for the time step
#   return(day)
# }

#' Extract vertical water movements at the boundary between shallow and deep compartments
#'
#' This function reads in a packet of daily netcdf files from the same month, and returns a single monthly dataframe
#' of vertical water movements.
#'
#' The function reads in a packet of monthly netcdf files and passes them to `interp_vertical`.
#'
#' `interp_vertical` interpolates to the depth boundary between the shallow and deep compartments. These daily dataframes represent a spatial plane.
#'
#' The daily dataframes are then passed back to `interp_month`, which averages again by month and variable, dropping the spatial information.
#'
#' @param data A list of paths to netcdf files which share the same month.
#' @param grid_points An SF object containing the grid points which fall within the domain polygons.
#' @return The function returns a dataframe of points averaged within a month and interpolated to the boundary depth.
#' @family NEMO-MEDUSA variable extractors
#' @export
interp_month <- function(data, grid_points) {

  month <- data %>%                                             # Take the month
    dplyr::mutate(data = purrr::map2(Path, File, .f = interp_vertical, Points = grid_points)) %>%  # Extract and average data from each day
    tidyr::unnest(data) %>%                                     # unnext the column holding the extraction
    dplyr::group_by(Year, Month, Flow) %>%                      # Group by information we want to retain
    dplyr::summarise(Value = mean(Value)) %>%                   # Average the values within a month
    dplyr::ungroup()                                            # Ungroup for speed
  return(month)
}

#### Making transects at the boundaries of model domains ####

#' Break an SF linestring into segments at corners
#'
#' This function is wrapped within the boundaries function and is used to convert the perimeter of
#' model zones into seperate segments between each corner. These segments can then be used independently to sample data.
#'
#' @param line The geographic perimeter of the inshore and offshore zone as SF object.
#' @param crs The coordinate reference system for the project. This is inherited from `boundaries` and should be specified in
#' the project region file.
#' @return The function returns a list of individual SF lines connecting each corner of the perimeters of model compartments.
#' @family Boundary sampling functions
#' @export
to_segments <- function(line, crs) {

  g = sf::st_geometry(line) %>% sf::st_cast("POINT")                        # Grab the geometry and extract the points defining the line
  hd = head(g, -1)                                                          # Take all points except the last
  tl = tail(g, -1)                                                          # Take all points except the first
  purrr::map2(hd, tl, function(x,y) sf::st_combine(sf::st_sfc(x, y, crs = crs))) %>% # Pair the points across vectors
    purrr::map(function(x) sf::st_cast(x, "LINESTRING"))                    # Turn pairs of points into lines
}

#' Create a series of transects along the boundaries of model compartments
#'
#' This function takes the domain polygon file and returns transects for use when sampling
#' fluxes across boundaries between model compartments.
#'
#' The function exposes the geometry for the domain polygons to the `to_segments` function.
#' `to_segments` converts the geometry into a collection of transects and passes them back to this function.
#' A unique id is added for each transect along with the length of the transect for use when
#' calculating weighted averages.
#'
#' @param domain An SF object containing polygons for the geographic extent of the inshore and offshore zones for the project.
#' @param crs The coordinate reference system for the project. This should be specified in the project region file.
#' @return The function returns a dataframe containg SF boundary transects around the model domain, a length column is
#' included for weighting transects in future averaging functions.
#' @family Boundary sampling functions
#' @export
boundaries <- function(domain, crs) {

  passed <- crs                                                               # Pass CRS to nested custom function

  segments <- domain %>%
    dplyr::pull(geometry) %>%                                                 # Open the list column containing the geometry
    purrr::map(to_segments, crs = passed) %>%                                 # Apply function to break each line in turn at the corners
    unlist(recursive = FALSE) %>%                                             # Bring all the segments into the same level
    do.call(c, .)                                                             # Bind

  Length <- sf::st_length(segments) %>%                                       # Add in the length of each segment for weighting
    as.data.frame() %>%
    tibble::rowid_to_column(var = "Segment") %>%                              # Add an ID for each segment
    sf::st_sf(geometry = segments)
  return(Length)
}

#### Labelling transects ####

#' Check if a segment runs parallel to latitude or longitude
#'
#' This function checks whether the transects created to sample the boundaries of model compartments lie on
#' a lat-lon grid. If they do, then we know that when sampling meridional and zonal currents only one is needed
#' per transect.
#'
#' The function takes a boundary segment (transect) and reprojects it on a lat-lon grid. The coordinates are extracted,
#' and checked for identical latitudes or longitudes.
#'
#' @param segment An SF line object representing a transect along the perimeter of a model domain polygon.
#' @return The function returns a dataframe detailing whether laitudues are the same and longitudes are the same along a segment.
#' @family Boundary sampling functions
#' @export
check_grid <- function(segment, Weighted) {

  latlon <- Weighted[segment,] %>%
    sf::st_transform(4326) %>%                  # Convert to latitude and longitude
    sf::st_coordinates() %>%                    # Pull the coordinates
    round(digits = 3)                           # Round to drop conversion error

  Lon_same <- latlon[1,"X"] == latlon[2, "X"]   # Are the longitudes the same either end of the line?
  Lat_same <- latlon[1,"Y"] == latlon[2, "Y"]   # Are the latitudes the same either end of the line?

  check <- c(Lon_same, Lat_same) %>%            # Connect queries as a row for a dataframe
    t() %>%
    as.data.frame()

  return(check)
}

#' Test whether "positive" currents flow in or out of a model compartment at a segment
#'
#' This function takes a transect from the edge of a model compartment and works out where it is spatially. This is
#' neccessary for correctly averaging variables along the whole boundary of a compartment.
#'
#' The function determines which model compartment is next to the focal compartment. The focal compartment is either the
#' offshore or inshore polygon, with the offshore polygon selected ahead of the inshore polygon when these share a boundary.
#' This defines the link between model compartments.
#'
#' Direction (in/out) is also calculated relative to the focal polygon. At the western boundary of a polygon, positive zonal
#' currents \emph{(West to East)} indicate water flowing into the model compartment. However, at the Eastern boundary,
#' flows from West to East now leave the model compartment. A True or False column called flip is included. If a positive value
#' for a water velocity is leaving the polygon, or a negative flow value is entering the polygon, flip can be used to multiply
#' by -1 to correctly sum the water movements from the perspective of the focal polygon.

#' @param segment An SF line object representing a transect along the perimeter of a model domain polygon.
#' @return The function returns a dataframe detailing which two model compartments are either side of the transect,
#' and the direction of the exchange.
#' @family Boundary sampling functions
#' @export
direction <- function(segment) {

 #segment <- 1                                                   # Testing function
 #segment <- 24246

  midpoint <- sf::st_line_sample(Checked[segment,], n = 1)       # Grab the midpoint of a line segment

  if(Checked[segment,]$Shore == "Offshore") domain <- dplyr::filter(domains, Shore == "Offshore")# Change the domain polygon to match the boundary segment
  if(Checked[segment,]$Shore == "Inshore") domain <- dplyr::filter(domains, Shore == "Inshore")  # Change the domain polygon to match the boundary segment
  if(Checked[segment,]$current == "Meridional") {
    minus <- c(0, -0.001)                                        # Adjust for point below the segment
    plus <-  c(0, +0.001)                                        # Adjust for point above the segment
    flow_shift <- c(+100, 0)                                     # Adjust for current indicator
  #  flow_plot <- ggplot2::geom_segment(aes(xend= min(flow[,"X"]), y = min(flow[,"Y"]), # PLotting line for current indicator
  #                                x = max(flow[,"X"]), yend = max(flow[,"Y"])), arrow = arrow())
  }            # Change the shift used for test point relative midpoint
  if(Checked[segment,]$current == "Zonal")      {
    minus <- c(-0.001, 0)
    plus <-  c(+0.001, 0)
    flow_shift <- c(0, +100)
  #  flow_plot <- ggplot2::geom_segment(aes(x = min(flow[,"X"]), y = min(flow[,"Y"]),
  #                                xend = max(flow[,"X"]), yend = max(flow[,"Y"])), arrow = arrow())
  }            # Based on current of interest

  coords <- midpoint %>% sf::st_transform(4326)                  # Transform to lat-lon to ensure test points are perpendicular to the line

  mid_minus <- coords + minus                                    # Shift from the mid point
  sf::st_crs(mid_minus) <- 4326                                  # set crs
  mid_minus <- sf::st_transform(mid_minus, crs) %>%              # change crs for plotting
    sf::st_cast("POINT")

  mid_plus <- coords + plus                                      # Shift from the mid point
  sf::st_crs(mid_plus) <- 4326                                   # set crs
  mid_plus <- sf::st_transform(mid_plus, crs) %>%                # change crs for plotting
    sf::st_cast("POINT")

  #flow <- sf::st_union(x = mid_plus, y = mid_minus) %>%          # Link and shift points to make an arrow to illustrate flow
  #  sf::st_cast("LINESTRING") %>%
  #  + flow_shift
  #sf::st_crs(flow) <- crs
  #flow <- sf::st_coordinates(flow)

  test <- sf::st_sf(geometry = c(mid_plus, mid_minus), side = c("plus", "minus")) %>% # Create a full SF object to allow point in polygon analysis
    sf::st_join(domain, join = st_intersects) %>%                # Are the points either side of the boundary in a domain polygon?
    dplyr::mutate(Contained = dplyr::if_else(is.na(Shore), "out", "in")) %>%   # Label the points
    dplyr::mutate(Flip = dplyr::if_else(side == "plus" & Contained == "out" |  # Determine whether positive currents flow in or out of the domain and so if they need to be flipped
                            side == "minus" & Contained == "in", TRUE, FALSE))

  neighbour <- dplyr::filter(test, Contained == "out") %>%       # Grab the SF which is outside the focal polygon
    sf::st_intersects(domains) %>%                               # Which polygon DOES this point sit in? (a list of 1 number)
    as.numeric() %>%                                             # Drop unneccessary formatting
    domains$Shore[.]                                             # Pull the shore label for the polygon

  # window <- st_bbox(Checked[segment,])                         # Get a zoom which fits the segment of interest

  # ggplot() +
  #   geom_sf(data = domain, fill = "yellow") +
  #   geom_sf(data = Checked[segment,]) +
  #   geom_sf(data = midpoint, colour = "red") +
  #   geom_sf(data = test, aes(colour = Contained), show.legend = "point") +
  #   flow_plot +        # The arrow still doesn't plot perpendicularly for all locations on the projection, but the tests are working fine
  #   theme_minimal() +
  #   labs(x = NULL, y = NULL, subtitle = paste0("Should currents be flipped (for + numbers into polygon)? = ", test$Flip[1])) +
  #   coord_sf(xlim = c(window$xmin, window$xmax), ylim = c(window$ymin, window$ymax))

  summary <- dplyr::filter(test, Contained == "in") %>%          # Take the metadata associated with the point in the focal polygon
    dplyr::select(Flip) %>%                                      # Keep only interesting columns
    sf::st_set_geometry(NULL) %>%                                # Drop unneccessary geometry
    dplyr::mutate(Neighbour = as.character(neighbour)) %>%       # Attach the neighbouring polygon
    tidyr::replace_na(list(Neighbour = "Ocean"))                 # If there wasn't a neighbouring polygon, assign ocean

  return(summary)
}

#### Sampling Flows transects ####

#' Extract the values from a grid under transects along the external boundaries of the model domain
#'
#' A data object of water velocities is inherited from `Sample` before using a precalculated set of
#' indices of where transects intersect the grid for speedy extraction.
#'
#' Extracted values are attached to transects before returning to `Sample` for further averaging.
#'
#' @param var The component of water movement to extract, either "Zonal" or "Meridional" velocities, inherited from `Sample`.
#' @param Depth The depth layer to extract data from. Either "S" or "D" inherited from `Sample`.
#' @param Data The data object as passed from `Sample`.
#' @param transects A nested list containing 4 sets of transects. Each set is an SF dataframe for a combination of var and Depth.
#' @param intersections A nested list containing 4 sets of intersections. Each set indicates which grid cells a transect touches, for a set of transects.
#' @return The function returns a dataframe of transects and their average zonal \strong{OR} meridional water velocities at a depth.
#' @family Boundary sampling functions
#' @export
extract <- function (var, Depth, Data, transects, intersections) {

  Data <- Data[[Depth]]

  Samples <- purrr::map(intersections[[Depth]][[var]],                       # Select the relevant set of intersections between transects and the grid
                        function(x) mean(Data[x,var], na.rm = T)) %>%        # calulate the average current per transect
    unlist() %>%                                                             # Collapse list to vector
    dplyr::mutate(transects[[Depth]][[var]], Sample = .) %>%                 # Copy transect meta-data over to the samples
    tidyr::drop_na(Sample) %>%                                               # NM points on land return NA, remove these
    dplyr::mutate(Sample = ifelse(Flip == TRUE, -1* Sample, Sample),         # Flip current direction if required
           Depth = Depth)

  return(Samples)
}

#' Sample water flows across transects along the boundaries of model compartments
#'
#' This function calculates the water exchanges between model compartments and external boundaries.
#'
#' A datafile containing water velocities is read in and split by depth before calling `extract` to
#' pull the correct velocities.
#'
#' `Sample` then groups exchanges by depth, shore zone, neighbour and direction. Flows
#' are summed to return net water movements for the given time step.
#'
#' @param file Path to the .rds object containing both Zonal and Meridional water velocities.
#' @param ... Additional arguments to pass to `extract`.
#' @return The function returns a dataframe of net water movements between model compartments
#' and the external environment in a given month.
#' @family Boundary sampling functions
#' @export
Sample <- function(file, ...) {

  #file <- "./Objects/Months/NM.1.1980.rds"
  #file <- "./Objects/Months/NM.1.1981.rds"

  Data <- readRDS(file) %>%                                                  # Read in a current file
    split(.$Depth) %>%                                                       # Seperate shallow and deep data
    purrr::map(.x = ., .f = dplyr::left_join, x = cells)

    Summary <- purrr::map2(rep(c("Meridional", "Zonal"), each = 2),
                  c("S", "D", "S", "D"), extract, Data = Data, ...) %>%      # Extract for the combinations of depth and current
    data.table::rbindlist() %>%                                              # Bind
    dplyr::mutate(Flow = Sample * Weights,                                   # Weight the velocities by transect area
           Direction = ifelse(Sample > 0, "In", "Out")) %>%           # Create tag for flows in and out of the box
    dplyr::group_by(Shore, Depth, Direction, Neighbour) %>%                  # Data to retain when summing
    dplyr::summarise(Flow = sum(Flow)) %>%                                   # Sum flows
    dplyr::ungroup() %>%                                                     # Ungroup for speed
    dplyr::mutate(Year = Data[["S"]]$Year[1], Month = Data[["S"]]$Month[1])  # Add date from original file

  return(Summary)
}

#### Sampling Boundaries transects ####

#' Extract the values from a grid under transects along the external boundaries of the model domain
#'
#' This function is a variant of extract which only considers the external boundaries (Out Of Bounds -> OOB) of the model domain
#' for variables other than water movement.
#'
#' A data object is inherited from `Sample_OOB` before using a precalculated set of indices of where
#' transects intersect the grid for speedy extraction.
#'
#' Extracted values are attached to transects before returning to `Sample_OOB` for further averaging.
#'
#' @param Depth The depth layer to extract data from. Either "S" or "D" and inherited from `Sample_OOB`.
#' @param Data The data object as passed from `Sample_OOB`.
#' @param transects A nested list containing 2 sets of transects. Each set is an SF dataframe for sampling either the shallow or deep layer.
#' @param intersections A nested list containing 2 sets of intersections. Each set indicates which grid cells a transect touches, for a depth layer.
#' @param variables The variables to extract, inherited from `Sample_OOB`.
#' @return The function returns a dataframe of transects and their average DIN, chlorophyll, temperature,
#' and salinity values by depth.
#' @family Boundary sampling functions
#' @export
extract_OOB <- function (Depth, Data, transects, intersections, variables) {

  Data <- Data[[Depth]]

Samples <- purrr::map(intersections[[Depth]], function(x) colMeans(Data[x, variables], na.rm = T)) %>% # Select the cells and calulate the average current per transect
  do.call(rbind, .) %>%                                                              # Collapse list to vector
  cbind(transects[[Depth]], .) %>%                                          # Copy transect meta-data over to the samples
  tidyr::drop_na() %>%                                                                 # NM points on land return NA, remove these
  dplyr::mutate(Depth = Depth)

  return(Samples)
}

#' Sample along the external boundaries of the model domain
#'
#' This function is a variant of Sample which only considers the external boundaries (Out Of Bounds -> OOB) of the model domain
#' for variables other than water movement.
#'
#' A datafile is imported and split by depth before calling `extract_OOB` to sample the grid.
#'
#' After extraction, `Sample_OOB` calculates the average for variables of interest by depth and shore zone.
#'
#' @param file Path to the .rds object containing data.
#' @param variables The variables to extract, inherited from `Sample_OOB`.
#' @param ... Additional arguments to pass to `extract_OOB`.
#' @return The function returns a dataframe of average DIN, chlorophyll, temperature, and salinity at
#' the external model boundary by compartment for a month.
#' @family Boundary sampling functions
#' @export
Sample_OOB <- function(file, variables, ...) {

  #file <- "./Objects/Months/NM.1.1980.rds"
  #file <- "./Objects/Months/NM.1.1981.rds"
#   variables <- c("DIN", "Chlorophyll", "Temperature", "Salinity")

  Data <- readRDS(file) %>%                                              # Read in a current file
    split(.$Depth) %>%                                                   # Seperate shallow and deep data
    purrr::map(.x = ., .f = dplyr::left_join, x = cells)                 # Reorder data onto the grid

  Summary <- rbind(extract_OOB("S", Data, variables = variables, ...),    # Extract data from shallow OOB transects
                   extract_OOB("D", Data, variables =variables, ...)) %>% # Extract data from deep OOB transects
    dplyr::group_by(Shore, Depth) %>%                                    # Data to retain when summing
    dplyr::select(eval(variables)) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::ungroup() %>%                                                 # Ungroup for speed
    dplyr::mutate(Year = Data[["S"]]$Year[1],
           Month = Data[["S"]]$Month[1],
           Date = as.Date(paste("15", Month, Year, sep = "/"), format = "%d/%m/%Y"),
           Compartment = paste(Shore, Depth)) %>%
    tidyr::pivot_longer(c(DIN, Chlorophyll, Temperature, Salinity), names_to = "Variable", values_to= "Measured")

  return(Summary)
}
