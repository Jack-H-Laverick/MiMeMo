

#' Plot Temporal Summaries
#'
#' This function builds a time series, split by model compartment.
#'
#' @param var A "quoted" column name denoting how to colour the map.
#' @return The function creates a time series by model compartment of a summarised variable.
#' The plot is saved into the NEMO-MEDUSA folder.
#' @family NEMO-MEDUSA plots
#' @export
ts_plot <- function(var) {

  ts <- ggplot2::ggplot(TS, aes(x=date, y= get(var), colour = Compartment)) +
    ggplot2::geom_line(size = 0.2) +
    ggplot2::geom_smooth(span = 0.008, size = 0.2, se = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::labs(caption = paste("NM", var, "time series by compartment"), y = var) +
    NULL
  ggplot2::ggsave(paste0("./Figures/NEMO-MEDUSA/TS_", var, ".png"), plot = ts, width = 16, height = 10, units = "cm", dpi = 500)

}

#' Map Spatial Summaries
#'
#' This function builds a map of the passed variable, facetted by month.
#'
#' @param data A dataframe containing the columns "Zonal" and "Meridional" for currents.
#' @param var A "quoted" column name denoting how to colour the map.
#' @param zoom A call to `coord_sf` specifying the project plotting window.
#' @return The function creates a facetted map of a summarised variable. Each facet shows a month.
#' The plot is saved into the NEMO-MEDUSA/grids folder.
#' @family NEMO-MEDUSA plots
#' @export
point_plot <- function(data, var, zoom) {

  decade <- data$Decade[1]; depth <- data$Depth[1]                             # Find out what the data is

  if(depth == "D" & var %in% c("Ice", "Ice_conc", "Ice_Thickness", "Snow_Thickness")) {
    print("Skipped deep ice plot") } else {

      print(paste("plotting", var, "for", decade, depth))                          # Show things are working

      map <- ggplot2::ggplot() +                                                            # Create the base
        #zoom +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Decade:", decade),
             subtitle = paste("Water layer:", depth), x = NULL, y = NULL) +
        ggplot2::geom_raster(data = data, aes(x=x, y=y, fill = get(var))) +
        viridis::scale_colour_viridis(option = "viridis", name = var, na.value = "red") +
        ggplot2::facet_wrap(vars(Month)) +
        #ggplot2::geom_sf(data = world) +
        # bathymetry
        #ggnewscale::new_scale_colour() +
        #ggplot2::geom_sf(data = lines, aes(colour = level), stroke = 0, size = 0.2, show.legend = "line") +
        #ggplot2::scale_colour_manual(name = 'Depth (m)', values = c("-1000" = "white", "-200" = "grey40", "-30" = "black")) +
        #zoom +
        NULL
      ggplot2::ggsave(paste0("./Figures/NEMO-MEDUSA/grids/map ", var, " ", depth, " ", decade, ".png"),
             plot = map, scale = 1, width = 32, height = 20, units = "cm", dpi = 500)
    }
}

#' Map Spatial Summaries
#'
#' This function builds a map of the passed variable, facetted by month.
#'
#' @param data A dataframe containing the columns "Zonal" and "Meridional" for currents.
#' @param var A "quoted" column name denoting how to colour the map.
#' @param zoom A call to `coord_sf` specifying the project plotting window.
#' @return The function creates a facetted map of a summarised variable. Each facet shows a month.
#' The plot is saved into the NEMO-MEDUSA/grids folder.
#' @family NEMO-MEDUSA plots
#' @export
point_plot_voronoi <- function(data, var, zoom) {

  decade <- data$Decade[1]; depth <- data$Depth[1]                             # Find out what the data is

  if(depth == "D" & var %in% c("Ice", "Ice_conc", "Ice_Thickness", "Snow_Thickness")) {
    print("Skipped deep ice plot") } else {

      print(paste("plotting", var, "for", decade, depth))                          # Show things are working

      map <- ggplot2::ggplot() +                                                            # Create the base
        zoom +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Decade:", decade),
                      subtitle = paste("Water layer:", depth), x = NULL, y = NULL) +
        ggvoronoi::geom_voronoi(data = data, aes(x=x, y=y, fill = get(var))) +
        viridis::scale_colour_viridis(option = "viridis", name = var, na.value = "red") +
        ggplot2::facet_wrap(vars(Month)) +
        ggplot2::geom_sf(data = world) +
        # bathymetry
        ggnewscale::new_scale_colour() +
        ggplot2::geom_sf(data = lines, aes(colour = level), stroke = 0, size = 0.2, show.legend = "line") +
        ggplot2::scale_colour_manual(name = 'Depth (m)', values = c("-1000" = "white", "-200" = "grey40", "-30" = "black")) +
        zoom +
        NULL
      ggplot2::ggsave(paste0("./Figures/NEMO-MEDUSA/grids/map ", var, " ", depth, " ", decade, "2.png"),
                      plot = map, scale = 1, width = 32, height = 20, units = "cm", dpi = 500)
    }
}

#' Create Stick Plots to Show Currents
#'
#' This function builds a stick plot to map water velocities in the project window.
#'
#' @param data A dataframe containing the columns "Zonal" and "Meridional" for currents.
#' @param zoom A call to `coord_sf` specifying the project plotting window.
#' @param pre A list of predefined preferences for saving maps in the project (scale, width, height, units, dpi)
#' @return The function creates a map of Water velocities. Colour is used to show speed, sticks are included to show direction.
#' The plot is saved into the NEMO-MEDUSA/Currents folder.
#' @family NEMO-MEDUSA plots
#' @export
stick_plot <- function(data, zoom, pre) {

  #data <- SP[[1]]

  current_uv_scalar <- 500000      # Establish the vector scalar for the currents

  direction <- dplyr::slice(data, seq(1,nrow(data), by = 50)) %>%
    tidyr::drop_na()  # Take a subset of arrows to plot for legibility, dropping arrows with NAs

  sticks <- ggplot(data) +
    ggplot2::geom_point(aes(x=x, y=y, colour = log(abs(Zonal)+abs(Meridional))+1), size = 0.65, stroke = 0, shape = 16) +
    viridis::scale_colour_viridis(option = "viridis", name = "Current\nspeed (m/s)", na.value = NA) +
    # Bathymetry
    ggnewscale::new_scale_colour() +
    ggplot2::geom_sf(data = lines, aes(colour = level), stroke = 0, size = 0.1, show.legend = "line") +
    ggplot2::scale_colour_manual(name = 'Depth (m)', values = c("-1000" = "red", "-200" = "pink" , "-30" = "white")) +
    # Sticks
    ggplot2::geom_segment(data = direction, aes(x=x, xend = x + (Zonal * current_uv_scalar),
                                       y=y, yend = y + (Meridional * current_uv_scalar)), colour = "black", size = 0.1) +
    ggplot2::geom_point(data = direction, aes(x=x, y=y), colour = "white", size = 0.2, stroke = 0, shape = 16) +
    # Land
    ggplot2::geom_sf(data = world, fill = "black", size = 0.2) +
    zoom +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.background = element_rect(fill="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggplot2::labs(caption = paste("Decade:", unique(data$Decade), "    Water layer:", unique(data$Depth)), x = NULL, y = NULL) +
    NULL

  ggplot2::ggsave(paste("./Figures/NEMO-MEDUSA/currents/currents ", unique(data$Depth), unique(data$Decade), ".png"), plot = sticks,
                  scale = pre$scale, width = pre$width, height = pre$height, units = pre$units, dpi = pre$dpi)
}
