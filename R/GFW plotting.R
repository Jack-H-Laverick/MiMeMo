
#### Global Fishing Watch ####

#' Map GFW data for a geartype, above a minimum level of fishing activity
#'
#' This function saves a static map of fishing activity.
#'
#' The function saves a map of fishing activity, as recorded by global fishing watch, with FAO fishing regions overlayed.
#' The user can specify which fishing gear type to visualise, and can also remove points below a certain amount
#' of fishing effort. This is useful for keeping the plot clean, but also allows you to remove any fishing signals
#' you might not trust from the global fishing watch model GFW model.
#'
#' @param gear The gear type used by boats to be mapped.
#' @param limit The minimum number of hours of fishing activity required for a pixel to be included in the map.
#' @return The function saves a map of fishing activity for a particular gear type above a minimum level of activity.
#' @family Global Fishing Watch functions
#' @export
#' @keyword internal
# jacks_lazy_avg <- function (gear, limit) {
#   plot <- ggplot() +
#     geom_point(data = filter(Seasonal, geartype == gear & fishing >= limit),
#                aes(x = long, y = lat), colour = "red", alpha = 0.1, size = 0.1, shape=16, stroke = 0) + # shape and stroke get the smallest point possible
#     geom_polygon(data = world, aes(x=long, y = lat, group = group), fill = "black") +
#     geom_polygon(data = FAO_arctic, aes(long, lat, group = group), colour = "black", fill = "NA", size = 0.1) +                                                        # Plots polygons with holes
#     coord_fixed(1.3, xlim = c(min(-50), max(75)),
#                 ylim = c(min(50), max(90))) +
#     theme_minimal() +
#     labs(x = 'Longitude (W)', y = 'Latitude (N)',
#          title = paste(gear, "activity longer than", limit, "hours in the Arctic (2012-2016)", sep = " ")) +
#     theme(legend.position = "None") +
#     NULL
#   plot
#
#   ggsave(paste("./Figures/GFW/", gear,".png"), plot = plot, scale = 1, width = 16, height = 10, units = "cm", dpi = 500)
#
#   return(plot)
#}

#' Animate a map of GFW data for a geartype, above a minimum level of fishing activity
#'
#' This functions saves an animated map of fishing activity by month.
#'
#' The function saves a map of fishing activity, as recorded by global fishing watch, with FAO fishing regions overlayed.
#' The user can specify which fishing gear type to visualise, and can also remove points below a certain amount
#' of fishing effort. This is useful for keeping the plot clean, but also allows you to remove any fishing signals
#' you might not trust from the global fishing watch model GFW model.
#'
#' @param gear The gear type used by boats to be mapped.
#' @param limit The minimum number of hours of fishing activity required for a pixel to be included in the map.
#' @return The function saves a map of fishing activity for a particular gear type above a minimum level of activity.
#' @family Global Fishing Watch functions
#' @export
#' @keyword internal
# jacks_lazy_animation <- function (gear, limit) {
#   plot <- ggplot() +
#     geom_polygon(data = world, aes(x=long, y = lat, group = group), fill = "black", colour = NA) +
#     geom_point(data = filter(Seasonal, geartype == gear & fishing >= limit),
#                aes(x = long, y = lat),  colour = "red", size = 0.1, shape=16, stroke = 0) + # shape and stroke get the smallest point possible
#     geom_polygon(data = FAO_arctic, aes(long, lat, group = group), colour = "black", fill = "NA", size = 0.1) +                                                        # Plots polygons with holes
#     coord_fixed(1.3, xlim = c(min(-50), max(75)),
#                 ylim = c(min(50), max(90))) +
#     theme_minimal() +
#     labs(x = 'Longitude (W)', y = 'Latitude (N)',
#          title = paste(gear, "activity longer than", limit, "hours in the Arctic (2012-2016)", sep = " ")) +
#     theme(legend.position = "None") +
#     transition_events(start = Month-1, enter_length = 1, exit_length = 1) +          # It fails if your just provide month. I think it's because the frame count drops off the end of the cycle
#     enter_fade() +
#     exit_fade() +
#     labs(subtitle = "Month: {round(frame_time)+1}") +                                         # Comment out to make static
#     NULL
#
#   gganimate::animate(plot, width = 16, height = 10, res = 300, units = "cm", fps = 10, nframes = 100)                                   # Run animation - user specified
#   anim_save(paste("./Figures/GFW/GIF2 ", gear), animation = last_animation())    # Save animation
#   # return(plot) uncomment this and comment out two lines above if you only want to return a single animation
# }

#' Map GFW data for a geartype, above a minimum level of fishing activity (polar projection)
#'
#' This function saves a static map of fishing activity, with a CRS appropriate for the arctic.
#'
#' The function saves a map of fishing activity, as recorded by global fishing watch, with FAO fishing regions overlayed.
#' The user can specify which fishing gear type to visualise, and can also remove points below a certain amount
#' of fishing effort. This is useful for keeping the plot clean, but also allows you to remove any fishing signals
#' you might not trust from the global fishing watch model GFW model.
#'
#' @param gear The gear type used by boats to be mapped.
#' @param limit The minimum number of hours of fishing activity required for a pixel to be included in the map.
#' @return The function saves a map of fishing activity for a particular gear type above a minimum level of activity.
#' @family Global Fishing Watch functions
#' @export
#' @keyword internal
# jacks_polar_static <- function (gear, limit)   {
#   plot <- ggplot() +
#     # bathymetry
#     geom_sf(data = lines, aes(colour = level), stroke = 0, size = 0.2, show.legend = "line") +
#     scale_colour_viridis(name = 'Depth (m)', discrete = TRUE,
#                          guide = guide_legend(override.aes = list(linetype = "solid", shape = NA))) +
#     # boats
#     geom_sf(data = filter(Seasonal_p, geartype == gear & fishing >= limit),
#             colour = "red", size = 0.1, shape=16, stroke = 0) + # shape and stroke get the smallest point possible
#     # FAO
#     geom_sf(data = FAO_arctic_p, colour = "black", fill = "NA", size = 0.1) +                                                        # Plots polygons with holes
#     # world
#     geom_sf(data = world_p) +
#
#     coord_sf(xlim = c(7000000, 000000), ylim = c(7500000, 3500000)) +
#     theme_minimal() +
#     labs(title = paste(gear, "activity longer than", limit, "hours in the Arctic (2012-2016)", sep = " ")) +
#     theme(legend.position = "None") +
#     #    transition_events(start = Month-1, enter_length = 1, exit_length = 1) +          # It fails if your just provide month. I think it's because the frame count drops off the end of the cycle
#     #    enter_fade() +
#     #    exit_fade() +
#     #    labs(subtitle = "Month: {round(frame_time)+1}") +                                         # Comment out to make static
#     NULL
#   ggsave(paste("./Figures/GFW/3", gear,".png"), plot = plot, scale = 1, width = 16, height = 10, units = "cm", dpi = 500)
#   return(plot)
# }     # Plot polar maps

#' Animate a map of GFW data for a geartype, above a minimum level of fishing activity (polar projection)
#'
#' This functions saves an animated map of fishing activity by month. with a CRS appropriate for the arctic.
#'
#' The function saves a map of fishing activity, as recorded by global fishing watch, with FAO fishing regions overlayed.
#' The user can specify which fishing gear type to visualise, and can also remove points below a certain amount
#' of fishing effort. This is useful for keeping the plot clean, but also allows you to remove any fishing signals
#' you might not trust from the global fishing watch model GFW model.
#'
#' @param gear The gear type used by boats to be mapped.
#' @param limit The minimum number of hours of fishing activity required for a pixel to be included in the map.
#' @return The function saves a map of fishing activity for a particular gear type above a minimum level of activity.
#' @family Global Fishing Watch functions
#' @export
#' @keyword internal
# jacks_polar_animate <- function (gear, limit)  {
#   plot <- ggplot() +
#     # bathymetry
#     #geom_sf(data = lines, aes(colour = level), stroke = 0, size = 0.2, show.legend = "line") +
#     #scale_colour_viridis(name = 'Depth (m)', discrete = TRUE,
#     #                     guide = guide_legend(override.aes = list(linetype = "solid", shape = NA))) +
#     # ice
#     geom_sf(data = filter(Ice_m, Sea_ice_avg > 0), aes(alpha = Sea_ice_avg),
#             colour = "turquoise", stroke = 0, size = 0.75, show.legend = "point") +
#     scale_alpha(range = c(0,1), name = 'Sea Ice\n Concentration', guide = guide_legend(override.aes = list(linetype = "blank", shape = 16, size = 4)))+
#     # boats
#     geom_sf(data = filter(Seasonal_p, geartype == gear & fishing >= limit),
#             colour = "red", size = 0.1, shape=16, stroke = 0) + # shape and stroke get the smallest point possible
#     # FAO
#     #geom_sf(data = FAO_arctic_p, colour = "black", fill = "NA", size = 0.1) +
#     # Domains
#     geom_sf(data = Greenland_off, aes(geometry = geometry), colour = "Yellow3", fill = NA, size=0.1) +
#     geom_sf(data = Greenland_in, aes(geometry = geometry), colour = "Yellow", fill = NA, size=0.1) +
#     geom_sf(data = Barents_in, aes(geometry = geometry), colour = "Yellow", fill = NA, size=0.1) +
#     geom_sf(data = Barents_off, aes(geometry = geometry), colour = "Yellow3", fill = NA, size=0.1) +
#     # world
#     geom_sf(data = world_p) +
#
#     coord_sf(xlim = c(6100000, 2800000), ylim = c(7500000, 4400000)) +
#     theme_minimal() +
#     labs(title = paste(gear, "activity longer than", limit, "hours in the Arctic (2012-2016)", sep = " ")) +
#     theme(legend.position = "None") +
#     enter_fade() +
#     exit_fade() +
#     transition_manual(Month) +
#     labs(subtitle = "Month: {current_frame}") +                                         # Comment out to make static
#     NULL
#
#   gganimate::animate(plot, width = 16, height = 10, res = 300, units = "cm", fps = 10, nframes = 100)                                   # Run animation - user specified
#   anim_save(paste("./Figures/GFW/GIF3_domains", gear), animation = last_animation())    # Save animation
#   # return(plot) uncomment this and comment out two lines above if you only want to return a single animation
# }     # Animate polar maps
#
