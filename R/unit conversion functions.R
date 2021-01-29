
#' Helper functions for common unit conversions
#'
#' These functions are named to make it intuitive to switch between scales, with pairs existing for
#' different directions of conversion.
#'
#' @details
#' Though some of these functions do the same thing numerically, having multiple names which match the desired
#' conversion should reduce friction when developing code. This should also help avoid typos when it comes to
#' 0s, and also stop you needing to constantly look up conversion factors!
#'
#' The `shortwave_to_einsteins` function is not exact and is based on a couple of assumptions to get from watts
#' to einsteins per m^2 per day. It expects an input of total Watts per day per m^2.
#' \itemize{
#'  \item{The proportion of PAR from shortwave downward radiation is 0.43.}
#'  \item{The energy content of a mole of photons is 217.503 kJ, taken from midband radiation (550nm).}}
#'
#' @param data A numeric vector.
#' @return A numeric vector following unit conversion.
#' @name Convert-units
NULL

#' @rdname Convert-units
#' @export
micro_to_full <- function(data) data / 1e6
#' @rdname Convert-units
#' @export
full_to_micro <- function(data) data * 1e6

#' @rdname Convert-units
#' @export
micro_to_milli <- function(data) data / 1e3
#' @rdname Convert-units
#' @export
milli_to_micro <- function(data) data * 1e3

#' @rdname Convert-units
#' @export
milli_to_full <- function(data) data / 1e3
#' @rdname Convert-units
#' @export
full_to_milli <- function(data) data * 1e3

#' @rdname Convert-units
#' @export
l_to_m3 <- function(data) data / 1e3
#' @rdname Convert-units
#' @export
m3_to_l <- function(data) data * 1e3

#' @rdname Convert-units
#' @export
sec_to_day <- function(data) data / 86400
#' @rdname Convert-units
#' @export
day_to_sec <- function(data) data * 86400

#' @rdname Convert-units
#' @export
shortwave_to_einstein <- function(data) {

  # 1 watt is 1 joule per second
  MJ <- (data * 86400)/1e6               # Watts Integrated per day to MegaJoules per day
  MJ_PAR  <- MJ * 0.43                    # PAR fraction of shortwave radiation
  einsteins <- MJ_PAR/(217.503/1000)      # 217 kJ per mol, so convert kJ to MJ and divide

  ## So we scaled shortwave radiation by the proportion as PAR and by the ratio
  ## of moles to energy for mid-band photons (550nm)
  ## This is not exact. 0.43 is a rule of thumb
  ## Also the energy content of a mole of photons is an approximation as it varies
  ## with latitude and across wavelengths

}

#' Convert a U-V velocity field to speed and direction in degrees
#'
#' This function takes a vector of u and v velocities and calculates the direction and speed of the combined movement.
#'
#'This function was lifted from the `Rsenal` package, where it was originally used to calculate wind speeds. All I've done
#'is built a wrapper which accounts for different conventions when describing wind and flow directions.
#'
#' @param u A vector of Zonal currents (from West to East).
#' @param v A vector of Meridional currents (from South to North).
#' @return a dataframe of two columns is returned. Speed contains the composite speed of both velocities on the same scale.
#' Direction is the resolved direction of the flow in degrees, 0 heads north, 90 East, 180 South, 270 West.
#' @family NEMO-MEDUSA spatial tools
#' @export
vectors_2_direction <- function (u, v) {
  u <- -u                                        # This function was built to use wind direction
  v <- -v                                        # Winds  are "opposite", people care about where wind comes from, not where it goes

  # Lovingly lifted from the "Rsenal" package

  degrees <- function(radians) 180 * radians/pi
  mathdegs <- degrees(atan2(v, u))
  wdcalc <- ifelse(mathdegs > 0, mathdegs, mathdegs + 360)
  uvDirection <- ifelse(wdcalc < 270, 270 - wdcalc, 270 - wdcalc + 360)
  uvSpeed <- sqrt(u^2 + v^2)
  return(cbind(uvDirection, uvSpeed))
}
