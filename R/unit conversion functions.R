
#' Helper functions for common unit conversions
#'
#' These functions are named to make it intuitive to switch between scales, with pairs existing for
#' different directions of conversion.
#'
#' Though some of these functions do the same thing numerically, having multiple names which match the desired
#' conversion should reduce friction when developing code. This should also help avoid typos when it comes to
#' 0s, and also stop you needing to constantly look up conversion factors!
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
