#' Change detection using the PVts-\eqn{\beta} approach
#'
#' This algorithm will allow to detect disturbances in the forests using
#' all the available Landsat set. In fact, it can also be run with sensors
#' such as MODIS.
#'
#' @author Yonatan Tarazona Coronel
#'
#' @section References:
#' Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving tropical
#' deforestation detection through using photosynthetic vegetation time
#' series (PVts-\eqn{\beta}). Ecological Indicators, 94, 367 379.
#'
#' @section Note:
#' In order to optimise the detections, it is advisable to make a smoothing
#' through the \link[ForesToolboxRS]{smootH} function before detecting changes. The smoothing will
#' allow to eliminate outliers that were not eliminated during the masking. If the object is a
#' time series, it must have the structure of a "ts", that is, the serie must have a start, end
#' and frequency. See \link[stats]{ts} for more details. In addition, in case the input is a
#' matrix, the first dimension must be rows*columns of the image, and the second dimension the
#' number of images.
#'
#' @section Details:
#' Regarding the example data, if the idea is to detect changes in 2008 (position 19),
#' then we will smooth the data only up to that position (i.e. \code{ndfi[1:19]}). This is to avoid
#' a smoothing of the monitor position value.
#'
#' @importFrom stats sd time ts
#' @importFrom graphics points abline polygon text grid legend plot
#' @importFrom raster values
#' @importFrom forecast na.interp
#'
#' @param x Vector, univariate time series (class "ts") or a matrix without NA's.
#' @param startm The start of the monitoring time.
#' @param endm The end of the monitoring time.
#' @param threshold The default threshold is 5 for photosynthetic vegetation,
#' while for indices such as NDVI and EVI the threshold is 3.
#' Please see Tarazona et al. (2018) for more details.
#' @param verbose This paramater is Logical. It Prints progress messages during execution.
#'
#' @examples
#' library(ForesToolboxRS)
#' library(raster)
#'
#' # Example
#' # Normalized Difference Fraction Index from 1990 to 2017, one index for each year.
#' ndfi <- c(0.86,0.93,0.97,0.91,0.95,0.96,0.91,0.88,0.92,0.89,
#'          0.90,0.89,0.91,0.92,0.89,0.90,0.92,0.84,0.46,0.20,
#'          0.27,0.22,0.52,0.63,0.61,0.67,0.64,0.86)
#'
#' # Detecting changes in 2008 (position 19)
#' # First, let´s apply a smoothing
#' ndfi_smooth <- ndfi
#' ndfi_smooth[1:19] <- smootH(ndfi[1:19])
#'
#' # Now, detect changes in 2008 (position 19)
#' cd <- pvts(x = ndfi_smooth, startm = 19, endm = 19, threshold = 5)
#' plot(cd)
#'
#' # Detect changes in 2008 (time serie)
#' ndfi_smt_ts <- ts(ndfi_smooth, start = 1990, end = 2017, frequency = 1)
#' cd <- pvts(x = ndfi_smt_ts, startm = 2008, endm = 2008,  threshold = 5)
#' plot(cd)
#'
#' @export
#'
pvts <- function(x, startm, endm, threshold = 5, verbose = FALSE) {

  if(!inherits(x, "numeric") | !inherits(x, "ts")) stop("x must be numeric or time serie (ts)", call. = TRUE)

  if (any(is.na(x))){
    stop("The object cannot contain NA. Please use the smootH function of this package to fill missing data and then smooth it.", call. = TRUE)
  }

  if (is(x, "vector")) {

    if(verbose){
      message(paste0(paste0(rep("*",10), collapse = ""), " Calculating the mean, standard deviation and lower limit" , paste0(rep("*",10), collapse = "")))
      print(model_algo)
    }

    mean.pvts <- mean(x[1:(startm-1)])

    std.pvts <- sd(x[1:(startm-1)])

    li <- mean.pvts - threshold*std.pvts

    value <- x[endm]

  } else if (is(x, "ts")) {

    if(verbose){
      message(paste0(paste0(rep("*",10), collapse = ""), " Calculating the mean, standard deviation and lower limit " , paste0(rep("*",10), collapse = "")))
      print(model_algo)
    }

    startm.pvts <- which(time(x) == startm)
    endm.pvts <- which(time(x) == endm)

    mean.pvts <- mean(x[1:(startm.pvts-1)])

    std.pvts <- sd(x[1:(startm.pvts-1)])

    li <- mean.pvts - threshold*std.pvts

    value <- x[endm.pvts]

  } else {

    stop(class(x), ' class is not supported', call. = TRUE)
  }

  if (value < li) {
    Monitoring_period = c(start = startm, end = endm)
    Breakpoint = c(Year_index = endm, value = value)
    Threshold = c(Threshold = threshold, Lower_limit = li)

    return(structure(list(Ts = x,
                          Breakpoint = Breakpoint,
                          Monitoring_period = Monitoring_period,
                          Threshold = Threshold), class = "pvts"))

  } else {
    Monitoring_period = c(start = startm, end = endm)
    Breakpoint = c(Year_index = NA, value = NA)
    Threshold = c(Threshold = threshold, Lower_limit = li)

    return(structure(list(Ts = x,
                          Breakpoint = Breakpoint,
                          Monitoring_period = Monitoring_period,
                          Threshold = Threshold), class = "pvts"))
  }

}
