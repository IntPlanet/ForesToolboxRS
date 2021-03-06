% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvtsRaster.R
\name{pvtsRaster}
\alias{pvtsRaster}
\title{Change detection using the PVts-\eqn{\beta} approach (raster version)}
\usage{
pvtsRaster(x, startm, endm, threshold = 5, img, vf = FALSE, verbose = FALSE)
}
\arguments{
\item{x}{Matrix, RasterStack or Rasterbrick without NA's.}

\item{startm}{The start of the monitoring time.}

\item{endm}{The end of the monitoring time.}

\item{threshold}{The default threshold is 5 for photosynthetic vegetation,
while for indices such as NDVI and EVI the threshold is 3.
Please see Tarazona et al. (2018) for more details.}

\item{img}{The image of the position immediately before the monitoring start,
i.e. the "start-1" position (in case "x" is a matrix). The matrix \code{x} and
the image \code{img} must have the same dimension (i.e. rows and cols).}

\item{vf}{If the monitoring is with Photosynthetic Vegetation series,
then switch to \code{TRUE}.}

\item{verbose}{This paramater is Logical. It Prints progress messages during execution.}
}
\description{
This algorithm will allow to detect disturbances in the forests using
all the available Landsat set. In fact, it can also be run with sensors
such as MODIS.
}
\section{References}{

Tarazona, Y., Mantas, V.M., Pereira, A.J.S.C. (2018). Improving tropical
deforestation detection through using photosynthetic vegetation time
series (PVts-\eqn{\beta}). Ecological Indicators, 94, 367 379.
}

\section{Note}{

In order to optimise the detections, it is advisable to make a smoothing before
detecting changes. So the matrix or the stack must go through a smoothing first.
The smoothing will allow to eliminate outliers that were not eliminated during the masking.
See \link[ForesToolboxRS]{smootH} for more details. In addition, in case the input is a
matrix, the first dimension must be rows*columns of the image, and the second dimension the
number of images.
}

\examples{
library(ForesToolboxRS)
library(raster)
data(FTdata)

# Detect changes in 2008 (position 19) using a raster(RasterStack)
cd <- pvtsRaster(x=serie_pv, startm=19, endm=19,  threshold= 5)

}
\author{
Yonatan Tarazona Coronel
}
