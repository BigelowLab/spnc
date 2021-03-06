% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L3SMI.R
\name{L3SMI_get_raster}
\alias{L3SMI_get_raster}
\title{Get a raster for L3SMIRefClass}
\usage{
L3SMI_get_raster(X, what = X$VARS[1], layer = 1, bb = X$BB,
  crs = get_projection("longlat"), flip = FALSE,
  time_fmt = "D\%Y\%j")
}
\arguments{
\item{X}{L3SMIRefClass object}

\item{what}{character one or more variable names or variable indices}

\item{layer}{numeric vector either a 1-based indices or POSIXct timestamps}

\item{bb}{a 4 element bounding box vector [left, right, bottom, top], defaults
to [-180, 180, -90, 90]}

\item{crs}{character, the coordinate reference system to apply}

\item{flip}{logical if TRUE then flip the raster in the y direction}

\item{time_fmt}{if multiple time layers are returned, this controls the layer names}
}
\value{
a \code{raster::brick} or \code{raster::layer} object or NULL
}
\description{
Get a raster for L3SMIRefClass
}
