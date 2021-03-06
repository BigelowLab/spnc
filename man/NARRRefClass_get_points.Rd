% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NARR.R
\name{NARRRefClass_get_points}
\alias{NARRRefClass_get_points}
\title{Retrieve points of data}
\arguments{
\item{x}{numeric vector of x locations (longitude), or a matrix, data frame or list
of x and y locations}

\item{y}{numeric vector of y locations - if NULL then these must be part of x}

\item{time}{numeric index or POSIXct}

\item{...}{further arguments for \code{raster::extract()}}
}
\value{
as specified by \code{raster::extract()}
}
\description{
Retrieve points of data
}
