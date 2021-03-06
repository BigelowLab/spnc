% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPNC.R
\name{SPNCRefClass_subset_coords}
\alias{SPNCRefClass_subset_coords}
\title{Retrieve the subset coordinates}
\arguments{
\item{lon}{numeric vector of lons (ascending order please!) to select from}

\item{lat}{numeric vector of lats (ditto)}
}
\value{
a list of \code{start} indices in x and y, \code{counts} in x and y and
   a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
}
\description{
Retrieve the subset coordinates
}
