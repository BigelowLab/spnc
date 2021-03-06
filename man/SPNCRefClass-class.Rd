% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPNC.R
\docType{class}
\name{SPNCRefClass-class}
\alias{SPNCRefClass-class}
\alias{SPNCRefClass}
\title{An object that connects to NetCDF (file or OpeNDAP) using ncdf4}
\description{
An object that connects to NetCDF (file or OpeNDAP) using ncdf4
}
\section{Fields}{

\describe{
\item{\code{flavor}}{a character vector of at least 'source' and 'type'}

\item{\code{NC}}{the ncdf4 class object}

\item{\code{BB}}{the 4 element bounding box}

\item{\code{DIM}}{the dimensions}

\item{\code{VAR}}{the variable names}

\item{\code{STEP}}{2 element vector of x and y step size (or NULL)}

\item{\code{Z}}{possibly a vector of z-values}

\item{\code{TIME}}{possibly a vector of times}
}}


