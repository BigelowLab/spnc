#' Retrieve a named list of global attributes
#'
#' @export
#' @param NC a ncdf4 object
#' @param rm_pattern character a pattern of characters to remove from the 
#'    attribute names.  By default 'NC_GLOBAL.'.  Set to "" or NA to skip
#' @param fixed logical by default TRUE but see \code{grepl}
#' @return named vector of global attributes
ncglobal_atts <- function(NC, rm_pattern = 'NC_GLOBAL.', fixed = TRUE){
   d <- if (!is.null(NC)) ncdf4::ncatt_get(NC, varid = 0) else NULL
   if (!is.null(d)){
      if (!is.na(rm_pattern) && (nchar(rm_pattern) > 0)){
         names(d) <- gsub("NC_GLOBAL.", "", names(d), fixed = fixed)
      }
   }
   return(d)
}

#' Retrieve a vector of dimensions
#'
#' @export
#' @param NC a ncdf4 object
#' @return a named vector of dimesnions
ncdim_get <- function(NC){
   stopifnot('dim' %in% names(NC))
   sapply(NC[['dim']], '[[', 'len')
}

#' Retrieve a vector of variable names
#'
#' @export
#' @param NC a ncdf4 object
#' @return a named vector variable names
ncvarname_get <- function(NC){
   stopifnot('dim' %in% names(NC))
   names(NC[['var']])
}


#' Retrieve a list of dimension vectors, one for each variable
#'
#' @export
#' @param NC a ncdf4 object
#' @return a named list of variable dimension vectors
ncvardim_get <- function(NC){
   stopifnot('dim' %in% names(NC))
   vn <- names(NC[['var']])
   names(vn) <- vn
   get_vardim <- function(nm, NC = NULL){
      d <- NC[['var']][[nm]]
      get_vardim_one <- function(x){
         len <- x[['len']]
         names(len) <- x[['name']]
         return(len)
      }
      
      dims <- d[['dim']]
      sapply(dims, get_vardim_one)
   }
   lapply(vn, get_vardim, NC = NC)
}

#' Retrieve a vector of timestamps for a multilayer NC object or NULL otherwise
#'
#' @export
#' @param NC a ncdf4 class object
#' @param name the name of the time-associated variable, by default 'time'
#' @param as_POSIXct logical, if TRUE then convert to POSIXct
#' @return a numeric vector of timestamps (possibly POSIXct) or NULL if none
nctime_get <- function(NC, name = 'time', as_POSIXct = TRUE){
   
   d <- ncdim_get(NC)
   if (!("time" %in% names(d)) ) return(NULL)
   
   v <- NC[["dim"]][['time']][['vals']]
   
   if (as_POSIXct){
      u <- NC[["dim"]][['time']][['units']]
      
      if (grepl(" since ", u, fixed = TRUE)){
         secsperhour <- 60 * 60
         secsperday <- 24 * secsperhour
         spaces <- gregexpr(" since ", u, fixed = TRUE)[[1]]
         incr <- substring(u, 1, spaces[1]-1)
         dt <- substring(u, spaces[1] + attr(spaces, 'match.length'))
         t0 <- as.POSIXct(dt, tz = "UTC")
         v <- switch(incr,
            "days" =  t0 + (v * secsperday),
            "seconds" = t0 + v,
            "hours" = t0 + (v * secsperhour),
            t0 + v)
      } else {
         cat("nctime_get: unknown time format for conversion to POSIXct\n")
      }  
   }
   invisible(v)
}

#' Compute the mean step.  If the step size is provided in the global attributes
#' then it is advised you use that information.  See \code{ncglobal_atts()}
#'
#' @export
#' @param NC a ncdf4 class object
#' @param dim_names a two element character vector providing the names of the 
#'   x and y coordinates.  By default \code{c('lon', 'lat')}
#' @return a 2 element numeric of the mean step size - always positive
nc_step <- function(NC,
    dim_names = c('lon','lat')){
    
    xn <- dim_names[1]
    yn <- dim_names[2]
    x <- abs(c(mean(diff(NC$dim[[xn]]$vals)), mean(diff(NC$dim[[yn]]$vals))))
    names(x) <- dim_names
    x
}


#' Compute subset coordinates given a bounding box - suitable for use with 
#' \code{ncvar_get()}.  This is sort of a generic function that may not work for 
#' all NetCDF files. Coordinates within the files are assumed to be pixel centers.
#'
#' Pixel center coordinates are typically specified within the 'dim' attribute in the 
#' ncdf4 object.  Step size between pixel centers is taken to be the mean for 
#' all pixels by default, but this can be overridden.
#' 
#' @export
#' @param NC a ncdf4 class object
#' @param bb the 4 element bounding box [left, right, bottom, top]
#' @param step one or two element numeric of the step size. 
#' @param dim_names a two element character vector providing the names of the 
#'   x and y coordinates.  By default \code{c('lon', 'lat')}
#' @return a 5 element list of 
#'  \itemize{
#'      \item{start a 2 element vector of indices to start}
#'      \item{count a 2 element vector number of columns/rows}
#'      \item{ext a 4 element extent that accounts for the pixels sizes - the envelope}
#'      \item{y_ascends logical TRUE if y values are stored low to high}
#'      \item{step two element step size (aka resolution)}
#'   }
nc_subset <- function(NC,
    bb = c(-77, -51.5, 37.9, 56.7),
    step = nc_step(NC),
    dim_names = c('lon','lat')){

    if (length(step) < 2) step <- c(step, step)
    s <- step/2.0

    xn <- dim_names[1]
    yn <- dim_names[2]
    
    # some ncdf files have lat from south to north, other north to south
    ascends <- diff(NC$dim[[yn]]$vals[1:2]) > 0

    # note we use inclusive boundaries - and select just the first and last
    ix <- which(NC$dim[[xn]]$vals >= bb[1] & NC$dim[[xn]]$vals <= bb[2])
    ix <- ix[c(1, length(ix))]

    iy <- which(NC$dim[[yn]]$vals >= bb[3] & NC$dim[[yn]]$vals <= bb[4])
    iy <- iy[c(1, length(iy))]

    # these are 'real' coordinates of pixel centers
    x <- NC$dim[[xn]]$vals[ix]
    nx <- length(x)
    y <- NC$dim[[yn]]$vals[iy]
    ny <- length(y)
    if (ascends){
        ext <- c( x[c(1, nx)],     y[c(1, ny)]  ) + c(-s[1], s[1], -s[2], s[2])
    } else {
        ext <- c( x[c(1, nx)], rev(y[c(1, ny)]) ) + c(-s[1], s[1], -s[2], s[2])
    }
    list(
        start       = c(ix[1], iy[1]),
        count       = c(ix[2] - ix[1] + 1,
                        iy[2] - iy[1] + 1),
        ext      = ext,
        y_ascends   = ascends,
        step        = step)
}
