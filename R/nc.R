#' Determine if cell locations are stored in ascending or descending order.
#'
#' @export
#' @param   x   netcd4 or raster object
#' @param   dir character either 'lon' or 'lat' (default)  d
#' @return      logical, TRUE if the specified cell locations ascend in value
ascends <- function(x, dir = 'lat'){
    if (inherits(x, 'ncdf4')){
        asc = diff(x$dim[[dir]]$vals[1:2]) > 0
    } else {
        asc = switch(dir,
            'lon' = xFromCol(x, 1:2),
            'lat' = yFromRow(x, 1:2))
        asc = diff(asc) > 0
    }
    asc
}

#' Get a 2 element vector of [dx, dy] resolution
#'
#' @export
#' @param   x   netcd4 or raster object
#' @param  flavor character vector from \code{spnc_flavor(x)}
#' @return      2 element numeric of resolution in [x,y]
get_res <- function(x, flavor = spnc_flavor(x)){
    if (inherits(x, 'ncdf4')){
        if (flavor$source == 'MURSST'){
            r = c(0.01, 0.01)
        } else {
            shape = get_shape(x)
            aa = spnc::ncglobal_atts(x)
            bb = unname(unlist(aa[c('westernmost_longitude', 'easternmost_longitude',
                       'southernmost_latitude', 'northernmost_latitude') ]))
            r = c((bb[2]-bb[1])/shape[1], (bb[4]-bb[3])/shape[2])
        }
    } else {
        r = res(x)
    }
    r
}

#' Get a 2 element vector of [nx, ny] dimensions
#'
#' @export
#' @param   x   netcd4 or raster object
#' @return      a two element vector of [nx, ny]
get_shape <- function(x){
    if (inherits(x, 'ncdf4')){
        r = c(x$dim$lon$len, x$dim$lat$len)
    } else {
        r = dim(x)
    }
    r
}

#' Get a vector of lon or lat cell locations
#'
#' @export
#' @param x     netcd4 or raster object
#' @param what  character either 'lon' or 'lat' (default)
#' @param  flavor character vector from \code{spnc_flavor(x)}
#' @return      a vector of the cell locations
get_loc <- function(x, what = 'lon', flavor = spnc_flavor(x)){
    if (inherits(x, 'ncdf4')){
        if (flavor$source == 'MURSST'){
            s = get_shape(x)
            r = get_res(x)
            if (what == 'lon'){
                v  = seq(from = -180 + r[1]/2, by = r[1], length = s[1])
            } else {
                v  = seq(from = -90 + r[2]/2, by = r[2], length = s[2])
            }
        } else {
            #
            aa = spnc::ncglobal_atts(x)
            bb = unname(unlist(aa[c('westernmost_longitude', 'easternmost_longitude',
                          'southernmost_latitude', 'northernmost_latitude') ]))
            res = get_res(x)
            v = switch(what[1],
                'lon'  = seq(from = bb[1] + 0.5 * res[1],
                             to = bb[2] - 0.5 * res[1],
                             by = res[1]),
                'lat'  = if(ascends(x))
                            seq(from = bb[3] + 0.5 * res[2],
                                to = bb[4] - 0.5 * res[2],
                                by = res[2])
                         else
                            seq(from = bb[4] - 0.5 * res[2],
                                to = bb[3] + 0.5 * res[2],
                                by = -res[2])
                )
        }
    } else {
        v  = switch(what[1],
                    'lon' = raster::xFromCol(x, 1:ncol(x)),
                    'lat' = raster::yFromRow(x, 1:nrow(x)) )
    }
    v
}

#' Generate the navigation elements required to subset an ncdf4 object.
#' Returns a list of start indices and count and an adjusted bounding box
#'
#' @export
#' @param r     raster object
#' @param bb    a 4 element vector describing the bounding box [left, right, bottom, top]
#' @param flip  character either 'x', 'y' or 'xy' to flip the orientation
#' @return      a 3 element vector of start [ix,iy], count, [nx, ny], and
#'    an adjusted bb [left, right, bottom, top] and ext (extent) same as bb
nav_from_bb <- function(r, bb, flip = ''){

    cr      = raster::crop(r, bb)
    cx      = get_loc(cr, 'lon')
    cy      = get_loc(cr, 'lat')
    rx      = get_loc(r, 'lon')
    ry      = get_loc(r, 'lat')
    if (grepl('x', flip[1], fixed = TRUE)){
        ix = which.min(abs(rev(rx) - cx[1]))
    } else {
        ix = which.min(abs(rx - cx[1]))
    }
    if (grepl('y', flip[1], fixed = TRUE)) {
        iy = which.min(abs(rev(ry) - cy[length(cy)]))
    } else {
        iy = which.min(abs(ry - cy[1]))
    }

    ext    = as.vector(extent(cr))
    list(
        start = c(ix,iy),
        count = c(ncol(cr), nrow(cr)),
        bb   = ext,
        ext  = ext)
}

#' Create a raster template from the contents of a ncdf4 object
#'
#' @export
#' @param   x ncdf4 object
#' @param   proj CRS projection string
#' @return  raster
nc_template <- function(x,
    proj = get_projection("longlat")){

    xres    = get_res(x)
    xshape  = get_shape(x)
    xlon    =  get_loc(x, 'lon')
    xlat    = get_loc(x, 'lat')

    raster::raster(
        xmn = min(xlon) - xres[1]/2,
        xmx = max(xlon) + xres[1]/2,
        ymn = min(xlat) - xres[2]/2,
        ymx = max(xlat) + xres[2]/2,
        ncol = xshape[1],
        nrow = xshape[2],
        crs = proj)
}

#' Crop a ncdf4 object
#'
#' @export
#' @param   x ncdf4 object
#' @param varname   character variable name 'sst' by default
#' @param bb        4 element bounding box [left, right, bottom, top]
#' @param template  raster template
#' @return raster
nc_crop <- function(x, varname = 'sst', bb = c(-180, 180, -90, 90),
    template = nc_template(x)){
    # ## S4 method for signature 'missing'
    # raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90,
    #       crs, ext, resolution, vals=NULL)

    nv = nav_from_bb(template, bb)
    raster::raster(t(ncvar_get(NC, varname, start = nv$start, count = nv$count)),
        crs = raster::projection(template),
        xmn = nv$bb[1], xmx = nv$bb[2],
        ymn = nv$bb[3], ymx = nv$bb[4])
}




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
