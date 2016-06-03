# CPCUGBDP


#' Test if an NCDF contains CPCUGBDP data.
#' 
#' @export
#' @param x ncdf4 object or SPNCRefClass
#' @return logical
is_CPCUGBDP <- function(x){
   ok <- FALSE
   if (inherits(x, "SPNCRefClass")){
      atts <- try(ncglobal_atts(x$NC))
   } else if(inherits(x, "ncdf4")){
      atts <- try(ncglobal_atts(x))
   } else {
      warning("input must be either SPNCRefClass or ncdf4 class object")
      return(ok)
   }
   natts <- names(atts) <- tolower(names(atts))
   if ('title' %in% natts)
      ok <- mgrepl("Unified Gauge-Based Analysis of Daily Precipitation", 
      atts[['title']], fixed = TRUE)  
   ok 
}


#' A subclass of SPNCRefClass for CPC .25x.25 Daily US Unified Gauge-Based Analysis of Precipitation \url{http://www.esrl.noaa.gov/psd/data/gridded/data.unified.daily.conus.html}
#' 
#' Lat is an ascending order [-90, 90] and while Lon in mapped to [0,360] 
#' bounding box requests must follow the [0,360] form.
#'
#' @include SPNC.R
#' @export
CPCUGBDPRefClass <- setRefClass("CPCUGBDPRefClass",
    contains = "SPNCRefClass",
    methods = list(
      init = function(nc){
         callSuper(nc)
         # hours since 1800-01-01 00:00:0.0
         hours <- nctime_get(nc, as_POSIXct = FALSE)
         time <- as.POSIXct("1800-01-01 00:00:0.0", tz = 'UTC') + (hours * 3600)
         .self$field("TIME", time)
      })
   )
#' Get the step size (resolution) as a guess - the space between the first 
#' on and lat values.
#'
#' @name CPCUGBDPRefClass_step
#' @return two element numeric vector of step size in x and y or NULL
#' These may be signed for descending values of lon and or lat
NULL
CPCUGBDPRefClass$methods(
   step = function(){
      c(0.25, 0.25)
   })


#' Craft subset indices into a ncdf array
#'
#' @name CPCUGBDPRefClass_subset_bbox
#'
#' @seealso \url{http://r.789695.n4.nabble.com/Lat-Lon-NetCDF-subset-td3394326.html}
#' @param bb numeric, four element bounding box [left, right, bottom, top]
#' @param time numeric time to retrieve
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
NULL
CPCUGBDPRefClass$methods(
   subset_bbox = function(bb = NULL, time = 1){
      
      if(!all(diff(time)==1)) stop("subset_bbox: time indicies must be ascending contiguous")
      
      llon <- .self$lon("leading")
      llat <- .self$lat("leading")
      if (is.null(bb)){
         s <- list(start = c(lon=1,lat=1, time=time[1]), 
            count = c(lon=length(llon), lat=length(llat), time=length(time)),
            bb = .self$get_extent() ) 
      } else {
         ix <- spnc::find_interval(bb[0:2], llon)
         iy <- spnc::find_interval(bb[3:4], llat)
         # get these in order
         ix <- range(ix)
         iy <- range(iy)
         # make sure they fit within the dims of the data
         ix <- spnc::coerce_within(ix, c(1, length(llon)))
         iy <- spnc::coerce_within(iy, c(1, length(llat)))
         s <- .self$STEP
         xx <- llon[ix] + if (s[1] < 0) c(s[1],0) else c(0, s[1])
         yy <- llat[iy] + if (s[2] < 0) c(s[2],0) else c(0, s[2]) 
         
         s <- list(start = c(lon=ix[1], lat=iy[1], time=time[1]), 
            count = c(lon=ix[2]-ix[1]+1, lat=iy[2]-iy[1]+1, time=length(time)),
            bb = c(range(xx), range(yy)) )
      }
      
      vnames <- c("lon", "lat", "time")
      s[['start']] <- s[['start']][vnames]
      s[['count']] <- s[['count']][vnames]
      s
   }) # subset_bbox



#' Compute indicies (start and count) for individual [lon,lat] points
#' 
#' @name CPCUGBDPRefClass_subset_points
#' @param x either a vector of longitude or a 2d matrix [lon,lat] or a data.frame
#'    with lon,lat columns
#' @param y numeric vector, if x is a vector, otherwise NULL
#' @return a list of start/count elements, one per xy input.  For inputs that
#'    fall beyond the bounds of the data then NA is returned.
NULL
CPCUGBDPRefClass$methods(
   subset_points = function(x, y = NULL, layer = 1){
   
      # first verify x and y
      if (is.null(y)){
         # this covers list(x=..., y=...) or data.frame(x=...,y=...)
         if (is.list(x)){
            if (!all(names(x) %in% c("x", "y"))){
               stop("if y is NULL, then x must have elements names x and y")
            } else {
               y <- x$y
               x <- x$x
            }
         } else if (is.matrix(x)){
            if (!all(colnames(x) %in% c("x", "y"))){
               stop("if x is a matrix, then x must have columns named x and y")
            } else {
               y <- x[,'y']
               x <- x[,'x']
            }
         } else {
            stop("if y is NULL, then x must be list, data.frame or matrix with x and y elements")
         }
      } # is y provided or NULL
   
      stopifnot(length(x) == length(y))
      
      # first we get extent, leading cell edges and indices
      e <- .self$get_extent()
      llon <- .self$lon("leading")
      llat <- .self$lat("leading")   
      ix <- spnc::find_interval(x, llon)
      iy <- spnc::find_interval(y, llat)
      # now convert to list of (start = [ix,iy], count = [nx,ny])
      # because we have points count is always 1
      iz <- mapply(function(x, y, layer = 1){         
            list(start = c(x, y, layer), count = c(1,1,1))
         },
         ix,iy, layer = layer, SIMPLIFY = FALSE)
      # now we have to flag any beyond extent as NA
      # we could use the indices for the leading edge (where ix or iy < 1)
      # but we would still need to dub around with the trailing edge
      # so it is just as easy to use the input x and y instead of ix and iy
      ixna <- (x < e[1]) | (x > e[2])
      iyna <- (y < e[3]) | (y > e[4])
      iz[ixna | iyna] <- NA
      return(iz)
   }) #subset_points
   
   
#' Get a raster
#' 
#' @name CPCUGBDPRefClass_get_raster
#' @param what character one or more variable names
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param time numeric vector either a 1-based indices or POSIXct timestamps.
#'   If POSIXct then the layer at or just before the specified times are returned.
#'   See \code{\link{findInterval}} for details. Must resolve to contiguous ascending 
#'   indices.
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NULL
CPCUGBDPRefClass$methods(
   get_raster = function(what = 'precip', time = 1, bb = .self$BB,
      crs = "+proj=longlat +datum=WGS84", 
      flip = TRUE, time_fmt = "D%Y%j"){
      R <- CPCUGBDP_get_raster(.self, what = what, time = time, bb = bb,
         crs = crs, flip = flip, time_fmt = time_fmt)
      return(R)
   }) # CPCUGBDP_get_raster


#' Get points
#' 
#' @name CPCUGBDPRefClass_get_points
#' @param x numeric, vector of longitude points or, if y is NULL, a matrix [x,y] 
#'     or a data.frame or list with x,y elements
#' @param y numeric vector of lattitude points or NULL  
#' @param what character one or more variable names or variable indices
#' @return numeric vector of values
NULL
CPCUGBDPRefClass$methods(
   get_points = function(x, y = NULL, what = .self$VARS[1], layer = 1){
      
      if (!all(what %in% .self$VARS))
         stop("one or more requested variable(s) not in data:", paste(what, collapse = "\n"))
      CPCUGBDP_get_points(.self, x=x, y=y, what=what, layer = layer)
   }) # L3SMI_get_points
   

##### Methods above
##### functions below

#' Get a raster
#' 
#' @export
#' @param X CPCUGBDPRefClass object
#' @param what character one or more variable names
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param time numeric vector either a 1-based indices or POSIXct timestamps.
#'   If POSIXct then the layer at or just before the specified times are returned.
#'   See \code{\link{findInterval}} for details.  Must resolve to contiguous ascending 
#'   indices.
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
CPCUGBDP_get_raster <- function(X, what = 'precip', time = 1, bb = X$BB,
      crs = "+proj=longlat +datum=WGS84", 
      flip = TRUE, time_fmt = "D%Y%j"){
  
      stopifnot(inherits(X, "CPCUGBDPRefClass"))
      if (inherits(time, "POSIXct") && (length(X$TIME) > 1) ){
         time <- findInterval(time, X$TIME, all.inside = TRUE)
      }
      subnav <- X$subset_bbox(bb, time = time)
      ext <- X$get_extent(bb = subnav[['bb']])
      
      R <- raster::raster(nrow = subnav[['count']]['lat'], 
         ncol = subnav[['count']]['lon'], ext = ext, crs = crs)
      x <- ncdf4::ncvar_get(X$NC, varid = what, 
            start = c(subnav[['start']]), 
            count = c(subnav[['count']]) )     
      
      dims <- dim(x)
      if (length(dims) > 2){
         r <- raster::raster(nrow = subnav[['count']]['lat'], 
            ncol = subnav[['count']]['lon'], ext = ext, crs = crs)
         for (itime in seq_len(dims[3])) 
            r <- raster::addLayer(r,raster::raster(t(x[,,itime]), template = R))
         R <- r
      } else {
         R <- raster::raster(t(x), template = R)
      }
         
      if (length(X$TIME) > 1){
         names(R) <- format(X$TIME[time], time_fmt)
      } else {
         names(R) <- paste("time", time, sep = "_")
      } 
      
      if (flip) R <-  raster::flip(R,"y")
      return(R)
   }


#' Get points for CPCUGBDPRefClass
#' 
#' @export
#' @param X CPCUGBDPRefClass object
#' @param x numeric, vector of longitude points or, if y is NULL, a matrix [x,y] 
#'     or a data.frame or list with x,y elements
#' @param y numeric vector of lattitude points or NULL
#' @return numeric vector of values
CPCUGBDP_get_points <- function(X, x, y = NULL, what = X$VARS[1], layer = 1){
   
   stopifnot(inherits(NC, "CPCUGBDPRefClass"))
   
   p <- X$subset_points(x,y=y, layer = layer)
      
      
   sapply(p, 
         function(x, nc = null, what = 1){
            if (!is.list(x) && is.na(x)){
               r <- NA
            } else {
               r <- ncdf4::ncvar_get(nc, what, 
                  start = x[['start']], 
                  count = x[['count']] )
            }
            r   
         },
         nc = X$NC, what = what[1])
}


