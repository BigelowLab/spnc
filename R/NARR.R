# NARRRefClass.R

#' A subclass of SPNCRefClass for [NARR](http://www.esrl.noaa.gov/psd/data/gridded/data.narr.html) NCEP North American Regional Reanalysis
#' 
#' This NetCDF data is projected and can be ready by the \code{raster} package. Unlike
#' other members of the SPNCRefClass, this class stores the data as raster.  The get
#' \code{get_raster()} and \code{get_points()} functions are still exposed for consistency.
#'
#' @seealso \url{http://www.esrl.noaa.gov/psd/data/narr/format.html}
#' @seealso \url{http://www.esrl.noaa.gov/psd/data/gridded/data.narr.html}
#' @include SPNC.R
#' @export  
NARRRefClass <- setRefClass("NARRRefClass",
   contains = "SPNCRefClass",
   fields = list(
      R = 'ANY',
      R_x = 'ANY', 
      R_y = 'ANY',
      R_lon = 'ANY',
      R_lat = 'ANY'
   ),  # fields
   methods = list(
      init = function(nc){
         .self$field('R', raster::raster(nc$filename))
         #.self$BB <- as.vector(raster::extent(.self$R))
         xy <- raster::xyFromCell(.self$R, 1:raster::ncell(.self$R))
         nx <- ncol(.self$R)
         ny <- nrow(.self$R)
         .self$R_x <- matrix(xy[,1], nrow = ny, ncol = nx, byrow = TRUE)
         .self$R_y <- matrix(xy[,2], nrow = ny, ncol = nx, byrow = TRUE)
         ll <- rgdal::project(xy, raster::projection(.self$R), inv = TRUE) 
         .self$R_lon <- matrix(ll[,1], nrow = ny, ncol = nx, byrow = TRUE)
         .self$R_lat <- matrix(ll[,2], nrow = ny, ncol = nx, byrow = TRUE)
         .self$field("DIMS", ncdim_get(nc))
         .self$field("VARS", names(nc[['var']]))
         .self$field("STEP", .self$step()) 
         if ("zlev" %in% names(.self$DIMS)) .self$field("Z", nc[["dim"]][['zlev']][['vals']])
         hours <- spnc::nctime_get(nc, as_POSIXct = FALSE)
         time <- as.POSIXct("1800-01-01 00:00:0.0", tz = 'UTC') + (hours * 3600)
         .self$field("TIME", time)
         
      }) # methods
   ) # setRefClass


#' Compute extent given a bounding box
#' 
#' @name NARRRefClass_get_extent
#' @param bb the bounding box needed, if missing the current one is used
#' @param as_longlat logical if TRUE inverse project the exent to long lat values
#' @return a raster::Extent object
NULL
NARRRefClass$methods(
   get_extent = function(bb, as_longlat = FALSE){
      if(missing(bb)){
         # this doesn't work in the method - I think there is a namespace smashup
         # ext <- base::as.vector(raster::extent(.self$R))s
         # so we do it manually
         ext <- raster::extent(.self$R)
         ext <- c(ext@xmin, ext@xmax, ext@ymin, ext@ymax)
      } else {
         ext <- bb
      }
      if(as_longlat) {
         ext <- base::as.vector(rgdal::project(matrix(ext, ncol = 2, nrow = 2),
            proj = raster::projection(.self$R), inv = TRUE))
      }
      raster::extent(ext)
   })

   
#' Get the nominal step size (resolution) in km
#'
#' @seealso http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID221
#' @name NARRRefClass_step
#' @return two element numeric vector of step size in x and y or NULL
NULL
NARRRefClass$methods(
   step = function(){
      c(32.46341, 32.46341)
   })
 
#' Get the longitude values of the grid.  Note the grid is projected, so 
#  return a matrix - lon per cell.
#'
#' @name NARRRefClass_lon
#' @param what the desired location - 'center' or 'native' (default) which are the same.
#' @return numeric vector or NULL
NARRRefClass$methods(
   lon = function(what = c('native', 'center')[1]){
      xy <- xyFromCell(.self$R, 1:ncell(.self$R))
      ll <- project(xy, projection(.self$R), inv = TRUE)
      ll[,1]
   })
   

#' Get the longitude values of the grid
#'
#' @name NARRRefClass_lon
#' @param what the desired location - 'center' or 'native' (default) which are the same.
#' @return numeric vector or NULL
NARRRefClass$methods(
   lon = function(what = c('native', 'center')[1]){
     .self$R_lon
   })

#' Get the latitude values of the grid
#'
#' @name NARRRefClass_lat
#' @param what the desired location - 'center' or 'native' (default) which are the same.
#' @return numeric vector or NULL
NARRRefClass$methods(
   lat = function(what = c('native', 'center')[1]){
     .self$R_lat
   })
   

#' Convert bb to longlat
#'
#' @name NARRREfClass_bbox_project
#' @return 

#' Retrieve a raster, possibly a cropped
#' 
#' @name NARRRefClass_get_raster
#' @param what ignored
#' @param time numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top]
#'   defaults to the extent
#' @param crs character, the coordinate reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NARRRefClass$methods(
   get_raster = function(what = NULL, time = 1, bb = .self$BB,
      crs = raster::projection(.self$R), flip = FALSE){     
      NARR_get_raster(.self, what = what, bb = bb, time = time, flip = flip, crs = crs)
   }) # get_raster

#' Retrieve points of data
#' 
#' @name NARRRefClass_get_points
#' @param x numeric vector of x locations (longitude), or a matrix, data frame or list
#'    of x and y locations
#' @param y numeric vector of y locations - if NULL then these must be part of x
#' @param time numeric index or POSIXct 
#' @param ... further arguments for \code{raster::extract()}
#' @param return as specified by \code{raster::extract()}
NARRRefClass$methods(
   get_points = function(x, y = NULL, time = 1, ...){
      NARR_get_points(.self, x, y = y, time = time, ...)
   })

#' Retrieve a raster, possibly a cropped
#' 
#' @name NARRRefClass_get_raster
#' @param X NARRRefClass object
#' @param what ignored
#' @param time numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordinate reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @return a \code{raster::RasterBrick} or \code{raster::RasterLayer} object or NULL
NARR_get_raster <- function(X, bb = NULL, time = 1, what = NULL,
   crs = raster::projection(X$R), flip = FALSE){
   
   if (is.null(bb)) stop("bounding box, bb, must be provided")
   tix <- X$time_index(time)
   spbb <- spnc::bbox_to_polygon(bb, projstring = crs)
   r <- raster::crop(X$R[[time]], spbb)
   if (flip) r <- raster::flip(r, 'y')
   r
}


#' Get data values by specifying long/lat locations
#' 
#' @export 
#' @param X NARRRefClass object
#' @param x numeric vector of x locations (longitude), or a matrix, data frame or list
#'    of x and y locations
#' @param y numeric vector of y locations - if NULL then these must be part of x
#' @param time numeric index or POSIXct 
#' @param ... further arguments for \code{raster::extract()}
#' @param return as specified by \code{raster::extract()}
NARR_get_points <- function(X, x, y = NULL, time = 1, ...){
   tix <- X$time_index(time)
   xy <- cbind(grDevices::xy.coords(x, y=y)[c('x', 'y')])
   xy <- rgdal::project(xy, raster::projection(X$R))
   raster::extract(X$R[[tix]], xy, ...)
}
