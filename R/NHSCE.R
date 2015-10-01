# NHSCE

#' A subclass of SPNCRefClass for NOAA/NCDC Climate Data Record of snow cover 
#' extent 
#' 
#' Reference \url{https://climatedataguide.ucar.edu/climate-data/snow-cover-extent-northern-hemisphere-climate-data-record-rutgers}
#' 
#' Source \url{http://www.ncdc.noaa.gov/thredds/dodsC/cdr/snowcover/nhsce_v01r01_19661004_latest.nc}
#' 
#' See \url{http://climate.rutgers.edu/snowcover/docs.php?target=vis}
#' 
#' See \url{http://www1.ncdc.noaa.gov/pub/data/sds/cdr/docs/NH-SCE-CDR-CATBD-Rev2.pdf}
#'
#' @include SPNC.R
#' @export
NHSCERefClass <- setRefClass("NHSCERefClass",
    contains = "SPNCRefClass",
    fields = list(
      MASK = "ANY",
      AREA = "ANY"),
    methods = list(
          # init is a bit more specific and may be overwriten by subclasses
      # in here we deal with extracting from the ncdf resource
      init = function(nc){
         .self$field("DIMS", ncdim_get(nc))
         .self$field("VARS", names(nc[['var']])) 
         if ("longitude" %in% names(.self$VARS)) .self$field("LON", ncvar_get(nc, 'longitude'))
         if ("latitude" %in% names(.self$VARS)) .self$field("LAT", ncvar_get(nc, 'latitude'))
         if ("land" %in% names(.self$VARS)) .self$field("MASK", ncvar_get(nc, 'land'))
         if ("area" %in% names(.self$VARS)) .self$field("AREA", ncvar_get(nc, 'area'))
         if ("zlev" %in% names(.self$DIMS)) .self$field("Z", nc[["dim"]][['zlev']][['vals']])
         if ("time" %in% names(.self$DIMS)) .self$field("TIME", nctime_get(nc))
         return(TRUE)
      })
      
)


#' Get a raster
#' 
#' @name MURSSTRefClass_get_raster
#' @param what character one or more variable names
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps.
#'   If POSIXct then the layer at or just before the specified times are returned.
#'   See \code{\link{findInterval}} for details.
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NULL
NHSCERefClass$methods(
   get_raster = function(what = 'snow_cover_extent', layer = 1, bb = .self$BB,
      crs = "+proj=longlat +datum=WGS84", 
      flip = TRUE, time_fmt = "D%Y%j"){
      R <- NHSCE_get_raster(.self, what = what, layer = layer, bb = bb,
         crs = crs, flip = flip, time_fmt = time_fmt)
      return(R)
   }) # NHSCE_get_raster



##### Methods above
##### functions below

#' Get a raster
#' 
#' @export
#' @param NC NHSCERefClass object
#' @param what character one or more variable names
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps.
#'   If POSIXct then the layer at or just before the specified times are returned.
#'   See \code{\link{findInterval}} for details.
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NHSCE_get_raster <- function(NC, what = 'snow_cover_extent', layer = 1, bb = NC$BB,
      crs = "+proj=longlat +datum=WGS84", 
      flip = TRUE, time_fmt = "D%Y%j"){
  
      subnav <- NC$subset_coords(bb, NC$LON, NC$LAT)
      ext <- raster::extent(subnav[['bb']])
      
      if (inherits(layer, "POSIXct") && (length(NC$TIME) > 1) ){
         layer <- findInterval(layer, NC$TIME, all.inside = TRUE)
      }
      
      getOneVar <- function(vname, NC, subnav, layer = 1,
         crs = "+proj=longlat +datum=WGS84"){
         ncdf4::ncvar_get(NC, vname, 
            start = c(subnav[['start']], layer[1]), 
            count = c(subnav[['count']], 1) )
      }
      getOneLayer <- function(layer, NC, subnav, what = names(NC[['var']])[[1]],
         crs = "+proj=longlat +datum=WGS84"){
         ncdf4::ncvar_get(NC, what, 
            start = c(subnav[['start']], layer[1]), 
            count = c(subnav[['count']], 1) )
      }
      
      
      R <- raster::raster(nrow = subnav[['count']][2], 
         ncol = subnav[['count']][1], ext = ext, crs = crs)
        
      if (length(what) > 1){   
         X <- lapply(what, getOneVar, NC$NC, subnav, crs = crs, layer = layer[1] )
         for (w in names(X)) R <-  raster::addLayer(R, 
             raster::raster(t(X[[w]]), template = R))
         names(R) <- what
      } else {
         X <- lapply(layer, getOneLayer, NC$NC, subnav, crs = crs, what = what[1]) 
         for (r in X) R <-  raster::addLayer(R, 
             raster::raster(t(r), template = R))
         if (length(NC$TIME) > 1){
            names(R) <- format(NC$TIME[layer], time_fmt)
         } else {
            names(R) <- paste("layer", layer, sep = "_")
         } 
      }
      if (flip) R <-  raster::flip(R,"y")
      return(R)
   }


