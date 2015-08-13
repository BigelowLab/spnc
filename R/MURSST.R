# MURSST

#' A subclass of SPNCRefClass for Multi-scale Ultra-high Resolution Sea Surface Temperature \url{http://mur.jpl.nasa.gov/}
#' 
#' @include SPNC.R
#' @export
MURSSTRefClass <- setRefClass("MURSSTRefClass",
    contains = "SPNCRefClass")


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
MURSSTRefClass$methods(
   get_raster = function(what = 'analysed_sst', layer = 1, bb = .self$BB,
      crs = "+proj=longlat +datum=WGS84", 
      flip = TRUE, time_fmt = "D%Y%j"){
      R <- MURSST_get_raster(.self, what = what, layer = layer, bb = bb,
         crs = crs, flip = flip, time_fmt = time_fmt)
      return(R)
   }) # MURSST_get_raster



##### Methods above
##### functions below

#' Get a raster
#' 
#' @export
#' @param NC MURSSTRefClass object
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
MURSST_get_raster <- function(NC, what = 'analysed_sst', layer = 1, bb = NC$BB,
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


