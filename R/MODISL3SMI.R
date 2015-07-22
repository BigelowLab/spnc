# MODISL3SMI

#' A subclass of SPNCRefClass for Modis Aqua \url{http://oceancolor.gsfc.nasa.gov/cms/}
#' 
#' @include SPNC.R
#' @export 
MODISL3SMIRefClass <- setRefClass("MODISL3SMIRefClass",
    contains = "SPNCRefClass",
    methods = list(
      init = function(...){
         callSuper(...)
         atts <- ncdf4::ncatt_get(.self$NC, varid = 0)
         nm <- strsplit(atts[['product_name']], ".", fixed = TRUE)[[1]][1]
         .self$TIME <- as.POSIXct(nm, format = "A%Y%j", tz = "UTC")   
      })
   )
   


#' Get a raster
#' 
#' @name MODISL3SMIRefClass_get_raster
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NULL
MODISL3SMIRefClass$methods(
   get_raster = function(what = .self$VARS, layer = 1,bb = .self$BB,
      crs = "+proj=longlat +datum=WGS84", flip = TRUE, time_fmt = time_fmt){

      subnav <- .self$subset_coords(bb, .self$LON, .self$LAT)
      ext <- raster::extent(subnav[['bb']])
      
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
         ncol = subnav[['count']][1],  ext = ext, crs = crs)
        
      if (length(what) > 1){   
         X <- lapply(what, getOneVar, .self$NC, subnav, crs = crs, layer = layer[1] )
         for (w in names(X)) R <- raster::addLayer(R, 
            raster::raster(t(X[[w]]), template = R))
         names(R) <- what
      } else {
         X <- lapply(layer, getOneLayer, .self$NC, subnav, crs = crs, what = what[1]) 
         for (r in X) R <- raster::addLayer(R, 
            raster::raster(t(r), template = R))
         if (length(x$TIME) > 1){
            names(R) <- format(.self$TIME[layer], time_fmt)
         } else {
            names(R) <- paste("layer", layer, sep = "_")
         } 
      }
      if (flip) R <- raster::flip(R,"y")
      return(R)
   }) # MODISL3SMI_get_raster


