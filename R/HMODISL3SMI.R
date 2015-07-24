# HMODISL3SMI

#' A subclass of SPNCRefClass for Modis Aqua via \url{http://thredds.jpl.nasa.gov}
#'
#' Warning - metadata looks corrupted - see get_global_atts() method
#' 
#' @include SPNC.R
#' @export
HMODISL3SMIRefClass <- setRefClass("HMODISL3SMIRefClass",
    contains = "SPNCRefClass",
    methods = list(
      init = function(...){
         callSuper(...)
         if (.self$is_local()) {
            atts <- .self$get_global_atts()
            nm <- strsplit(atts[['Product Name']], ".", fixed = TRUE)[[1]][1]
            .self$TIME <- as.POSIXct(nm, format = "A%Y%j", tz = "UTC") 
         }   
      })
   )
   
#' Get global attributes
#'
#' @name HMODISL3SMIRefClass_get_global_atts
#' @return a named list of global attributes or NULL
NULL
HMODISL3SMIRefClass$methods(
   get_global_atts = function(){
      d <- if (!is.null(.self$NC)) ncdf4::ncatt_get(.self$NC, varid = 0) else NULL
      names(d) <- gsub('%2520', ' ', names(d), fixed = TRUE)
      return(d)
   })
   


#' Get a raster
#' 
#' @name L3SMIRefClass_get_raster
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NULL
HMODISL3SMIRefClass$methods(
   get_raster = function(what = .self$VARS[1], layer = 1,bb = .self$BB,
      crs = "+proj=longlat +datum=WGS84", flip = FALSE, time_fmt = "D%Y%j"){

      
      subnav <- .self$subset_coords(bb, lon=.self$LON, lat=.self$LAT)
      ext <- raster::extent(subnav[['bb']])
      
      R <- raster::raster(nrow = subnav[['count']][2], 
         ncol = subnav[['count']][1],  ext = ext, crs = crs) 
        
      # local means just one thing in the NC plus a palette
      # http means possibly multiple layers 
      if (.self$flavor[['local']] == TRUE){
         
         x <- ncdf4::ncvar_get(.self$NC, what[1], 
               start = subnav[['start']], 
               count = subnav[['count']] )
               
         R <- raster::raster(t(x), template = R)
         
      } else {
      
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
         
         
         if (length(what) > 1){   
            X <- lapply(what, getOneVar, .self$NC, subnav, crs = crs, layer = layer[1] )
            for (w in names(X)) R <- raster::addLayer(R, 
               raster::raster(t(X[[w]]), template = R))
            names(R) <- what
         } else {
            X <- lapply(layer, getOneLayer, .self$NC, subnav, crs = crs, what = what[1]) 
            for (r in X) R <- raster::addLayer(R, 
               raster::raster(t(r), template = R))
            if (length(.self$TIME) > 1){
               names(R) <- format(.self$TIME[layer], time_fmt)
            } else {
               names(R) <- paste("layer", layer, sep = "_")
            } 
         } # multiple variable?
      
      } # local or OpenDAP?
      
      if (flip) R <- raster::flip(R,"y")
      
      return(R)
   }) # L3SMI_get_raster
