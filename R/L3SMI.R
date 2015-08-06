# L3SMI.R

#' A subclass of SPNCRefClass for OBPG L3 Standard Mapped Image
#' 
#' @include SPNC.R
#' @export 
L3SMIRefClass <- setRefClass("L3SMIRefClass",
    contains = "SPNCRefClass",
    methods = list(
      init = function(...){
         callSuper(...)
         #if (.self$is_local()) {
            atts <- .self$get_global_atts()
            nm <- strsplit(atts[['product_name']], ".", fixed = TRUE)[[1]][1]
            .self$TIME <- as.POSIXct(nm, format = "A%Y%j", tz = "UTC") 
         #}   
      })
   )
   
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
L3SMIRefClass$methods(
   get_raster = function(what = .self$VARS[1], layer = 1,bb = .self$BB,
      crs = "+proj=longlat +datum=WGS84", flip = FALSE, time_fmt = "D%Y%j"){

      R <- L3SMI_get_raster(.self, what=what, layer = layer, bb = bb,
         crs = crs, flip = flip, time_fmt = time_fmt)
      return(R)
   }) # L3SMI_get_raster

   
##### methods above
##### functions below


#' Get a raster for L3SMIRefClass
#' 
#' @param NC L3SMIRefClass object
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordinate reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
L3SMI_get_raster <- function(NC, what = NC$VARS[1], layer = 1,bb = NC$BB,
      crs = "+proj=longlat +datum=WGS84", flip = FALSE, time_fmt = "D%Y%j"){
   
   stopifnot(inherits(NC, "L3SMIRefClass"))
   
   subnav <- NC$subset_coords(bb, lon=NC$LON, lat=NC$LAT)
   ext <- raster::extent(subnav[['bb']])
   
   R <- raster::raster(nrow = subnav[['count']][2], 
      ncol = subnav[['count']][1],  ext = ext, crs = crs) 
     
   if (NC$flavor[['local']] == TRUE){
      
      x <- ncdf4::ncvar_get(NC$NC, what[1], 
            start = subnav[['start']], 
            count = subnav[['count']] )
            
      R <- raster::raster(t(x), template = R)
      
   } else {
   
      getOneVar <- function(vname, NC, subnav,
         crs = "+proj=longlat +datum=WGS84"){
         ncdf4::ncvar_get(NC, vname, 
            start = c(subnav[['start']], layer[1]), 
            count = c(subnav[['count']], 1) )
      }
      getOneLayer <- function(layer, NC, subnav, what = names(NC[['var']])[[1]],
         crs = "+proj=longlat +datum=WGS84"){
         ncdf4::ncvar_get(NC, what, 
            start = c(subnav[['start']]), 
            count = c(subnav[['count']]) )
      }
      
      
      if (length(what) > 1){   
         X <- lapply(what, getOneVar, NC$NC, subnav, crs = crs, layer = layer[1] )
         for (w in names(X)) R <- raster::addLayer(R, 
            raster::raster(t(X[[w]]), template = R))
         names(R) <- what
      } else {
         X <- lapply(layer, getOneLayer, NC$NC, subnav, crs = crs, what = what[1]) 
         for (r in X) R <- raster::addLayer(R, 
            raster::raster(t(r), template = R))
         if (length(NC$TIME) > 1){
            names(R) <- format(NC$TIME[layer], time_fmt)
         } else {
            names(R) <- paste("layer", layer, sep = "_")
         } 
      } # multiple variable?
   
   } # local or OpenDAP?
   
   if (flip) R <- raster::flip(R,"y")
   
   return(R)
}
