# HMODISL3SMI

HMODISL3SMIRefClass <- setRefClass("HMODISL3SMIRefClass",
    contains = "SPNCRefClass",
    methods = list(
      init = function(...){
         callSuper(...)
      })
   )
   
#' Get a raster
#' 
#' @name HMODISL3SMIRefClass_get_raster
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NULL
HMODISL3SMIRefClass$methods(
   get_raster = function(what = .self$VARS, layer = 1,
      crs = "+proj=longlat +datum=WGS84", flip = TRUE, 
      time_fmt = "D%Y%j"){

      subnav <- .self$subset_coords()
      ext <- raster::extent(subnav[['bb']])
      
      getOneVar <- function(vname, NC, subnav, layer = 1,
         crs = "+proj=longlat +datum=WGS84"){
         ncvar_get(NC, vname, 
            start = c(subnav[['start']], layer[1]), 
            count = c(subnav[['count']], 1) )
      }
      getOneLayer <- function(layer, NC, subnav, what = names(NC[['var']])[[1]],
         crs = "+proj=longlat +datum=WGS84"){
         ncvar_get(NC, what, 
            start = c(subnav[['start']], layer[1]), 
            count = c(subnav[['count']], 1) )
      }
      
      
      R <- raster(nrow = subnav[['count']][2], ncol = subnav[['count']][1], 
         ext = ext, crs = crs)
        
      if (length(what) > 1){   
         X <- lapply(what, getOneVar, x$NC, subnav, crs = crs, layer = layer[1] )
         for (w in names(X)) R <- addLayer(R, raster(t(X[[w]]), template = R))
         names(R) <- what
      } else {
         X <- lapply(layer, getOneLayer, x$NC, subnav, crs = crs, what = what[1]) 
         for (r in X) R <- addLayer(R, raster(t(r), template = R))
         if (length(x$TIME) > 1){
            names(R) <- format(x$TIME[layer], time_fmt)
         } else {
            names(R) <- paste("layer", layer, sep = "_")
         } 
      }
      if (flip) R <- flip(R,"y")
      return(R)
   }) # HMODISL3SMI_get_raster


#' Retrieve the subset coordinates 
#' 
#' @name HMODISL3SMIRefClass_subset_coords
#' @param bb 4 element vector in standard form [left, right, bottom, top] [-180,180]
#' @param lon numeric vector of lons (ascending order please!) to select from
#' @param lat numeric vector of lats (ditto)
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
NULL
SPNCRefClass$methods(
   subset_coords = function(bb = .self$BB, 
      lon = .self$LON, 
      lat = .self$LAT){
   
      if (is.null(bb)){
         return( list(start = c(1,1), count = c(length(lon), length(lat)),
            bb = c( range(lon), range(lat) ) ) )
      }
      
      # [-180,180] lon so we are fine here
      #> head(HMO$LON)
      #[1] -179.9792 -179.9375 -179.8958 -179.8542 -179.8125 -179.7708
      # [90, -90] so we have descending lat
      #> head(HMO$LAT)
      #[1] 89.97917 89.93750 89.89583 89.85417 89.81250 89.77083
      
      yDescends <- (lat[2] - lat[1]) < 0
      
      ix <- findInterval(bb[0:2], lon, all.inside = TRUE)
      if (ix[2] < length(lon)) ix[2] <- ix[2] + 1
      
      if (yDescends){
         iy <- findInterval(bb[3:4], rev(lat), all.inside = TRUE)
         if (iy[2] < length(lat)) iy[2] <- iy[2] + 1
      } else {
         iy <- findInterval(bb[3:4], lat, all.inside = TRUE)
         if (iy[2] < length(lat)) iy[2] <- iy[2] + 1
      }
      
      
      list(start = c(ix[1], iy[1]), 
         count = c(ix[2]-ix[1]+1, iy[2]-iy[1]+1),
         bb = c(lon[ix], lat[iy])  )
   })