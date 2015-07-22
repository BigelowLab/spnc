# OISST

#' A subclass of SPNCRefClass for NOAA Optimum Interpolation (OI) Sea Surface Temperature \url{http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html}
#' 
#' @include SPNC.R
#' @export
OISSTRefClass <- setRefClass("OISSTRefClass",
    contains = "SPNCRefClass")



#' Get a raster
#' 
#' @name OISSTRefClass_get_raster
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NULL
OISSTRefClass$methods(
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
         ncol = subnav[['count']][1], ext = ext, crs = crs)
        
      if (length(what) > 1){   
         X <- lapply(what, getOneVar, .self$NC, subnav, crs = crs, layer = layer[1] )
         for (w in names(X)) R <- raster::addLayer(R, raster::
            raster(t(X[[w]]), template = R))
         names(R) <- what
      } else {
         X <- lapply(layer, getOneLayer, .self$NC, subnav, crs = crs, what = what[1]) 
         for (r in X) R <- raster::addLayer(R, 
            raster::raster(t(r), template = R))
         if (length(.self$TIME) > 1){
            names(R) <- format(x$TIME[layer], time_fmt)
         } else {
            names(R) <- paste("layer", layer, sep = "_")
         } 
      }
      if (flip) R <- raster::flip(R,"y")
      return(R)
   }) # OISST_get_raster


#' Retrieve the subset coordinates 
#' 
#' @name SPNCRefClass_subset_coords
#' @param lon numeric vector of lons (ascending order please!) to select from
#' @param lat numeric vector of lats (ditto)
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
NULL
SPNCRefClass$methods(
   subset_coords = function(bb = .self$BB, lon = c(-180,180), lat = c(-90,90)){
   
      if (is.null(bb)){
         return( list(start = c(1,1), count = c(length(lon), length(lat)),
            bb = c( range(lon), range(lat) ) ) )
      }
      
      yDescends <- (lat[1] - lat[2]) < 0
      
      ix <- findInterval(bb[0:2], lon, all.inside = TRUE)
      if (ix[2] < length(lon)) ix[2] <- ix[2] + 1
      
      if (yDescends){
         iy <- findInterval(bb[3:4], lat, all.inside = TRUE)
         if (iy[2] < length(lat)) iy[2] <- iy[2] + 1
      } else {
         iy <- findInterval(bb[3:4], rev(lat), all.inside = TRUE)
         if (iy[2] < length(lat)) iy[2] <- iy[2] + 1
      }
      
      
      list(start = c(ix[1], iy[1]), 
         count = c(ix[2]-ix[1]+1, iy[2]-iy[1]+1),
         bb = c(lon[ix], lat[iy])  )
   })