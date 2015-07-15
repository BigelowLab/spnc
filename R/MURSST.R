# MURSST

MURSSTRefClass <- setRefClass("MURSSTRefClass",
    contains = "SPNCRefClass")


#' Get a raster
#' 
#' @name MURSSTRefClass_get_raster
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NULL
MURSSTRefClass$methods(
   get_raster = function(what = .self$VARS, layer = 1,
      crs = "+proj=longlat +datum=WGS84", flip = TRUE, time_fmt = time_fmt){
  
      subnav <- .self$subset_coords(.self$BB, .self$LON, .self$LAT)
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
   }) # MURSST_get_raster
