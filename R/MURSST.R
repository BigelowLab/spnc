# MURSST


#' Test if an NCDF contains MURSST data.
#' 
#' @export
#' @param x ncdf4 object or SPNCRefClass
#' @return logical
is_MURSST <- function(x){
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
      ok <- mgrepl('MUR',atts[['title']], fixed = TRUE) 
   ok  
}


#' A subclass of SPNCRefClass for Multi-scale Ultra-high Resolution Sea Surface Temperature \url{http://mur.jpl.nasa.gov/}
#' 
#' @include SPNC.R
#' @export
MURSSTRefClass <- setRefClass("MURSSTRefClass",
    contains = "SPNCRefClass")


#' Get the step size (resolution) as a guess - the space between the first 
#' lon and lat values.  Subclasses with unambiguous steps sizes, like L3SMI and  
#' MURSST, should override this method.
#'
#' @name MURSSTClass_step
#' @return two element numeric vector of step size in x and y or NULL
#' These may be signed for descending values of lon and or lat
NULL
MURSSTRefClass$methods(
   step = function(){
      c(0.01, 0.01)
   })

  
#' Compute extent given a bounding box
#' 
#' @name MURSSTRefClass_get_extent
#' @param bb the bounding box needed, if missing the current one is used
#' @return a raster::Extent object
NULL
MURSSTRefClass$methods(
   get_extent = function(bb){
      if (missing(bb)){
         x <- range(.self$NC$dim$lon$vals)
         y <- range(.self$NC$dim$lat$vals)
         s <- .self$step()
         #bb <- c(range(.self$lon()[c(1,nx)]), range(.self$lat()[c(1,ny)]))
         bb <- c( x + c(-s[1], s[1]), y + c(-s[2], s[2]))
         return(raster::extent(bb))
      }
      if (identical(bb, .self$BB)) return(raster::extent(.self$BB))
      
      llon <- .self$lon("leading")
      llat <- .self$lat("leading")
      ix <- find_interval(bb[1:2], llon)
      ix[ix < 1] <- 1
      iy <- find_interval(bb[3:4], llat)
      iy[ix < 1] <- 1
      
      s <- .self$STEP
      xx <- llon[ix] + if (s[1] < 0) c(s[1],0) else c(0, s[1])
      yy <- llat[iy] + if (s[2] < 0) c(s[2],0) else c(0, s[2])       
      raster::extent(c(range(xx), range(yy)) )
   })   

#' Compute indicies (start and count) for individual [lon,lat] points
#' 
#' @name MURSSTRefClass_subset_points
#' @param x either a vector of longitude or a 2d matrix [lon,lat] or a data.frame
#'    with lon,lat columns
#' @param y numeric vector, if x is a vector, otherwise NULL
#' @return a list of start/count elements, one per xy input.  For inputs that
#'    fall beyond the bounds of the data then NA is returned.
NULL
MURSSTRefClass$methods(
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


#' Get points
#' 
#' @name MURSSTRefClass_get_points
#' @param x numeric, vector of longitude points or, if y is NULL, a matrix [x,y] 
#'     or a data.frame or list with x,y elements
#' @param y numeric vector of lattitude points or NULL  
#' @param what character one or more variable names or variable indices
#' @param layer numeric layer index
#' @return numeric vector of values
NULL
MURSSTRefClass$methods(
   get_points = function(x, y = NULL, what = .self$VARS[1], layer = 1){
      
      if (!all(what %in% .self$VARS))
         stop("one or more requested variable(s) not in data:", paste(what, collapse = "\n"))
      MURSST_get_points(.self, x=x, y=y, what=what, layer = layer)
   }) # L3SMI_get_points
   

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
  
      stopifnot(inherits(NC, "MURSSTRefClass"))
      subnav <- NC$subset_bbox(bb)
      ext <- NC$get_extent(bb = subnav[['bb']])
      
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
         ncdf4::ncvar_get(NC, varid = what, 
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


#' Get points for MURSSTRefClass
#' 
#' @export
#' @param NC MURSSTRefClass object
#' @param x numeric, vector of longitude points or, if y is NULL, a matrix [x,y] 
#'     or a data.frame or list with x,y elements
#' @param y numeric vector of lattitude points or NULL
#' @param what character one or more variable names or variable indices
#' @param layer the layer index
#' @return numeric vector of values
MURSST_get_points <- function(NC, x, y = NULL, what = NC$VARS[1], layer = 1){
   
   stopifnot(inherits(NC, "MURSSTRefClass"))
   
   p <- NC$subset_points(x,y=y, layer = layer)
      
      
   sapply(p, 
         function(x, nc = NULL, what = 1){
            if (!is.list(x) && is.na(x)){
               r <- NA
            } else {
               r <- ncdf4::ncvar_get(nc, what, 
                  start = x[['start']], 
                  count = x[['count']] )
            }
            r   
         },
         nc = NC$NC, what = what[1])
}


