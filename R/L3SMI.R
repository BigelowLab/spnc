# L3SMI.R




#' Test if an NCDF contains L3SMI data.
#' 
#' @export
#' @param x ncdf4 object or SPNCRefClass
#' @return logical
is_L3SMI <- function(x){
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
      ok <- mgrepl('Level-3 Standard Mapped Image', 
         atts[['title']], fixed = TRUE)   
   ok
}


#' A subclass of SPNCRefClass for OBPG L3 Standard Mapped Image
#' 
#' @include SPNC.R
#' @export 
L3SMIRefClass <- setRefClass("L3SMIRefClass",
    contains = "SPNCRefClass",
    methods = list(
      init = function(...){
         callSuper(...)
         atts <- .self$get_global_atts()
         nm <- strsplit(atts[['product_name']], ".", fixed = TRUE)[[1]][1]
         .self$TIME <- as.POSIXct(nm, format = "A%Y%j", tz = "UTC")   
      })
   )



#' Get the step size (resolution) as reported in the global attributes.
#'
#' @name L3SMIRefClass_step
#' @return two element numeric vector of step size in x and y or NULL.
#' These may be signed for descending values of lon and or lat
NULL
L3SMIRefClass$methods(
   step = function(){
      atts <- .self$get_global_atts()
      lat <- .self$lat()
      onelat <- if (lat[2]-lat[1] > 0) ? 1 else -1
      c(atts[['longitude_step']], onelat*atts[['latitude_step']])
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
L3SMIRefClass$methods(
   get_raster = function(what = .self$VARS[1], layer = 1,bb = .self$BB,
      crs = "+proj=longlat +datum=WGS84", flip = FALSE, time_fmt = "D%Y%j"){

      if (!all(what %in% .self$VARS))
         stop("one or more requested variable(s) not in data:", paste(what, collapse = "\n"))
         
      R <- L3SMI_get_raster(.self, what=what, layer = layer, bb = bb,
         crs = crs, flip = flip, time_fmt = time_fmt)
      return(R)
   }) # L3SMI_get_raster

   
#' Get points
#' 
#' @name L3SMIRefClass_get_points
#' @param x numeric, vector of longitude points or, if y is NULL, a matrix [x,y] 
#'     or a data.frame or list with x,y elements
#' @param y numeric vector of lattitude points or NULL  
#' @param what character one or more variable names or variable indices
#' @return numeric vector of values
NULL
L3SMIRefClass$methods(
   get_points = function(x, y = NULL, what = .self$VARS[1]){
      
      if (!all(what %in% .self$VARS))
         stop("one or more requested variable(s) not in data:", paste(what, collapse = "\n"))
      L3SMI_get_points(.self, x=x, y=y, what=what)
   }) # L3SMI_get_points

##### methods above
##### functions below


#' Get a raster for L3SMIRefClass
#' 
#' @export
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
   subnav <- NC$subset_bbox(bb)
   ext <- NC$get_extent(bb = subnav[['bb']])
   
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



#' Get points for L3SMIRefClass
#' 
#' @export
#' @param NC L3SMIRefClass object
#' @param x numeric, vector of longitude points or, if y is NULL, a matrix [x,y] 
#'     or a data.frame or list with x,y elements
#' @param y numeric vector of lattitude points or NULL
#' @param what character one or more variable names or variable indices
#' @return numeric vector of values
L3SMI_get_points <- function(NC, x, y = NULL, what = NC$VARS[1]){
   
   stopifnot(inherits(NC, "L3SMIRefClass"))
   
   p <- spnc::subset_points(x,y=y,lon=NC$lon("leading"), lat=NC$lat("leading"))
   
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

