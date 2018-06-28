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
         .self$rastertemplate = nc_template(.self$NC, proj = get_projection('longlat'))
      })
   )


#' Get the longitude locations
#'
#' @name L3SMIRefClass_lon
#' @param what the desired location - 'leading', 'trailing', 'center' or 'native' (default)
#' Values 'trailing' and 'leading' are about the direction lon/lat are stored.  For example,
#' if lon/lat are each stored in ascending order then the leading edge is bottom and left.
#' On the other hand, is lon descends (east-to-west) then the leading edge is
#' bottom and right.  Values 'native and 'center' are assumed the same - if that
#' is not a correct for a class, then a method override is required
#' @return numeric vector or NULL
NULL
L3SMIRefClass$methods(
   lon = function(what = c('native', 'leading', 'trailing', 'center')[1]){

      # native is assumed to the center
      #x <- if ('lon' %in% names(.self$NC$dim)) .self$NC$dim$lon$vals else  NULL
      #if (is.null(x)) return(NULL)
      if (!inherits(.self$NC, "ncdf4")) return(NULL)
      x <- get_loc(.self$NC, "lon")
      s <- .self$step/2
      switch(tolower(what[1]),
         'leading' = x + if (s[1] > 0) -s[1] else s[1],
         'trailing' = x + if(s[1] > 0) s[1] else -s[1],
         x)
   })

#' Get the latitude locations
#'
#' @name L3SMIRefClass_lat
#' @param what the desired location -  'leading', 'trailing', 'center' or 'native' (default)
#' Values 'trailing' and 'leading' are about the direction lon/lat are stored.  For example,
#' if lon/lat are each stored in ascending order then the leading edge is bottom and left.
#' On the other hand, is lat descends (north-to-south) then the leading edge is
#' top and left.  Values 'native and 'center' are likely the same.
#' @return numeric vector or NULL
NULL
L3SMIRefClass$methods(
   lat = function(what = c('native', 'leading', 'trailing', 'center')[1]){
      # native is assumed to the center
      #y <- if ('lat' %in% names(.self$NC$dim)) .self$NC$dim$lat$vals else NULL
      #if (is.null(y)) return(NULL)
      if (!inherits(.self$NC, "ncdf4")) return(NULL)
      y <- get_loc(.self$NC, "lat")
      s <- .self$step()/2
      switch(tolower(what[1]),
         'leading' = y + if (s[2] > 0) -s[2] else s[2],
         'trailing' = y + if (s[2] > 0) s[2] else -s[2],
         y)
   })

#' Get the step size (resolution) as reported in the global attributes.
#'
#' @name L3SMIRefClass_step
#' @return two element numeric vector of step size in x and y or NULL.
NULL
L3SMIRefClass$methods(
   step = function(){
      #c(.self$GATTS$longitude_step, .self$GATTS$latitude_step)
      #c((.self$GATTS[['easternmost_longitude']] - .self$GATTS[['westernmost_longitude']])/.self$DIMS[['lon']],
      #  (.self$GATTS[['northernmost_latitude']] - .self$GATTS[['southernmost_latitude']])/.self$DIMS[['lat']])
      get_res(.self$NC)
})



#' Compute extent given a bounding box
#'
#' @name L3SMIRefClass_get_extent
#' @param bb the bounding box needed, if missing the current one is used
#' @return a raster::Extent object
NULL
L3SMIRefClass$methods(
   get_extent = function(bb){
      s <- .self$step()
      if (missing(bb)){
         bb <- as.vector(raster::extent(.self$rastertemplate))
      }
      if (identical(bb, .self$BB)) return(raster::extent(bb))
      raster::extent(raster::crop(.self$rastertemplate, bb))

   })


#' Craft subset indices into a ncdf array
#'
#' @name L3SMIRefClass_subset_bbox
#'
#' @seealso \url{https://stat.ethz.ch/pipermail/r-help/2011-March/272641.html}
#' @param bb numeric, four element bounding box [left, right, bottom, top]
#' @return a list of
#' \itemize{
#'  \item{start indices in x and y}
#'  \item{counts in x and y}
#'  \item{bb vector of [left, right, bottom, top] pixel centers}
#'  \item{ext extent vector [left, right, bottom, top] pixel outer edges}
#'  }
NULL
L3SMIRefClass$methods(
   subset_bbox = function(bb = NULL){


      if (is.null(bb)) bb <- .self$get_extent()
      nav_from_bb(.self$rastertemplate, bb)
   }) # subset_bbox



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
      crs = get_projection("longlat"),
      flip = FALSE, time_fmt = "D%Y%j"){

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
#' @param X L3SMIRefClass object
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordinate reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
L3SMI_get_raster <- function(X, what = X$VARS[1], layer = 1, bb = X$BB,
      crs = get_projection("longlat"),
      flip = FALSE, time_fmt = "D%Y%j"){


    if (FALSE){
        layer = 1
        crs = get_projection("longlat")
        flip = FALSE
        time_fmt = "D%Y%j"
    }

    stopifnot(inherits(X, "L3SMIRefClass"))
    subnav <- X$subset_bbox(bb)
    template <- raster::crop(X$rastertemplate, bb)


    if (X$flavor[['local']] == TRUE){

       x <- ncdf4::ncvar_get(X$NC, what[1],
             start = subnav[['start']],
             count = subnav[['count']] )

       R <- raster::raster(t(x), template = template)

    } else {

       getOneVar <- function(vname, NC, subnav, crs = crs){
          ncdf4::ncvar_get(NC, vname,
             start = c(subnav[['start']], layer[1]),
             count = c(subnav[['count']], 1) )
       }
       getOneLayer <- function(layer, NC, subnav, what = names(NC[['var']])[[1]], crs = crs){
          ncdf4::ncvar_get(NC, what,
             start = c(subnav[['start']]),
             count = c(subnav[['count']]) )
       }


       if (length(what) > 1){
          RR <- lapply(what, getOneVar, X$NC, subnav, crs = crs, layer = layer[1] )
          for (w in names(RR)) R <- raster::addLayer(R, raster::raster(t(RR[[w]]), template = template))
          names(R) <- what
       } else {
          RR <- lapply(layer, getOneLayer, X$NC, subnav, crs = crs, what = what[1])
          R <- template
          for (r in RR) R <- raster::addLayer(R, raster::raster(t(r), template = template))
          if (length(X$TIME) > 1){
             names(R) <- format(X$TIME[layer], time_fmt)
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
#' @param X L3SMIRefClass object
#' @param x numeric, vector of longitude points or, if y is NULL, a matrix [x,y]
#'     or a data.frame or list with x,y elements
#' @param y numeric vector of lattitude points or NULL
#' @param what character one or more variable names or variable indices
#' @return numeric vector of values
L3SMI_get_points <- function(X, x, y = NULL, what = X$VARS[1]){

   stopifnot(inherits(X, "L3SMIRefClass"))

   p <- spnc::subset_points(x,y=y,lon=X$lon("leading"), lat=X$lat("leading"))

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
         nc = X$NC, what = what[1])
}

