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
    contains = "SPNCRefClass",
    methods  = list(
        init = function(...){
            # this extent is from downloading an example via FTP and then reading
            # with x = raster::raster(ncdf_file)
            # then dput(as.vector(raster::extent(x)))
            # this builds an empty raster identical to that read from file as far
            # as extent and res goes.
            e = c(-179.99500549324, 180.005000000076, -89.9949978636508, 89.9949978636508)
            d = c(17999, 36000)
            .self$rastertemplate = raster::raster(
                ext = raster::extent(e),
                ncol = d[2], nrow = d[1],
                crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

            }
        )
    )


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


#' Craft subset indices into a ncdf array
#'
#' @name MURSSTRefClass_subset_bbox
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
MURSSTRefClass$methods(
   subset_bbox = function(bb = NULL){
      if (is.null(bb)) bb <- .self$get_extent()
      nav_from_bb(.self$rastertemplate, bb, flip = 'y')
   }) # subset_bbox




#' Compute extent given a bounding box
#'
#' @name MURSSTRefClass_get_extent
#' @param bb the bounding box needed, if missing the current one is used
#' @return a raster::Extent object
NULL
MURSSTRefClass$methods(
   get_extent = function(bb){
      s <- .self$step()
      if (missing(bb)){
         bb <- as.vector(raster::extent(.self$rastertemplate))
      }
      if (identical(bb, .self$BB)) return(raster::extent(.self$BB))
      raster::extent(raster::crop(.self$rastertemplate, bb))
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
   }) # MURSST_get_points


##### Methods above
##### functions below

#' Get a raster
#'
#' @export
#' @param X MURSSTRefClass object
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
MURSST_get_raster <- function(X, what = 'analysed_sst', layer = 1, bb = X$BB,
      crs = get_projection("longlat"), #"+proj=longlat +datum=WGS84",
      flip = TRUE, time_fmt = "D%Y%j"){

      stopifnot(inherits(X, "MURSSTRefClass"))
      subnav <- X$subset_bbox(bb)
      R <- raster::crop(X$rastertemplate, bb)
      ext <- X$get_extent(bb = subnav[['bb']])

      if (inherits(layer, "POSIXct") && (length(X$TIME) > 1) ){
         layer <- findInterval(layer, X$TIME, all.inside = TRUE)
      }

      getOneVar <- function(vname, X, subnav, layer = 1,
         crs = "+proj=longlat +datum=WGS84"){
         ncdf4::ncvar_get(X, vname,
            start = c(subnav[['start']], layer[1]),
            count = c(subnav[['count']], 1) )
      }
      getOneLayer <- function(layer, X, subnav, what = names(X[['var']])[[1]],
         crs = "+proj=longlat +datum=WGS84"){
         ncdf4::ncvar_get(X, varid = what,
            start = c(subnav[['start']], layer[1]),
            count = c(subnav[['count']], 1) )
      }


      #R <- raster::raster(nrow = subnav[['count']][2],
      #  ncol = subnav[['count']][1], ext = ext, crs = crs)

      if (length(what) > 1){
         xx <- lapply(what, getOneVar, X$NC, subnav, crs = crs, layer = layer[1] )
         for (w in names(X)) R <-  raster::addLayer(R,
             raster::raster(t(xx[[w]]), template = R))
         names(R) <- what
      } else {
         xx <- lapply(layer, getOneLayer, X$NC, subnav, crs = crs, what = what[1])
         for (r in xx) R <-  raster::addLayer(R,
             raster::raster(t(r), template = R))
         if (length(X$TIME) > 1){
            names(R) <- format(X$TIME[layer], time_fmt)
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


