# OISST

#' A subclass of SPNCRefClass for NOAA Optimum Interpolation (OI) Sea Surface Temperature \url{http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html}
#' 
#' Lon and Lat appear to be cell centers with 0.25 x 0.25 degree resolution
#' Lat is an asceinfing order
#'
#' @include SPNC.R
#' @export
OISSTRefClass <- setRefClass("OISSTRefClass",
    contains = "SPNCRefClass")


#' Craft subset indices into a ncdf array
#'
#' @name OISSTRefClass_subset_bbox
#'
#' @seealso \url{http://r.789695.n4.nabble.com/Lat-Lon-NetCDF-subset-td3394326.html}
#' @param bb numeric, four element bounding box [left, right, bottom, top]
#' @param zlev numeric zlevel to retrieve
#' @param time numeric time to retrieve
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
NULL
OISSTRefClass$methods(
   subset_bbox = function(bb = NULL, zlev = 1, time = 1){
      llon <- .self$lon("leading")
      llat <- .self$lat("leading")
      if (is.null(bb)){
         s <- list(start = c(lon=1,lat=1, zlev=zlev, time=time), 
            count = c(lon=length(llon), lat=length(llat), zlev=length(zlev), time=length(time)),
            bb = .self$get_extent() ) 
      } else {
         bb[1:2] <- to360(bb[1:2])
         ix <- spnc::find_interval(bb[0:2], llon)
         iy <- spnc::find_interval(bb[3:4], llat)
         # get these in order
         ix <- range(ix)
         iy <- range(iy)
         # make sure they fit within the dims of the data
         ix <- spnc::coerce_within(ix, c(1, length(llon)))
         iy <- spnc::coerce_within(iy, c(1, length(llat)))
         s <- .self$STEP
         xx <- llon[ix] + if (s[1] < 0) c(s[1],0) else c(0, s[1])
         yy <- llat[iy] + if (s[2] < 0) c(s[2],0) else c(0, s[2]) 
         
         s <- list(start = c(lon=ix[1], lat=iy[1], zlev=zlev, time=time), 
            count = c(lon=ix[2]-ix[1]+1, lat=iy[2]-iy[1]+1, zlev=length(zlev), time=length(time)),
            bb = c(range(xx), range(yy)) )
      }
      #.self$order_subset(s)
      s
   }) # subset_bbox

#' Get the step size (resolution) as a guess - the space between the first 
#' on and lat values.  Subclasses with unambiguous steps sizes, like L3SMI and  
#' MURSST, should override this method.
#'
#' @name OISSTRefClass_step
#' @return two element numeric vector of step size in x and y or NULL
#' These may be signed for descending values of lon and or lat
NULL
OISSTRefClass$methods(
   step = function(){
      c(0.250, 0.250)
   })


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
      crs = "+proj=longlat +datum=WGS84", flip = TRUE, time_fmt = "D%Y%j"){

      if (!all(what %in% .self$VARS))
         stop("one or more requested variable(s) not in data:", paste(what, collapse = "\n"))
      
      R <- OISST_get_raster(.self, what=what, layer = layer, bb = bb,
         crs = crs, flip = flip, time_fmt = time_fmt)
      return(R)
   }) # OISST_get_raster
   
######################### methods above ########################################
######################### functions below ######################################
#' Get a raster for OISSTRefClass
#' 
#' @export
#' @param NC OISSTRefClass object
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordinate reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
OISST_get_raster <- function(NC, what = NC$VARS[1], layer = 1,bb = NC$BB,
      crs = "+proj=longlat +datum=WGS84", flip = TRUE, time_fmt = "D%Y%j"){
   
   stopifnot(inherits(NC, "OISSTRefClass"))
   subnav <- NC$subset_bbox(bb)
   ext <- NC$get_extent(bb = subnav[['bb']])
   
   R <- raster::raster(
      nrow = subnav[['count']]['lat'], 
      ncol = subnav[['count']]['lon'],
      ext = ext,
      crs = crs) 
     
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


