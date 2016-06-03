# OISST

#' Test if an NCDF contains OISST data.
#' 
#' @export
#' @param x ncdf4 object or SPNCRefClass
#' @return logical
is_OISST <- function(x){
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
      ok <- mgrepl('Daily-OI-V2, Final, Data', atts[['title']], fixed = TRUE)   
   ok
}


#' A subclass of SPNCRefClass for NOAA Optimum Interpolation (OI) Sea Surface Temperature \url{http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html}
#' 
#' Lon and Lat appear to be cell centers with 0.25 x 0.25 degree resolution
#' Lat is an ascending order [-90, 90] and while Lon in mapped to [0,360] 
#' bounding box requests must follow the [0,360] form.
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
      
      # not every OISST is made the same - some have zlev and others don't
      # so we match our subnav to the actual dims
      VNAMES <- names(.self$NC[['dim']])
      if ('zlev' %in% VNAMES) {
         vnames <- c("lon", "lat", "zlev", "time") 
      } else {
         vnames <- c("lon", "lat", "time")
      }
      s[['start']] <- s[['start']][vnames]
      s[['count']] <- s[['count']][vnames]
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
#' Requests can be made for multiple variables at multiple times.  In such cases
#' it is impotant to know that the order of rasters is each time for each variable
#' requested.  Thus requesting variables 'v1', 'v2', 'v3' at times t1 and t2
#' results in a stack with v1-t1, v1-t2, v2-t1, v2-t2, v3-t1, v3-t2 order.
#'
#' @name OISSTRefClass_get_raster
#' @param what character one or more variable names or variable indices
#' @param time numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NULL
OISSTRefClass$methods(
   get_raster = function(what = .self$VARS[1], time = 1,bb = .self$BB,
      crs = "+proj=longlat +datum=WGS84", flip = TRUE){

      if (!all(what %in% .self$VARS))
         stop("one or more requested variable(s) not in data:", paste(what, collapse = "\n"))
      
      R <- OISST_get_raster(.self, what=what, time = time, bb = bb,
         crs = crs, flip = flip)
      return(R)
   }) # OISST_get_raster
   
######################### methods above ########################################
######################### functions below ######################################
#' Get a raster for OISSTRefClass
#' 
#' @export
#' @param NC OISSTRefClass object
#' @param what character one or more variable names or variable indices
#' @param time numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordinate reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
OISST_get_raster <- function(NC, what = NC$VARS[1], time = 1, bb = NC$BB,
      crs = "+proj=longlat +datum=WGS84", flip = TRUE){
   
   stopifnot(inherits(NC, "OISSTRefClass"))
   subnav <- NC$subset_bbox(bb)
   ext <- NC$get_extent(bb = subnav[['bb']])
   
   R <- raster::raster(
      nrow = subnav[['count']]['lat'], 
      ncol = subnav[['count']]['lon'],
      ext = ext,
      crs = crs) 
      
   # for each what
   # if (is_contiguous(time))
   #    get_block_of_what_when
   # else  
   #     for each time get a block
   
   tix <- NC$time_index(time)
   is_contiguous <- all(diff(tix) == 1)
   for (w in what){
      if (is_contiguous && (length(tix) > 1)){
         snav <- subnav
         snav[['start']][['time']] <- tix[1]
         snav[['count']][['time']] <- length(tix)
         x <- ncdf4::ncvar_get(NC$NC, w, start = c(snav[['start']]), count = c(snav[['count']]) )
         if (length(dim(x)) == 3){
            for(i in seq_len(dim(x)[3])) R <- addLayer(R, raster::raster(t(x[,,i]), template = R))
         } else {
            R <- raster::addLayer(R, raster::raster(t(x), template = R))
         }
      } else {
         for (ix in tix){
            snav <- subnav
            snav[['start']][['time']] <- ix
            snav[['count']][['time']] <- 1
            x <- ncdf4::ncvar_get(NC$NC, w, start = c(snav[['start']]), count = c(snav[['count']]) )
            R <- addLayer(R, raster::raster(t(x), template = R))
         }
      }
   }   
   
   if (flip) R <- raster::flip(R,"y")
   
   return(R)
}


