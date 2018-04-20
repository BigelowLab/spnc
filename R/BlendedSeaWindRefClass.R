#' Test if an NCDF contains Blended Sea Winds data.
#' 
#' @export
#' @param x ncdf4 object or SPNCRefClass
#' @return logical
is_BSW <- function(x){
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
      ok <- mgrepl(glob2rx('*Blended*Sea Surface Winds*'), 
         atts[['title']], fixed = FALSE)   
   ok
}

#' A subclass of SPNCRefClass for Blended Sea Winds
#' 
#' @include SPNC.R
#' @export 
BlendedSeaWindsRefClass <- setRefClass("BlendedSeaWindsRefClass",
    contains = "SPNCRefClass",
    methods = list(
      init = function(...){
         callSuper(...)
         #atts <- .self$get_global_atts()
         #nm <- strsplit(atts[['product_name']], ".", fixed = TRUE)[[1]][1]
         #.self$TIME <- as.POSIXct(nm, format = "A%Y%j", tz = "UTC")   
      })
   )


#' Craft subset indices into a ncdf array
#'
#' @name BlendedSeaWindsBlendedSeaWindsRefClass_subset_bbox
#'
#' @seealso \url{https://stat.ethz.ch/pipermail/r-help/2011-March/272641.html}
#' @param bb numeric, four element bounding box [left, right, bottom, top]
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
NULL
BlendedSeaWindsRefClass$methods(
   subset_bbox = function(bb = NULL){
      llon <- .self$lon("leading")
      llat <- .self$lat("leading")
      if (is.null(bb)){
         return(
            list(start = c(lon=1,lat=1, zlev = 1, time = 1), 
                count = c(lon=length(llon), lat=length(llat),
                    zlev = 1, time = 1),
                bb = .self$get_extent() ) 
            )
      }
      
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
      
      list(start = c(lon=ix[1], lat=iy[1], zlev = 1, time = 1), 
         count = c(lon=ix[2]-ix[1]+1, lat=iy[2]-iy[1]+1, zlev = 1, time = 1),
         bb = c(range(xx), range(yy)) )
   }) # subset_bbox
   


#' Get a raster
#' 
#' @name BlendedSeaWindsRefClass_get_raster
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [0, 360, -90, 90]
#' @param crs character, the coordinate reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
BlendedSeaWindsRefClass$methods(
    get_raster = function(what = .self$VARS[1], layer = 1, bb = .self$BB,
        crs = "+proj=longlat +datum=WGS84", flip = TRUE, time_fmt = "D%Y%j"){
    
        if (!all(what %in% .self$VARS))
            stop("one or more requested variable(s) not in data:", 
                paste(what, collapse = "\n"))
        
        BSW_get_raster(.self, what=what, layer = layer, bb = bb,
         crs = crs, flip = flip, time_fmt = time_fmt)
    })
    
#' Get a raster for BlendedSeaWindsRefClass
#' 
#' @export
#' @param X BlendedSeaWindsRefClass object
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [0, 360, -90, 90]
#' @param crs character, the coordinate reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
BSW_get_raster <- function(X, what = X$VARS[1], layer = 1, bb = X$BB,
      crs = "+proj=longlat +datum=WGS84", flip = TRUE, time_fmt = "D%Y%j"){
   
    
    getOneLayer <- function(layer, NC, subnav, what = names(NC[['var']])[[1]],
       crs = "+proj=longlat +datum=WGS84"){
       #subnav[['start']][['time']] <- layer[1]
       ncdf4::ncvar_get(NC, what, 
          start = subnav[['start']], 
          count = subnav[['count']] )
        }

    if (TRUE){
        what = X$VARS[1]
        layer = as.POSIXct("2015-08-15 10:00:00", tz = 'UTC')
        bb = X$BB
        crs = "+proj=longlat +datum=WGS84"
        flip = TRUE
        time_fmt = "D%Y%j"
    }
   
    stopifnot(inherits(X, "BlendedSeaWindsRefClass"))
    
    
    if (inherits(layer, "POSIXct") ){
        layer <- findInterval(layer, X$TIME, all.inside = TRUE)
    }
      
    subnav <- X$subset_bbox(bb)
    subnav[['start']][['time']] <- layer[1]
    
    ext <- X$get_extent(bb = subnav[['bb']])
    
    R <- raster::raster(nrow = subnav[['count']][2], 
        ncol = subnav[['count']][1],  ext = ext, crs = crs) 
      
    if (X$flavor[['local']] == TRUE){
       x <- ncdf4::ncvar_get(X$NC, what[1], 
             start = subnav[['start']], 
             count = subnav[['count']] )
       R <- raster::raster(t(x), template = R)
    } else {       
        xx <- lapply(layer, getOneLayer, X$NC, subnav, crs = crs, what = what[1]) 
        for (r in xx) R <- raster::addLayer(R, raster::raster(t(r), template = R))
        if (length(X$TIME) > 1){
           names(R) <- format(X$TIME[layer], time_fmt)
        } else {
           names(R) <- paste("layer", layer, sep = "_")
        } 
    
    } # local or OpenDAP?
    
    if (flip) R <- raster::flip(R,"y")
    
    return(R)
}
