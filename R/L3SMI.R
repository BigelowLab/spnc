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
NULL
L3SMIRefClass$methods(
   step = function(){
      c(.self$GATTS$longitude_step, .self$GATTS$latitude_step) 
      #c((.self$GATTS[['easternmost_longitude']] - .self$GATTS[['westernmost_longitude']])/.self$DIMS[['lon']],
      #  (.self$GATTS[['northernmost_latitude']] - .self$GATTS[['southernmost_latitude']])/.self$DIMS[['lat']])
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
         rx <- range(.self$lon('center'))
         ry <- range(.self$lat('center'))
         bb <- c(
            rx[1] - s[1]/2, 
            rx[2] + s[1]/2,
            ry[1] - s[2]/2, 
            ry[2] + s[2]/2)
         
      }
      if (identical(bb, .self$BB)) return(raster::extent(bb))
      
      llon <- .self$lon("center")
      llat <- .self$lat("center")
      ix <- find_interval(bb[1:2], llon)
      ix[ix < 1] <- 1
      iy <- find_interval(bb[3:4], llat)
      iy[iy < 1] <- 1
      
      xx <- llon[ix] + c(-s[1], s[1])/2
      yy <- llat[iy] + c(-s[2], s[2])/2       
      raster::extent(c(range(xx), range(yy)) )
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
      return(nc_subset(.self$NC, bb = bb, step = .self$STEP))
       
      # llon <- .self$lon("leading")
      # llat <- .self$lat("leading")
      # if (is.null(bb)){
      #    return(
      #       list(start = c(lon=1,lat=1), count = c(lon=length(llon), lat=length(llat)),
      #          bb = .self$get_extent() ) 
      #       )
      # }
      # 
      # ix <- range(spnc::find_interval(bb[0:2], llon))
      # iy <- range(spnc::find_interval(bb[3:4], llat))
      # 
      # # make sure they fit within the dims of the data
      # ix <- spnc::coerce_within(ix, c(1, length(llon)))
      # iy <- spnc::coerce_within(iy, c(1, length(llat)))
      # s <- .self$STEP
      # #xx <- llon[ix] + if (s[1] < 0) c(s[1],0) else c(0, s[1])
      # #yy <- llat[iy] + if (s[2] < 0) c(s[2],0) else c(0, s[2]) 
      # 
      # xx <- range(.self$lon("center")[ix])
      # yy <- range(.self$lat("center")[iy])
      # 
      # list(start = c(lon=ix[1], lat=iy[1]), 
      #    count = c(lon=ix[2]-ix[1]+1, lat=iy[2]-iy[1]+1),
      #    bb = c(xx, yy),
      #    ext = c(xx + c(-0.5,0.5)*s[1], yy + c(-0.5, 0.5)*s[2]) )
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
      crs = "+proj=longlat +datum=WGS84", flip = FALSE, time_fmt = "D%Y%j"){
   

    if (FALSE){
        layer = 1
        crs = "+proj=longlat +datum=WGS84"
        flip = FALSE
        time_fmt = "D%Y%j"
    }
    
    stopifnot(inherits(X, "L3SMIRefClass"))
    subnav <- X$subset_bbox(bb)
        
    template <- raster::raster(
        #nrow = subnav[['count']][2], 
        #ncol = subnav[['count']][1],  
        resolution = X$step(),
        ext = raster::extent(subnav[['ext']]),
        crs = crs) 
      
    if (X$flavor[['local']] == TRUE){
       
       x <- ncdf4::ncvar_get(X$NC, what[1], 
             start = subnav[['start']], 
             count = subnav[['count']] )
             
       R <- raster::raster(t(x), template = template)
       
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
          RR <- lapply(what, getOneVar, X$NC, subnav, crs = crs, layer = layer[1] )
          for (w in names(RR)) R <- raster::addLayer(R, 
             raster::raster(t(RR[[w]]), template = template))
          names(R) <- what
       } else {
          RR <- lapply(layer, getOneLayer, X$NC, subnav, crs = crs, what = what[1])
          R <- template
          for (r in RR) R <- raster::addLayer(R, 
             raster::raster(t(r), template = template))
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

