#' An object that connects to NetCDF (file or OpeNDAP) using ncdf4
#' 
#' @export
#' @field flavor a character vector of at least 'source' and 'type'
#' @field NC the ncdf4 class object
#' @field BB the 4 element bounding box
#' @field DIM the dimensions
#' @field VAR the variable names
#' @field STEP 2 element vector of x and y step size (or NULL)
#' @field Z possibly a vector of z-values
#' @field TIME possibly a vector of times
SPNCRefClass <- setRefClass("SPNCRefClass",
   fields = list(
      flavor = 'list',# source=path, type=raster|point, local=logical
      NC = 'ANY',          # the ncdf4 class object
      BB = 'numeric',      # the 4 element bounding box
      DIMS = 'numeric',    # the dimensions
      VARS = 'character',  # the variable names
      STEP = "ANY",
      Z = "ANY",           # possibly a vector of z-values
      TIME = 'ANY'),       # possibly a vector of times
   methods = list(
      # initialize if a general kickoff
      initialize = function(nc = NULL, bb = NULL, ...){
         .self$field("flavor", list(source = "", type = "", local = NA))
         .self$field("NC", NULL)
         .self$field("BB",  c(-180, 180, -90, 90))
         .self$field("DIMS", numeric())
         .self$field("VARS", character())
         .self$field("STEP", NULL)
         .self$field("Z", NULL)
         .self$field("TIME", NULL)
         
         if (inherits(nc, "ncdf4")) {
            .self$NC <- nc
            .self$flavor <- spnc_flavor(nc)
            .self$init(nc)
         }
         if (!is.null(bb)) .self$BB <- bb
         
      },
      finalize = function(){
         ok <- .self$close()
      }, 
      # init is a bit more specific and may be overwriten by subclasses
      # in here we deal with extracting from the ncdf resource
      init = function(nc){
         .self$field("DIMS", ncdim_get(nc))
         .self$field("VARS", names(nc[['var']]))
         .self$field("STEP", .self$step()) 
         e <- .self$get_extent()
         #.self$BB <- raster::as.vector(e) # c(e@xmin, e@xmax, e@ymin, e@ymax)
         .self$BB <- c(e@xmin, e@xmax, e@ymin, e@ymax)
         if ("zlev" %in% names(.self$DIMS)) .self$field("Z", nc[["dim"]][['zlev']][['vals']])
         if ("time" %in% names(.self$DIMS)) .self$field("TIME", nctime_get(nc))
         return(TRUE)
      })
   )
   
#' Show the contents in a tidy fashion
#' 
#' @name SPNCRefClass_show
NULL
SPNCRefClass$methods(
   show = function(){
      state <- if( .self$is_open() ) 'opened' else 'closed'
      cat("Reference Class:", classLabel(class(.self)), "\n")
      
      cat("  flavor:", paste(names(.self$flavor), .self$flavor, sep = "=", collapse = " "), "\n") 
      cat("  state:", state, "\n")
      cat("  bounding box:", sprintf("%0.4f %0.4f %0.4f %0.4f", 
         .self$BB[1],.self$BB[2], .self$BB[3], .self$BB[4]),  "\n")
      cat("  VARS:", paste(.self$VARS, collapse = " "), "\n")
      cat("  DIMS:", paste(paste(names(.self$DIMS), .self$DIMS, sep = "="), collapse = " "), "\n")
      ext <- .self$get_extent()
      if (!is.null(ext)){
         cat(sprintf("  LON: [ %0.2f, %0.2f]", ext[1], ext[2]), "\n")
         cat(sprintf("  LAT: [ %0.2f, %0.2f]", ext[3], ext[4]), "\n")
      }
      if (!is.null(.self$Z))
         cat(sprintf("  Z: [ %f, %f]", .self$Z[1], .self$Z[.self$DIMS[['zlev']]]), "\n")
       if (!is.null(.self$TIME))
         cat(sprintf("  TIME: [ %s, %s]", .self$TIME[1], .self$TIME[length(.self$TIME)]), "\n")   
      invisible(NULL)
   }) # show     


#' Test if this is a local file 
#'
#' @name SPNCRefClass_is_local
#' @return logical
NULL
SPNCRefClass$methods(
   is_local = function(){
   
      # flavor[['local']] could be NA, TRUE or FALSE.  If NA then we return FALSE
      if (is.na(.self$flavor[['local']])) {
         return(FALSE)
      } else {
         return(.self$flavor[['local']])
      }
   })

#' Test if the ncdf4 path is open
#'
#' @name SPNCRefClass_is_open
#' @return logical
NULL
SPNCRefClass$methods(
      is_open = function(){
         !is.null(.self$NC) && inherits(.self$NC, "ncdf4")
   }) # is_open
   
#' Open the ncdf4 object
#'
#' @name SPNCRefClass_open
#' @param path the path to the connection, if not present then try the 
#'     object's path
#' @param ... furtehr arguments for ncfd4::nc_open
#' @return logical
NULL
SPNCRefClass$methods(
   open = function(path, ...){
      ok <- .self$is_open()
      if (ok == TRUE){
         cat("the connection is already open!\nPlease close before reopening\n")
         return(FALSE)
      }
      nc <- try(ncdf4::nc_open(path,...))
      
      if (inherits(nc, "try-error")){
         cat("unable to nc_open() the path:", path, "\n")
         .self$field("NC", NULL)
         return(FALSE)
      }
      
      .self$NC <- nc
      .self$flavor <- spnc_flavor(nc)
      .self$init(nc)
      return(r)
   }) # open



#' Close the ncdf4 object if it is not closed already
#'
#' @name SPNCRefClass_close
#' @return logical
NULL
SPNCRefClass$methods(
   close = function(){
      ok <- .self$is_open()
      if (ok) ncdf4::nc_close(.self$NC)
      .self$NC <- NULL
      return(TRUE)
   }) # close

#' Get the longitude locations
#'
#' @name SPNCRefClass_lon
#' @param what the desired location - 'leading', 'trailing', 'center' or 'native' (default)
#' Values 'trailing' and 'leading' are about the direction lon/lat are stored.  For example,
#' if lon/lat are each stored in ascending order then the leading edge is bottom and left. 
#' On the other hand, is lon descends (east-to-west) then the leading edge is
#' bottom and right.  Values 'native and 'center' are assumed the same - if that
#' is not a correct for a class, then a method override is required
#' @return numeric vector or NULL
NULL
SPNCRefClass$methods(
   lon = function(what = c('native', 'leading', 'trailing', 'center')[1]){
      # native is assumed to the center
      x <- if ('lon' %in% names(.self$NC$dim)) .self$NC$dim$lon$vals else  NULL
      if (is.null(x)) return(NULL)
      s <- .self$STEP/2
      switch(tolower(what[1]),
         'leading' = x + if (s[1] > 0) -s[1] else s[1],
         'trailing' = x + if(s[1] > 0) s[1] else -s[1],
         x)
   })

#' Get the latitude locations
#'
#' @name SPNCRefClass_lat
#' @param what the desired location -  'leading', 'trailing', 'center' or 'native' (default)
#' Values 'trailing' and 'leading' are about the direction lon/lat are stored.  For example,
#' if lon/lat are each stored in ascending order then the leading edge is bottom and left. 
#' On the other hand, is lat descends (north-to-south) then the leading edge is
#' top and left.  Values 'native and 'center' are likely the same.
#' @return numeric vector or NULL
NULL
SPNCRefClass$methods(
   lat = function(what = c('native', 'leading', 'trailing', 'center')[1]){
      # native is assumed to the center
      y <- if ('lat' %in% names(.self$NC$dim)) .self$NC$dim$lat$vals else NULL
      if (is.null(y)) return(NULL)
      s <- .self$STEP/2
      switch(tolower(what[1]),
         'leading' = y + if (s[2] > 0) -s[2] else s[2],
         'trailing' = y + if (s[2] > 0) s[2] else -s[2], 
         y)
   })

#' Get the step size (resolution) as a guess - the space between the first 
#' lon and lat values.  Subclasses with unambiguous steps sizes, like L3SMI and  
#' MURSST, should override this method.
#'
#' @name SPNCRefClass_step
#' @return two element numeric vector of step size in x and y or NULL
#' These may be signed for descending values of lon and or lat
NULL
SPNCRefClass$methods(
   step = function(){
      x <- NULL
      y <- NULL
      if (all(c('lat','lon') %in% names(.self$NC$dim))){
            xr <- range(.self$NC$dim$lon$vals)
            d <- .self$get_dims()
            x <- (xr[2]-xr[1])/d[['lon']]
            yr <- range(.self$NC$dim$lat$vals)
            y <- (yr[2]-yr[1])/d[['lat']]
      }
      if (is.null(x) || is.null(y)) return(NULL)
      c(x, y)
   })

#' Compute extent given a bounding box
#' 
#' @name SPNCRefClass_get_extent
#' @param bb the bounding box needed, if missing the current one is used
#' @return a raster::Extent object
NULL
SPNCRefClass$methods(
   get_extent = function(bb){
      if (missing(bb)){
         nx <- .self$DIMS[['lon']]
         ny <- .self$DIMS[['lat']]
         s <- .self$step()
         bb <- c(range(.self$lon()[c(1,nx)]), range(.self$lat()[c(1,ny)]))
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


#' Compute time indices (which must be contiguous) from POSIXt, Date or index.
#' An warning is thrown if the dates requested are not contiguous.  If the object
#' does not have TIME element then 1 is returned.
#'
#' @name SPNCRefClass_time_index
#'
#' @param when numeric, POSIXct or Date times
#' @param no_zero logical, any times requested before the first time are mapped
#'    to the first time
#' @return one or more indices 
NULL
SPNCRefClass$methods(
   time_index = function(when = 1, no_zero = TRUE){
      if (is.null(.self$TIME)) return(1)
      if (length(.self$TIME) <= 1) return(1)
      if (inherits(when, 'POSIXt') || inherits(when, 'Date')) {    
         ix <- find_interval(when, .self$TIME)
      } else {
         ix <- find_interval(when, seq_len(.self$NC[["dim"]][['time']][['len']]))
      }
      if (no_zero) ix[ix <= 0] <- 1
      if (!all(diff(ix) == 1)) warning("time indices are not contiguous")
      ix
   })


#' Get global attributes
#'
#' @name SPNCRefClass_get_global_atts
#' @param ... further arguments for \code{ncglobal_atts} 
#' @return a named list of global attributes or NULL
NULL
SPNCRefClass$methods(
   get_global_atts = function(...){
      return(ncglobal_atts(.self$NC, ...))
   })
   
#' Retrieve the variable names
#'
#' @name SPNCRefClass_get_varnames
#' @return character vector or NULL
NULL
SPNCRefClass$methods(
   get_varnames = function(){
      d <- if (!is.null(.self$NC)) ncvarname_get(.self$NC) else NULL
      return(d)
   })


#' Retrieve a named vector of dimensions
#'
#' @name SPNCRefClass_get_dims
#' @return a named numeric vector of dimensions or NULL
NULL
SPNCRefClass$methods(
   get_dims = function(){
      d <- if (!is.null(.self$NC)) ncdim_get(.self$NC) else NULL
      return(d)
   })


#' Retrieve the subset coordinates 
#' 
#' @name SPNCRefClass_subset_coords
#' @param lon numeric vector of lons (ascending order please!) to select from
#' @param lat numeric vector of lats (ditto)
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
NULL
SPNCRefClass$methods(
   subset_coords = function(bb = .self$BB, lon = .self$lon("leading"), 
      lat = .self$lat("leading")){
      subset_nav(bb=bb,lon=lon,lat=lat)  
   })


#' Order a subset coordinate system into the order matching that stored in
#' the ncdf file.  For ncdf4 the order should be X,Y,Z,T (which might be 
#' lon, lat, zlev, time, but not necessarily those names).  Honestly, this
#' probably should not even be a method - but use at your own risk/pain.
#'
#' @name SPNCRefClass_order_subset
#' @param s a subset list comprised of at least start and count
#' @return the same list with elements of start and count properly ordered
NULL
SPNCRefClass$methods(
   order_subset = function(s){
     stopifnot(all(c("start", "count") %in% names(s)))
     vnames <- names(.self$NC$dim)
     s[['start']] <- s[['start']][vnames]
     s[['count']] <- s[['count']][vnames]
     return(s)
   })

#' Craft subset indices into a ncdf array
#'
#' @name SPNCRefClass_subset_bbox
#'
#' @seealso \url{https://stat.ethz.ch/pipermail/r-help/2011-March/272641.html}
#' @param bb numeric, four element bounding box [left, right, bottom, top]
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
NULL
SPNCRefClass$methods(
   subset_bbox = function(bb = NULL){
      llon <- .self$lon("leading")
      llat <- .self$lat("leading")
      if (is.null(bb)){
         return(
            list(start = c(lon=1,lat=1), count = c(lon=length(llon), lat=length(llat)),
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
      
      list(start = c(lon=ix[1], lat=iy[1]), 
         count = c(lon=ix[2]-ix[1]+1, lat=iy[2]-iy[1]+1),
         bb = c(range(xx), range(yy)) )
   }) # subset_bbox
   

#' Compute indicies (start and count) for individual [lon,lat] points
#' 
#' @name SPNCRefClass_subset_points
#' @param x either a vector of longitude or a 2d matrix [lon,lat] or a data.frame
#'    with lon,lat columns
#' @param y numeric vector, if x is a vector, otherwise NULL
#' @return a list of start/count elements, one per xy input.  For inputs that
#'    fall beyond the bounds of the data then NA is returned.
NULL
SPNCRefClass$methods(
   subset_points = function(x, y = NULL){
   
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
      iz <- mapply(function(x, y){list(start = c(lon=x,lat=y), count = c(lon=1,lat=1))},
         ix,iy, SIMPLIFY = FALSE)
      # now we have to flag any beyond extent as NA
      # we could use the indices for the leading edge (where ix or iy < 1)
      # but we would still need to dub around with the trailing edge
      # so it is just as easy to use the input x and y instead of ix and iy
      ixna <- (x < e[1]) | (x > e[2])
      iyna <- (y < e[3]) | (y > e[4])
      iz[ixna | iyna] <- NA
      return(iz)
   }) #subset_points


#' Get point data
#'
#' @name SPNCRefClass_get_points
#' @param what character one or more variable names or variable indices
#' @param crs character, the coordiante reference system to apply
#' @return a \code{SpatialPointsDataFrame} or NULL
NULL
SPNCRefClass$methods(
   get_points = function(what = .self$VARS, 
      layer = 1, 
      crs = "+proj=longlat +datum=WGS84"){
       cat(paste0(classLabel(class(.self)),"$get_points: not implemented\n"))
      return(NULL)  
   })


#' Get path (line) data
#'
#' @name SPNCRefClass_get_path
#' @param what character one or more variable names or variable indices
#' @param crs character, the coordiante reference system to apply
#' @return a \code{SpatialLinesDataFrame} or NULL
NULL
SPNCRefClass$methods(
   get_path = function(x, y, time, what = .self$VARS, 
      crs = "+proj=longlat +datum=WGS84"){
       cat(paste0(classLabel(class(.self)),"$get_path: not implemented\n"))
      return(NULL)  
   })



   
#' Get one or more variables of one or more layers (time) By default just the 
#' first layer (time) is returned.
#'
#' @name SPNCRefClass_get_raster
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NULL
SPNCRefClass$methods(
   get_raster = function(what = .self$VARS, layer = 1, bb = .self$BB,
      crs = "+proj=longlat +datum=WGS84", flip = TRUE,
      time_fmt = "D%Y%m%d"){
      cat(paste0(classLabel(class(.self)),"$get_raster: not implemented\n"))
      return(NULL)
 
   }) # get_raster


##### methods above
##### functions below

#' A function to create an SPNCRefClass or subclass reference
#' 
#' @export
#' @param nc 'ncdf4' class object or path to one
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to NULL
#' @param nc_verbose logical, passed to ncdf4::nc_open
#' @param n_tries numeric, when nc is a path we try to open upt to n_tries before failing
#'    This helps accomodate occasional network/server issues.
#' @param ... futher arguments
#' @return a SPNCRefClass object or subclass or NULL
SPNC <- function(nc, bb = NULL, nc_verbose = FALSE, n_tries = 3, ...){
   
   # if (grepl("windows", .Platform$OS.type, fixed = TRUE)){
   #    cat("Windows platform detected\n")
   #    cat("Be advised that OpeNDAP access with ncdf4 are problematic.\n")
   #    cat("Local NetCDF file access is trouble free.\n")
   # }
   
   if (!inherits(nc, "ncdf4")){
      path <- nc[1]
      i <- 1
      nc <- NULL
      while(i <= n_tries){
         nc <- try(ncdf4::nc_open(path, verbose = nc_verbose))
         if (inherits(nc, 'try-error')){
            if (i < n_tries) {
               utils::flush.console()
               cat(sprintf("  attempt %i failed, trying again\n", i))
            } else {
               utils::flush.console()
               cat("  ***\nexhausted permitted tries, returning NULL\n  ***\n")
            }
            nc <- NULL
            i <- i + 1
         } else {
            if (i > 1) cat(sprintf("  whew!  attempt %i successful\n", i))
            break
         }
      }
      if (is.null(nc) || inherits(nc, 'try-error')) return(NULL)
      
   }
   flvr <- spnc_flavor(nc)
   X <- switch(tolower(flvr[['source']]),
      'l3smi' = L3SMIRefClass$new(nc, bb = bb, ...),
      'oisst' = OISSTRefClass$new(nc, bb = bb, ...),
      'mursst' = MURSSTRefClass$new(nc, bb = bb, ...),
      'cpcugbdp' = CPCUGBDPRefClass$new(nc, bb = bb, ...),
      'narr' = NARRRefClass$new(nc, bb = bb, ...),
      'namanl' = NAMANLRefClass$new(nc, bb = bb, ...),
      'ncep' = NCEPRefClass$new(nc, bb = bb, ...),
      'bsw' = BlendedSeaWindsRefClass$new(nc, bb, ...),
      SPNCRefClass$new(nc, bb = bb, ...))
   
   invisible(X)
}  # SPNC


###### methods above
###### functions below
