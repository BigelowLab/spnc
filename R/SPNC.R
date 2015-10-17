#' An object that connects to NetCDF (file or OpeNDAP) using ncdf4
#' 
#' @export
#' @field flavor a character vector of at least 'source' and 'type'
#' @field the ncdf4 class object
#' @field the 4 element bounding box
#' @field the dimensions
#' @field the variable names
#' @field possibly a vector of lons
#' @field possibly a vector of lats
#' @field possibly a vector of z-values
#' @field possibly a vector of times
SPNCRefClass <- setRefClass("SPNCRefClass",
   fields = list(
      flavor = 'list',# source=path, type=raster|point, local=logical
      NC = 'ANY',          # the ncdf4 class object
      BB = 'numeric',      # the 4 element bounding box
      DIMS = 'numeric',    # the dimensions
      VARS = 'character',  # the variable names
      LON = "ANY",         # possibly a vector of lons
      LAT = "ANY",         # possibly a vector of lats
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
         .self$field("LON", NULL)
         .self$field("LAT", NULL)
         .self$field("Z", NULL)
         .self$field("TIME", NULL)
         
         if (inherits(nc, "ncdf4")) {
            .self$NC <- nc
            .self$flavor <- spnc_flavor(nc)
            .self$init(nc)
         }
         if (!is.null(bb)) .self$BB <- bb
         callSuper(...)
      },
      # init is a bit more specific and may be overwriten by subclasses
      # in here we deal with extracting from the ncdf resource
      init = function(nc){
         .self$field("DIMS", ncdim_get(nc))
         .self$field("VARS", names(nc[['var']])) 
         if ("lon" %in% names(.self$DIMS)) .self$field("LON", nc[["dim"]][['lon']][['vals']])
         if ("lat" %in% names(.self$DIMS)) .self$field("LAT", nc[["dim"]][['lat']][['vals']])
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
      cat("  bounding box:", paste(.self$BB, collapse = " "), "\n")
      cat("  VARS:", paste(.self$VARS, collapse = " "), "\n")
      cat("  DIMS:", paste(paste(names(.self$DIMS), .self$DIMS, sep = "="), collapse = " "), "\n")
      if (!is.null(.self$LON))
         cat(sprintf("  LON: [ %f, %f]", .self$LON[1], .self$LON[.self$DIMS[['lon']]]), "\n")
      if (!is.null(.self$LAT))
         cat(sprintf("  LAT: [ %f, %f]", .self$LAT[1], .self$LAT[.self$DIMS[['lat']]]), "\n")
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
         !is.null(.self$NC)
   }) # is_open
   
#' Open the ncdf4 object
#'
#' #' @name SPNCRefClass_open
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


#' Get global attributes
#'
#' @name SPNCRefClass_get_global_atts
#' @return a named list of global attributes or NULL
NULL
SPNCRefClass$methods(
   get_global_atts = function(){
      d <- if (!is.null(.self$NC)) ncdf4::ncatt_get(.self$NC, varid = 0) else NULL
      return(d)
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
   subset_coords = function(bb = .self$BB, lon = .self$LON, lat = .self$LAT){
      subset_nav(bb=bb,lon=lon,lat=lat)  
   })



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
      cat("SPNCRefClass$get_points: not implemented\n")
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
      cat("SPNCRefClass$get_points: not implemented\n")
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
      
      cat("SPNCRefClass$get_raster: not implemented\n")
      return(NULL)
 
   }) # get_raster


##### methods above
##### functions below

#' A function to create an SPNCRefClass or subclass reference
#' 
#' @export
#' @param nc 'ncdf4' class object or path to one
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param nc_verbose logical, passed to ncdf4::nc_open
#' @param n_tries numeric, when nc is a path we try to open upt to n_tries before failing
#'    This helps accomodate occasional network/server issues.
#' @param ... futher arguments
#' @return a SPNCRefClass object or subclass or NULL
SPNC <- function(nc, bb = c(-180, 180, -90, 90), nc_verbose = FALSE, n_tries = 3, ...){
   
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
               flush.console()
               cat(sprintf("  attempt %i failed, trying again\n", i))
            } else {
               flush.console()
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
      'modisl3smi' = MODISL3SMIRefClass$new(nc, bb = bb, ...),
      'hmodisl3smi' = HMODISL3SMIRefClass$new(nc, bb = bb, ...),
      'nhsce' = NHSCERefClass$new(nc, bb = bb, ...),
      SPNCRefClass$new(nc, bb = bb, ...))
   
   invisible(X)
}  # SPNC


###### methods above
###### functions below
