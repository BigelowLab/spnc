# An object that connects to NetCDF (file or OpeNDAP) using ncdf4
SPNCRefClass <- setRefClass("SPNCRefClass",
   fields = list(
      flavor = 'character',# a vector of at least 'source' and 'type'
      NC = 'ANY',          # the ncdf4 class object
      BB = 'numeric',      # the 4 element bounding box
      DIMS = 'numeric',    # the dimensions
      VARS = 'character',  # the variable names
      LON = "ANY",         # possibly a vector of lons
      LAT = "ANY",         # possibly a vector of lats
      Z = "ANY",           # possibly a vector of z-values
      TIME = 'ANY'),       # possibly a vector of times
   methods = list(
         initialize = function(nc = NULL, bb = NULL, ...){
         
         .self$field("flavor", c(source = "", type = ""))
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
      cat("  flavor:", paste(.self$flavor, collapse = ":"), "\n") 
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
#' @return logical
NULL
SPNCRefClass$methods(
   open = function(path){
      ok <- .self$is_open()
      if (ok == TRUE){
         cat("the connection is already open!\nPlease close before reopening\n")
         return(FALSE)
      }
      nc <- try(nc_open(path))
      
      if (inherits(NC, "try-error")){
         cat("unable to nc_open() the path:", path, "\n")
         .self$field("NC", NULL)
         return(FALSE)
      }
      
      .self$NC <- nc
      .self$flavor <- spnc_flavor(nc)
      .self$init(nc)
      return(r)
   }) # open


#' Init for GENERIC data
#'
#' @name SPNCRRefClass_init
#' @return logical
NULL
SPNCRefClass$methods(
   init = function(nc){
      .self$field("DIMS", ncdim_get(nc))
      .self$field("VARS", names(nc[['var']])) 
      if ("lon" %in% names(.self$DIMS)) .self$field("LON", nc[["dim"]][['lon']][['vals']])
      if ("lat" %in% names(.self$DIMS)) .self$field("LAT", nc[["dim"]][['lat']][['vals']])
      if ("zlev" %in% names(.self$DIMS)) .self$field("Z", nc[["dim"]][['zlev']][['vals']])
      if ("time" %in% names(.self$DIMS)) .self$field("TIME", nctime_get(nc))
      return(TRUE)
   })
   

#' Close the ncdf4 object if it is not closed already
#'
#' @name SPNCRefClass_close
#' @return logical
NULL
SPNCRefClass$methods(
   close = function(){
      ok <- .self$is_open()
      if (ok) nc_close(.self$NC)
      .self$NC <- NULL
      return(TRUE)
   }) # close


#' Retrieve the subset coordinates 
#' 
#' @name SPNCRefClass_subset_coords
#' @param lon numeric vector of lons (ascending order please!) to select from
#' @param lat numeric vector of lats (ditto)
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
NULL
SPNCRefClass$methods(
   subset_coords = function(bb = .self$BB, lon = c(-180,180), lat = c(-90,90)){
   
      if (is.null(bb)){
         return( list(start = c(1,1), count = c(length(lon), length(lat)),
            bb = c( range(lon), range(lat) ) ) )
      }
      
      yDescends <- (lat[1] - lat[2]) < 0
      
      ix <- findInterval(bb[0:2], lon, all.inside = TRUE)
      if (ix[2] < length(lon)) ix[2] <- ix[2] + 1
      
      if (yDescends){
         iy <- findInterval(bb[3:4], lat, all.inside = TRUE)
         if (iy[2] < length(lat)) iy[2] <- iy[2] + 1
      } else {
         iy <- findInterval(bb[3:4], rev(lat), all.inside = TRUE)
         if (iy[2] < length(lat)) iy[2] <- iy[2] + 1
      }
      
      
      list(start = c(ix[1], iy[1]), 
         count = c(ix[2]-ix[1]+1, iy[2]-iy[1]+1),
         bb = c(lon[ix], lat[iy])  )
   })


#' Get one or more variables of one or more layers (time) By default just the 
#' first layer (time) is returned.
#'
#' @name SPNCRefClass_get_raster
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NULL
SPNCRefClass$methods(
   get_raster = function(what = .self$VARS, layer = 1, 
      crs = "+proj=longlat +datum=WGS84", flip = TRUE,
      time_fmt = "D%Y%m%d"){
      
      cat("SPNCRefClass$get_raster: not implemented\n")
      return(NULL)
      # if (.self$flavor[['type']] != 'raster'){
      #    cat("SPNC$get_raster: flavor is not of type raster\n")
      #    return(NULL)
      # }
      # 
      # if (!nzchar(.self$flavor[['source']])){
      #    cat("SPNC$get_raster: flavor source is not specified\n")
      #    return(NULL)
      # }
      # 
      # if (!.self$is_open()) {
      #    cat("SPNC$get_raster: Please open the file first\n")
      #    return(NULL)
      # }
      # # convert varid index to varid character
      # if (is.numeric(what)) what <- unname(.self$VARS[what])
      # # check the varid
      # if ( !all(what %in% .self$VARS) ) {
      #    cat("one or more varid not found:", paste(varid, collapse = " "), "\n")
      #    return(NULL)
      # }
      # if (inherits(layer, "POSIXct") && (length(.self$TIME) > 1) ){
      #    ilayer <- findInterval(layer, .self$TIME, all.inside = TRUE)
      # } else {
      #    ilayer = as.numeric(layer)
      # }
      # 
      # nwhat <- length(what)
      # nlayer <- length(ilayer)
      # if ( (nlayer > 1) && (nwhat > 1)){ 
      #    cat("Either length or what must be length 1\n")
      #    return(NULL)
      # }
      # 
      # 
      # subnav <- switch(.self$flavor[['source']],
      #    "OISST" = subset_coords(.self$BB, to180(.self$LON), .self$LAT),
      #    subset_coords(.self$BB, .self$LON, .self$LAT))
      # 
      # names(what) <- what
      # 
      # R <- switch(.self$flavor[['source']],
      #    'MURSST' = MUR_get_raster(.self, what = what, layer = ilayer, crs = crs, subnav = subnav, time_fmt = time_fmt),
      #    'OISST' = OISST_get_raster(.self, what = what, layer = ilayer, crs = crs, subnav = subnav),
      #    'MODISL3SMI' = MODISL3SMI_get_raster(.self, what = what, ilayer = layer, crs = crs, subnav = subnav),
      #    NULL)
      # return(R)
   }) # get_raster



##### methods above
##### functions below

#' A function to create an SPNCRefClass or subclass reference
#' 
#' @param nc 'ncdf4' class object or path to one
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param ... futher arguments
#' @return a SPNCRefClass object or subclass or NULL
SPNC <- function(nc,
   bb = c(-180, 180, -90, 90), ...){
   
   if (!inherits(nc, "ncdf4")){
      nc <- try(nc_open(nc[1]))
      if ((inherits(nc, "try-error"))) {
         cat("SPNC: unable to open", nc[1], "\n")
         cat("SPNC: nc argument must be path or ncdf4 class object\n")
         return(NULL)
      }
   }
   flvr <- spnc_flavor(nc)
   X <- switch(tolower(flvr[['source']]),
      'oisst' = OISSTRefClass$new(nc, bb = bb, ...),
      'mursst' = MURSSTRefClass$new(nc, bb = bb, ...),
      'modisl3smi' = MODISL3SMIRefClass$new(nc, bb = bb, ...),
      'hmodisl3smi' = HMODISL3SMIRefClass$new(nc, bb = bb, ...),
      SPNCRefClass$new(nc, bb = bb, ...))
   
   invisible(X)
}  # SPNC


###### methods above
###### functions below
