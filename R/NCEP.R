# NCEP.R

#' Test if an NCDF contains NCEP/NMC data.
#' 
#' @export
#' @param x ncdf4 object or SPNCRefClass
#' @return logical
is_NCEP <- function(x){
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
      ok <- mgrepl(c("NMC", "NCEP"), atts[['title']], fixed = TRUE)   
   ok
}


#' A subclass of SPNCRefClass for NCEP/NCM
#' 
#' This NetCDF data is projected and can be ready by the \code{raster} package. Unlike
#' other members of the SPNCRefClass, this class stores the data as raster.  The get
#' \code{get_raster()} and \code{get_points()} functions are still exposed for consistency.
#'
#' @include SPNC.R
#' @export 
#' @field R a Raster* class object
NCEPRefClass <- setRefClass("NCEPRefClass",
    contains = "SPNCRefClass",
    fields = list(
      R = 'ANY'
      ),  # fields
    methods = list(
      init = function(nc){
         if (file.exists(nc$filename)) {
            .self$R <- raster::brick(nc$filename)
         } else {
            .self$R <- raster::raster(nc$filename)
         }
         callSuper(nc)
         # hours since 1800-01-01 00:00:0.0
         hours <- nctime_get(nc, as_POSIXct = FALSE)
         time <- as.POSIXct("1800-01-01 00:00:0.0", tz = 'UTC') + (hours * 3600)
         .self$field("TIME", time)
         if (inherits(.self, 'RasterBrick')) .self$close()
      })
      
   )
   
#' Compute extent given a bounding box
#' 
#' @name NCEPRefClass_get_extent
#' @param bb the bounding box needed, if missing the current one is used
#' @return a raster::Extent object
NULL
NCEPRefClass$methods(
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


#' Get the latitude locations
#'
#' @name NCEPRefClass_lat
#' @param what the desired location -  'leading', 'trailing', 'center' or 'native' (default)
#' Values 'trailing' and 'leading' are about the direction lon/lat are stored.  For example,
#' if lon/lat are each stored in ascending order then the leading edge is bottom and left. 
#' On the other hand, is lat descends (north-to-south) then the leading edge is
#' top and left.  Values 'native and 'center' are likely the same.
#' @return numeric vector or NULL
NULL
NCEPRefClass$methods(
   lat = function(what = c('native', 'leading', 'trailing', 'center')[1]){
      # native is assumed to the top/leading
      y <- raster::yFromRow(.self$R, 1:nrow(.self$R))
      s <- .self$STEP/2
      switch(tolower(what[1]),
         'center' = y + if (s[2] > 0) -s[2] else s[2],
         'trailing' = y + if (s[2] > 0) s[2] else -s[2], 
         y)
   })
   
   
#' Get the longitude locations
#'
#' @name NCEPRefClass_lon
#' @param what the desired location - 'leading', 'trailing', 'center' or 'native' (default)
#' Values 'trailing' and 'leading' are about the direction lon/lat are stored.  For example,
#' if lon/lat are each stored in ascending order then the leading edge is bottom and left. 
#' On the other hand, is lon descends (east-to-west) then the leading edge is
#' bottom and right.  Values 'native and 'center' are assumed the same - if that
#' is not a correct for a class, then a method override is required
#' @return numeric vector or NULL
NULL
NCEPRefClass$methods(
   lon = function(what = c('native', 'leading', 'trailing', 'center')[1]){
      # native is assumed to the leading
      x <- raster::xFromCol(.self$R, 1:ncol(.self$R))
      s <- .self$STEP/2
      switch(tolower(what[1]),
         'center' = x + if (s[1] > 0) s[1] else -s[1],
         'trailing' = x + if(s[1] > 0) s[1] else -s[1],
         x)
   })

#' Get the step size (resolution) as a guess - the space between the first 
#' lon and lat values.  Subclasses with unambiguous steps sizes, like L3SMI and  
#' MURSST, should override this method.
#'
#' @name NCEPRefClass_step
#' @return two element numeric vector of step size in x and y or NULL
#' These may be signed for descending values of lon and or lat
NULL
NCEPRefClass$methods(
   step = function(){
      x <- NULL
      y <- NULL
      if(!is.null(.self$NC)){
         if (all(c('lat','lon') %in% names(.self$NC$dim))){
               xr <- range(.self$NC$dim$lon$vals)
               d <- .self$get_dims()
               x <- (xr[2]-xr[1])/d[['lon']]
               yr <- range(.self$NC$dim$lat$vals)
               y <- (yr[2]-yr[1])/d[['lat']]
         }
      } else {
         x <- raster::xres(.self$R)
         y <- raster::yres(.self$R)
      }
      if (is.null(x) || is.null(y)) return(NULL)
      c(x, y)
   })   
   
# > ff <- "/Users/Shared/data/ecocast/wx/ncep-surface/air.sig995.2002.nc"
# > nc <- nc_open(ff)
# > nc
# File /Users/Shared/data/ecocast/wx/ncep-surface/air.sig995.2002.nc (NC_FORMAT_NETCDF4_CLASSIC):
# 
#      1 variables (excluding dimension variables):
#         float air[lon,lat,time]   
#             long_name: mean Daily Air temperature at sigma level 995
#             units: degK
#             precision: 2
#             least_significant_digit: 1
#             GRIB_id: 11
#             GRIB_name: TMP
#             var_desc: Air temperature
#             dataset: NCEP Reanalysis Daily Averages
#             level_desc: Surface
#             statistic: Mean
#             parent_stat: Individual Obs
#             missing_value: -9.96920996838687e+36
#             actual_range: 195.839996337891
#              actual_range: 317.700012207031
#             valid_range: 185.160003662109
#              valid_range: 331.160003662109
# 
#      3 dimensions:
#         lon  Size:144
#             units: degrees_east
#             long_name: Longitude
#             actual_range: 0
#              actual_range: 357.5
#             standard_name: longitude
#             axis: X
#         lat  Size:73
#             units: degrees_north
#             actual_range: 90
#              actual_range: -90
#             long_name: Latitude
#             standard_name: latitude
#             axis: Y
#         time  Size:365   *** is unlimited ***
#             long_name: Time
#             delta_t: 0000-00-01 00:00:00
#             avg_period: 0000-00-01 00:00:00
#             standard_name: time
#             axis: T
#             units: hours since 1800-01-01 00:00:0.0
#             actual_range: 1770696
#              actual_range: 1779432
# 
#     6 global attributes:
#         Conventions: COARDS
#         title: mean daily NMC reanalysis (2002)
#         description: Data is from NMC initialized reanalysis
# (4x/day).  These are the 0.9950 sigma level values.
#         platform: Model
#         history: created 03/08/18 by Hoop (netCDF2.3)
# Converted to chunked, deflated non-packed NetCDF4 2014/09
#         References: http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html