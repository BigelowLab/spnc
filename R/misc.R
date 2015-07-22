# misc.R

#' Determine the type and source of the netcdf data
#'
#' @export
#' @param x the netcdf4 or SPNCRefClass object
#' @param return character vector of source and type.  If each is "" then no flavor determined.
spnc_flavor <- function(x){
   
   flvr <- c(source = "", type = "")
      
   lut = c(
      OISST = 'Daily-OI-V2, Final, Data (Ship, Buoy, AVHRR: NOAA19, METOP, NCEP-ice)',
      MODISL3SMI = 'MODIS Level-3 Standard Mapped Image',
      HMODISL3SMI = 'HMODISA Level-3 Standard Mapped Image',
      MURSST = 'Multi-scale Ultra-high Resolution (MUR) SST analysis'
      )

   if (inherits(x, "SPNCRefClass")){
      atts <- try(ncdf4::ncatt_get(x$NC, varid = 0))
   } else {
      if (!inherits(x, "ncdf4")) {
         cat("spnc_flavor: input must be SPNCRefClass or ncdf4 class\n")
         return(flvr)
      }
      atts <- try(ncdf4::ncatt_get(x, varid = 0))
   }
    
   if (inherits(atts, "try-error")){
     cat("spnc_flavor: error accessing global variables\n")
     return(flvr)
   }
   
   natts <- names(atts)

   # try by title
   if ("title" %in% natts){
      if (grepl(lut[['OISST']], atts[['title']], fixed = TRUE)){
         flvr <- c(source = "OISST", type = "raster")
      } else if (grepl(lut[['MODISL3SMI']], atts[['title']], fixed = TRUE)){
         flvr <- c(source = "MODISL3SMI", type = "raster")
      } else if (grepl(lut[['HMODISL3SMI']], atts[['title']], fixed = TRUE)){
         flvr <- c(source = "HMODISL3SMI", type = "raster")
      } else if (grepl(lut[['MURSST']], atts[['title']], fixed = TRUE)){
         flvr <- c(source = "MURSST", type = "raster")
      }
      if ( nzchar(flvr[['source']]) ) return(flvr)
   }
   
   return(flvr)
}


#' Convert [0,360] longitudes to [-180, 180]
#' 
#' @export
#' @param x numeric vector, no check is done for being withing 0,360 range
#' @return numeric vector
to180 <- function(x) {ix <- x > 180 ; x[ix] <- x[ix]-360; x}

#' Convert [-180,180] longitudes to [0,360]
#' 
#' @export
#' @param x numeric vector, no check is done for being withing 0,360 range
#' @return numeric vector
to360 <- function(x) {ix <- x < 0 ; x[ix] <- x[ix]+ 360; x}

#' Craft subset indices into a ncdf array
#' @seealso \url{http://r.789695.n4.nabble.com/Lat-Lon-NetCDF-subset-td3394326.html}
#' @param bb numeric, four element bounding box [left, right, bottom, top]
#' @param lon numeric vector of lons (ascending order please!) to select from
#' @param lat numeric vector of lats (ditto)
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
subset_coords <- function(bb = NULL, lon = c(-180,180), lat = c(-90,90)){
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
}