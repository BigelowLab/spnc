# misc.R

#' Determine the type and source of the netcdf data
#'
#' @export
#' @param x the netcdf4 or SPNCRefClass object
#' @param return list of 
#'  \itemize{
#'    \item source character of 'MURSST', 'OISST', etc. "" means unknown
#'    \item type character of 'raster', 'points', etc. "" means unknown
#'    \item local logical (or NA) flag to indicate local or http source. NA mean unknown
#'  }
spnc_flavor <- function(x){
   
   flvr <- list(source = "", type = "raster", local = NA)
   filename <- NULL 
   lut = c(
      OISST = 'Daily-OI-V2, Final, Data (Ship, Buoy, AVHRR: NOAA19, METOP, NCEP-ice)',
      MODISL3SMI = 'MODIS Level-3 Standard Mapped Image',
      HMODISL3SMI = 'HMODISA Level-3 Standard Mapped Image',
      MURSST = 'Multi-scale Ultra-high Resolution (MUR) SST analysis',
      VIIRS = 'VIIRS Level-3 Standard Mapped Image',
      L3SMI = 'Level-3 Standard Mapped Image',
      NHSCE = "Climate Data Record (CDR) of Northern Hemisphere (NH) Snow Cover Extent (SCE) (CDR Name: Snow_Cover_Extent_NH_IMS_Robinson)"
      )

   if (inherits(x, "SPNCRefClass")){
      atts <- try(ncdf4::ncatt_get(x$NC, varid = 0))
      filename <- x$NC[['filename']]
   } else {
      if (!inherits(x, "ncdf4")) {
         cat("spnc_flavor: input must be SPNCRefClass or ncdf4 class\n")
         return(flvr)
      }
      atts <- try(ncdf4::ncatt_get(x, varid = 0))
      filename <- x[['filename']]
   }
    
   if (inherits(atts, "try-error")){
     cat("spnc_flavor: error accessing global variables\n")
     return(flvr)
   }
   
   # assign the local flag 
   if (!is.null(filename)) flvr[['local']] <- !grepl("^http", filename[1])
   
   natts <- names(atts) <- tolower(names(atts))


   # try by title
   # note that the default type is 'raster' but others, as needed,
   # must be assigned with flvr[['type']] <- 'points' or whatever
   if ("title" %in% natts){
      
      #if ( nzchar(flvr[['source']]) ) return(flvr)
      if (grepl(lut[['HMODISL3SMI']], atts[['title']], fixed = TRUE)){
         flvr[['source']] <- "HMODISL3SMI"
      } else if (grepl(lut[['L3SMI']], atts[['title']], fixed = TRUE)){
         flvr[['source']] <- "L3SMI"  # OBPG
      } else if (grepl(lut[['OISST']], atts[['title']], fixed = TRUE)){
         flvr[['source']] <- "OISST"
      } else if (grepl(lut[['MODISL3SMI']], atts[['title']], fixed = TRUE)){
         flvr[['source']] <- "MODISL3SMI"
      } else if (grepl(lut[['MURSST']], atts[['title']], fixed = TRUE)){
         flvr[['source']] <- "MURSST"
      } else if (grepl(lut[['VIIRS']], atts[['title']], fixed = TRUE)){
         flvr[['source']] <- "VIIRS"
      } else if (grepl(lut[['NHSCE']], atts[['title']], fixed = TRUE)){
         flvr[['source']] <- "NHSCE"
      }
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
#'
#' @export
#' @seealso \url{http://r.789695.n4.nabble.com/Lat-Lon-NetCDF-subset-td3394326.html}
#' @param bb numeric, four element bounding box [left, right, bottom, top]
#' @param lon numeric vector of lons to select from - these come from the LON attribute in nc
#' @param lat numeric vector of lats to select from - these come from the LAT attribute in nc
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
subset_nav <- function(bb = NULL, lon = c(-180,180), lat = c(-90,90)){
   if (is.null(bb)){
      return( list(start = c(1,1), count = c(length(lon), length(lat)),
         bb = c( range(lon), range(lat) ) ) )
   }
   
   # when data is stored north to south
   yDescends <- (lat[2] - lat[1]) < 0
   
   ix <- findInterval(bb[0:2], lon, all.inside = TRUE)

   # when we pop out of this structure we need iy to be the indices of [start,stop]
   # into LAT (thus the indices into the NC data   
   if (yDescends){
      iy <- which(findInterval(lat, bb[3:4], rightmost.closed = TRUE) == 1)
      iy <- range(iy)
   } else {
      iy <- findInterval(bb[3:4], lat, all.inside = TRUE)
   }
   
   newbb <- c(lon[ix], range(lat[iy]))
   
   list(start = c(ix[1], iy[1]), 
      count = c(ix[2]-ix[1]+1, iy[2]-iy[1]+1),
      bb = newbb  )
}


#' Extract a matrix of data from a SPNCRefClass of raster flavor
#' 
#' @export
#' @param x the SPNCRefClass or subclass object
#' @param var character, the name of the variable
#' @param layer index of the layer (by default the first)
#' @param ... further arguments for SPNC$get_raster() method
#' @return a numeric matrix or NULL
get_matrix <- function(x, var, layer = 1,...){
   
   if (missing(x) || missing(var)) {
      cat("get_matrix: x and var are required\n")
      return(NULL)
   }
   
   if (x$flavor[['type']] != 'raster'){
      cat("get_matrix: SPNC flavor must be of type raster\n")
      return(NULL)
   }
   
   r <- try(x$get_raster(var, layer = layer,...))
   if (inherits(r, "try-error")){
      cat("get_matrix: error in get_raster method\n")
      return(NULL)
   }
   
   m <- as.matrix(r[[1]])
   rownames(m) <- raster::yFromRow(r)
   colnames(m) <- raster::xFromCol(r)
   invisible(m)  
}
