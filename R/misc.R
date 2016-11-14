# misc.R

#' Determine the type and source of the netcdf data
#'
#' @export
#' @param x the netcdf4 or SPNCRefClass object
#' @return return list of 
#'  \itemize{
#'    \item source character of 'MURSST', 'OISST', etc. "" means unknown
#'    \item type character of 'raster', 'points', etc. "" means unknown
#'    \item local logical (or NA) flag to indicate local or http source. NA mean unknown
#'  }
spnc_flavor <- function(x){
   
   flvr <- list(source = "", type = "raster", local = NA)
   filename <- NULL 

   if (inherits(x, "SPNCRefClass")){
      atts <- try(x$get_global_atts())
      #atts <- try(ncdf4::ncatt_get(x$NC, varid = 0))
      filename <- x$NC[['filename']]
   } else {
      if (!inherits(x, "ncdf4")) {
         cat("spnc_flavor: input must be SPNCRefClass or ncdf4 class\n")
         return(flvr)
      }
      atts <- try(ncglobal_atts(x))
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
      
   if (is_L3SMI(x)){
      flvr[['source']] <- "L3SMI"
   } else if (is_OISST(x)){
      flvr[['source']] <- "OISST"
   } else if (is_MURSST(x)){
      flvr[['source']] <- "MURSST"
   } else if (is_CPCUGBDP(x)){
      flvr[['source']] <- "CPCUGBDP"
   } else if (is_NARR(x)){
      flvr[['source']] <- "NARR"
   } else if (is_NAMANL(x)){
      flvr[['source']] <- "NAMANL"
   } else if (is_NCEP(x)){
      flvr[['source']] <- "NCEP"
   }
   
   
   return(flvr)
}


#' Perform grepl on multiple patterns; it's like OR-ing successive grepl statements.
#' 
#' @export
#' @param pattern character vector of patterns
#' @param x the character vector to search
#' @param ... further arguments for \code{grepl}
#' @return logical vector of matches where at least one element of \code{pattern} 
#'    is found in \code{x}
mgrepl <- function(pattern, x, ...){
   ix <- do.call(rbind, lapply(pattern, 
      function(p, x = "", ...){ grepl(p, x, ...) },
      x = x, ...))  
   ok <- apply(ix, 2, any)
}



#' Convert bounding box [0,360] longitudes to [-180, 180]
#' 
#' Bounding boxes are 4 element vectors of [left, right, bottom, top]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to180BB <- function(x) {x[1:2] <- to180(x[1:2]) ; x}

#' Convert [-180,180] bounding box longitudes to [0,360]
#' 
#' Bounding boxes are 4 element vectors of [left, right, bottom, top]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to360BB <- function(x) {x[1:2] <- to360(x[1:2]) ; x}

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


#' A wrapper around base::findInterval() that allows decreasing values in the 
#' value of the vector within which we wish to place values of x.
#'
#' When \code{vec} is in ascending order we use \code{base::findInterval()}, but
#' when \code{vec} is in descending order we implement an adaptation of the 
#' \code{locate()} function from Numerical Recipes for C \url{http://apps.nrbook.com/c/index.html}
#' 
#' @export
#' @param x numeric values we wish to located within \code{vec}
#' @param vec numeric vector of sorted values (ascending or descending order) 
#'    within which we wish to find the placement of \code{x}
#' @param rightmost.closed see \link{findInterval}
#' @param all.inside see \link{findInterval}
#' @return see \link{findInterval}
find_interval <- function(x, vec, rightmost.closed = FALSE, all.inside = FALSE){

   # locate one value of x within v
   # @param v ordered numeric vector
   # @param x one numeric to locate within v
   # @return index into v
   locate_one <- function(v, x){
      n <- length(v)
      ascnd <- v[n] >= v[1]
      iL <- 1
      iU <- n
      while((iU-iL) > 1){
         iM <- bitwShiftR((iU+iL),1)
         if (ascnd){
            if (x >= v[iM]){
               iL <- iM
            } else {
               iU <- iM
            }
         } else {
            if (x <= v[iM]){
               iL <- iM
            } else {
               iU <- iM
            }
         } 
      }
      
      if (ascnd) {
			if ( x < v[1]) {
				index <- 0
			} else if (x >= v[n]) {
				index <- n
			} else {
				index <- iL
			}
		} else {
			if ( x > v[1]) {
				index <- 0
			} else if (x <= vec[n]) {
				index <- n
			} else {
				index <- iL
			}
		}
		return(index)
	}  # locate_one

   ascending <- vec[length(vec)] >= vec[1]
   
   if (!ascending) {
      # here we do our own implementation (with a performance hit)
      j <- sapply(x, function(x, v=NULL) locate_one(v,x), v = vec)   
      nv <- length(vec)
      if (all.inside){
         j[j < 1] <- 1
         j[j >= nv] <- nv - 1
      } 
      if (rightmost.closed){
         j[x <= vec[nv]] <- nv - 1
      }
   } else {
      # this is plain vanilla stuff we pass to findInterval
      j <- base::findInterval(x, vec, 
         rightmost.closed = rightmost.closed, all.inside = all.inside)
   }
   j
}  # find_interval

#' Coerce values of x two be within the range of the bounds b.
#' 
#' @export
#' @param x numeric vector
#' @param b two element numeric vector of [min,max] defining bounds
#' @return a numeric vector same length as x with all values within bounds b
coerce_within <- function(x, b = c(0.0, 1.0)){
   x[x < b[1]] <- b[1]
   x[x > b[2]] <- b[2]
   x
} 


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


#' Craft subset indices into a ncdf array
#'
#' @export
#' @seealso \url{http://r.789695.n4.nabble.com/Lat-Lon-NetCDF-subset-td3394326.html}
#' @param bb numeric, four element bounding box [left, right, bottom, top]
#' @param lon numeric vector of lons to select from - these come from the LON attribute in nc
#' @param lat numeric vector of lats to select from - these come from the LAT attribute in nc
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
subset_bbox <- function(bb = NULL, lon = c(-180,180), lat = c(-90,90)){
   if (is.null(bb)){
      return( list(start = c(1,1), count = c(length(lon), length(lat)),
         bb = c( range(lon), range(lat) ) ) )
   }
   
   # when data is stored north to south
   #yDescends <- (lat[2] - lat[1]) < 0
   
  
   
   # # when we pop out of this structure we need iy to be the indices of [start,stop]
   # # into LAT (thus the indices into the NC data   
   # if (yDescends){
   #    iy <- which(findInterval(lat, bb[3:4], rightmost.closed = TRUE) == 1)
   #    iy <- range(iy)
   # } else {
   #    iy <- findInterval(bb[3:4], lat, all.inside = TRUE)
   # }
   
   ix <- find_interval(bb[0:2], lon)
   iy <- find_interval(bb[3:4], lat)
   # get these in order
   ix <- range(ix)
   iy <- range(iy)
   # make sure they fit within the dims of the data
   ix <- coerce_within(ix, c(1, length(lon)))
   iy <- coerce_within(iy, c(1, length(lat)))
   
   newbb <- c(range(lon[ix]), range(lat[iy]))
   list(start = c(ix[1], iy[1]), 
      count = c(ix[2]-ix[1]+1, iy[2]-iy[1]+1),
      bb = newbb  )
}


#' Converts Lon/Lat to start and count suitable for extracting points
#'
#' @export
#' @param x either a vector of longitude or a 2d array [lon,lat] or a data.frame
#' with lon,lat columns
#' @param y numeric vector, if x is a vector, otherwise NULL
#' @param lon numeric vector of lons to select from - these come from the LON attribute in nc
#' @param lat numeric vector of lats to select from - these come from the LAT attribute in nc
#' @return a list of start/count elements, one per xy input.  For inputs that
#'    fall beyond the bounds of the data then NA is returned.
subset_points <- function(x, y = NULL, lon = c(-180,180), lat = c(-90, 90)){
   
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
   }
   
            
   # when data is stored north to south
   yDescends <- (lat[2] - lat[1]) < 0
   
   dx <- lon[2] - lon[1]
   dy <- abs(lat[1]-lat[2])
   
   ix <- find_interval(x, lon)
   ix[ix < 1] <- 1

   iy <- find_interval(y, lat)
   iy[iy < 1] <- 1 
   
   iz <- mapply(function(x, y){
         list(start = c(x,y), count = c(1,1))
      },
      ix,iy, SIMPLIFY = FALSE)
   
   ix <- (x < lon[1]) | (x > (lon[length(lon)] + dx))
   iy <- if (yDescends){
         (y > lat[1]) | (y < (lat[length(lat)] - dy))
      } else {
         (y < lat[1]) | (y > (lat[length(lat)] + dy))
      } 
   iz[ix | iy] <- NA
   iz
} # subset_points


#' Convert a bbox to sp::polygon, sp::Polygons, or sp::SpatialPolygons object
#' 
#' @export
#' @param bb numeric, a 4-element bbox vector [left, right, bottom, top]
#' @param projstring character sutiable to pass to \code{sp::CRS}, by default "+proj=longlat +datum=WGS84"
#' @param id character, the polygon ID, by default 'bbox'
#' @param output_class character, either "SpatialPolygons", "Polygons" or "Polygon"
#' @return sp::SpatialPolygons object or NULL
bbox_to_polygon <- function(bb, 
   projstring = "+proj=longlat +datum=WGS84",
   id = 'bbox',
   output_class = c("SpatialPolygons", "Polygons",  "Polygon")[1]){

   if(missing(bb)) stop("bounding box, bb, is required")
   stopifnot(length(bb) == 4)
   
   # clockwise to form an 'island'
   bb <- cbind(
      x = c(bb[1], bb[1], bb[2], bb[2], bb[1]),
      x = c(bb[3], bb[4], bb[4], bb[3], bb[3]) )
   # project
   bbp <- rgdal::project(bb, projstring)
   # make a Polygon
   Poly <- sp::Polygon(bbp)
   if (output_class == 'Polygon') return(Poly)
   # make into Polygons
   Polys <- sp::Polygons(list(Poly), id)
   if (output_class == 'Polygons') return(Polys)
   sp::SpatialPolygons(list(Polys), proj4string = CRS(projstring))
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
