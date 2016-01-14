# nc.R


#' Retrieve a named list of global attributes
#'
#' @export
#' @param NC a ncdf4 object
#' @param rm_pattern character a pattern of characters to remove from the 
#'    attribute names.  By default 'NC_GLOBAL.'.  Set to "" or NA to skip
#' @param fixed logical by default TRUE but see \code{grepl}
#' @return named vector of global attributes
ncglobal_atts <- function(NC, rm_pattern = 'NC_GLOBAL.', fixed = TRUE){
   d <- if (!is.null(NC)) ncdf4::ncatt_get(NC, varid = 0) else NULL
   if (!is.null(d)){
      if (!is.na(rm_pattern) && (nchar(rm_pattern) > 0)){
         names(d) <- gsub("NC_GLOBAL.", "", names(d), fixed = TRUE)
      }
   }
   return(d)
}

#' Retrieve a vector of dimensions
#'
#' @export
#' @param NC a ncdf4 object
#' @return a named vector of dimesnions
ncdim_get <- function(NC){
   stopifnot('dim' %in% names(NC))
   sapply(NC[['dim']], '[[', 'len')
}

#' Retrieve a vector of variable names
#'
#' @export
#' @param NC a ncdf4 object
#' @return a named vector variable names
ncvarname_get <- function(NC){
   stopifnot('dim' %in% names(NC))
   names(NC[['var']])
}


#' Retrieve a list of dimension vectors, one for each variable
#'
#' @export
#' @param NC a ncdf4 object
#' @return a named list of variable dimension vectors
ncvardim_get <- function(NC){
   stopifnot('dim' %in% names(NC))
   vn <- names(NC[['var']])
   names(vn) <- vn
   get_vardim <- function(nm, NC = NULL){
      d <- NC[['var']][[nm]]
      get_vardim_one <- function(x){
         len <- x[['len']]
         names(len) <- x[['name']]
         return(len)
      }
      
      dims <- d[['dim']]
      sapply(dims, get_vardim_one)
   }
   lapply(vn, get_vardim, NC = NC)
}

#' Retrieve a vector of timestamps for a multilayer NC object or NULL otherwise
#'
#' @export
#' @param NC a ncdf4 class object
#' @param name the name of the time-associated variable, by default 'time'
#' @param as_POSIXct logical, if TRUE then convert to POSIXct
#' @return a numeric vector of timestamps (possibly POSIXct) or NULL if none
nctime_get <- function(NC, name = 'time', as_POSIXct = TRUE){
   
   d <- ncdim_get(NC)
   if (!("time" %in% names(d)) ) return(NULL)
   
   v <- NC[["dim"]][['time']][['vals']]
   
   if (as_POSIXct){
      u <- NC[["dim"]][['time']][['units']]
      
      if (grepl(" since ", u, fixed = TRUE)){
         secsperday <- 24 * 60 * 60
         spaces <- gregexpr(" since ", u, fixed = TRUE)[[1]]
         incr <- substring(u, 1, spaces[1]-1)
         dt <- substring(u, spaces[1] + attr(spaces, 'match.length'))
         t0 <- as.POSIXct(dt, tz = "UTC")
         v <- switch(incr,
            "days" =  t0 + (v * secsperday),
            "seconds" = t0 + v,
            t0 + v)
      } else {
         cat("nctime_get: unknown time format for conversion to POSIXct\n")
      }  
   }
   invisible(v)
}
