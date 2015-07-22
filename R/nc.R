# nc.R

#' Retrieve a vector of dimensions
#'
#' @export
#' @param NC a ncdf4 object
#' @return a named vector of dimesnions
ncdim_get <- function(NC){
   sapply(NC[['dim']], '[[', 'len')
}

#' Retrieve a vector of variable names
#'
#' @export
#' @param NC a ncdf4 object
#' @return a named vector variable names
ncvarname_get <- function(NC){
   names(NC[['var']])
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
