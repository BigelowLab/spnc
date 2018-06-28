#' Test if an NCDF contains NAM-ANL data.
#'
#' @export
#' @param x ncdf4 object or NAMANLRefClass
#' @return logical
is_NAMANL <- function(x){
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
      ok <- grepl("MESO NAM Model Forecast", atts[['title']], fixed = TRUE)
    ok
}


#' A subclass of SPNCRefClass for NAM-ANL
#'
#'
#' @include SPNC.R
#' @export
#' @field lambert a Raster* class object
#' @field mercator  a Raster* class object
NAMANLRefClass <- setRefClass("NAMANLRefClass",
    contains = "SPNCRefClass",
    fields = list(
        lccR = 'ANY',
        longlatR = 'ANY'
        ),  # fields
    methods = list(
        init = function(nc){

            lccProj <- "+proj=lcc +lat_1=25 +lat_2=25 +lat_0=25 +lon_0=-95 +x_0=0 +y_0=0 +a=6371200 +b=6371200 +units=m +no_defs"
            longlatProj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

            .self$lccR <- raster::raster(
                nrows=428,
                ncols=614,
                crs = lccProj,
                ext = raster::extent(c(
                    xmin = -4232183.26089321,
                    xmax = 3253090.73910679,
                    ymin = -838789.970814831,
                    ymax = 4378958.02918517)),
                resolution = c(12191, 12191),
                vals=1)

            .self$longlatR <- raster::projectRaster(from = lccR,
                crs = longlatProj)

            #"hours since 2006-06-01T00:00:00Z"
            s <- strsplit(nc[['var']][['time1_bounds']][['units']], " ")[[1]]
            dt <- ncvar_get(nc, 'time1_bounds')
            dt <- switch(tolower(s[1]),
               'hours' = dt * 3600,
               "days" = dt * 3600 * 24)
            t0 <- as.POSIXct(s[length(s)], format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
            .self$field("TIME", t0 + dt)
            .self$field("DIMS", ncdim_get(nc))
            .self$field("VARS", names(nc[['var']]))
            .self$field("STEP", .self$step())
            e <- .self$get_extent()
            .self$BB <- raster::as.vector(e)
            #callSuper(nc)

        }) # methods

   )


#' Get the step size (resolution) as a guess - the space between the first
#' lon and lat values.  Subclasses with unambiguous steps sizes, like L3SMI and
#' MURSST, should override this method.
#'
#' @name SNAMANLRefClass_step
#' @param proj character string either "lcc" or "longlat" (default) to specify the
#'  projection of the returned
#' @return two element numeric vector of step size in x and y or NULL
#' These may be signed for descending values of lon and or lat
NULL
NAMANLRefClass$methods(
   step = function(proj = c("lcc", "longlat")[2]){
      switch(tolower(proj[1]),
        'lcc' = raster::res(.self$lccR),
        raster::res(.self$longlatR))
   })

#' Compute extent given a bounding box
#'
#' @name NAMANLRefClass_get_extent
#' @param bb the bounding box needed, if missing the current one is used, always
#'  assumed to be longlat
#' @param proj character string either "lcc" or "longlat" (default) to specify the
#'  projection of the returned
#' @return a raster::Extent object
NULL
NAMANLRefClass$methods(
   get_extent = function(bb, proj = c("lcc", "longlat")[2]){
      if (missing(bb)){
         ext <- switch(tolower(proj[1]),
            'lcc' = raster::extent(.self$lccR),
            raster::extent(.self$longlatR) )
         return(raster::as.vector(ext))
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
      if (tolower(proj) == "lcc"){
        xy <- sp::SpatialPoints(cbind(xx,yy), proj4string=sp::CRS(raster::projection(.self$longlatR)))
        xy <- sp::spTransform(xy, CRS(raster::projection(.self$lccR)))
        r <- raster::extent(xy)
      } else {
        r <- raster::extent(c(range(xx), range(yy)) )
      }
      r
   })

#' Get the longitude locations
#'
#' @name NAMANLRefClass_lon
#' @param what the desired location - 'leading', 'trailing', 'center' or 'native' (default)
#' Values 'trailing' and 'leading' are about the direction lon/lat are stored.  For example,
#' if lon/lat are each stored in ascending order then the leading edge is bottom and left.
#' On the other hand, is lon descends (east-to-west) then the leading edge is
#' bottom and right.  Values 'native and 'center' are assumed the same - if that
#' is not a correct for a class, then a method override is required
#' @return numeric vector or NULL
NULL
NAMANLRefClass$methods(
    lon = function(what = c('native', 'leading', 'trailing', 'center')[1]){
        # native is assumed to the center
        x <- raster::xFromCol(.self$longlatR)
        if (is.null(x)) return(NULL)
        s <- .self$STEP/2
        x <- switch(tolower(what[1]),
            'leading' = x + if (s[1] > 0) -s[1] else s[1],
            'trailing' = x + if(s[1] > 0) s[1] else -s[1],
            x)
   })

#' Get the latitude locations
#'
#' @name NAMANLRefClass_lat
#' @param what the desired location -  'leading', 'trailing', 'center' or 'native' (default)
#' Values 'trailing' and 'leading' are about the direction lon/lat are stored.  For example,
#' if lon/lat are each stored in ascending order then the leading edge is bottom and left.
#' On the other hand, is lat descends (north-to-south) then the leading edge is
#' top and left.  Values 'native and 'center' are likely the same.
#' @return numeric vector or NULL
NULL
NAMANLRefClass$methods(
   lat = function(what = c('native', 'leading', 'trailing', 'center')[1]){
      # native is assumed to the center
      y <- yFromRow(.self$longlatR)
      if (is.null(y)) return(NULL)
      s <- .self$STEP/2
      switch(tolower(what[1]),
         'leading' = y + if (s[2] > 0) -s[2] else s[2],
         'trailing' = y + if (s[2] > 0) s[2] else -s[2],
         y)
   })



#' Get a raster
#'
#' @name NAMANLRefClass_get_raster
#' @param what character one or more variable names or variable indices
#' @param layer numeric vector either a 1-based indices or POSIXct timestamps
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordiante reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NULL
NAMANLRefClass$methods(
   get_raster = function(what = .self$VARS[1], bb = .self$BB,
      crs = "+proj=longlat +datum=WGS84", flip = FALSE, time_fmt = "D%Y%j", ...){

      if (!all(what %in% .self$VARS))
         stop("one or more requested variable(s) not in data:", paste(what, collapse = "\n"))

      R <- NAMANL_get_raster(.self, what=what, bb = bb,
         crs = crs, flip = flip, time_fmt = time_fmt)
      return(R)
   }) # NAMANL_get_raster

#' Craft subset indices into a ncdf array
#'
#' @name NAMANLRefClass_subset_bbox
#'
#' @seealso \url{https://stat.ethz.ch/pipermail/r-help/2011-March/272641.html}
#' @param bb numeric, four element bounding box [left, right, bottom, top]
#' @return a list of \code{start} indices in x and y, \code{counts} in x and y and
#'    a possibly updated copy of \code{bb} vector of [left, right, bottom, top]
NULL
NAMANLRefClass$methods(
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


#' Get a raster for NAMANLRefClass
#'
#' @export
#' @param NC NAMANLRefClass object
#' @param what character one or more variable names or variable indices
#' @param bb a 4 element bounding box vector [left, right, bottom, top], defaults
#'    to [-180, 180, -90, 90]
#' @param crs character, the coordinate reference system to apply
#' @param flip logical if TRUE then flip the raster in the y direction
#' @param time_fmt if multiple time layers are returned, this controls the layer names
#' @param ... further arguments to select the parameter by name
#' @return a \code{raster::brick} or \code{raster::layer} object or NULL
NAMANL_get_raster <- function(NC, what = NC$VARS[1], bb = NC$BB,
    crs = "+proj=longlat +datum=WGS84", flip = FALSE, time_fmt = "D%Y%j", ...){

    stopifnot(inherits(NC, "NAMANLRefClass"))

    subnav <- NC$subset_bbox(bb)
    ext <- NC$get_extent(bb = subnav[['bb']])

    vd <- ncvardim_get(NC$NC)[[what]]
    ix <- names(vd) %in% c('x','y')
    vd <- vd[!ix]
    p <- list(...)

    # determine which are missing
    ip <- names(vd) %in% names(p)
    names(ip) <- names(vd)

    # for each required element in p named n
    # if (ip[n] == TRUE)
    #    if (TIME)
    #        convert to index via POSIXct or index
    #    else
    #        covert to index via indexing only
    #    start = coerce_within(index)[[1]]
    #    count = length(index)

    P <- vector(mode = 'list', length = length(ip))
    for (n in names(ip)){
        if (ip[[n]]){
            if(grepl('time', n, fixed = TRUE)){
                if(inherits(p[[n]], 'POSIXct')){
                    ix <- find_interval(p[[n]], NC$TIME)
                } else {
                    ix <- p[[n]]
                }
            } else {

            }# time?

        }

    }


    if(TRUE) return(list(vd=vd, p=p, subnav = subnav))

    x <- ncdf4::ncvar_get(NC$NC, what[1])
    R <- raster::raster(t(x), template = NC$lccR)
    R <- raster::projectRaster(from = R, to = NC$longlatR)
    R <- raster::crop(R, bb)

    if (flip) R <- raster::flip(R,"y")

    return(R)
}

