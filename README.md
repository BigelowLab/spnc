### spnc: Spatial and NetCDF

`spnc` provides a uniform interface to [NetCDF](http://www.unidata.ucar.edu/software/netcdf) files for accessing data in the `sp` framework of classes.  Using Unidata's netcdf_4 library permits access to both locally stored data (as .nc files) as well as network served data (as OpeNDAP/DODS/Thredds/etc/).

Under the hood we use David Pierce's [ncdf4](http://cran.r-project.org/web/packages/ncdf4/index.html) R package to interface with the Unidata netcdf_4 library.  We try to keep the underlying resources opaque to the user so that we may switch to another interface to Unidata's netcdf_4 library as the need arises.

#### Requirements

+ [R](http://www.r-project.org/)
+ [NetCDF](http://www.unidata.ucar.edu/software/netcdf)  C-library

R packages... may require substantial effort to install on some platforms

+ [ncdf4](http://cran.r-project.org/web/packages/ncdf4/index.html)
+ [sp](http://cran.r-project.org/web/packages/sp/)
+ [raster](http://cran.r-project.org/web/packages/raster/)
    
#### Installation

The repository is set up for installation using the [devtools](http://cran.r-project.org/web/packages/devtools/) package.

```R
library(devtools)
# for your own library
install_github("btupper/spnc")
# or for system wide installation (requires admin-level permissions)
with_lib(.Library, install_github("btupper/spnc"))
```

If you are developing new data source readers and have the repo on your local platform then use the following.

```R
library(devtools)
pkg_path <- '/path/to/spnc'
install(pkg_path)
```

If you edit the code (with documentation as needed!) then you can easily re-install - sort of on the fly.

```R
library(spnc)
x <- SPNC(some_resource)

# edit, document and then redo the document/install process
library(devtools)
pkg_path <- '/path/to/spnc'

document(pkg_path)
install(pkg_path)

# you may need to reinstantiate your reference *if* the definition or methods have changed in your edits
x <- SPNC(some_resource)
```

### Usage

We try to follow the KISS principle by minimizing the exposure of details.  See the [wiki](https://github.com/btupper/spnc/wiki) for examples.

#### Instantiate

+ `SPNC(ncdf_resource_or_filename, bb)` create an instance of SPNCRefClass or subclass with the provided bounding box.

`bb` is a 4 element vector that follows the pattern of the `Extent` class in the `sp` package: [xmin, xmax, ymin, ymax].  If you set the bounding box when you instantiate the SPNC object, then it is used as the default for subsequent calls to the methods shown below, unless you specifically override the value. 

#### Properties
#### Methods

+ `open()` open the NetCDF connection
+ `close()` close the NetCDF connection
+ `flavor()` retrieve the flavor of the NetCDF data (*e.g.* [source='OISST',type='raster', local = TRUE/FALSE/NA])
+ `get_raster(bb, t)` retrieve a n-d array as SpatialRaster within the bounding box for the given times
+ `get_path(bb, t)` retrieve a path SpatialLinesDataFrame within the bounding box for the given times (not sure what this is just yet!)
+ `get_points(bb, t)` retrieve set of points as SpatialPointsDataFrame within the bounding box for the given times

#### Known sources of data

+ `L3SMI`  [Ocean Color Level 3 Standard Mapped Image](http://oceancolor.gsfc.nasa.gov/cms)

The Ocean Color group provides Level3 Standard Mapped Images (L3SMI) for a variety of products. See more in the [wiki](https://github.com/btupper/spnc/wiki/MODIS-and-OBPG) on how to programmatically search their database.  Tested on MODISA CHL and SST.  Note, a companion package, [obpgcrawler](https://github.com/btupper/obpgcrawler), is available to programmatically find [Ocean Color](http://oceancolor.gsfc.nasa.gov/cms) data.


+ `MURSST` [Multi-scale Ultra-high Resolution Sea Surface Temperature](http://mur.jpl.nasa.gov/)

JPL's PO.DAAC provides this data as an NCML (filename.ncml) - an XML file that the NetCDF-4 library happily reads as if it were a NetCDF file.   THis provides easy access to data at different times, so not only does one collect data from geographic subregions but for different times.  See the [wiki](https://github.com/btupper/spnc/wiki/MUR-and-NCML) for details.


+ `OISST` [NOAA Optimum Interpolation (OI) Sea Surface Temperature ](http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html)

Currently this will access (a) locally downloaded files and (b) OpeNDAP files for a single time (dailies).  Access to the aggregate NCML OpeNDAP resources is still under development. Note that longitude is stored in the 0-360 format which requires some juggling internally.  You can still specify bounding boxes in [-180, 180] longitude form.  Here we access the same data two ways... 

+ (a) download and unpack this file 

ftp://eclipse.ncdc.noaa.gov/pub/OI-daily-v2/NetCDF/2004/AVHRR-AMSR/amsr-avhrr-v2.20040101.nc.gz

+ (b) connect via OpeNDAP to this resource 

http://www.ncdc.noaa.gov/thredds/dodsC/oisst/NetCDF/AVHRR-AMSR/2004/AVHRR-AMSR/amsr-avhrr-v2.20040101.nc

... and visually compare the two (which are the same!)

```R
library(raster)
library(spnc)

BB <- c(-72,-63,39,46)
localFile <- '/Users/ben/Downloads/amsr-avhrr-v2.20040101.nc'
opendapFile <- 'http://www.ncdc.noaa.gov/thredds/dodsC/oisst/NetCDF/AVHRR-AMSR/2004/AVHRR-AMSR/amsr-avhrr-v2.20040101.nc'

X <- SPNC(localFile)
x <- X$get_raster(what = 'sst', bb = BB)

Y <- SPNC(opendapFile)
y <- Y$get_raster(what = 'sst', bb = BB)

spplot(stack(x,y))

# not yet for NCML
# OISST_file = 'http://www.ncdc.noaa.gov/thredds/dodsC/OISST-V2-AVHRR_agg'
# OI <- SPNC(OISST_file, bb = bb)
```

+ `NHSCE` [Northern Hemisphere Snow Cover Extent](https://climatedataguide.ucar.edu/climate-data/snow-cover-extent-northern-hemisphere-climate-data-record-rutgers) 

```R
snow <- SPNC('http://www.ncdc.noaa.gov/thredds/dodsC/cdr/snowcover/nhsce_v01r01_19661004_latest.nc')
```

Under development. Provides mask and areal extent using a modified grid.
