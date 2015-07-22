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

# you may need to reinstantiate your reference *if* the definition or methods have changed in you edits
x <- SPNC(some_resource)
```

### Usage

We try to follow the KISS principle by minimizing the exposure of details.

#### Instantiate
+ `SPNC(ncdf_resource_or_filename, bb)` create an instance of SPNCRefClass or subclass with the provided bounding box.

#### Properties
#### Methods

+ `open()` open the NetCDF connection
+ `close()` close the NetCDF connection
+ `flavor()` retrieve the flavor of the NetCDF data (*e.g.* [source='OISST',type='raster', local = TRUE/FALSE/NA])
+ `get_raster(bb, t)` retrieve a n-d array as SpatialRaster within the bounding box for the given times
+ `get_path(bb, t)` retrieve a path SpatialLinesDataFrame within the bounding box for the given times (not sure what this is just yet!)
+ `get_points(bb, t)` retrieve set of points as SpatialPointsDataFrame within the bounding box for the given times

#### Known sources of data

Establish the basics...

```R
librayr(raster)
library(spnc)
bb <- c(-72,-63,39,46)
```

+ `MURSST` [Multi-scale Ultra-high Resolution Sea Surface Temperature](http://mur.jpl.nasa.gov/)

```R
MUR_file <- 'http://thredds.jpl.nasa.gov/thredds/dodsC/ncml_aggregation/OceanTemperature/ghrsst/aggregate__ghrsst_JPL-L4UHfnd-GLOB-MUR.ncml'
MUR <- SPNC(MUR_file, bb = bb)
# Reference Class: "MURSSTRefClass" (from the global environment) 
#   flavor: MURSST:raster 
#   state: opened 
#   bounding box: -72 -63 39 46 
#   VARS: analysed_sst analysis_error mask sea_ice_fraction 
#   DIMS: lat=16384 lon=32768 time=3963 
#   LON: [ -179.994995, 179.994995] 
#   LAT: [ -89.994499, 89.994499] 
#   TIME: [ 2002-06-01 09:00:00, 2013-04-06 09:00:00] 

r <- MUR$get_raster(what = 'analysed_sst', layer = 1)
raster::spplot(r)
```

+ `OISST` [NOAA Optimum Interpolation (OI) Sea Surface Temperature ](http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html)

```R
OISST_file = 'http://www.ncdc.noaa.gov/thredds/dodsC/OISST-V2-AVHRR_agg'
OI <- SPNC(OISST_file, bb = bb)
```

+ `L3SMI`  [Ocean Color Level 3 Standard Mapped Image](http://oceancolor.gsfc.nasa.gov/cms)

```R
filename <- 'V20151932015200.L3m_8D_NPP_CHL_chlor_a_4km.nc'
L3 <- SPNC(filename, bb = bb)
raster::spplot(log10(L3))
```

+ `MODISA`  [Modis Aqua](http://oceancolor.gsfc.nasa.gov/cms/)

```R
MODISA_file <- 'http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2008/001/A2008001.L3m_DAY_CHL_chlor_a_4km.nc'
MO <- SPNC(MODISA_file, bb = bb)
```


