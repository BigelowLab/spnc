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
install_github("btupper/spnc")
```

### Usage

We try to follow the KISS principle by minimizing the exposure of details.

#### Instantiate
+ `SPNC(ncdf_resource_or_filename, bb)` create an instance of SPNCRefClass or subclass with the provided bounding box.

#### Properties
#### Methods

+ `open()` open the NetCDF connection
+ `close()` close the NetCDF connection
+ `flavor()` retrieve the flavor of the NetCDF data (*e.g.* [source='OISST',type='raster'])
+ `get_raster(bb, t)` retrieve a n-d array as SpatialRaster
+ `get_path(x, y, t)` retrieve a path SpatialLinesDataFrame
+ `get_points(x, y, t)` retrieve set of points as SpatialPointsDataFrame 

#### Known sources of data

+ `MUR`
+ `OISST`
+ `MODISA`
