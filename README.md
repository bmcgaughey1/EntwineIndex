
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EntwineIndex

<!-- badges: start -->

<!-- badges: end -->

EntwineIndex is a simple code repository for code that merges the
Entwine lidar data
[index](https://raw.githubusercontent.com/hobu/usgs-lidar/master/boundaries/resources.geojson)
created by Howard Butler with the WESM index
(<ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/WESM.gpkg>)
for the USGS 3DEP data collection. The goal is to add lidar project
information to the Entwine index to facilitate querying the index for
data for specific dates.

Ultimately, this code will run every night/morning to maintain
synchronization with the Entwine index.

The direct link to the
[index](https://raw.githubusercontent.com/bmcgaughey1/EntwineIndex/main/Index/ENTWINEBoundaries.gpkg)
is:
`https://raw.githubusercontent.com/bmcgaughey1/EntwineIndex/main/Index/ENTWINEBoundaries.gpkg`

When using the index, it is probably best to grab a copy of the index
and store it locally. However, you can read the index directly from
GitHub.

``` r
url <- "https://raw.githubusercontent.com/bmcgaughey1/EntwineIndex/main/Index/ENTWINEBoundaries.gpkg"

# use rgdal to read the geopackage
library(rgdal)

# read directly
projects <- rgdal::readOGR(dsn = url, verbose = F, stringsAsFactors = FALSE)

# download to local file, then read using utils library
library(utils)
if (!utils::download.file(url, "ENTWINEBoundaries.gpkg", mode = "wb",)) {
  projects <- rgdal::readOGR(dsn = "ENTWINEBoundaries.gpkg", verbose = F, stringsAsFactors = FALSE)
}
```
