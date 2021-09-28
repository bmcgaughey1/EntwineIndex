
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EntwineIndex

<!-- badges: start -->

<!-- badges: end -->

EntwineIndex is a simple code repository for code that merges the
Entwine lidar data
[index](https://raw.githubusercontent.com/hobu/usgs-lidar/master/boundaries/resources.geojson)
created by Howard Butler with the WESM index
(<https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/WESM.gpkg>)
for the USGS 3DEP data collection. The goal is to add lidar project
information to the Entwine index to facilitate querying the index for
data for specific dates.

Ultimately, this code will run every night/morning to maintain
synchronization with the Entwine index.

The direct link to the
[index](https://raw.githubusercontent.com/bmcgaughey1/EntwineIndex/main/Index/ENTWINEBoundaries.gpkg)
is:
`https://raw.githubusercontent.com/bmcgaughey1/EntwineIndex/main/Index/ENTWINEBoundaries.gpkg`

When using the index, it may be best to grab a copy of the index and
store it locally. However, you can read the index directly from GitHub.
It is fairly small so the download only takes a few seconds.

``` r
url <- "https://raw.githubusercontent.com/bmcgaughey1/EntwineIndex/main/Index/ENTWINEBoundaries.gpkg"

# use rgdal to read the geopackage
library(sf)

# read directly
projects <- sf::st_read(url, stringsAsFactors = FALSE)

# download to local file using utils library, then read
library(utils)
if (!utils::download.file(url, "ENTWINEBoundaries.gpkg", mode = "wb",)) {
  projects <- sf::st_read("ENTWINEBoundaries.gpkg", stringsAsFactors = FALSE)
}
```
