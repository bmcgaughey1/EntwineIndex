
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EntwineIndex

<!-- badges: start -->

[![UpdateEntwinePlusIndex](https://github.com/bmcgaughey1/EntwineIndex/actions/workflows/main.yml/badge.svg)](https://github.com/bmcgaughey1/EntwineIndex/actions/workflows/main.yml)
<!-- badges: end -->

EntwineIndex is a simple code repository for code that merges the
Entwine lidar data
[index](https://raw.githubusercontent.com/hobu/usgs-lidar/master/boundaries/resources.geojson)
created by Howard Butler (GitHub repository:
[usgs-lidar](https://github.com/hobu/usgs-lidar)) with the [WESM
index](https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/WESM.gpkg)
for the USGS 3DEP data collection. The goal is to add lidar project
information to the Entwine index to facilitate querying the index for
data covering specific locations and dates.

This code is run daily to maintain synchronization with the Entwine and
WESM index files.

URL for the index is:
`https://raw.githubusercontent.com/bmcgaughey1/EntwineIndex/main/Index/ENTWINEBoundaries.gpkg`

You can grab a copy of the
[index](https://raw.githubusercontent.com/bmcgaughey1/EntwineIndex/main/Index/ENTWINEBoundaries.gpkg)
and store it locally. However, the index is fairly small so direct reads
take only a few seconds. You can read the index directly from GitHub
using *st\_read* from the *sf* package as shown in the example below.

The index is used with my
[USGSlidar](https://github.com/bmcgaughey1/USGSlidar) R package to help
discover lidar data for specific locations and dates. The index can be
accessed in the package using *fetchUSGSProjectIndex(type =
“entwineplus”)*.

``` r
url <- "https://raw.githubusercontent.com/bmcgaughey1/EntwineIndex/main/Index/ENTWINEBoundaries.gpkg"

# download to local file using utils library
library(utils)
if (!utils::download.file(url, "ENTWINEBoundaries.gpkg", mode = "wb",)) {
  projects <- sf::st_read("ENTWINEBoundaries.gpkg", stringsAsFactors = FALSE)
}

# use the sf package to read the index directly from GitHub
library(sf)
projects <- sf::st_read(url, stringsAsFactors = FALSE)
```
