
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EntwineIndex

<!-- badges: start -->

[![UpdateEntwinePlusIndex](https://github.com/bmcgaughey1/EntwineIndex/actions/workflows/main.yml/badge.svg)](https://github.com/bmcgaughey1/EntwineIndex/actions/workflows/main.yml)
<!-- badges: end -->

EntwineIndex is a code repository for R code that merges the Entwine
lidar data index
[index](https://raw.githubusercontent.com/hobu/usgs-lidar/master/boundaries/resources.geojson)
created by Howard Butler (GitHub repository:
[usgs-lidar](https://github.com/hobu/usgs-lidar)) with the [WESM
index](https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/WESM.gpkg)
for the USGS 3DEP data collection. The goal is to add lidar project
information to the Entwine index to facilitate querying the index for
data covering specific locations and dates.

Some projects are removed from the Entwine index. These include
FullState projects for IA, KY, and MN and two counties in Indiana that
are incorrectly located in Kansas. For these Indiana projects duplicates
exist in the Entwine collection that are correctly located.

This code is run every few days to maintain synchronization with the
Entwine and WESM index files.

URL for the index is:
`https://raw.githubusercontent.com/bmcgaughey1/EntwineIndex/main/Index/ENTWINEBoundaries.gpkg`

You can grab a copy of the
[index](https://raw.githubusercontent.com/bmcgaughey1/EntwineIndex/main/Index/ENTWINEBoundaries.gpkg)
and store it locally. However, I do not recommend maintaining a local
copy as the index changes frequently. The index is fairly small so
direct reads take only a few seconds. You can read the index directly
from GitHub using *st_read* from the *sf* package as shown in the
example code below.

The index is used with my
[USGSlidar](https://github.com/bmcgaughey1/USGSlidar) R package to help
discover lidar data for specific locations and dates. The index can be
accessed in the package using *fetchUSGSProjectIndex(type =
“entwineplus”)*.

The enhanced index can also be downloaded directly from this GitHub
repository and read using the following R code:

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

The original version of the code downloaded the USGS WESM index with
geometries. Over time, this index has grown quite large and can take a
long time to download (depends on the load on the rockyweb server). In
August 2024, i added a new way to build the enhanced index that uses
only the attributes from the WESM index (stored on rockyweb as .csv
file). The attribute file is much smaller than the index with geometries
so it is more likely to be successfully downloaded when updating the
index. However, because the code is not using the geometries (boundary
polygons for lidar projects), it can’t use the logic that searches for
USGS projects that contain the centroid of entwine projects. To maintain
the information for these projects, the existing enhanced index is read
and the projects with matchMethod = 3, 4, or 5 are copied into the new
index. All of the projects matched using the centroid logic are old and
have incomplete attributes so I expect this approach will continue to
work as new projects are added to both collections. It may be necessary
to use the older version of the code once in a while to update all
projects. However, downloading the WESM index hasn’t worked consistently
enough to rely on this download to update the index.
