# work on USGS data to get details for entwine collection
#
# Definitions:
#   Entwine collection: The collection of data that has been organized using entwine. Resides in amazon S3 public bucket
#   USGS collection: Entire USGS LPC collection that resides on the rockyftp server
#
# Both data collections are USGS-collected products. However, the entwine collection doesn't contain all of the data found
# in the USGS collection. In addition, you can't pull data directly from the rockyftp server. Instead you have to download
# tiles covering the target location and then pull data from the local copies of the tiles. The PDAL pipeline will be
# similar for both collections (readers will change).
#
library(rgdal)
library(sf)
library(PROJ)
library(raster)
library(jsonlite)
library(dplyr)
library(mapview)

# -------------------------------------------------------------------------------------------------
#                                         Configuration...
# -------------------------------------------------------------------------------------------------

# links to various project index files
# if (tolower(type) == "fesm") {
#   url <- "ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Elevation/LPC/FullExtentSpatialMetadata/FESM_LPC_Proj.gpkg"
# } else if (tolower(type) == "wesm") {
#   url <- "ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/WESM.gpkg"
# } else if (tolower(type) == "entwine") {
#   url <- "https://raw.githubusercontent.com/hobu/usgs-lidar/master/boundaries/resources.geojson"
# }
# as of 7/16/2021 the rockyftp server may be dead. You can access the WESM index on amazon at this link:
# http://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/metadata/WESM.gpkg

# ---------->folder and filenames
Folder <- "G:\\R_Stuff\\PlotClipping\\"
EntwinePolygonFile <- "resources.geojson"   # only for local file...https link is hard-coded for Howard's github location
EntwinePolygonLayer <- "resources"
USGSPolygonFile <- "FESM_LPC_Proj.gpkg"   # only for local file...https link is hard-coded for rockyftp
USGSPolygonLayer <- "FESM_LPC_PROJ"
# this is the file that is now being used (after 5/7/2020). As of 6/11/2020 it is not fully populated
# USGSPolygonFile <- "WESM.gpkg"   # only for local file...link is hard-coded for rockyftp
# USGSPolygonLayer <- "main.WESM"

# ---------->parameters
commonProjection <- CRS(SRS_string="EPSG:3857")
#commonProjection <- proj_create("EPSG:3857", format = 1)

# *****Field names are different depending on the source of the USGS polygons
# FESM_LPC_Proj.gpkg
USGSProjectIDField <- "project_id"
USGSDataURLField <- "url"
# WESM.gpkg
# USGSProjectIDField <- "WorkUnit"
# USGSDataURLField <- "LPC_Link"

# ---------->flags to control program flow
UseLocalEntwinePolygonFile <- FALSE
UseLocalUSGSPolygonFile <- TRUE
UseUSGS_WESM <- TRUE

# flag to display final project boundaries
showMaps <- FALSE

# -------------------------------------------------------------------------------------------------
#                                   Start of important bits...
# -------------------------------------------------------------------------------------------------
if (UseUSGS_WESM) {
  USGSPolygonFile <- "WESM_9_21_2021.gpkg"   # only for local file...link is hard-coded for rockyweb
  USGSPolygonLayer <- "WESM"

  # USGSProjectIDField <- "WorkUnit"
  # USGSDataURLField <- "LPC_Link"
  USGSProjectIDField <- "workunit"
  USGSDataURLField <- "lpc_link"
}

# verify that we have the geoJSON reader...not needed
#"GeoJSON" %in% ogrDrivers()$name

# read USGS-entwine boundaries...in LatLong coords WGS84.
# In theory, you should be able to grab this file directly from GitHub.
# However, R starts to get really sluggish after I use getURL() so it may
# be better to download the resources.geojson file and then read the local copy.
if (UseLocalEntwinePolygonFile) {
  File <- paste(Folder, EntwinePolygonFile, sep = "")

  # look at the file properites...not necessary but may be useful
  # ogrInfo(File)
  # ogrListLayers(File)

#  boundaries <- readOGR(File, EntwinePolygonLayer, stringsAsFactors = FALSE)
} else {
  # get the URL by displaying the map and then using the download button
  # this is slower than reading a local file but should get you the latest information
  # site URL: https://github.com/hobu/usgs-lidar/tree/master/boundaries
#  fileURL <- c("https://raw.githubusercontent.com/hobu/usgs-lidar/master/boundaries/resources.geojson")
  File <- "https://raw.githubusercontent.com/hobu/usgs-lidar/master/boundaries/resources.geojson"

  #  ogrInfo(File)

#  boundaries <- readOGR(File, EntwinePolygonLayer, stringsAsFactors = FALSE)
}
boundaries <- readOGR(File, EntwinePolygonLayer, stringsAsFactors = FALSE)

# reproject project boundaries to web mercator
EntwineboundariesWebMerc <- spTransform(boundaries, CRS = commonProjection)

# clean up
rm(boundaries)

if (UseLocalUSGSPolygonFile) {
  # read geopackage for entire USGS collection
  File <- paste(Folder, USGSPolygonFile, sep = "")
} else {
  # this is MUCH slower than reading a local file but it does work and will get you the latest information
  # file is 45Mb+ so you have to download the entire file to open...takes at least 15 minutes
  #File <- "ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/LPC/FullExtentSpatialMetadata/FESM_LPC_Proj.gpkg"
  #File <- "ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/metadata/WESM.gpkg"
  File <- "https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/WESM.gpkg"
}
boundaries <- readOGR(File, USGSPolygonLayer, stringsAsFactors = FALSE)
USGSboundariesWebMerc <- spTransform(boundaries, CRS = commonProjection)

rm(boundaries)

# cleanup USGS project name...contain "-", " ", and leading/trailing spaces. Nothing very magical about
# the cleanup process...just necessary to get things to match
# replace "-" with "_"
USGSboundariesWebMerc@data[, USGSProjectIDField] <- chartr(old = "-", new = "_", USGSboundariesWebMerc@data[, USGSProjectIDField])
#USGSboundariesWebMerc@data$project_id <- chartr(old = "-", new = "_", USGSboundariesWebMerc@data$project_id)

# get rid of leading and trailing spaces
USGSboundariesWebMerc@data[, USGSProjectIDField] <- trimws(USGSboundariesWebMerc@data[, USGSProjectIDField])

# replace any remaining spaces with "_"
USGSboundariesWebMerc@data[, USGSProjectIDField] <- chartr(old = " ", new = "_", USGSboundariesWebMerc@data[, USGSProjectIDField])

# get rid of any features with NA in the project_id field
USGSboundariesWebMerc <- USGSboundariesWebMerc[!is.na(USGSboundariesWebMerc@data[, USGSProjectIDField]), ]

# add a sequential field
USGSboundariesWebMerc@data$SEQ <- seq_len(nrow(USGSboundariesWebMerc))

# add project name to entwine collection by parsing url
EntwineboundariesWebMerc@data$ENTpid <- basename(dirname(EntwineboundariesWebMerc@data[, "url"]))

# do clean-up of entwine project identifiers...probably not necessary but done anyway to be safe
EntwineboundariesWebMerc@data$ENTpid <- chartr(old = "-", new = "_", EntwineboundariesWebMerc@data$ENTpid)

# get rid of leading and trailing spaces
EntwineboundariesWebMerc@data$ENTpid <- trimws(EntwineboundariesWebMerc@data$ENTpid)

# replace any remaining spaces with "_"
EntwineboundariesWebMerc@data$ENTpid <- chartr(old = " ", new = "_", EntwineboundariesWebMerc@data$ENTpid)

# now we have the sets of boundaries in web mercator projection
# we want to match up projects in the entwine collection to their counterpart in the USGS collection so
# we can populate the entwine polygons with attributes (like the collection date)

# I have spent way too long trying to get things to work. I started by matching entwine entries to those
# in USGS collection to populate additional information to the entwine polygons. Then I decided it was
# better to populate the USGS polygons with the entwine project info. However, I failed to recognize that
# there are multiple projects in the entwine collection that are actually covered by a single polygon
# in the USGS collection. Now I am back to populating the entwine polygons...basically want to join all
# fields from the USGS polygon to the record for the entwine polygon. We will do this once using names
# to match and then try a different approach for those entwine projects that don't have a match.
#
# This logic seems appropriate when working with the entwine collection. If you want to work with the
# full USGS collection, it may make sense to simply clean up the project identifier in the USGS attributes
# and then work directly with the USGS boundaries.

# additional field to add to entwine collection that points to USGS record number
EntwineboundariesWebMerc@data$eSEQ <- -1

# additional field to indicate how match was done
# meaning:
# 0 = no match
# 1 = single match...USGS collection name is found in 1 Entwine collection name
# 2 = multiple matches...USGS collection name is found in several Entwine collection names
# 3 = centroid match...centroid of Entwine polygon (largest polygon in boundary) is within USGS polygon (also largest in boundary)
# 4 = single match...Entwine collection name is found in 1 USGS collection name
# 5 = multiple matches...Entwine collection name is found in several USGS collection names
EntwineboundariesWebMerc@data$MatchMethod <- 0

# search for project name in the USGSProjectIDField field in the USGS collection for matching
# entries in the entwine collection. In most cases, the entwine labels are longer.
#
# I am noticing that the Entwine collection often has more than 1 project that matches the USGS name. It
# looks like this happens when a portion of the project is either delivered late or reflown. The Entwine
# collection keeps the area with later delivery as a separate area but the USGS collection does not. As the
# logic below works, it only considers the first Entwine project as a match, leaving subsequent matches
# labeled as unmatched. Given that I want to match plot measurement dates with lidar acquisition dates,
# it may not make sense to drop the information for the later delivery that may be "found" using the
# centroid search logic. It may make sense to do the name search after the centroid search to see if
# there are areas that will otherwise be unmatched...kind of a final alternative assignment that is better
# than nothing. However, I suspect that the date information for the possible later acquisition will be lost.
#
#
# In a run on 2/17/2021, I only have 3 projects in the Entwine collection that don't get a match. One is the
# full state data for Kentucky that was assembled for a special project, one is NM_Albuquerque_2010, and the
# last is USGS_LPC_SD_MORiver_Woolpert_B3_2016_LAS_2018. The last is actually matched by SD_MORiver_Woolpert_B3_2016
# but there were 2 matches involving this USGS project and only the first match was used to populate the
# Entwine collection.
matches <- 0
multiplematches <- 0
nomatches <- 0
for (thePoly in 1:nrow(USGSboundariesWebMerc)) {
#  if (!is.na(USGSboundariesWebMerc@data[thePoly, USGSProjectIDField])) {
    result <- grep(tolower(USGSboundariesWebMerc@data[thePoly, USGSProjectIDField]), tolower(EntwineboundariesWebMerc@data[, "ENTpid"]), fixed = TRUE)

    if (length(result) == 0) {
#      cat("No match for polygon", thePoly, "USGS:", USGSboundariesWebMerc@data[thePoly, USGSProjectIDField], "\n")
      nomatches <- nomatches + 1
    } else if (length(result) > 1) {
#      cat("***Multiple matches for polygon", thePoly, "USGS:", USGSboundariesWebMerc@data[thePoly, USGSProjectIDField], "Entwine:", EntwineboundariesWebMerc@data[result, "ENTpid"], "\n")

      EntwineboundariesWebMerc@data[result[1], "eSEQ"] <- USGSboundariesWebMerc@data[thePoly, "SEQ"]

      EntwineboundariesWebMerc@data[result[1], "MatchMethod"] <- 2

      multiplematches <- multiplematches + 1
    } else {
#      cat("Match for polygon", thePoly, "USGS:", USGSboundariesWebMerc@data[thePoly, USGSProjectIDField], "Entwine:", EntwineboundariesWebMerc@data[result, "ENTpid"], "\n")

      EntwineboundariesWebMerc@data[result, "eSEQ"] <- USGSboundariesWebMerc@data[thePoly, "SEQ"]

      EntwineboundariesWebMerc@data[result, "MatchMethod"] <- 1

      matches <- matches + 1
    }
#  }
}

cat("USGS collection name is contained in Entwine collection name:\n",
  "  No matches ", nomatches, "\n",
  "  Multiple matches ", multiplematches, "\n",
  "  Single matches ", matches, "\n"
)

# The above approach gets us a match for most projects in the entwine collection. There
# could also be additional projects due to the multiple matches (USGS project ID is found in several
# entwine project labels). Looking at the polygon sets, it looks like there are matches for just about
# all entwine projects on the USGS FTP server (this makes total sense because the entwine collection
# should be a subset of all USGS projects). Looking at the database for the entwine collection after
# running through the matching logic, there are a few projects that don't have a match. Cleaning up the
# names in the USGS names was essential for the matching process.
#
# Another approach...compute the centroid (or label point) for the entwine project polygons and intersect
# these points with the USGS polygons. This might give us matches for the cases where names are different.
# However, there are some non-closed polygons so we could end up with centroids that are outside the
# project polygon so we still might have missing cases. Try this to find matches for those without a name
# match and maybe for those with multiple name matches. This approach works but some centroids are
# covered by more than one collection. Project name match should take the highest priority, then
# use the centroid approach to find possible matches for entwine projects without a matching name. In addition,
# there is the problem where data for a USGS project was delivered in multiple batches. I don't know if
# this is due to reflights or just reprocessing for some tiles but several projects have "LAS_date"
# (with "date" containing a year) added to the workunit name. The "problem" with this is that
# there could be different dates associated with the data so if the logic makes a bad match, we
# will end up with a mismatch.
#
# Not really part of this code but a solution might be to use the GPS time stamp in the LAS data when
# doing plot clips to get a definitive acquisition date for the point data. Then the date would be
# available when matching plot clips to FIA plot data.
#
# compute centroids for entwine polygons...use of_largest_polygon = TRUE to give centroid of the
# largest polygon in multi-polygon features
missing <- EntwineboundariesWebMerc[EntwineboundariesWebMerc@data$eSEQ == -1, ]

cat(length(missing), "Entwine polygons with no match\n")

# This throws a warning regarding attributes are constant over geometries of x but I
# can't figure out why. Function seems to work fine and returns good points.
# the conversion to sf object and use of st_agr() is to prevent the warning. Previous
# code did the conversion in the call to st_intersect (commented line)
missing_sf <- st_as_sf(missing)
st_agr(missing_sf) <- "constant"
cent <- st_centroid(missing_sf, of_largest_polygon = TRUE)
#cent <- st_centroid(st_as_sf(missing), of_largest_polygon = TRUE)

USGSboundariesWebMerc <- rgeos::gBuffer(USGSboundariesWebMerc,
               byid = TRUE,
               width = 0
)

# intersect centroids with USGS polygons...this throws a warning regarding proj4 strings but I can't figure
# out what to do to make the warning go away...both datasets are in the same projection
# The only difference is that the proj4 string for cent contains "+wktext"
t <- raster::intersect(as(cent, "Spatial"), USGSboundariesWebMerc)

cat(length(t), "intersections with entwine polygon centroids\n")

# set the USGS record in the entwine data...entwine id starts at 0
NewEntwineboundariesWebMerc <- EntwineboundariesWebMerc
NewEntwineboundariesWebMerc@data[t@data$id + 1, "eSEQ"] <- t@data$SEQ

NewEntwineboundariesWebMerc@data[t@data$id + 1, "MatchMethod"] <- 3

names(NewEntwineboundariesWebMerc@data)[names(NewEntwineboundariesWebMerc@data) == "eSEQ"] <- "SEQ"

# merge columns from USGS polygons into entwine polygons
final <- NewEntwineboundariesWebMerc@data %>% left_join(USGSboundariesWebMerc@data, by = "SEQ")
NewEntwineboundariesWebMerc@data <- final

# count the number of entwine polygons without a match
missing <- NewEntwineboundariesWebMerc[NewEntwineboundariesWebMerc@data$SEQ == -1, ]
cat("Still have", length(missing), "entwine polygons with no matching USGS polygon\n")
missing@data$name

# write off new entwine polygons...includes polygons that don't have a match in the USGS collection.
# for geopackage, you provide the full file name in dsn and the layer name in layer
# for shapefiles, you provide the folder name without the trailing / in dsn and the file name without .shp in the layer
writeOGR(NewEntwineboundariesWebMerc,
  dsn = paste(dirname(Folder), "/", basename(Folder), "/", format(Sys.Date(), format = "%Y_%m_%d"), "_", "ENTWINEBoundaries.gpkg", sep = ""),
  layer = "ENTWINEBoundaries",
  overwrite_layer = TRUE,
  driver = "GPKG")

# write to Index folder so the updated boundaries will be uploaded to github
writeOGR(NewEntwineboundariesWebMerc,
         dsn = "G:/R_Stuff/EntwineIndex/Index/ENTWINEBoundaries.gpkg",
         layer = "ENTWINEBoundaries",
         overwrite_layer = TRUE,
         driver = "GPKG")

# write copy with date to Index folder so the boundaries will be uploaded to github
IndexFolder <- "G:\\R_Stuff\\EntwineIndex\\Index\\"
writeOGR(NewEntwineboundariesWebMerc,
         dsn = paste(dirname(IndexFolder), "/", basename(IndexFolder), "/", format(Sys.Date(), format = "%Y_%m_%d"), "_", "ENTWINEBoundaries.gpkg", sep = ""),
         layer = "ENTWINEBoundaries",
         overwrite_layer = TRUE,
         driver = "GPKG")

if (showMaps) mapview(list(NewEntwineboundariesWebMerc, USGSboundariesWebMerc))
if (showMaps) mapview(NewEntwineboundariesWebMerc)

