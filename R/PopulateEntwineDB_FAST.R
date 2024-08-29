#
# PopulateEntwineDB
#
# Definitions:
#   Entwine collection: The collection of data that has been organized using entwine. 
#   Resides in amazon S3 public bucket
#   USGS collection: Entire USGS LPC collection that resides on the rockyweb server
#
# Both data collections are USGS-collected products. However, the entwine collection doesn't contain
# all of the data found in the USGS collection. In addition, you can't pull data directly from the
# rockyweb server for specific areas of interest. Instead you have to download tiles covering the
# target location and then pull data from the local copies of the tiles.
#  
# changes:
#   11/2/2021: Fixed the way I was handling multiple matches for the USGS identifier. Previously,
#   only the first match was being assigned to an Entwine polygon so I was ending up with Entwine
#   projects with no match when there actually was a match.
#
#   Also improved the logic used when matching projects using entwine project centroids. There are
#   still a few projects that are being matched to the wrong USGS project. For most of these, the 
#   match is coming from data collected prior to 2010. However, there are two projects in Indiana
#   that are placed in Kansas (Entwine names: USGS_LPC_IN_Johnson_2011_LAS_2016 and
#   USGS_LPC_IN_Marion_2011_LAS_2016). These projects also show up in Kansas in the Entwine data 
#   collection so the problem is not on my end. I will try to drop these from the enhanced index
#   but they appear to be duplicates and I don't know if there is a clean way to identify them.
#
#   Added a flag and logic to control removal of FullState projects from the Entwine collection.
#
#   Added a flag to do some final cleaning of the Entwine collection to remove projects that
#   don't have a match in the USGS collection and projects that are incorrectly located.
#
#   8/28/2024: The USGS WESM index is about 2Gb and usually takes a long time to download. Occasionally,
#   the download is fast (when running on GitHub) taking 10-30 minutes. In addition, the request for 
#   the WESM index often fails so the overall script fails. USGS has a .csv file with all the attribute
#   information that is small (~2Mb). Potentially, I could do all of the matching using the .csv
#   file and then try downloading the WESM index for the centroid search. There are about 30 polygons
#   matched using the centroid logic so these would not be updated if the WESM download fails.
#   Alternatively, I could save a file with the results of the centroid match that could be used 
#   when the WESM download fails (or maybe always?). All of the projects matched using the centroid
#   logic are older projects so the match is unlikely to change as newer data are added.
#
#   Implemented a new version that downloads the .csv file to get project attributes and reads
#   the existing enhanced index for the projects matched using the centroid search logic. This
#   new logic will work until there are projects that can't be matched using the attributes. 
#   I expect this will be rare. All of the projects matched using the centroid search are older
#   projects so the attribute data is not complete (or has errors).
library(sf)
library(dplyr)

# -------------------------------------------------------------------------------------------------
#                                         Configuration...
# -------------------------------------------------------------------------------------------------

# show current working directory
message(getwd())

# ---------->folder and filenames
# this is only needed when running from local instance of rstudio
# it will fail when run on GitHub so we catch the error
tryCatch(setwd("G:/R_Stuff/EntwineIndex"), error = function(e) NULL)

Folder <- ""
EntwinePolygonFile <- "resources.geojson"   # only for local file...https link is hard-coded for Howard's github location
EntwinePolygonLayer <- "resources"

# ---------->parameters
commonProjection <- 3857

# ---------->flags to control program flow
#UseLocalEntwinePolygonFile <- TRUE     # local runs
UseLocalEntwinePolygonFile <- FALSE    # github runs
#UseLocalUSGSPolygonFile <- TRUE      # for local runs
UseLocalUSGSPolygonFile <- FALSE     # for GitHub runs
UseUSGS_WESM <- TRUE

# flag to remove "FullState" projects. There are 3 of these as of 11/2/2021 and they
# cause some problems when searching for data coverage. Only the IA area has a corresponding
# project in the USGS collection.
# as of 11/2/2021, set this to TRUE to remove the FUllState areas.
#
# *****See additional comments near the code that drops these projects.
removeFullState <- TRUE

# flag to do some "manual" cleaning to remove Entwine projects without a matching USGS project
# and to remove projects that are incorrectly located in the Entwine collection. As of 11/2/2021,
# there are 2 duplicated projects from IN that have their polygons in KS.
manualClean <- TRUE

# flag to display final project boundaries
# we always want this off for GitHub runs
showMaps <- FALSE
if (showMaps)
  library(mapview)

# flag to save an extra copy of the enhanced index identified by date...not really needed and will add to amount of 
# storage on GitHub
saveDatedIndex <- FALSE

# -------------------------------------------------------------------------------------------------
#                                   Start of important bits...
# -------------------------------------------------------------------------------------------------
if (UseUSGS_WESM) {
  USGSPolygonFile <- "WESM.gpkg"   # only for local file...link is hard-coded for rockyweb
  USGSPolygonLayer <- "WESM"

  USGSProjectIDField <- "workunit"
  USGSDataURLField <- "lpc_link"
} else {
  USGSPolygonFile <- "FESM_LPC_Proj.gpkg"   # only for local file...https link is hard-coded for rockyftp
  USGSPolygonLayer <- "FESM_LPC_PROJ"
  
  USGSProjectIDField <- "project_id"
  USGSDataURLField <- "url"
}

# read USGS-entwine boundaries...in LatLong coords WGS84.
if (UseLocalEntwinePolygonFile) {
  File <- paste0(Folder, EntwinePolygonFile)
} else {
  File <- "https://raw.githubusercontent.com/hobu/usgs-lidar/master/boundaries/resources.geojson"
}
boundaries <- tryCatch(
  st_read(File, EntwinePolygonLayer, stringsAsFactors = FALSE),
  error = function(cond) {
    message("Entwine boundaries not available from HOBU GitHub repository")
    NULL
  }
)
if (is.null(boundaries)) stop()

#boundaries <- st_read(File, EntwinePolygonLayer, stringsAsFactors = FALSE)

# drop "FullState" projects. These are aggregations of other projects so, essentially, duplicate
# data. I think there is a corresponding "project" for the IA_FullState area in the USGS index
# but the project details are not consistent with other areas. Prior to the modifications to the
# centroid matching logic, attributes from a single project were assigned to the "FullState" areas
# but these attributes are not correct over the entire area. The FullState projects in MN is
# assigned attributes for a single project that may not be correct for the data. The KY FullState
# is not matched.
#
# The statewide or fullstate projects are problematic. For some states, the fullstate project duplicates
# all the data in the state but for others this is not the case. Some projects also have "statewide"
# in project names but the project may only cover a single county (e.g., Indiana). For Iowa, the fullstate
# project has older data covering the entire state but about 2/3 of the state has 2020 data not included
# in the full state project. For Minnesota, the fullstate project has data for areas not covered by any other 
# project.
#
# Overall, it may be better to keep these projects and drop clips based on some criteria that finds duplicate
# clips. However, some of these aggregated projects may not have collection date information, etc because
# the matching logic below fails to find a match.
#
# 4/12/2024 modified the code so it finds any project with "fullstate" in its name and drops them.
#
if (removeFullState) {
  t <- grepl("fullstate", boundaries$name, ignore.case = TRUE)
  
  boundaries <- subset(boundaries, !t)

  message("Removed FullState projects from Entwine collection")
}

# reproject project boundaries to web mercator
EntwineboundariesWebMerc <- st_transform(boundaries, crs = commonProjection)

# clean up
rm(boundaries)

if (UseLocalUSGSPolygonFile) {
  # read geopackage for entire USGS collection
  File <- paste(Folder, USGSPolygonFile, sep = "")
} else {
  File <- "https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/WESM.csv"
  #File <- "https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/metadata/WESM.gpkg"
}
attributes <- tryCatch(
  read.csv(File, stringsAsFactors = FALSE),
  error = function(cond) {
    message("USGS WESM attributes not available from rockyweb server")
    NULL
  }
)
if (is.null(attributes)) stop()

#attributes <- read.csv(File, stringsAsFactors = FALSE)
#system.time({boundaries <- st_read(File, USGSPolygonLayer, stringsAsFactors = FALSE)})

USGSboundariesWebMerc <- attributes
#USGSboundariesWebMerc <- st_transform(boundaries, crs = commonProjection)

rm(attributes)

# We need to manipulate values in columns (replace characters, remove whitespace, etc). This is 
# problematic if we want to select columns using [] because the return for an sf object is a
# data frame instead of a vector. We could use $ and the column name but we want our column
# name to be a variable so we can use different index files easily. The work-around is to copy 
# the geometry to a temporary object, manipulate the columns, then re-attach the geometry.
# However, when manipulating the columns, you can't remove or add any rows or else the geometry
# won't match.

# work on USGS index
# copy the geometry
#g <- st_geometry(USGSboundariesWebMerc)

# drop the geometry
#USGSboundariesWebMerc <- st_drop_geometry(USGSboundariesWebMerc)

# manipulate the columns...don't delete or rearrange rows
# cleanup USGS project name...contain "-", " ", and leading/trailing spaces. Nothing very magical about
# the cleanup process...just necessary to get things to match
# replace "-" with "_"
USGSboundariesWebMerc[, USGSProjectIDField] <- gsub("-", "_", USGSboundariesWebMerc[, USGSProjectIDField])

# get rid of leading and trailing spaces
USGSboundariesWebMerc[, USGSProjectIDField] <- trimws(USGSboundariesWebMerc[, USGSProjectIDField])

# replace any remaining spaces with "_"
USGSboundariesWebMerc[, USGSProjectIDField] <- gsub(" ", "_", USGSboundariesWebMerc[, USGSProjectIDField])

# add a sequential field
USGSboundariesWebMerc$SEQ <- seq_len(nrow(USGSboundariesWebMerc))

# add a field indicating that the project_id or workunit is NA
USGSboundariesWebMerc$drop <- is.na(USGSboundariesWebMerc[, USGSProjectIDField])

# re-attach the geometry
#USGSboundariesWebMerc <- st_set_geometry(USGSboundariesWebMerc, g)

# get rid of any features with NA in the project_id (or workunit) field
# this MUST be done after re-attaching the geometry
USGSboundariesWebMerc <- USGSboundariesWebMerc[!USGSboundariesWebMerc$drop, ]

# drop projects that are 1/9 Arc Second
# not a good idea...you end up with many older projects in the Entwine collection being
# matched to newer projects in the USGS collection. These matches will result in bad
# attributes (collection dates, crs info, etc) for the older projects
#USGSboundariesWebMerc <- subset(USGSboundariesWebMerc, !grepl("1/9", USGSboundariesWebMerc$project))

#rm(g)

# work on Entwine index...don't need to deal with geometry because we are using $ to
# select columns
# add project name to entwine collection by parsing url
EntwineboundariesWebMerc$ENTpid <- basename(dirname(EntwineboundariesWebMerc$url))

# do clean-up of entwine project identifiers...probably not necessary but done anyway to be safe
EntwineboundariesWebMerc$ENTpid <- gsub("-", "_", EntwineboundariesWebMerc$ENTpid)

# get rid of leading and trailing spaces
EntwineboundariesWebMerc$ENTpid <- trimws(EntwineboundariesWebMerc$ENTpid)

# replace any remaining spaces with "_"
EntwineboundariesWebMerc$ENTpid <- gsub(" ", "_", EntwineboundariesWebMerc$ENTpid)

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
EntwineboundariesWebMerc$eSEQ <- -1

# additional field to indicate how match was done
# meaning:
# 0 = no match
# 1 = single match...USGS collection name is found in 1 Entwine collection name
# 2 = multiple matches...USGS collection name is found in several Entwine collection names
# 3 = centroid match...centroid of Entwine polygon (largest polygon in boundary) is within a single USGS polygon (also largest polygon in boundary)
# 4 = centroid match...centroid of Entwine polygon (largest polygon in boundary) is within a multiple USGS polygon (also largest polygon in boundary)
#                       best match is determined by lpc_pub_date and last 4 digits of entwine project name
# 5 = manual match
EntwineboundariesWebMerc$MatchMethod <- 0

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
  if (USGSProjectIDField == "project_id") {
    result <- grep(tolower(USGSboundariesWebMerc$project_id[thePoly]), tolower(EntwineboundariesWebMerc$ENTpid), fixed = TRUE)
  } else if (USGSProjectIDField == "workunit") {
    result <- grep(tolower(USGSboundariesWebMerc$workunit[thePoly]), tolower(EntwineboundariesWebMerc$ENTpid), fixed = TRUE)
  }
  
  if (length(result) == 0) {
    nomatches <- nomatches + 1
  } else if (length(result) > 1) {
    for (m in 1:length(result)) {
      EntwineboundariesWebMerc$eSEQ[result[m]] <- USGSboundariesWebMerc$SEQ[thePoly]
        
      EntwineboundariesWebMerc$MatchMethod[result[m]] <- 2
        
      multiplematches <- multiplematches + 1
    }
  } else {
    EntwineboundariesWebMerc$eSEQ[result] <- USGSboundariesWebMerc$SEQ[thePoly]
      
    EntwineboundariesWebMerc$MatchMethod[result] <- 1
      
    matches <- matches + 1
  }
}

message("Using ", nrow(USGSboundariesWebMerc), " items in the USGS collection and ",
    nrow(EntwineboundariesWebMerc), " items in the Entwine collection:\n",
    "  USGS ", USGSProjectIDField, " found in Entwine ENTpid:\n",
    "    No matches: ", nomatches, "\n",
    "    Multiple matches: ", multiplematches, "\n",
    "    Single matches: ", matches
)

missing <- EntwineboundariesWebMerc[EntwineboundariesWebMerc$eSEQ == -1, ]

message("After project identifier matching, there are ", nrow(missing), " Entwine polygons with no match")
rm(missing)

# read the existing index and get the entries with MatchMethod = 3 or 4. In general,
# these are older projects that were matched using centroid logic so the matches
# shouldn't change as new projects are added to the USGS collection.
oldIndex <- tryCatch(
  st_read(paste0(Folder, "Index/ENTWINEBoundaries.gpkg"),
          layer = "ENTWINEBoundaries", stringsAsFactors = FALSE),
  error = function(cond) {
    message("Enhanced entwine index is not available from the EntwineIndex repository")
    NULL
  }
)
if (is.null(oldIndex)) stop()

#oldIndex <- st_read(paste0(Folder, "Index/ENTWINEBoundaries.gpkg"),
#                    layer = "ENTWINEBoundaries", stringsAsFactors = FALSE)

# drop entries matched using attributes
oldIndex <- oldIndex[oldIndex$MatchMethod >= 3,]

message("Read ", nrow(oldIndex), " entries from existing enhanced index to get results of centroid matching")

# get rid of entwine entries that have no match
NewEntwineboundariesWebMerc <- EntwineboundariesWebMerc[EntwineboundariesWebMerc$MatchMethod > 0, ]

# there is 1 area in IN that has a typo in the USGS project name and it has no lpc_pub_date value.
# this results in an incorrect match that we can correct manually
# entwine ENTpid = USGS_LPC_IN_Hendricks_2011_LAS_2016
# USFS workunit = IN_HendricksCo_2011
NewEntwineboundariesWebMerc$eSEQ[NewEntwineboundariesWebMerc$ENTpid == "USGS_LPC_IN_Hendricks_2011_LAS_2016"] <- USGSboundariesWebMerc$SEQ[USGSboundariesWebMerc$workunit == "IN_HendricksCo_2011"]
NewEntwineboundariesWebMerc$MatchMethod[NewEntwineboundariesWebMerc$ENTpid == "USGS_LPC_IN_Hendricks_2011_LAS_2016"] <- 5


names(NewEntwineboundariesWebMerc)[names(NewEntwineboundariesWebMerc) == "eSEQ"] <- "SEQ"

# merge columns from USGS polygons into entwine polygons
final <- NewEntwineboundariesWebMerc %>% left_join(USGSboundariesWebMerc, by = "SEQ")

# change name of geometry column...not sure why it is "geom"
names(oldIndex)[names(oldIndex) == "geom"] <- "geometry"
st_geometry(oldIndex) <- "geometry"

# append entries read from old index to the matched entries
NewEntwineboundariesWebMerc <- rbind(final, oldIndex)

# count the number of entwine polygons without a match
missing <- NewEntwineboundariesWebMerc[NewEntwineboundariesWebMerc$SEQ == -1, ]
if (nrow(missing)) {
  message("After all matching, there are still ", nrow(missing), " Entwine polygons with no matching USGS polygon:")
  for (name in missing$name) message(name)
}

if (manualClean) {
  # remove projects without a match
  NewEntwineboundariesWebMerc <- subset(NewEntwineboundariesWebMerc, !(SEQ == -1))
  
  # remove projects that are in the wrong place
  # 2 projects for IN are showing up in KS. There are duplicate projects in the correct locations
  NewEntwineboundariesWebMerc <- subset(NewEntwineboundariesWebMerc, !(ENTpid == "USGS_LPC_IN_Johnson_2011_LAS_2016" & workunit == "KS_FTRILEY_2006"))
  NewEntwineboundariesWebMerc <- subset(NewEntwineboundariesWebMerc, !(ENTpid == "USGS_LPC_IN_Marion_2011_LAS_2016" & workunit == "KS_RILEY_2010"))
  
  message("Manual cleaning complete")
}

# write off new entwine polygons...includes polygons that don't have a match in the USGS collection.
if (saveDatedIndex) {
  st_write(NewEntwineboundariesWebMerc,
    paste0(Folder, "Index/", format(Sys.Date(), format = "%Y_%m_%d"), "_", "ENTWINEBoundaries.gpkg"),
    layer = "ENTWINEBoundaries",
    delete_dsn = TRUE
  )
}

# write to Index folder so the updated boundaries will be uploaded to github
st_write(NewEntwineboundariesWebMerc,
         paste0(Folder, "Index/ENTWINEBoundaries.gpkg"),
         layer = "ENTWINEBoundaries",
         delete_dsn = TRUE
)

# # write copy with date to Index folder so the boundaries will be uploaded to github
# IndexFolder <- "G:\\R_Stuff\\EntwineIndex\\Index\\"
# write_sf(NewEntwineboundariesWebMerc,
#          paste(dirname(IndexFolder), "/", basename(IndexFolder), "/", format(Sys.Date(), format = "%Y_%m_%d"), "_", "ENTWINEBoundaries.gpkg", sep = ""),
#          layer = "ENTWINEBoundaries"
# )

if (showMaps) mapviewOptions(fgb = FALSE)
if (showMaps) mapview(list(NewEntwineboundariesWebMerc, USGSboundariesWebMerc))
if (showMaps) mapview(as(NewEntwineboundariesWebMerc, "Spatial"))
