# This script is for ensuring global lake and river data (hydrology) are formatted appropriately.
# These files will then be used for determining a measure of habitat complexity across space.
# Due to the size of the datafiles - particularly the river data - this script and therefore its associated file directories are for working on the BioTime server.
# To keep server clear all necessary raw files need loading to the server and removing when finished with.

# Date last modified: May 2022
# Primary Author: James Cant
# ------------------------------------------------------------------------------------------------------------

#Clear the workspace
rm(list=ls(all=TRUE))

# load required packages
library(raster)
library(terra)
library(sf)
library(sp)
library(maptools)
library(tidyr)
library(rgdal)
library(gdalUtilities)
library(readr)

# Define load directory
raw_file_dir <- "~/James/Raw geodiversity data/"

# Save directory
fileDest <- '~/James/Geodiversity rasters/'

# temporary file directory
TempDest <- '~/James/TempFile/'

#################################################
# STEP 1: Lake Diversity
#################################################

# Load the lakes shapefiles (automatically loads with a WGS84 projection)
LAKES <- st_read(paste0(raw_file_dir, "HydroLAKES_polys_v10.shp"))

# firstly the lakes shapefile needs converting to a raster file (this will just be a 'simple' raster showing whether areas of the globe are lake or not)
# Load a raster of the world - this will serve as a template over which to re-project the lakes shapefile
# generate a world raster template with NAs for ocean and zeros for land.
temp_raster <- raster(paste0(raw_file_dir, 'worldRaster.tif'), resolution = c(0.008333333, 0.008333333), ext = extent(c(-180, 180, -90, 90)))
# set non-ocean pixels in the template raster to 0. Because this is a large raster this will need to be done in raster blocks
# determine for the HWSD raster what the appropriate number and size of each block should be for this
bSize <- blockSize(temp_raster)

# initiate the use of raster blocking
worldRaster <- writeStart(temp_raster, filename = paste0(raw_file_dir, 'worldRaster.grd'))

# run the value replacement
for(ii in seq_along(bSize$row)){
  # progress read out
  print(ii)
  # extract first raster chunk
  rasterSegment <- getValues(temp_raster, row = bSize$row[ii], nrows = bSize$nrows[ii])
  # convert all non-NA values to zero
  rasterSegment[!is.na(rasterSegment)] <- 0
  # replace raster values with corresponding numerical soil classification 
  writeValues(worldRaster, rasterSegment, start = bSize$row[ii])
}

# close memory pathway used for raster blocking
worldRaster <- writeStop(worldRaster)
# remove temporary raster file using in conversion
rm(temp_raster)

# Now rasterize the lakes polygon file.
Lake_raster <- fasterize(LAKES, worldRaster, field = 'Lake_type', fun = 'sum')
# Now were not actually interested in lake type - only whether a lake is present so the raster values need reassigning

# again initiate the use of raster blocking
LakeRaster <- writeStart(Lake_raster, filename = paste0(raw_file_dir, 'LakesRaster.grd'))

# run the value replacement
for(ii in seq_along(bSize$row)){
  # progress read out
  print(ii)
  # extract first raster chunk
  rasterSegment <- getValues(Lake_raster, row = bSize$row[ii], nrows = bSize$nrows[ii])
  # convert all non-NA values to zero
  rasterSegment[!is.na(rasterSegment)] <- 1
  # replace raster values with corresponding numerical soil classification 
  writeValues(LakeRaster, rasterSegment, start = bSize$row[ii])
}

# close memory pathway used for raster blocking
LakeRaster <- writeStop(LakeRaster)
# and merge the raster with the base global template raster (the raster now displays 0 for no lake and 1 for lake)
LakeRaster <- merge(LakeRaster, worldRaster)

# although the lake raster will require merging with the river raster it is worth saving the file here.
raster::writeRaster(LakeRaster, paste0(fileDest, 'LakeRaster.grd'))

#################################################
# STEP 2: River Diversity
#################################################

# free up some memory space
rm(rasterSegment,LAKES,LakeRaster,Lake_raster,ii)

# define pathway for loading the rivers shapefiles (automatically loads with a WGS84 projection)
RIVERS <- st_read(paste0(raw_file_dir, "HydroRIVERS_v10.shp"))

# This is a large shapefile, and needs to be broken down for processing.
setwd(TempDest) # processing will be done into a temporary file location
# determine how many sections to split the RIVERS data into for processing
linestrings <- 1:nrow(RIVERS) # number of entries
# Split entries into n parts
n <- 250
parts <- split(linestrings, cut(linestrings, n))

# now define a function that will work through chunks of the RIVERS data to implement rasterisation - storing the output as a temporary file
manualRastTrim <- function(ii, sfData, tmplt, parts){
  # progress read out
  print(ii)
  # extract chunk of shapefile
  try(sfVect <- sfData[parts[[ii]], c('ORD_STRA','geometry')])
  # convert to a spatVector (for use in terra to speed things up)
  try(SptVect <- vect(sfVect))
  # run rasterising using terra (without saving to the output to the server)
  try(terra::writeRaster(trim(terra::rasterize(SptVect, tmplt, field = 'ORD_STRA', fun = 'sum')), filename = paste0('TrimmedRaster_', ii, '.tif'), overwrite = FALSE))
}  
# and run rasterization (first converting the worldRaster template into a spatRaster)
worldSpat <- rast(worldRaster)
lapply(1:n, manualRastTrim, sfData = RIVERS, tmplt = worldSpat, parts = parts)

# randomly some of those run-through fail so they need troubleshooting here
# determine which iterations failed
non_error_files <- list.files("/home/jic2/James/TempFile", '.tif')
non_error_files <- parse_number(non_error_files)
error_files <- c(1:n)
error_files <- error_files[!(error_files %in% non_error_files)]
# a couple of manual run-throughs suggest that these should have worked (so the failed entries will be re-run here)
lapply(error_files, manualRastTrim, sfData = RIVERS, tmplt = worldSpat, parts = parts)

# Identify all available trimmed rasters 
trimmedfiles <- list.files(pattern = '^TrimmedRaster')
# collate together these rasters
gdalbuildvrt(trimmedfiles, 'River.vrt')
gdalwarp('River.vrt', 'RiverRast.tif', co=c("BIGTIFF=YES", "COMPRESS=DEFLATE", "TILED=TRUE"))

# and save raster file (works as a checkpoint)
terra::writeRaster(rast('RiverRast.tif'), filename = paste0(fileDest,'RiverRast.tif'), overwrite = FALSE)

#################################################
# STEP 3: Merge Rasters
#################################################

# free up some memory space
rm(error_files,parts,RIVERS, linestrings, n, non_error_files, trimmedfiles, manualRastTrim)

# reassign working directory
setwd(fileDest)

# Load lake and river rasters
LakeRaster <- rast('LakeRaster.grd')
River_Raster <- raster('RiverRast.tif') # needed in this format for applying chunking

# Before merging the lake and river rasters it is necessary to +1 onto all river order classifications.
# This way the hydrology raster will simultaneously differentiate between No water (0), Lakes (1) & Rivers (2+), and different rivers orders (2-9)

# Given the size of the raster it is necessary to apply this adjustment across raster chunks.
bSize <- blockSize(River_Raster)

# Initiate the use of raster blocking
RiverRast <- writeStart(River_Raster, filename = paste0(TempDest, 'TmpRiverRast.tif'), overwrite = TRUE)

# run the value replacement
for(ii in seq_along(bSize$row)){
  # progress read out
  print(ii)
  # extract first raster chunk
  rasterSegment <- getValues(River_Raster, row = bSize$row[ii], nrows = bSize$nrows[ii])
  # a little fix to remove NaNs
  rasterSegment[is.nan(rasterSegment)] <- NA
  # +1 to all non-NA values
  rasterSegment[!is.na(rasterSegment)] <- rasterSegment[!is.na(rasterSegment)] + 1
  # replace raster values with corresponding numerical soil classification 
  writeValues(RiverRast, rasterSegment, start = bSize$row[ii])
}

# close memory pathway used for raster blocking
RiverRast <- writeStop(RiverRast)

# Convert river raster to spatraster
RiverRast <- rast(RiverRast)

# and merge the two hydrology rasters (the raster now displays 0 for no water, 1 for lake, and 2+ for river)
HydroRast <- merge(RiverRast, LakeRaster)

# And save to complete Hydroraster processing
writeRaster(HydroRast, filename = paste0(fileDest, 'HydroRaster.tif'))

#--------------------------------------------------------------------- End of Code ----------------------------------