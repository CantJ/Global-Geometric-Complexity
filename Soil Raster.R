# This script is nessecary for ensuring the lithology (soil diversity) raster is formatted appropriately.
# This file will then be used for determining a measure of habitat complexity across space.

# Date last modified: April 2022
# Primary Author: James Cant
# ------------------------------------------------------------------------------------------------------------

#Clear the workspace
rm(list=ls(all=TRUE))

# load required packages
library(raster)
library(terra)
library(sf)
library(sp)

# Define load directory
raw_file_dir <- "/Volumes/Pocillopora/Geodiversity Data/Raw/Soil Database/"

# Save directory
fileDest <- '~/James/Geodiversity rasters/'

# temporary file directory
TempDest <- '/Volumes/Pocillopora/Geodiversity Data/Processed/'

# Import the Harmonized World Soil database raster file
hwsd <- raster(paste0(raw_file_dir, "hwsd.bil"))
# What details does this file contain
ncol(hwsd); nrow(hwsd); res(hwsd); extent(hwsd); projection(hwsd) # lots of data with a global extent, but no projection format
# add projection format
crs(hwsd) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# Take a look at the raster file.
plot(hwsd, col=bpy.colors(length(unique(hwsd)))) # currently the colours are associated with unique indexing values.
# the next step is to associate these unique index values with actual soil type categories.

# The harmonized world soil database has a soil attribute file that details the different soil categories and their associated indexing value.
# read in the attribute data
HWSD_cat90 <- read.csv(paste0(raw_file_dir, "HWSD Categories90.csv"), row.names = 1)
HWSD_cat85 <- read.csv(paste0(raw_file_dir, "HWSD Categories85.csv"), row.names = 1)
HWSD_cat74 <- read.csv(paste0(raw_file_dir, "HWSD Categories74.csv"), row.names = 1)
HWSD_data <- read.csv(paste0(raw_file_dir, "HWSD Data.csv"), row.names = 1)
# currently multiple soil types are associated with each indexing value (due to soil mixes) therefore for each index value with multiple representations I need to remove any duplicates
HWSD_data <- HWSD_data[which(HWSD_data$SEQ == 1),]

# Assign actual category names to the entries in the main attribute data file
HWSD_data$NAME <- NA
for(ii in 1:dim(HWSD_data)[1]){
  if(is.na(HWSD_data[ii,"SU_CODE90"])){ # if there is no classification data from 1990 then use the classification provided in 1985
    if(is.na(HWSD_data[ii,"SU_CODE85"])){ # if there is no classification data from neither 1990 nor 1985 then use the classification data provided in 1974
      HWSD_data[ii,'NAME'] <- HWSD_cat74[which(HWSD_cat74$CODE == HWSD_data[ii,'SU_CODE74']),]$VALUE
    } else {
      HWSD_data[ii,'NAME'] <- HWSD_cat85[which(HWSD_cat85$CODE == HWSD_data[ii,'SU_CODE85']),]$VALUE
    }
  } else {
    HWSD_data[ii,"NAME"] <- HWSD_cat90[which(HWSD_cat90$CODE == HWSD_data[ii,'SU_CODE90']),]$VALUE
  }
}

# For the rest of this section I am only interested in the Name and index identifier variables within the HWSD_data file
HWSD_data <- HWSD_data[,c('MU_GLOBAL','NAME')]

# I am also not interested in non-litholigical categories - particularly those that will be represented in one of the other geodiversity categories (ie. ocean and other water bodies)
# These categories can be omitted here.
HWSD_data[which(HWSD_data$NAME %in% c('???', 'Fishpond', 'Glaciers', 'Humanly disturbed', 'Island', 'No Data', 'Urban, Mining, etc.', 'Water bodies')), 'NAME'] <- NA

# Now rasters work on numerical values so I need to assign each category of soil classification with a unqiue indentifier (not the same as the index value already provided - this is for spatial referencing)
HWSD_data$CAT <- NA
soil.class <- unique(HWSD_data$NAME); soil.class <- soil.class[!is.na(soil.class)]
for(ii in 1:dim(HWSD_data)[1]){
  # progress read out
  print(ii)
  # find location of the matching soil classification from category list
  catNum <- match(HWSD_data[ii, 'NAME'], soil.class)
  # assign numerical category
  HWSD_data[ii, 'CAT'] <- catNum
}

# Now link (and replace) the raster index values with their corresponding soil classification number (inserting NA for index values for which there is now no category to be assigned)
# considering the size of the raster this is a big job - luckily there are functions that can be used to breakdown large raster before applying the functions over each smaller chunck
# determine for the HWSD raster what the appropriate number and size of each block should be for this
bSize <- blockSize(hwsd)
# create blank raster
HWSD_raster <- hwsd; values(HWSD_raster) <- NA

# initiate the use of raster blocking
HWSD_raster <- writeStart(HWSD_raster, filename = paste0(fileDest, 'SoilRaster'))

# run the value replacement
for(ii in seq_along(bSize$row)){
  # progress read out
  print(ii)
  # extract first raster chunk
  rasterSegment <- getValues(hwsd, row = bSize$row[ii], nrows = bSize$nrows[ii])
  # determine matching indexing values and there associated soil classification category
  indexCat <- HWSD_data[match(rasterSegment, HWSD_data$MU_GLOBAL),'CAT']
  # replace raster values with corresponding numerical soil classification 
  writeValues(HWSD_raster, indexCat, start = bSize$row[ii])
}

# close memory pathway used for raster blocking
HWSD_raster <- writeStop(HWSD_raster)

# Now for the moment of truth - what does the new soil raster look like
plot(HWSD_raster, col = bpy.colors(length(unique(HWSD_raster)))) # Now the colours are associated with identified soil categories.

# Handily the block processing of the hwsd raster has already saved this new raster file.

# ------------------------------------------------------------ End of Code -----------------------------------