# This script is for processing bathymetry and topography data from across the globe into a single raster file.
# This file will then be used to estimate fractal dimension and surface rugosity at varying spatial scales.
# N.B. A secondary output of this script will be a (terrestrially focused) global raster of high resolution topography from which various measures of geomorphometry can be obtained.

# Date last modified: May 2022
# Primary Author: James Cant
# ------------------------------------------------------------------------------------------------------------

#Clear the workspace
rm(list=ls(all=TRUE))

# load required packages
library(raster)
library(terra)
library(gdalUtilities)
library(sp)
library(fasterize)
library(rgdal)

# Define file directories
bathy_dir <- "/Volumes/Pocillopora/Geodiversity_Data/Raw/Bathymetry/"
DEM_dir <- "/Volumes/Acropora/NASADEM/" # NASADEM hgt files
NASADEM_dir <- '/Volumes/Pocillopora/Geodiversity_Data/Raw/NASADEM/' # Processed NASADEM tif files.

# Define temporary file directories
TmpDir <- '/Users/james/Documents/R/R_TempFiles/' # (to be cleared at the end of each session)
VRTDir <- '/Volumes/Pocillopora/Geodiversity_Data/Raw/VRTFiles/' # for storing vrt and NASADEM tif tiles

# Save directory
fileDest <- '/Volumes/Pocillopora/Geodiversity_Data/Processed/'
TmpfileDest <- '/Volumes/Pocillopora/Geodiversity_Data/Raw/'

#################################################
# STEP 1: Generate Bathymetry Raster
#################################################

# Load in rasters - these files have been downloaded from GEBCO (gebco.net).
# They are global gridded 15 arc second (450 meter) resolution topography and bathymetry files. 

# identify the raster files (.tif) available.
setwd(bathy_dir)
tifs <- list.files(pattern = "\\.tif$")
# define file paths for each raster file
filepaths <- paste0(bathy_dir, tifs)

# generate raster files (working with terra not raster)
bathy_rast_list <- src(lapply(filepaths, rast))

# and merge into a singular file
BathRaster <- terra::mosaic(bathy_rast_list)

# take a look
BathRaster <- raster(BathRaster); plot(BathRaster, col=bpy.colors(20))

# Save file as a check point
raster::writeRaster(BathRaster, paste0(fileDest, 'BathRaster.grd'))

# Note this raster contains both bathymetry and topographic data.
# However the NASADEM data below provides topographical elevation data at a higher resolution.
# Therefore, the NASADEM data will be used to re-project the Bathymetry raster (to ensure matching resolution), before being merged with it.

#################################################
# STEP 2: Generate NASADEM Raster
#################################################

# reset the working directory
setwd(DEM_dir)

# determine the zip files present
zipfiles <- list.files()
# identify hgt files names
hgtfiles <- substr(zipfiles, 13, 19)

# work through each NASADEM zip file to extract raster information and save as tif file (for easier processing)
# define processing function
RastExtract <- function(ii){
  # Progress read out
  print(ii)
  # unzip file
  unzip(zipfiles[ii], ex = TmpDir)
  # generate and save Spatraster (faster processing)
  tempRast <- rast(paste0(TmpDir, hgtfiles[ii], '.hgt'))
  terra::writeRaster(tempRast, filename = paste0(NASADEM_dir, hgtfiles[ii], '.tif'))
  # remove zip file contents from temp directory
  file.remove(list.files(TmpDir, full.names = TRUE))
}

# run NASADEM Raster extraction
lapply(1:length(zipfiles), RastExtract)

# and collate together the raster files as a Spatraster collections
# reassign the working directory
setwd(NASADEM_dir)
rastfiles <- list.files(pattern = '.tif') # Identify available files

# there is a lot of files here and processing them all at once overloads GDAL so it needs to be done in batches.
# determine how many sections to split the data into for processing
DataCount <- 1:length(rastfiles) # number of entries
# Split entries into n parts
n <- 145
parts <- split(DataCount, cut(DataCount, n))

# now define a function that will work through chunks of the data to implement raster mosaicing
VRTtoTIF <- function(ii){
  # progress read out
  print(ii)
  # collate selected raster files together within a virtual raster tile 
  gdalbuildvrt(rastfiles[parts[[ii]]], paste0(VRTDir, 'NASADEM', ii, '.vrt'))
  # reformat the virtual raster tile into a .tif file.
  gdal_translate(src_dataset = paste0(VRTDir, 'NASADEM', ii, '.vrt'), dst_dataset = paste0(VRTDir, 'NASADEM', ii, '.tif'), co=c("BIGTIFF=YES", "COMPRESS=DEFLATE", "TILED=TRUE"))
}  

# Run VRT tiling
lapply(1:n, VRTtoTIF)

# This produces 145 large raster tiles that will again be pieced together in a VRT and translated into one large singular raster file
setwd(VRTDir)
# How many large raster tiles?
tifFiles <- list.files(pattern = '.tif')
tifFiles <- tifFiles[!grepl('.tif.aux.xml', tifFiles)]
# collate these large tile files into a single virtual raster file
gdalbuildvrt(tifFiles, paste0(TmpfileDest, 'NASADEM_FULL.vrt'))
# and reformat into a complete global (terrestrial) raster
gdal_translate(src_dataset = paste0(TmpfileDest, 'NASADEM_FULL.vrt'), dst_dataset = paste0(fileDest, 'NASADEM_FULL.tif'), co=c("BIGTIFF=YES", "COMPRESS=DEFLATE", "TILED=TRUE"))

# open raster (to check all has worked)
TopoRaster <- rast(paste0(fileDest, 'NASADEM_FULL.tif'))
plot(TopoRaster)

#################################################
# STEP 3: Unify Rasters
#################################################

# Load in required rasters (as Spatrasters as these process faster)
BathRaster <- raster(paste0(fileDest, 'BathRaster.grd'))
TopoRaster <- rast('NASADEM.tif')

# save file raster file
GlobalDEM

#################################################
# STEP 4: Determine measures of complexity
#################################################

# Define initial parameters for initiating the fractal dimension assessment.
# Starting coordinates
xmn <- ext(GlobalDEM)[1]; xmx <- ext(GlobalDEM)[2]; ymn <- ext(GlobalDEM)[3]; ymx <- ext(GlobalDEM)[4]
# Specify dimensions of grid cells for which fractal dimension should be estimated.
grid_dim <- seq(x) 
L0 <- res(GlobalDEM)[1] # minimum scale that height range can be estimated for

test <- projectRaster(from = GlobalDEM, crs = "+proj=laea +datum=WGS84 +units=m +no_defs")

# Add grid cells across the raster
D = 3 - S

S = 


floor(((long + 180)/6)+1)


terra::project()


## -------------------------------------------------------- End of Code -----------------------------------








