# This script is for processing bathymetry and topography data from across the global into a single raster file.
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
library(terrainr)
library(sp)
library(fasterize)

# Define file directories
bathy_dir <- "/Volumes/Pocillopora/Geodiversity_Data/Raw/Bathymetry/"
DEM_dir <- "/Volumes/Pocillopora/Geodiversity_Data/Raw/NASADEM/"

# Save directory
fileDest <- '/Volumes/Pocillopora/Geodiversity_Data/Processed/'

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


