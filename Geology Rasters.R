# This script is for ensuring the global geological diversity data (lithologyl) are formatted appropriately.
# The raster file produced will then be used for determining a measure of habitat complexity across space.

# Date last modified: May 2022
# Primary Author: James Cant
# ------------------------------------------------------------------------------------------------------------

#Clear the workspace
rm(list=ls(all=TRUE))

# load required packages
library(raster)
library(sf)
library(sp)
library(fasterize)

# Define load directory
raw_file_dir <- "~/James/Raw geodiversity data/"

# Save directory
fileDest <- '~/James/Geodiversity rasters/'

#################################################
# STEP 1: Load and format raw data
#################################################

# load in Global Lithology data (data originally downloaded from the Global Lithological Map Database v1.0 (gridded to 0.5° spatial resolution) [PANGAEA, https://doi.org/10.1594/PANGAEA.788537])
# Open geodatabase file.
gdb <- path.expand(paste0(raw_file_dir, "LiMW_GIS 2015.gdb"))
# Identify available layers
ogrListLayers(gdb)
# read in datafile (opens as a spatial polygon file)
GLiM <- readOGR(gdb, "GLiM_export")

# Load in a global raster template (consistent with the template used across the other geodiversity variables)
worldRaster <- raster(paste0(raw_file_dir, 'worldRaster.grd'))
# store the crs associated with this template
crs_use <- crs(worldRaster)

# reformat the data ready for rasterising
# convert file to simple features
GLiM <- as(GLiM, 'sf')
# convert lithology categories into a factor 
GLiM$Litho <- as.factor(GLiM$Litho)
# re-project Geology data
GLiM <- st_transform(GLiM, crs = crs_use)

#################################################
# STEP 2: Generate raster file
#################################################

# Rasterise spatial data
RockRaster <- fasterize(GLiM, worldRaster, field = 'Litho')
# check it has worked
plot(RockRaster, col = bpy.colors(50))
# and save
writeRaster(RockRaster, paste0(fileDest, 'RockRaster.tif'))

# ----------------------------------------------------------------------- End of Code -----------------------------------