# This script is for generating outputs rasters detailing global patterns in fractal dimension, surface rugosity, and height range at a resolution of ~1km (1240m)
# This script deals with a large quantity of data and therefore is slow to complete and requires a lot of RAM.

# Date last modified: July 2023
# Primary Author: James Cant (jic2@st-andrews.ac.uk)
# ------------------------------------------------------------------------------------------------------------

# Clear workspace
rm(list=ls(all=TRUE))

# Load required packages
library(arrow)
library(raster)
library(dplyr)
library(sf)
library(terra)

######################################
# STEP 1: Housekeeping
######################################

# define raw data file pathway
filePath <- '~/James/Raw geodiversity data/'

# define file save pathway
fileSave <- '~/James/FinalisedRasters/'

# Load datafile using memory effective approach
ComplexValues <- read_parquet(paste0(filePath, 'GlobalComplexity.parquet'), as_data_frame = F)

# define resolution
L <- 1240
L0 <- L/5 # map resolution is equal to 5 times the resolution of the original DEM.

# Define spatial boundaries
options(digits = 22) # Allow R to cope with more digits to prevent values being rounded down.
maxY <- as.numeric(ComplexValues %>% summarise(across(Lat, max)) %>% collect()) + L
maxX <- as.numeric(ComplexValues %>% summarise(across(Lon, max)) %>% collect()) + L
minY <- as.numeric(ComplexValues %>% summarise(across(Lat, min)) %>% collect())
minX <- as.numeric(ComplexValues %>% summarise(across(Lon, min)) %>% collect())

# Define global projection
Molle <- crs('ESRI:54009')

######################################
# STEP 2: Generate initial global rasters
######################################

# Extract spatial coordinates
coords <- cbind(as.numeric(ComplexValues[['Lon']]), as.numeric(ComplexValues[['Lat']]))

# Convert pixel values of interest into a spatial data frame (be patient, this is a slow running process)
# Complexity values
HRange <- suppressWarnings(SpatialPixelsDataFrame(coords, data = data.frame(as.numeric(ComplexValues[['H']])), proj4string = CRS(Molle)))
Rugosity <- suppressWarnings(SpatialPixelsDataFrame(coords, data = data.frame(as.numeric(ComplexValues[['R']])), proj4string = CRS(Molle)))
FracD <- suppressWarnings(SpatialPixelsDataFrame(coords, data = data.frame(as.numeric(ComplexValues[['D']])), proj4string = CRS(Molle)))
# Terrestrial mask
LandMask <- suppressWarnings(SpatialPixelsDataFrame(coords, data = data.frame(as.numeric(ComplexValues[['Terrestrial']])), proj4string = CRS(Molle)))
# Assign raster value IDs
names(HRange) <- 'H'
names(Rugosity) <- 'R'
names(FracD) <- 'D'
names(LandMask) <- 'Land'

# And convert to raster files
# Height range
HRast <- rast(HRange, layer = 'H')
# Rugosity
RRast <- rast(Rugosity, layer = 'R')
# Fractal Dimension
DRast <- rast(FracD, layer = 'D')
# Land Template
MaskRast <- rast(LandMask, layer = 'Land')

# view
plot(HRast)
plot(RRast)
plot(DRast)
plot(MaskRast) # is being projected upside down.
# reverse raster display
MaskRast <- flip(MaskRast, direction = 'v')
plot(MaskRast)

######################################
# STEP 3: Final Quality control
######################################

# Evaluate data ranges
# Fractal Dimension
range(values(DRast), na.rm = TRUE) # should be between 2 and 3.
hist(values(DRast)) # not quite a normal distribution - but this could be due to the mix of terrestrial and marine environments

# Rugosity 
range(values(RRast), na.rm = TRUE) # should be greater than or equal to 1
hist(values(RRast)) # Poisson distribution
# Following Torres-Pulliza et al (2020) Rugosity can be standardized to (R^2)-1
R_Rast <- app(RRast, fun = function(x) { (x^2)-1 })
hist(log10(values(R_Rast))) # transformed to normal

# Height Range
range(values(HRast), na.rm = TRUE) # must be non-negative
hist(values(HRast)) # Poisson distribution
# Following Torres-Pulliza et al (2020) Height range can be standardized to HR/sqrt(2)*L0
HR_Rast <- app(HRast, fun = function(x) { x/(sqrt(2)*L0) })
hist(log10(values(HR_Rast))) # transformed to normal

# Ensure the extents and resolutions match across rasters
res(HR_Rast) == res(R_Rast)
res(R_Rast) == res(DRast)
res(DRast) == res(MaskRast)
ext(HR_Rast) == ext(R_Rast)
ext(R_Rast) == ext(DRast)
ext(DRast) == ext(MaskRast)

######################################
# STEP 4: Save finalized rasters
######################################

writeRaster(HR_Rast, paste0(fileSave, 'GlobalHeightRange.tif'))
writeRaster(R_Rast, paste0(fileSave, 'GlobalRugosity.tif'))
writeRaster(DRast, paste0(fileSave, 'GlobalFractalDimension.tif'))
writeRaster(MaskRast, paste0(fileSave, 'LandMask.tif'))

## ------------------------------------------------------------- End of Code ---------------------------------------