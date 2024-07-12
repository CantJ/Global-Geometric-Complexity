# This script is for reprojecting the global climate rasters downloaded from worldClim before estimating a measure of monthly thermal variability 
# for the global terrestrial landscape.

# Primary Author: James Cant
# Date last modified: July 2024
# ---------------------------------------------------------------

# Clear working directory
rm(list = ls())

# Load required packages
library(terra)

# set up file directory pathways
LoadPath <- 'C:/Users/jic2/OneDrive - University of St Andrews/Documents/Complexity Analyses/WorldClim/Raw/'
RastPath <- 'C:/Users/jic2/OneDrive - University of St Andrews/Documents/Complexity Analyses/Final Data Products/'
SavePath <- 'C:/Users/jic2/OneDrive - University of St Andrews/Documents/Complexity Analyses/WorldClim/Reprojected/'

# Load Fractal Dimension raster as reprojection template.
DRast <- rast(paste0(RastPath, "GlobalFractalDimension.tif"))

# List raster files to be reprojected
climFiles <- list.files(LoadPath, full.names = FALSE)
# Generate vector of names for new rasters produced by reprojection
FileNames <- substr(climFiles, 12, 23) # extract variable, year and month.

# work through each raw file to reproject and save
for(ii in 1:length(climFiles)){
  # Progress read out
  print(ii)
  # open selected raster
  tmpRast <- rast(paste0(LoadPath, climFiles[ii]))
  # reproject using selected template
  tmpRast <- project(tmpRast, DRast, method = 'bilinear')
  # save new raster
  writeRaster(tmpRast, paste0(SavePath, FileNames[ii],'.tif'))
} # depending on the number of climate rasters downloaded this loop can take a minute to run.

# Match up tmax and tmin pairs to estimate mean monthly temperature patterns
RastMonths <- unique(substr(FileNames, 6,12))

for(ii in 1:length(RastMonths)){  
  # Progress read out
  print(ii)
  # Identify corresponding tmax and tmin rasters
  RastSelect <- list.files(SavePath, pattern = RastMonths[ii], full.names = T)
  tmpRast <- rast(RastSelect)
  # take average
  tmpRast <- mean(tmpRast)
  # save new raster
  writeRaster(tmpRast, paste0(SavePath, 'MeanTemp', RastMonths[ii], '.tif'))
}

# Compute the variability in monthly temperatures values across the global to generate a single climate variability raster

# ------------------------------------- End of Code ------------------------------------------