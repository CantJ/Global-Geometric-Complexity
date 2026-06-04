# This script is for re-projecting the global climate rasters downloaded from worldClim before estimating a measure of monthly thermal variability 
# for the global terrestrial landscape.

# Primary Author: James Cant
# Date last modified: July 2024
# ---------------------------------------------------------------

# Load Fractal Dimension raster as re-projection template.
DRast <- rast(paste0(FilePath, "GlobalFractalDimension.tif"))

# List raster files to be reprojected
climFiles <- list.files(ClimatePath, full.names = FALSE)
# Generate vector of names for new rasters produced by re-projection
FileNames <- substr(climFiles, 12, 23) # extract variable, year and month.

# work through each raw file to re-project and save
for(ii in 1:length(climFiles)){
  # Progress read out
  print(ii)
  # open selected raster
  tmpRast <- rast(paste0(ClimatePath, climFiles[ii]))
  # re-project using selected template
  tmpRast <- project(tmpRast, DRast, method = 'bilinear')
  # save new raster
  writeRaster(tmpRast, paste0(tmpSave, FileNames[ii],'.tif'))
} # depending on the number of climate rasters downloaded this loop can take a minute to run.

# Match up tmax and tmin pairs to estimate mean monthly temperature patterns
RastMonths <- unique(substr(FileNames, 6,12))

for(ii in 1:length(RastMonths)){  
  # Progress read out
  print(ii)
  # Identify corresponding tmax and tmin rasters
  RastSelect <- list.files(tmpSave, pattern = RastMonths[ii], full.names = T)
  tmpRast <- rast(RastSelect)
  # take average
  tmpRast <- mean(tmpRast)
  # save new raster
  writeRaster(tmpRast, paste0(tmpSave,'MeanTemp_', RastMonths[ii], '.tif'))
}

# Compute the variability in monthly temperatures values across the global to generate a single climate variability raster
# Match up monthly mean temperatures from across different years to compute within month temperature variability.
meanFiles <- list.files(tmpSave, pattern = '.tif')
MonthList <- unique(substr(meanFiles, 15, 16))

for(ii in 1:length(MonthList)){  
  # Progress read out
  print(ii)
  # Identify matching monthly rasters
  RastSelect <- list.files(tmpSave, pattern = paste0('-',MonthList[ii]), full.names = T)
  tmpRast <- sds(RastSelect)
  # calculate standard deviation
  tmpRast <- app(tmpRast, sd, na.rm = TRUE)
  # save new raster
  writeRaster(tmpRast, paste0(tmpSave, 'MonthlySD_', MonthList[ii], '.tif'))
}

#Finally, estimate the mean monthly temperature variability across the period 2010-2021.
MMVRast <- sds(list.files(tmpSave, pattern = 'MonthlySD_', full.names = T)) 
MMVRast <- app(MMVRast, mean, na.rm = TRUE)

# Save Temperature variability raster
writeRaster(MMVRast, paste0(FilePath, 'ClimVar_Mollweide_1870m.tif'))

# ------------------------------------- End of Code ------------------------------------------