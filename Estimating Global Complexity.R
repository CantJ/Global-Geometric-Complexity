# This script is for estimating fractal dimension, surface rugosity, and height range from a global digital elevation model (Topography and Bathymetry combined), 
# at a resolution of 187*187m. 
# The outputs from this script are a large csv file (stored as a compressed parquet to allow for memory efficient storage) and
# rasters depicting global patterns in the complexity measures of height range, rugosity and fractal dimension.

# Date last modified: August 2023
# Primary Author: James Cant
# ------------------------------------------------------------------------------------------------------------

# Clear workspace
rm(list=ls(all=TRUE))

# Define directory pathways
filePath <- 'FileDirectory1/'
fileSave <- 'FileDirectory2/' # temporary file directory for saving .csv sub-files
FinalSave <- 'FileDirectory3/' # final output save location
# Define selected raster files
GlobalDEM <- 'GlobalDEMfile'# Global raster (187m resolution)
# This file is a global DEM comprising terrestrial topography and ocean bathymetry. 
# This file was produced by mosaicing tiles obtained from GEBCO.
# This DEM was then reprojected to an World Mollweide equal area projection.
MaskRast <- 'LandTemplatefile' # This raster will be used to differentiate between complexity values associated with terrestrial and marine environments (it has also been set to a resolution equal to L and reprojected to a Mollweide projection)
# This raster has been generated from the Land Cover 2020 data product from the Copernicus Climate Change service, by removing the coverage category associated with permenant water bodies.

# load necessary packages
library(tidyr)
library(sp)
library(raster)
library(parallel)
library(doParallel)
library(terra)
library(dplyr)
library(readr)
library(data.table)
library(arrow)
library(sf)

#################################################
# STEP 1: Define complexity function
#################################################
# *** This function has been adapted from Torres-Pulliza et al. (2020) Nature Ecol. & Evol 4. 1495-1501 ***

# Estimate complexity values for a selected Longitudinal Band
ExtractComplexity <- function(ii, DEMBand, L, L0, GridY, x0, progress = FALSE){
  
  # progress read out if requested (good for troubleshooting). Defaults to no progress display.
  if(progress == TRUE) { print(ii) }
  
  # Isolate selected pixel from longitudinal raster bands provided.
  DEMCell <- crop(DEMBand, extent(x0, x0+L, GridY[ii], GridY[ii]+L))
  
  # Extract mean cell elevation
  E <- cellStats(DEMCell, stat = 'mean', na.rm = TRUE)
  
  # Estimate Height range for selected DEM cell
  H <- suppressWarnings(max(DEMCell[!is.na(DEMCell)]) - min(DEMCell[!is.na(DEMCell)]))
  
  # Generate Spatial Grid data frame
  TmpDEMCell <- as(DEMCell, 'SpatialGridDataFrame')
  
  # Estimate Rugosity
  R <- surfaceArea(TmpDEMCell) / L^2
  
  # Estimate Fractal Dimension
  D <- suppressWarnings(3 - log10(H / (sqrt(2) * L0 * sqrt(R^2 - 1))) / log10(L / L0))
  
  # A little fix for grid cells that contain no height range
  # If the selected grid cell contains no height variation (i.e. H and H0 == 0) 
  # then it violates the theory of fractal dimension (because it is not fractal)
  # and thus cannot be processed.
  if(H == 0){
    R <- 1
    H <- 0
    D <- 2
  } # elevation can remain as calculated
  
  # A further little fix incase the selected Grid Cell contains NA's
  if(is.na(max(values(DEMCell))) & is.na(min(values(DEMCell)))){
    R <- NA
    H <- NA
    D <- NA
    E <- NA
  }
  
  # Return extracted complexity values
  return(list(Y = ii, Lat = GridY[ii], H = H, R = R, D = D, E = E))
}

#################################################
# STEP 2: Determine scope/extent, scales of variation and resolution parameters
#################################################

# load rasters
DEM <- rast(paste0(filePath, GlobalDEM))
LandMask <- rast(paste0(filePath, MaskRast))
# check projections
DEM 
LandMask # should all be Mollweide equal area
plot(DEM)
plot(LandMask)

# Store the crs projection
Molle <- crs(DEM)

# define L0 (minimum resolution)
L0 <- res(DEM)[1]

# define L
# Previous approaches have ensured the difference between L0 and L covers 2 orders of magnitude.
# However a sensitivity test performed previously indicates that estimates of complexity will remain consistent at ~0.8 orders of magnitude
# Therefore the scales selected will allow for calculating complexity at close to a 1km scale.
L <- L0*6

# define dimensions for sub-setting DEM in to manageable longitudinal bands equal to L in width to ease processing demands.
minX <- ext(DEM)[1]
maxX <- ext(DEM)[2] - L
minY <- ext(DEM)[3]
maxY <- ext(DEM)[4] - L

# define longitudinal band sequence
GridX <- seq(minX, maxX, by = L)
# and define latitudinal indexing sequence
# This index works by dividing each longitudinal band into latitudinal cells each L*L in dimension
# remembering that L is equal to 6 pixels.
GridY <- seq(minY, maxY, by = L)

#################################################
# STEP 3: Estimate complexity
#################################################

# Open multicore interface
nCores = detectCores()
cl <- makeCluster(nCores*0.5) # not to overload the computer system
registerDoParallel(cl)

# Run multicore processing
# loop through XY grid extracting complexity measures
initial_x <- length(list.files(fileSave, 'ComplexityValues'))+1 # allows for quick indexing following processing errors
# Run process across cluster
ComplexValues <- foreach(xx = initial_x:length(GridX), .inorder = FALSE) %dopar% { # cluster processing along longitudinal grid coordinates
  
  # load required packages across cluster
  library(sp)
  library(raster)
  library(terra)
  library(data.table)
  library(readr)
  
  # Allow all cluster nodes access to the necessary rasters (rasters cannot be called from the global environment)
  DEM <- rast(paste0(filePath, GlobalDEM)) # Digital Elevation Model
  LandMask <- rast(paste0(filePath, MaskRast)) # Land Mask raster
  
  # define origin coordinate for selected longitudinal band
  x0 <- GridX[xx]
  
  # Clip rasters to focus on selected longitudinal band
  tmpDEMRast <- raster(terra::crop(DEM, ext(x0, x0+L, minY, maxY+L)))
  tmpMask <- terra::crop(LandMask, ext(x0, x0+L, minY, maxY+L))
  
  # Extract complexity values
  CellComplex <- as.data.frame(t(sapply(1:length(GridY), ExtractComplexity, DEMBand = tmpDEMRast, L = L, L0 = L0, GridY = GridY, x0 = x0)))
  # transform out of list type
  CellComplex[c('Y','Lat','H','R','D','E')] <- as.numeric(unlist(CellComplex[c('Y','Lat','H','R','D','E')]))
  
  # Identify the realm associated with each pixel
  Realm <- numeric(dim(CellComplex)[1])
  Realm[!(is.na(rev(values(tmpMask))))] <- 1 # if pixel is terrestrial, value = 1
  
  # Attach Longitudinal details and realm identifier
  XValues <- data.frame(X = rep(xx, dim(CellComplex)[1]),
                        Lon = rep(x0, dim(CellComplex)[1]),
                        Terrestrial = Realm) # Marine = 0 & Terrestrial = 1
  CellComplex <- cbind(XValues, CellComplex)
  
  # Reorder columns to make viewing easier
  CellComplex <- CellComplex[, c('X','Y','Lon','Lat','Terrestrial','E','H','R','D')]
  
  # and save outputs (in-case of system meltdowns)
  fwrite(CellComplex, paste0(fileSave, 'ComplexityValues_X', xx, '.csv'), row.names = FALSE) # This will be one of 32157 csv subfiles.
  
  # return cluster output
  CellComplex
  # end of loop
}

# Close multicore cluster
stopCluster(cl)

# Save extracted complexities 
# In-case of loop crashes and the .csv file subsets need to be called.
# Compile the .csv sub file together
ComplexValues <- as.data.frame(list.files(fileSave, full.names = TRUE) %>% 
                                 lapply(read_csv, show_col_types = FALSE) %>% 
                                 bind_rows)
# and write to file (Checkpoint)
# using the arrow package this large csv file is saved in a memory efficient format
write_parquet(ComplexValues, paste0(filePath,'ComplexityData.parquet'))

#################################################
# STEP 4: Clean Complexity estimates
#################################################

# Reopen dataset using memory efficient format
ComplexValues <- read_parquet(paste0(filePath,'ComplexityData.parquet'), as_data_frame = F)

# How many NAs (due to mollweide projection)
ComplexValues %>%
  filter(is.na(H)) %>% collect() %>% dim() # 517052403 pixels of which 111057113 are NA

# check for inappropriate entries
ComplexValues %>%
  filter(is.infinite(H)) %>% collect() %>% dim()
ComplexValues %>%
  filter(is.infinite(R)) %>% collect() %>% dim() # no inf entries in rugosity & Height range
ComplexValues %>%
  filter(is.infinite(D)) %>% collect() %>% dim() # 6515 inf entries in Fractal dimension
ComplexValues %>%
  filter(!is.infinite(D) & D < 2 | D > 3) %>% collect() %>% dim() # a further 9827 Fractal dimension values are outside the range of 2 to 3
ComplexValues %>%
  filter(!is.infinite(D) & R < 1) %>% collect() %>% dim() # No rugosity values are lower than 1.
ComplexValues %>%
  filter(!is.infinite(D) & H < 0) %>% collect() %>% dim() # Height range values are all positive.
# carry out initial data cleaning
ComplexValues |>
  mutate(H = ifelse(is.infinite(D), NA, H),
         R = ifelse(is.infinite(D), NA, R),
         D = ifelse(is.infinite(D), NA, D)) |>
  mutate(D = ifelse(D<2 | D>3, NA, D)) |>
  collect() |>
  write_parquet(paste0(filePath, 'ComplexityData.parquet')) # overwrite output to memory efficient object

# Reload updated data file
ComplexValues <- read_parquet(paste0(filePath,'ComplexityData.parquet'), as_data_frame = F)

# Confirm that for all instances where D = 2, then R = 1 and H = 0 and visa-versa.
ComplexValues %>%
  filter(D == 2) %>% count(R, sort = TRUE) %>% collect() # When D = 2 then R = 1 in all cases.
ComplexValues %>%
  filter(D == 2) %>% count(H, sort = TRUE) %>% collect() # When D = 2 then H = 0 in all cases.
ComplexValues %>%
  filter(H == 0) %>% count(R, sort = TRUE) %>% collect() # When H = 0 then R = 1 in all cases.
ComplexValues %>%
  filter(H == 0) %>% count(D, sort = TRUE) %>% collect() # When H = 0 then D = 2 in all cases.
ComplexValues %>%
  filter(R == 1) %>% count(D, sort = TRUE) %>% collect()
ComplexValues %>%
  filter(R == 1) %>% count(H, sort = TRUE) %>% collect() # However there are 11364 instances where R = 1 but neither D = 2 nor H = 0.
# These are likely rounding errors and need to be corrected.
# Isolate offending cells. 
RugosityErrors <- ComplexValues %>%
  filter(R == 1 & H!=0) %>% 
  collect %>%
  tibble()
# Allow R to cope with more digits to prevent values very close to 1 being rounded down.
options(digits = 22) 
# Run correction loop
for (ii in 1:dim(RugosityErrors)[1]) {
  # progress read out
  print(ii)
  # Extract grid cell corresponding with error
  tmpRast <- crop(DEM, ext(RugosityErrors$Lon[ii], RugosityErrors$Lon[ii]+L, RugosityErrors$Lat[ii], RugosityErrors$Lat[ii]+L))
  # calculate height range
  H <- max(tmpRast[!is.na(tmpRast)]) - min(tmpRast[!is.na(tmpRast)])
  RugosityErrors$H2[ii] <- H
  # Calculate rugosity
  TmpGridCell <- as(raster(tmpRast), 'SpatialGridDataFrame')
  R <- surfaceArea(TmpGridCell) / L^2
  RugosityErrors$R2[ii] <- R
  # and calculate fractal dimension
  RugosityErrors$D2[ii] <- suppressWarnings(3 - log10(H / (sqrt(2) * L0 * sqrt(R^2 - 1))) / log10(L / L0))
}
# Drop repeated variables that aren't needed for merging dataframes
RugosityErrors <- RugosityErrors[,c('X','Y','Lon','Lat','H2','R2','D2')]
# Insert the corrected values back into the original data file
ComplexValues %>% 
  left_join(RugosityErrors,  
            by = c("X" = 'X', "Y" = 'Y', "Lon" = 'Lon', "Lat" = 'Lat')) %>%
  mutate(R = ifelse(!is.na(R2), R2, R)) %>% 
  select(Lon, Lat, Terrestrial, E, H, R, D) %>%
  write_parquet(paste0(FinalSave, 'ComplexityData.parquet')) # Save final data file output

# Reload updated data file to confirm successful data cleaning.
ComplexValues <- read_parquet(paste0(FinalSave,'ComplexityData.parquet'), as_data_frame = F)

# Check all rugosity values of 1 correspond with fractal dimension values of 2
ComplexValues %>%
  filter(R == 1) %>% count(D, sort = TRUE) %>% collect() # all good to proceed.

#################################################
# STEP 5: Generate Initial Global Rasters
#################################################

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

# view to confirm appropriate mollweide projection
plot(HRast)
plot(RRast)
plot(DRast)
plot(MaskRast)

#################################################
# STEP 6: Final Quality Control
#################################################

# Evaluate data ranges
# Fractal Dimension
range(values(DRast), na.rm = TRUE) # should be between 2 and 3.
hist(values(DRast)) # skewed normal distribution - could be due to the mix of estimates from across both terrestrial and marine environments?

# Rugosity 
range(values(RRast), na.rm = TRUE) # should be greater than or equal to 1
hist(values(RRast)) # Clear Poisson distribution
# Following Torres-Pulliza et al (2020) Rugosity can be standardized to (R^2)-1
R_Rast <- app(RRast, fun = function(x) { (x^2)-1 })
hist(log10(values(R_Rast))) # transformed to normal

# Height Range
range(values(HRast), na.rm = TRUE) # must be non-negative
hist(values(HRast)) # Again, a clear poisson distribution
# Following Torres-Pulliza et al (2020) Height range can be standardized to HR/sqrt(2)*L0
HR_Rast <- app(HRast, fun = function(x) { x/(sqrt(2)*L0) })
hist(log10(values(HR_Rast))) # transformed to normal

# Ensure the extents and resolutions match expectations and across rasters
res(HR_Rast) == c(L,L)
res(R_Rast) == c(L,L)
res(DRast) == c(L,L)
res(MaskRast) == c(L,L)
ext(HR_Rast) == ext(R_Rast)
ext(R_Rast) == ext(DRast)
ext(DRast) == ext(MaskRast)

#################################################
# STEP 7: Save Finalized Raster Outputs
#################################################

writeRaster(HR_Rast, paste0(FinalSave, 'GlobalHeightRange.tif'))
writeRaster(R_Rast, paste0(FinalSave, 'GlobalRugosity.tif'))
writeRaster(DRast, paste0(FinalSave, 'GlobalFractalDimension.tif'))
writeRaster(MaskRast, paste0(FinalSave, 'LandMask.tif'))

## ------------------------------------------------------------- End of Code ---------------------------------------