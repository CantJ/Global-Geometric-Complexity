# This script is for estimating fractal dimension, surface rugosity, and height range from a global digital elevation model (Topography and Bathymetry combined), 
# at a resolution of 187*187m. 
# The outputs from this script are a large csv file (stored as a compressed parquet to allow for memory efficient storage) and
# rasters depicting global patterns in the complexity measures of height range, rugosity and fractal dimension.

# Date last modified: May 2024
# Primary Author: James Cant
# ------------------------------------------------------------------------------------------------------------
# A word of warning. If running this script to completion for the first time it can take a very long time to complete.
# It is advised to use parallel processing and HPC computing were possible.

# Clear workspace
rm(list=ls(all=TRUE))

# Define directory pathways
filePath <- '/File_Directory_1/'
TilePath <- '/File_Directory_2/' # directory for storing DEM tile files
fileSave <- '/File_Directory_3/' # temporary file directory for saving .csv sub-files
FinalSave <- '/File_Directory_4/' # final output save location
# Define selected raster files
GlobalDEM <- 'GlobalDEM_Mollweide_187m.tif'# Global raster (187m resolution)
# This file is a global DEM comprising terrestrial topography and ocean bathymetry. 
# This file was produced by mosaicing tiles obtained from GEBCO.
# This DEM was then reprojected to an World Mollweide equal area projection.
LandMask <- 'LandCover2022_Mollweide_1870m.tif' # this is a raster created using the 2022 CCCI land cover survey product.

# Is this the first time this code is being implemented? I.e. is DEM tiling required
FirstRun <- FALSE

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
library(habtools)
library(xpectr)
library(gtools)

#################################################
# STEP 1: Define complexity extraction function
#################################################
# *** This function has been adapted from Torres-Pulliza et al. (2020) Nature Ecol. & Evol 4. 1495-1501 and the habtools package***

# Estimate complexity values for a selected Longitudinal Band
ExtractComplexity <- function(ii, DEMCells, L, L0, progress = FALSE){
  
  # progress read out if requested (good for troubleshooting). Defaults to no progress display.
  if(progress == TRUE) { print(ii) }
  
  # Isolate selected pixel from longitudinal raster bands provided.
  DEMCell <- raster(DEMCells[ii])
  
  # No need to process the selected Grid Cell if it only contains NA's
  if(is.na(max(values(DEMCell))) & is.na(min(values(DEMCell)))){
    R <- NA
    H <- NA
    D <- NA
    E <- NA
  } else {
    # Extract mean cell elevation
    E <- cellStats(DEMCell, stat = 'mean', na.rm = TRUE)
    
    # Estimate Height range for selected DEM cell
    H <- suppress_mw(hr(DEMCell))
    
    # If the selected grid cell contains no height variation (i.e. H and H0 == 0) 
    # then it violates the theory of fractal dimension (because it is not fractal)
    # and thus cannot be processed.
    if(H == 0){
      R <- 1
      D <- 2
    } else { # otherwise rugosity and fractal dimension can be estimated
      # Rugosity
      R <- suppress_mw(rg(DEMCell, method = 'area', L0 = L0))
      # Fractal Dimension
      D <- suppress_mw(rdh_theory(R = R, H = H, L = L, L0 = L0)$D)
    }
  }
  
  # Return extracted complexity values
  return(list(Y = ii, Lat = DEMCell@extent[4], H = H, R = R, D = D, E = E))
}

#################################################
# STEP 2: Determine scope/extent, scales of variation and resolution parameters
#################################################

# load rasters
DEM <- rast(paste0(filePath, GlobalDEM))
maskRast <- rast(paste0(filePath, LandMask))

# check projections
DEM 
maskRast # should all be Mollweide equal area

# plot rasters
plot(DEM)
plot(maskRast)

# Store the crs projection
Molle <- crs(DEM)

# define L0 (minimum resolution)
L0 <- res(DEM)[1]

# define L
# Previous approaches have ensured the difference between L0 and L covers 2 orders of magnitude.
# However a sensitivity test performed previously indicates that estimates of complexity will remain consistent at ~0.8 orders of magnitude
# Therefore the scales selected will allow for calculating complexity at close to a 1.5km scale (One order of magnitude).
L <- L0*10

# define dimensions for sub-setting DEM in to manageable longitudinal and latitudinal bands equal to L in width 
# This will help to set up a pixel grid of L*L pixels.
minX <- ext(DEM)[1]
maxX <- ext(DEM)[2] - L
minY <- ext(DEM)[3]
maxY <- ext(DEM)[4] - L

# define longitudinal band sequence
GridX <- seq(minX, maxX, by = L)
# and define latitudinal indexing sequence
# This index works by dividing each longitudinal band into latitudinal cells each L*L in dimension (this is to determine the possible number of L*L pixels that can be generated)
# remembering that L is equal to 10 pixels.
GridY <- seq(minY, maxY, by = L)

# Set up tiling sequence for parallelisation
if(FirstRun == TRUE) { setwd(TilePath) # set working directory for temporary saving of raster tiles
  # Implement tiling
  DEMTiles <- terra::makeTiles(x = DEM, y = c(nrow(DEM), L/L0), filename = 'Tile_.tif', overwrite = FALSE)
  gc() # clear available memory
} else {
  DEMTiles <- mixedsort(list.files(TilePath))
}
# N.B: This tiling is a very slow process.
# Therefore the tiling process only needs to be called if this is the first run.

#################################################
# STEP 3: Estimate complexity
#################################################

# Open multicore interface
nCores = detectCores()
cl <- makeCluster(floor(nCores*0.7)) # not to overload the computer system
registerDoParallel(cl)

# Run multicore processing
# loop through XY grid extracting complexity measures
initial_x <- length(list.files(fileSave, 'ComplexityValues'))+1 # allows for quick indexing following processing errors
# Run process across cluster 
ComplexValues <- foreach(xx = initial_x:length(GridX), .inorder = FALSE) %dopar% { # cluster processing along longitudinal DEM tiles (focusing on only tiles with a width equal to L)
  
  # load required packages across cluster
  library(sp)
  library(raster)
  library(terra)
  library(data.table)
  library(readr)
  library(habtools)
  library(xpectr)
  
  # Allow all cluster nodes access to the geodiversity and land Mask rasters (rasters cannot be called from the global environment)
  GDRast <- rast(paste0(filePath, GeoRast)) # Geodiversity
  maskRast <- rast(paste0(filePath, LandMask)) # Land Cover
  
  # define origin coordinate for selected longitudinal band
  x0 <- GridX[xx]
  
  # Clip rasters to focus on selected longitudinal band
  tmpGDRast <- terra::crop(GDRast, ext(x0, x0+L, minY, maxY+L))
  tmpMask <- terra::crop(maskRast, ext(x0, x0+L, minY, maxY+L))
  
  # Load DEM tile
  tmpDEMRast <- rast(paste0(TilePath, DEMTiles[xx]))
  
  # Break selected DEM tile into temporary cells
  filename <- paste0(tempfile(), "_.tif")
  TempCells <- terra::makeTiles(tmpDEMRast, y = c(L/L0, L/L0), filename, overwrite = TRUE, extend = T)
  
  # Extract complexity values
  CellComplex <- as.data.frame(t(sapply(1:length(GridY), ExtractComplexity, DEMCells = TempCells, L = L, L0 = L0)))
  # transform out of list type
  CellComplex[c('Y','Lat','H','R','D','E')] <- as.numeric(unlist(CellComplex[c('Y','Lat','H','R','D','E')]))
  
  # Extract the Land Cover values associated with each pixel
  LC <- as.numeric(values(tmpMask))
  
  # Identify the realm associated with each pixel
  Realm <- numeric(dim(CellComplex)[1])
  Realm[!(is.na(LC))] <- 1 # if pixel is terrestrial, value = 1
  
  # Attach Longitudinal details, geodiversity estimates, and realm identifier
  XValues <- data.frame(X = rep(xx, dim(CellComplex)[1]),
                        Lon = rep(x0, dim(CellComplex)[1]),
                        Realm = Realm) # Marine = 0 & Terrestrial = 1
  CellComplex <- cbind(XValues, CellComplex)
  
  # Reorder columns to make viewing easier
  CellComplex <- CellComplex[, c('X','Y','Lon','Lat','Realm','E','H','R','D')]
  
  # and save outputs (in-case of system meltdowns)
  fwrite(CellComplex, paste0(fileSave, 'ComplexityValues_X', xx, '.csv'), row.names = FALSE) # The number of resultant csv files will = length(GridX)
  
  # return cluster output
  CellComplex
  # end of loop
}

# Close multicore cluster
stopCluster(cl)

# Save extracted complexities 
# In-case of loop crashes and the .csv file subsets need to be called.
# Compile the .csv sub files together
ComplexValues <- as.data.frame(list.files(fileSave, full.names = TRUE) %>% 
                                 lapply(read_csv, show_col_types = FALSE) %>% 
                                 bind_rows)
# and write to file (Checkpoint)
# using the arrow package this large csv file is saved in a memory efficient format
write_parquet(ComplexValues, paste0(FinalSave,'GlobalComplexity.parquet'))

#################################################
# STEP 4: Clean Complexity estimates
#################################################

# Reopen dataset using memory efficient format
ComplexValues <- read_parquet(paste0(FinalSave,'GlobalComplexity.parquet'), as_data_frame = F)

# How many NAs (due to mollweide projection)
ComplexValues %>%
  filter(is.na(H)) %>% collect() %>% dim() # 186,129,218 pixels of which 39,976,460 are NA

# check for inappropriate entries
ComplexValues %>%
  filter(is.infinite(H)) %>% collect() %>% dim()
ComplexValues %>%
  filter(is.infinite(R)) %>% collect() %>% dim() # no inf entries in rugosity & Height range
ComplexValues %>%
  filter(is.infinite(D)) %>% collect() %>% dim() # 17 inf entries in Fractal dimension
ComplexValues %>%
  filter(!is.infinite(D) & D < 2 | D > 3) %>% collect() %>% dim() # a further 85046 Fractal dimension values are outside the range of 2 to 3
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
  write_parquet(paste0(FinalSave, 'GlobalComplexity.parquet')) # overwrite output to memory efficient object

# Reload updated data file
ComplexValues <- read_parquet(paste0(FinalSave,'GlobalComplexity.parquet'), as_data_frame = F)

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
  filter(R == 1) %>% count(H, sort = TRUE) %>% collect() # However there are 137 instances where R = 1 but neither D = 2 nor H = 0.
# These are likely rounding errors and need to be omitted
ComplexValues |>
  mutate(H = ifelse(R == 1 & H != 0, NA, H),
         D = ifelse(R == 1 & D != 2, NA, D),
         R = ifelse(is.na(D) & is.na(H), NA, R)) |>
  collect() |>
  write_parquet(paste0(FinalSave, 'GlobalComplexity.parquet'))

# Reload updated data file to confirm successful data cleaning.
ComplexValues <- read_parquet(paste0(FinalSave,'GlobalComplexity.parquet'), as_data_frame = F)

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
LandMask <- suppressWarnings(SpatialPixelsDataFrame(coords, data = data.frame(as.numeric(ComplexValues[['Realm']])), proj4string = CRS(Molle)))
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
plot(HRRast)
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
RRast <- app(RRast, fun = function(x) { (x^2)-1 })
hist(log10(values(RRast))) # transformed to normal

# Height Range
range(values(HRast), na.rm = TRUE) # must be non-negative
hist(values(HRast)) # Again, a clear poisson distribution
# Following Torres-Pulliza et al (2020) Height range can be standardized to HR/sqrt(2)*L0
HRRast <- app(HRast, fun = function(x) { x/(sqrt(2)*L0) })
hist(log10(values(HRRast))) # transformed to normal

# Ensure the extents and resolutions match expectations and across rasters
res(HRRast) == c(L,L)
res(RRast) == c(L,L)
res(DRast) == c(L,L)
res(MaskRast) == c(L,L)
ext(HRRast) == ext(RRast)
ext(RRast) == ext(DRast)
ext(DRast) == ext(MaskRast)

#################################################
# STEP 7: Save Finalized Raster Outputs
#################################################

writeRaster(HRRast, paste0(FinalSave, 'GlobalHeightRange.tif'))
writeRaster(RRast, paste0(FinalSave, 'GlobalRugosity.tif'))
writeRaster(DRast, paste0(FinalSave, 'GlobalFractalDimension.tif'))
writeRaster(MaskRast, paste0(FinalSave, 'LandMask.tif'))

## ------------------------------------------------------------- End of Code ---------------------------------------