# This script is for estimating fractal dimension, surface rugosity, and height range from a global digital elevation model (Topography and Bathymetry combined), 
# at a resolution of 248*248m. 
# The outcome of this script is a large csv file (stored as a compressed parquet to allow for memory efficient storage).

# Date last modified: July 2023
# Primary Author: James Cant
# ------------------------------------------------------------------------------------------------------------

# Clear workspace
rm(list=ls(all=TRUE))

# Define directory pathways
filePath <- 'C:/Users/jic2/OneDrive - University of St Andrews/Documents/Complexity Analyses/'
fileSave <- 'C:/Users/jic2/OneDrive - University of St Andrews/Documents/Complexity Analyses/Complexity Files/' # temporary file directory for saving .csv sub-files
# Define selected raster file
SelectDEM <- 'GlobalDEM_Mollweide_187m.tif'# Global raster (187m resolution)
# This file is a global DEM comprising terrestrial topography and ocean bathymetry. 
# This file was produced by mosaicing tiles obtained from GEBCO.
# This DEM was then reprojected to an World Mollweide equal area projection.
LandMask <- 'LandMask_mollweide.tif' # This is a raster template denoting pixels associated with terrestrial and marine environments.

# load necessary packages
library(tidyr)
library(rgdal)
library(sp)
library(raster)
library(parallel)
library(doParallel)
library(terra)
library(dplyr)
library(readr)
library(data.table)
library(arrow)


#################################################
# STEP 1: Define complexity function
#################################################
# *** This function has been adapted from Torres-Pulliza et al. (2020) Nature Ecol. & Evol 4. 1495-1501 ***

# Estimate complexity values for a selected Longitudinal Band
ExtractComplexity <- function(ii, DEMBand, L, L0, GridY, x0, progress = FALSE){
  
  # progress read out if requested (good for troubleshooting). Defaults to no progress display.
  if(progress == TRUE) { print(ii) }
  
  # Isolate Grid Cell from Longitudinal DEM band provided.
  GridCell <- crop(DEMBand, extent(x0, x0+L, GridY[ii], GridY[ii]+L))
  
  # Estimate Height range for selected DEM cell
  H <- suppressWarnings(max(GridCell[!is.na(GridCell)]) - min(GridCell[!is.na(GridCell)]))
  
  # Generate Spatial Grid data frame
  TmpGridCell <- as(GridCell, 'SpatialGridDataFrame')
  
  # Estimate Rugosity
  R <- surfaceArea(TmpGridCell) / L^2
  
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
  }
  
  # A further little fix incase the selected Grid Cell contains NA's
  if(is.na(max(values(GridCell))) & is.na(min(values(GridCell)))){
    R <- NA
    H <- NA
    D <- NA
  }
  
  # Return extracted complexity values
  return(list(Y = ii, Lat = GridY[ii], H = H, R = R, D = D))
}

#################################################
# STEP 2: Determine scope/extent, scales of variation and resolution parameters
#################################################

# load DEM raster
GlobalDEM <- rast(paste0(filePath, SelectDEM))
# check projection
GlobalDEM # should be Mollweide equal area
plot(GlobalDEM)

# Store the crs projection
Molle <- crs(GlobalDEM)

# define L0 (minimum resolution)
L0 <- res(GlobalDEM)[1]

# define L
# Previous approaches have ensured the difference between L0 and L covers 2 orders of magnitude.
# However a sensitivity test performed previously indicates that estimates of complexity will remain consistent at ~0.7 orders of magnitude
# Therefore the scales selected will allow for calculating complexity at close to a 1km scale.
L <- L0*6

# define dimensions for sub-setting DEM in to manageable longitudinal bands equal to L in width to ease processing demands.
minX <- ext(GlobalDEM)[1]
maxX <- ext(GlobalDEM)[2] - L
minY <- ext(GlobalDEM)[3]
maxY <- ext(GlobalDEM)[4] - L

# define longitudinal band sequence
GridX <- seq(minX, maxX, by = L)
# and define latitudinal indexing sequence
# This index works by dividing each longitudinal band into latitudinal cells each L*L in dimension
# remembering that L is equal to 5 pixels.
GridY <- seq(minY, maxY, by = L)

#################################################
# STEP 3: Estimate complexity
#################################################

# Open multicore interface
nCores = detectCores()
cl <- makeCluster(nCores*0.7) # not to overload the computer system
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
  DEM <- rast(paste0(filePath, SelectDEM))  
  MaskRast <- rast(paste0(filePath, LandMask)) # Load Mask raster
  
  # define origin coordinate for selected longitudinal band
  x0 <- GridX[xx]
  
  # Clip rasters to focus on selected longitudinal band
  tmpDEMRast <- raster(terra::crop(DEM, ext(x0, x0+L, minY, maxY+L)))
  tmpMask <- terra::crop(MaskRast, ext(x0, x0+L, minY, maxY+L))
  
  # Extract complexity values
  CellComplex <- as.data.frame(t(sapply(1:length(GridY), ExtractComplexity, DEMBand = tmpDEMRast, L = L, L0 = L0, GridY = GridY, x0 = x0)))
  # transform out of list type
  CellComplex[c('Y','Lat','H','R','D')] <- as.numeric(unlist(CellComplex[c('Y','Lat','H','R','D')]))
  
  # Identify the realm associated with each pixel
  Realm <- numeric(dim(CellComplex)[1])
  Realm[values(tmpMask) != 0] <- 1 # if pixel is terrestrial, value = 1
  # remove final Realm reading in order to match up vector sizes
  Realm <- Realm[1:dim(CellComplex)[1]]
  
  # Attach Longitudinal details and realm identifier
  XValues <- data.frame(X = rep(xx, dim(CellComplex)[1]),
                        Lon = rep(x0, dim(CellComplex)[1]),
                        Terrestrial = Realm) # Marine = 0 & Terrestrial = 1
  CellComplex <- cbind(XValues, CellComplex)
  
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
write_parquet(ComplexValues, paste0(filePath,'GlobalComplexity.parquet'))

#################################################
# STEP 4: Clean Complexity estimates
#################################################

# Reopen dataset using memory efficient format
ComplexValues <- read_parquet(paste0(filePath,'GlobalComplexity.parquet'), as_data_frame = F)

# check for inappropriate entries
ComplexValues %>%
  filter(is.infinite(H)) %>% collect() %>% dim()
ComplexValues %>%
  filter(is.infinite(R)) %>% collect() %>% dim() # no inf entries in rugosity & Height range
ComplexValues %>%
  filter(is.infinite(D)) %>% collect() %>% dim() # 7531 inf entries in Fractal dimension
ComplexValues %>%
  filter(!is.infinite(D) & D < 2 | D > 3) %>% collect() %>% dim() # a further 493 Fractal dimension values are outside the range of 2 to 3
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
  write_parquet(paste0(filePath, 'GlobalComplexity.parquet')) # overwrite output to memory efficient object

# Reload updated data file
ComplexValues <- read_parquet(paste0(filePath,'GlobalComplexity.parquet'), as_data_frame = F)

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
  filter(R == 1) %>% count(H, sort = TRUE) %>% collect() # However there are 14702 instances where R = 1 but neither D = 2 nor H = 0.
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
  tmpRast <- crop(GlobalDEM, ext(RugosityErrors$Lon[ii], RugosityErrors$Lon[ii]+L, RugosityErrors$Lat[ii], RugosityErrors$Lat[ii]+L))
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
RugosityErrors <- RugosityErrors[,-(c('H','D','R','Terrestrial'))]
# Insert the corrected values back into the original data file
ComplexValues %>% 
  left_join(RugosityErrors,  
            by = c("X" = 'X', "Y" = 'Y', "Lon" = 'Lon', "Lat" = 'Lat')) %>%
  mutate(R = ifelse(!is.na(R2), R2, R)) %>% 
  select(X, Y, Lon, Lat,Terrestrial, H, R, D) %>%
  write_parquet(paste0(filePath, 'GlobalComplexity.parquet')) # overwrite output to memory efficient object
# This is the final output of this script.

# Reload updated data file to confirm successful data cleaning.
ComplexValues <- read_parquet(paste0(filePath,'GlobalComplexity.parquet'), as_data_frame = T)

ComplexValues %>%
  filter(R == 1) %>% count(D, sort = TRUE) # all good to proceed.

# Save Data as a csv file.
write_rds(ComplexValues, paste0(filePath, 'GlobalComplexity.csv'))

##### ------------------------------------------------------------- End of Code --------------------------------##### ------------------------------------------------------------- End of Code --------------------------------