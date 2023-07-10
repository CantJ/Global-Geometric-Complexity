# This script is for testing the accuracy of increasing map resolution on estimates of fractal dimension and rugosity
# This script will also compare complexity estimate calculated using the planar equation and the surface area function
#
# Main author: James Cant (jic2@st-andrews.ac.uk)
# Date Last Modified: July 2023
# -----------------------------------------------------------

# Clear workspace
rm(list=ls(all=TRUE))

# Define directory pathways
filePath <- '~/James/Raw geodiversity data/'
# Define raster files
GlobalDEM <- 'GlobalDEM_Mollweide_187m.tif'# Global raster (187m resolution)
# This file is a global DEM comprising terrestrial topography and ocean bathymetry. 
# This file was produced by mosaicing tiles obtained from GEBCO.
# This DEM was then reprojected to an World Mollweide equal area projection.


# load necessary packages
library(rgdal)
library(sp)
library(raster)
library(terra)
library(dplyr)
library(SpaDES)
library(parallel)
library(data.table)
library(ggplot2)
library(viridis)

#################################################
# STEP 1: Define complexity functions
#################################################
# *** These functions have been adapted from Torres-Pulliza et al. (2020) Nature Ecol. & Evol 4. 1495-1501 ***
# Importantly these different functions will allow for testing the association between Rugosity and fractal dimension 
# and the accuracy of different approaches used to calculate them.

# Rugosity function
Rugosity <- function(H0, L0) {
  R <- sqrt((H0^2) / (2 * L0^2) + 1)
  # output
  return(R)
}

# Height Range function for calculating height variation across different spatial windows
estimateHR <- function(x, y, s, Rast, x0, y0) {
  bx <- extent(cbind(c(x0 + x, y0 + y), c(x0 + x + s, y0 + y + s)))
  # output
  return(diff(range(getValues(crop(Rast, bx)), na.rm = TRUE)))
}

# Convert DEM into moving height range windows
height_variation <- function(Rast, scl, L, y0, x0) {
  # define storage output
  temp <- data.frame()
  # loop through each scale
  for (s in scl) {
    # define scale grid
    inc <- seq(0, L - s, s)
    x <- rep(inc, L/s)
    y <- rep(inc, each = L/s)
    # estimate height range for each grid cell
    temp <- rbind(temp, data.frame(L0 = s, x = x, y = y, H0 = mcmapply(estimateHR, x = x, y = y, s = s, MoreArgs = list(Rast = Rast, x0 = x0, y0 = y0))))
  }
  # return output
  return(temp)
}

# Estimate complexity values for the selected grid cell
# This cannot be used for the Terrestrial data set due to the presence of 0 values along the margins of water bodies.
# As zero values cannot be log transformed without an adjustment factor. However introducing an adjustment factor influences estimates of Fractal Dimension.
ComplexityValues <- function(data, L, L0) {
  # log10 transform
  data$H0 <- log10(data$H0)
  data$L0 <- log10(data$L0)
  # A little fix for Height ranges of zero (if there are any)
  if(dim(data[is.infinite(data$H0),])[1] > 0){ data[is.infinite(data$H0),]$H0 <- 0 }
  
  # Estimate mean height range across each spatial scale
  data_m <- aggregate(H0 ~ L0, data, mean)
  
  # Extract height range at both L and L0
  H <- (10^data_m$H0[data_m$L0 == log10(L)])
  H0 <- (10^data_m$H0[data_m$L0 == log10(L0)])
  
  # Calculate slope coefficient (S)
  mod <- lm(H0 ~ L0, data_m)
  D <- 3 - coef(mod)[2] # D = 3 - S
  
  # Calculate rugosity
  R <- Rugosity(H0, L0)		
  
  # A little fix for grid cells that contain no height range
  # If the selected grid cell contains no height variation (i.e. H and H0 == 0) 
  # then it violates the theory of fractal dimension (because it is not fractal)
  # and thus cannot be processed.
  if(H == 0 & H0 ==0){
    R <- 1
    H <- 0
    D <- 2
  }
  
  # outputs
  return(list(D = D, R = R, H = H))
}

# Alternative method (Combining all complexity measures)
# This function is for when log transformation is not possible (due to the presence of zeros)
ExtractComplexity <- function(GridCell, L, L0){
  # Estimate Height range for selected grid cell
  H <- max(GridCell[!is.na(GridCell)]) - min(GridCell[!is.na(GridCell)])
  
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
  
  return(list(H = H, R = R, D = D))
}

#################################################
# STEP 2: Implement the effect of changing resolution on complexity estimates
#################################################

# Define function to test the varying effects of changing L.
# test both informs the function as to whether both the planar equation and surfaceArea function methods are to be compared.
VaryingL <- function(ii, test.both = TRUE){
  
  # set step counter for progress monitoring
  stepCount <- 1
  
  # Define scale over which to calculate height variation
  scl <- seq(from = L0, to = L[ii], by = L0)
  
  # define upper grid dimensions (for sub-setting DEM in to manageable chunks equal to L*L in dimension for processing)
  maxX <- Xlimit - L[ii]
  maxY <- Ylimit - L[ii]
  # and define Terrestrial and Marine grids
  GridX <- seq(minX, maxX, by = L[ii])
  GridY <- seq(minY, maxY, by = L[ii])
  
  # define storage
  PixelValues <- data.frame()
  
  # Run for loop
  for(xx in 1:length(GridX)){
    for(yy in 1:length(GridY)) { # work through each grid cell to estimate its corresponding measures of complexity.
      # update step count
      cat(ii, '.', stepCount, '\n')
      stepCount <- stepCount+1
      
      # define origin coordinates for selected grid cell
      x0 <- GridX[xx]
      y0 <- GridY[yy]
      
      # Clip raster to focus on selected grid cell (L*L)
      tmpRast <- raster(terra::crop(DEM, ext(x0, x0+L[ii], y0, y0+L[ii])))
      
      if(test.both == TRUE) {
        # For each grid cell calculate complexity using the two approaches
        # Approach 1: Calculate Height range and rugosity using the plane-equations.
        # This portion of the code is the processing time bottleneck.
        cellHR <- height_variation(Rast = tmpRast, scl = scl, L = L[ii], x0 = x0, y0 = y0)
        
        # Extract complexity values
        CellComplex1 <- ComplexityValues(data = cellHR, L = L[ii], L0 = L0)
        # Approach 2: Calculate Height range and rugosity using inbuilt R functions.
        CellComplex2 <- ExtractComplexity(GridCell = tmpRast, L = L[ii], L0 = L0)
        
        # Save outputs
        Outputs <- c(L[ii], # scale
                     xx, # Grid index values
                     yy, 
                     x0, # Longitude
                     y0, # Latitude
                     as.numeric(CellComplex1$H),
                     as.numeric(CellComplex1$R),
                     as.numeric(CellComplex1$D),
                     as.numeric(CellComplex2$H),
                     as.numeric(CellComplex2$R),
                     as.numeric(CellComplex2$D))
        
        # And store
        PixelValues <- rbind(PixelValues, Outputs)
        colnames(PixelValues) <- c('L','X','Y','Lon','Lat','H_theory','R_theory','D_theory','H','R','D')
      } else {
        # For each grid cell calculate complexity using only in the inbuilt R functions.
        CellComplex1 <- ExtractComplexity(GridCell = tmpRast, L = L[ii], L0 = L0)
        
        # Save outputs
        Outputs <- c(L[ii], # scale
                     xx, # Grid index values
                     yy, 
                     x0, # Longitude
                     y0, # Latitude
                     as.numeric(CellComplex1$H),
                     as.numeric(CellComplex1$R),
                     as.numeric(CellComplex1$D))
        
        # And store
        PixelValues <- rbind(PixelValues, Outputs)
        colnames(PixelValues) <- c('L','X','Y','Lon','Lat','H','R','D')
      }
    }      
  }    
  
  # and return
  return(PixelValues)
}

# load Global DEM raster
DEM <- rast(paste0(filePath, GlobalDEM))

# define L0 (minimum resolution)
# Minimum resolution requires 2*2 grid cell
L0 <- res(DEM)[1]*2

# The code below will work on only a subsection of the global DEM to test the effect of changing resolutions on complexity estimates
# This subsection is equal to ~9 cells at a spatial scale 2 orders of magnitude larger than L0.
# Define largest spatial scale
Lmax <- L0*(10^2) 
# Define boundaries of map subsample.
minX <- 3.5e+06; minY <- 0 
Xlimit <- minX+(Lmax*3); Ylimit <- minY+(Lmax*3)
# Outline scale over which to vary L
L <- c(L0*c(2,3,4,5,6,10^1,10^2))

# Test how varying L between two and less than one order of magnitude (37.4km and <1km) influences complexity estimates
LTest <- lapply(length(L):1, VaryingL)
# Convert outputs into a data frame
LTestdf <- rbindlist(LTest)
# and save data as a checkpoint
write.csv(LTestdf, paste0(filePath, 'ResolutionSens.csv'), row.names = FALSE)

#################################################
# STEP 3: Test the effect of changing resolution on complexity estimates
#################################################

# Load data from checkpoint (if needed)
# LTestdf <- read.csv(paste0(filePath, 'ResolutionSens.csv'))

# Firstly, it appears using the planar equation returns incorrect estimates of fractal dimension (i.e. D < 2)
# Visualize how the complexity estimates calculated using the two approaches compare to one another.
# Patterns in Rugosity
ggplot(LTestdf, aes(x=R_theory, y=R, col = L)) +
  geom_point(size = 2) +
  xlab('\nR (Equation)') +
  ylab('R (R function)\n') +
  scale_color_gradientn(colours = viridis(n = 4, option = 'cividis'),
                        guide = guide_colorbar(ticks = F, title = 'Scale (L)',
                                               reverse = F, label = T,
                                               na.value = "white")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.05,0.6), 
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# Patterns in Height Range
ggplot(LTestdf, aes(x=H_theory, y=H, col = L)) +
  geom_point(size = 2) +
  xlab('\nH (Equation)') +
  ylab('H (R function)\n') +
  scale_color_gradientn(colours = viridis(n = 4, option = 'cividis'),
                        guide = guide_colorbar(ticks = F, title = 'Scale (L)',
                                               reverse = F, label = T,
                                               na.value = "white")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.05,0.6), 
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# Patterns in Fractal Dimension
ggplot(LTestdf, aes(x=D_theory, y=D, col = L)) +
  geom_point(size = 2) +
  xlab('\nD (Equation)') +
  ylab('D (R function)\n') +
  scale_color_gradientn(colours = viridis(n = 4, option = 'cividis'),
                        guide = guide_colorbar(ticks = F, title = 'Scale (L)',
                                               reverse = F, label = T,
                                               na.value = "white")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.05,0.6), 
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))
# Patterns are largely consistent across scales

# Rugosity and Fractal Dimension are additive values and so regardless of the scale used the mean values of each should correspond for a given area.
aggregate(LTestdf[,c('H_theory','R_theory','D_theory','H','R','D')], list(LTestdf$L), mean)
# The approach using the plane equation does not appear to reflect this. However the approach using the inbuilt R functions does, with the R and D values at 0.5, 1, and 2 orders of magnitude reflecting similar scales of complexity

# ------------------------------------------
# Given that all further analyses will now focus on using the inbuilt R functions to calculate complexity the value of L0 can be set as equal to a single DEM pixel.
# It is therefore necessary to test how this effects the influence of resolution on estimates of complexity

# Re-define L0 (minimum resolution)
L0 <- res(DEM)[1]
# Re-define largest spatial scale
Lmax <- L0*(10^2) 
# Re-outline scale over which to vary L
L <- c(L0*c(2,3,4,5,6,10^1,10^2))

# Test how varying L between two and one order of magnitude (25km and 0.5km) influences complexity estimates
LTest2 <- lapply(length(L):1, VaryingL, test.both = FALSE) # this analysis will only focus on the approach using the inbuilt R surfaceArea function to calculate R and D.
# Convert outputs into a data frame
LTestdf2 <- rbindlist(LTest2)
# and save data as a checkpoint
write.csv(LTestdf2, paste0(filePath, 'ResolutionSens2.csv'), row.names = FALSE)

# Load data from checkpoint (if needed)
# LTestdf2 <- read.csv(paste0(filePath, 'ResolutionSens2.csv'))

# Rugosity and Fractal Dimension are additive values and so regardless of the scale used the mean values of each should correspond for a given area.
aggregate(LTestdf2[,c('H','R','D')], list(LTestdf2$L), mean)
# Additive estimates of rugosity and fractal dimension remain largely consistent down to a scale of 1240m (~km).

# Confirm variable distributions
range(LTestdf2$H)
range(LTestdf2$R)
range(LTestdf2$D) # A couple of pixels with values higher than expected.