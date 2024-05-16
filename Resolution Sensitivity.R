# This script is for testing the sensitivity of global estimates of fractal dimension and rugosity to increasing map resolution.
# Accordingly, this script also tests the sensitivity of fractal dimension and rugosity estimates to the changing order of magnitude between L0 and L.
#
# Main author: James Cant (jic2@st-andrews.ac.uk)
# Date Last Modified: May 2024
# -----------------------------------------------------------

# Clear workspace
rm(list=ls(all=TRUE))

# Define directory pathways
filePath <- '/FILE_DIRECTORY_1/'
# Define raster files
GlobalDEM <- 'GlobalDEM_Mollweide_187m.tif'# Global raster (187m resolution)
# This file is a global DEM comprising terrestrial topography and ocean bathymetry produced by mosaicing tiles obtained from GEBCO (see Build Global DEM.R)

# load necessary packages
library(sp)
library(raster)
library(terra)
library(dplyr)
library(parallel)
library(data.table)
library(ggplot2)
library(viridis)
library(habtools)
library(xpectr)

#################################################
# STEP 1: Define complexity extraction function 
#################################################

# *** This function is a compilation of functions obtained from the habtools packages, which itself is based on the analytical tools 
# presented in Torres-Pulliza et al. (2020) Nature Ecol. & Evol 4. 1495-1501 ***
ExtractComplexity <- function(ii, DEMCell, L, L0){
  
  # No need to process the selected Grid Cell if it only contains NA's
  if(is.na(max(values(DEMCell))) & is.na(min(values(DEMCell)))){
    R <- NA
    H <- NA
    D <- NA
  } else {
    
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
  return(list(H = H, R = R, D = D))
}

#################################################
# STEP 2: Test the effect of changing resolution on complexity estimates
#################################################

# Define function to test the varying effects of changing L on complexity estimates
VaryingL <- function(ii){
  
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
      
      # For each grid cell calculate complexity using only in the inbuilt R functions.
      CellComplex <- ExtractComplexity(DEMCell = tmpRast, L = L[ii], L0 = L0)
        
      # Save outputs
      Outputs <- c(L[ii], # scale
                   xx, # Grid index values
                   yy, 
                   x0, # Longitude
                   y0, # Latitude
                   as.numeric(CellComplex$H),
                   as.numeric(CellComplex$R),
                   as.numeric(CellComplex$D))
        
      # And store
      PixelValues <- rbind(PixelValues, Outputs)
      colnames(PixelValues) <- c('L','X','Y','Lon','Lat','H','R','D')
    }
  }      
  
  # and return
  return(PixelValues)
}

# load Global DEM raster
DEM <- rast(paste0(filePath, GlobalDEM))

# define L0 (minimum resolution)
L0 <- res(DEM)[1]

# The code below will work on only a subsection of the global DEM to test the effect of changing resolutions on complexity estimates
# The dimensions of this subsection is equal to ~9 pixels each at a resolution equal to 2 orders of magnitude larger than L0.
# Define largest spatial scale
Lmax <- L0*(10^2) 
# Define location of map subsample.
minX <- 2.5e+06; minY <- 0 # essentially this location can be random across the centre of the DEM.
Xlimit <- minX+(Lmax*3); Ylimit <- minY+(Lmax*3)
# Outline scale over which to vary L
L <- c(L0*c(2,3,4,5,6,10^1,10^2))

# Test how varying L between two and less than one order of magnitude (18.7km and <1km) influences complexity estimates
LTest <- lapply(length(L):1, VaryingL)
saveRDS(LTest, paste0(filePath, "LTest.rds")) # save as a checkpoint

# Convert outputs into a data frame
LTestdf <- rbindlist(LTest)

#################################################
# STEP 3: Test the effect of changing resolution on complexity estimates
#################################################

# Rugosity and Fractal Dimension are additive values and so regardless of the scale used the mean values of each should correspond for a given area.
meanDF <- aggregate(LTestdf[,c('H','R','D')], list(LTestdf$L), mean)
meanDF
# Additive estimates of Rugosity and fractal dimension remain largely consistent down to a scale of 748m (4*L0).

### Display scale sensitivity ----------------------
# calculate difference between Rugosity and Fractal Dimension scales
coeff <- max(meanDF$D)/max(meanDF$R)
# generate plot
ggplot(meanDF, aes(x=log10(Group.1))) +
  geom_line(aes(y=D), col = 'blue', linewidth = 1.5, linetype = 'dashed') +
  geom_line(aes(y=R*coeff), col = 'red', linewidth = 1.5, linetype = 'dashed') +
  geom_vline(xintercept = log10(1870), col = 'black', linetype = 'dotted', linewidth = 0.8) +
  xlab('\nScale (L)') +
  scale_x_continuous(labels = function(i) round(10^i)) +
  scale_y_continuous(name = 'D\n',
                     labels = function(i) round(i, digits = 2),
                     sec.axis = sec_axis(~./coeff, 
                                         name = 'R\n',
                                         labels = function(i) round(i, digits = 2))) +
  coord_cartesian(ylim = c(2,3)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.05,0.6), 
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))
# Following this sensitivity analysis it is appropriate to select a resolution of 1870m as a balance between enhancing resolution whilst ensuring the accuracy of complexity estimates.


# Confirm variable distributions
range(LTestdf$H)
range(LTestdf$R)
range(LTestdf$D)
# What percentage of D estimates exceed 3?
(length(LTestdf$D[LTestdf$D > 3])/length(LTestdf$D))*100

#--------------------------------------------------------------- End of Code ------------------------------