# This script is for producing schematic plots to showcase global patterns in geometric complexity 
# and how different ecosystem types and land use scenarios align with these measures
# The ecosystem classification used follows the typology outlined in Keith et al. (2022) Nature.
# The Land Use classifications have been derived from the Land Cover 2020 data product from the Copernicus Climate Change service

# Date last modified: May 2024
# Primary Author: James Cant
# -----------------------------------------------------------------------------------------

# clear working directory
rm(list = ls())

# load required files
library(terra)
library(raster)
library(ggplot2)
library(viridis)
library(stringr)
library(dplyr)
library(brms)
library(habtools)
library(fishualize)
library(sf)

# is this the first time running the below script? i.e. does the ecosystem typology zipped folder need unpacking
FirstRun <- FALSE

# Set random number seed
set.seed(458967)
# prevent rounding of small values
options(digits = 22)

###############################
# STEP 1: Plot global patterns
###############################

# Define directory pathway to raster files
ComplexRast <- '/FILE_DIRECTORY_1/'

# Load Complexity rasters
DRast <- rast(paste0(ComplexRast, 'GlobalFractalDimension.tif'))
RRast <- rast(paste0(ComplexRast, 'GlobalRugosity.tif'))
HRast <- rast(paste0(ComplexRast, 'GlobalHeightRange.tif'))

# plot maps
# Fractal Dimension
plot(DRast, col = fish(50, direction = -1, option = 'Variola_louti'),
     buffer = FALSE,
     plg = list(x = 'bottom', at = c(2,3), digits = 1, tic = 'none', size = c(1,2.5)),
     box = FALSE, 
     axes = FALSE)

# Rugosity
plot(log10(RRast), col = fish(50, direction = 1, option = 'Ostracion_whitleyi'),
     buffer = FALSE, 
     plg = list(x = 'bottom', at = c(-12,0), tic = 'none', size = c(1,2.5)),
     box = FALSE, 
     axes = FALSE)

# Height Range
plot(log10(HRast), col = fish(50, direction = 1, option = 'Acanthurus_leucosternon'), 
     buffer = FALSE, 
     plg = list(x = 'bottom', at = c(-5,0), tic = 'none', size = c(1,2.5)),
     box = FALSE, 
     axes = FALSE)

###############################
# STEP 2: Ecosystem typology and Land Use Data Download
###############################

# Define file pathways for processing the Land use and ecosystem typology files
EcoTypes <- '/FILE_DIRECTORY_2/' # These maps can be downloaded from the following repository: 
# Keith, D. et al. (2020). Indicative distribution maps for Ecosystem Functional Groups - Level 3 of IUCN Global Ecosystem Typology (Version 2.1.1) [Data set]. Zenodo. DOI: 10.5281/zenodo.3546513
# and unzipped using untar() [see below].
FilePath <- '/FILE_DIRECTORY_3/'

# Extract typology maps from downloaded typology .tar file
if(FirstRun == TRUE) { untar(paste0(FilePath, "all-maps-raster-geotiff.tar.bz2"), exdir = EcoTypes) }
# only needed once to unzip downloaded maps.
# Note: once these ecosystem typology maps have been extracted then need to be reprojected to match the mollweide projection, resolution, and extent of the complexity rasters.
# This process is faster if done in batch within QGIS, first using the warp function to change the desired resolution and CRS of the IUCN maps before clipping these reprojected maps to the same extent as the complexity maps.
# The reprojected rasters are then saved into the 'Ecotypes' folder ahead of the processing below.

###############################
# STEP 3: Identify Classification Categories (only needed if first time running script)
###############################

### Ecosystem types -----------------------------------

if(FirstRun == TRUE) {
  # List available ecosystem type maps
  fileNames <- list.files(EcoTypes, pattern = '.tif$')
  fileNames <- str_extract(fileNames, '[^.]*.[^.]*') # remove unnecessary file name details
  fileNames <- sub("_", " ", fileNames)
  fileNames <- sub("_", " ", fileNames) # remove second underscore
  fileNames <- word(fileNames, 3)
  # Extract necessary file details
  EcoTypesDF <- data.frame(Realm = gsub("[^a-zA-Z]", "", fileNames),
                         Biome = gsub(".*?([0-9]+).*", "\\1", fileNames),
                         EFG = sapply(strsplit(fileNames, '\\D+'),'[',3),
                         # Place holder for complexity details
                         Dmean =  rep(NA, length.out = length(fileNames)),
                         Dmin =  rep(NA, length.out = length(fileNames)),
                         Dmax = rep(NA, length.out = length(fileNames)),
                         Dsigma = rep(NA, length.out = length(fileNames)),
                         Rmean =  rep(NA, length.out = length(fileNames)),
                         Rmin =  rep(NA, length.out = length(fileNames)),
                         Rmax = rep(NA, length.out = length(fileNames)),
                         Rsigma = rep(NA, length.out = length(fileNames)))
  
  # Output storage for density plots
  EcoTypesDensityD <- EcoTypesDensityR <- list()
  
  # regenerate file list for extraction process
  fileNames <- list.files(EcoTypes, pattern = '.tif$', full.names = TRUE)

  ### Land Use types -----------------------------------

  # Open land use raster
  LandUse <- rast(paste0(ComplexRast, 'LandCover2022_Mollweide_1870m.tif'))
  # List unique scenario codes
  LU_cats <- unique((values(LandUse)))
  LU_cats <- LU_cats[!is.na(LU_cats)]
  # List unique scenario names
  CatNames <- c('Cropland', 'Mosaic Vegetation/Cropland', 'Tree Cover (Broadleaved)', 'Tree Cover (Needleleaved)',
                'Mixed Tree Cover', 'Mosaic Vegetation', 'Shrubland', 'Grassland', 'Lichens & Mosses', 'Sparse Vegetation',
                'Wetlands', 'Urban', 'Bare substrate', 'Permenant Snow & Ice')
  
  # Define land use scenario list
  LU_df <- data.frame(Cat_Name = CatNames)

  # Place holder for complexity details
  LU_df$Dmean =  rep(NA, length(CatNames))
  LU_df$Dmin =  rep(NA, length(CatNames))
  LU_df$Dmax = rep(NA, length(CatNames))
  LU_df$Dsigma = rep(NA, length(CatNames))
  LU_df$Rmean =  rep(NA, length(CatNames))
  LU_df$Rmin =  rep(NA, length(CatNames))
  LU_df$Rmax = rep(NA, length(CatNames))
  LU_df$Rsigma = rep(NA, length(CatNames))
  
  # Output storage for density plots
  LUDensityD <- LUDensityR <- list()

  ###############################
  # STEP 4: Extract Complexities (again only needed if first time running script)
  ###############################

  ### 1. Ecosystem Types
  # Working through each ecosystem type, overlay ecosystem map onto complexity map and extract corresponding complexity details 

  for(ii in 1:dim(EcoTypesDF)[1]) { # working through each ecosystem typology file
    # progress read out
    print(ii)

    # Open selected ecosystem raster (and re-project to match complexity raster)
    RastSelect <- rast(fileNames[ii])
    # Remove minor occurrences of selected ecosystem
    RastSelect[RastSelect != 1] <- NA
    # re-load required complexity rasters (multicore processing can't call the files from the global environment)
    tmpD <- rast(paste0(ComplexRast, 'GlobalFractalDimension.tif'))
    tmpR <- rast(paste0(ComplexRast, 'GlobalRugosity.tif'))
    # Ensure rasters align
    RastSelect <- resample(RastSelect, tmpR)
  
    # Isolate corresponding complexity values
    tmpD[is.na(RastSelect)] <- NA
    tmpR[is.na(RastSelect)] <- NA
  
    # Extract descriptive values
    EcoTypesDF$Dmax[ii] <- max(values(tmpD),na.rm = TRUE)
    EcoTypesDF$Dmin[ii] <- min(values(tmpD),na.rm = TRUE)
    EcoTypesDF$Dmean[ii] <- mean(values(tmpD),na.rm = TRUE)
    EcoTypesDF$Dsigma[ii] <- sd(values(tmpD),na.rm = TRUE)
    EcoTypesDF$Rmax[ii] <- max(values(tmpR),na.rm = TRUE)
    EcoTypesDF$Rmin[ii] <- min(values(tmpR),na.rm = TRUE)
    EcoTypesDF$Rmean[ii] <- mean(values(tmpR),na.rm = TRUE)
    EcoTypesDF$Rsigma[ii] <- sd(values(tmpR),na.rm = TRUE)
    
    # Extract distribution of complexity values
    EcoTypesDensityD[[ii]] <- density(na.omit(values(tmpD)))
    tmpR <- log10(tmpR) # Apply data transformation to rugosity values
    tmpR[is.infinite(tmpR)] <- NA
    EcoTypesDensityR[[ii]] <- density(na.omit(values(tmpR)))
  }
  
  # Add an additional ecosystem clustering variable (to group some of the realms)
  EcoTypesDF$Realm_2 <- NA
  EcoTypesDF[EcoTypesDF$Realm %in% c('F'),]$Realm_2 <- 'Freshwater'
  EcoTypesDF[EcoTypesDF$Realm %in% c('MFT','MT'),]$Realm_2 <- 'Coastal'
  EcoTypesDF[EcoTypesDF$Realm %in% c('S','SF','SM'),]$Realm_2 <- 'Subterranean'
  EcoTypesDF[EcoTypesDF$Realm %in% c('T'),]$Realm_2 <- 'Terrestrial'
  EcoTypesDF[EcoTypesDF$Realm %in% c('TF'),]$Realm_2 <- 'Wetland'
  EcoTypesDF[EcoTypesDF$Realm %in% c('M','FM'),]$Realm_2 <- 'Marine'
  
  # Clean data
  EcoTypesDF[sapply(EcoTypesDF, is.infinite)] <- NA # remove infinite's
  # remove categories with no data
  EcoTypesDensityD <- EcoTypesDensityD[complete.cases(EcoTypesDF)] 
  EcoTypesDensityR <- EcoTypesDensityR[complete.cases(EcoTypesDF)] 
  EcoTypesDF <- EcoTypesDF[complete.cases(EcoTypesDF),] 
  
  ### 2. Land Use types

  for(ii in 1:length(CatNames)) { # working through each land use category
    # progress read out
    print(ii)
    
    # Identify scenario codes associated with the selected category
    if(CatNames[ii] == 'Cropland'){ CatCode <- c(10,11,12,20) }
    if(CatNames[ii] == 'Mosaic Vegetation/Cropland'){ CatCode <-  c(30,40) }
    if(CatNames[ii] == 'Tree Cover (Broadleaved)'){ CatCode <-  c(50:62) }
    if(CatNames[ii] == 'Tree Cover (Needleleaved)'){ CatCode <-  c(70:82) }
    if(CatNames[ii] == 'Mixed Tree Cover'){ CatCode <-  c(90) }
    if(CatNames[ii] == 'Mosaic Vegetation'){ CatCode <-  c(100:110) }
    if(CatNames[ii] == 'Shrubland'){ CatCode <-  c(120:122) }
    if(CatNames[ii] == 'Grassland'){ CatCode <-  c(130) }
    if(CatNames[ii] == 'Lichens & Mosses'){ CatCode <-  c(140) }
    if(CatNames[ii] == 'Sparse Vegetation'){ CatCode <-  c(150:153) }
    if(CatNames[ii] == 'Wetlands'){ CatCode <-  c(160:180) }
    if(CatNames[ii] == 'Urban'){ CatCode <-  c(190) }
    if(CatNames[ii] == 'Bare substrate'){ CatCode <-  c(200:202) }
    if(CatNames[ii] == 'Permenant Snow & Ice'){ CatCode <-  c(220) }
    
    # Duplicate Land use raster
    RastSelect <- LandUse
    # Isolate selected Land Use scenario
    RastSelect[!(values(RastSelect) %in% CatCode)] <- NA 
    # re-load required complexity rasters (multicore processing can't call the files from the global environment)
    tmpD <- rast(paste0(ComplexRast, 'GlobalFractalDimension.tif'))
    tmpR <- rast(paste0(ComplexRast, 'GlobalRugosity.tif'))
    # Ensure rasters align
    ext(RastSelect) <- ext(tmpD)
  
    # Isolate corresponding complexity values
    tmpD[is.na(RastSelect)] <- NA
    tmpR[is.na(RastSelect)] <- NA
    
    # Extract values
    LU_df$Dmax[ii] <- max(values(tmpD),na.rm = TRUE)
    LU_df$Dmin[ii] <- min(values(tmpD),na.rm = TRUE)
    LU_df$Dmean[ii] <- mean(values(tmpD),na.rm = TRUE)
    LU_df$Dsigma[ii] <- sd(values(tmpD),na.rm = TRUE)
    LU_df$Rmax[ii] <- max(values(tmpR),na.rm = TRUE)
    LU_df$Rmin[ii] <- min(values(tmpR),na.rm = TRUE)
    LU_df$Rmean[ii] <- mean(values(tmpR),na.rm = TRUE)
    LU_df$Rsigma[ii] <- sd(values(tmpR),na.rm = TRUE)
    
    # Extract distribution of complexity values
    LUDensityD[[ii]] <- density(na.omit(values(tmpD)))
    tmpR <- log10(tmpR) # Apply data transformation
    tmpR[is.infinite(tmpR)] <- NA
    LUDensityR[[ii]] <- density(na.omit(values(tmpR)))
  }

  # Save files (data checkpoint)
  # Ecosystem types
  write.csv(EcoTypesDF, paste0(FilePath, 'EcosystemComplexity.csv'), row.names = FALSE)
  saveRDS(EcoTypesDensityD, paste0(FilePath, 'EcoTypesDensityD.RData'))
  saveRDS(EcoTypesDensityR, paste0(FilePath, 'EcoTypesDensityR.RData'))
  # Land Use types
  write.csv(LandUseDF, paste0(FilePath, 'LandUseComplexity.csv'), row.names = FALSE)
  saveRDS(LUDensityD, paste0(FilePath, 'LandUseDensityD.RData'))
  saveRDS(LUDensityR, paste0(FilePath, 'LandUseDensityR.RData'))
}

###############################
# STEP 5: Visualize complexity patterns
###############################

# Load files
EcoTypesDF <- read.csv(paste0(FilePath, 'EcosystemComplexity.csv'))
LandUseDF <- read.csv(paste0(FilePath, 'LandUseComplexity.csv'))

### Ecosystem Type Plot -------------------------------

# Select representative environments
Location1 <- EcoTypesDF[EcoTypesDF$Dmean == min(EcoTypesDF$Dmean),] # Seamounts
Location2 <- EcoTypesDF[EcoTypesDF$Rmean == min(EcoTypesDF$Rmean),] # Sub-glacial Lakes
Location3 <- EcoTypesDF[EcoTypesDF$Rmean == max(EcoTypesDF$Rmean),] # Oceanic temperate Rainforests
Location4 <- EcoTypesDF[which(EcoTypesDF$Realm == 'T' & EcoTypesDF$Biome == 6 & EcoTypesDF$EFG == 3),] # Polar Tundra
Location5 <- EcoTypesDF[EcoTypesDF$Dmean == max(EcoTypesDF$Dmean),] # Coastal River Deltas
Location6 <- EcoTypesDF[which(EcoTypesDF$Realm == 'F' & EcoTypesDF$Biome == 2 & EcoTypesDF$EFG == 6),] # Permanent salt lakes
Location7 <- EcoTypesDF[which(EcoTypesDF$Realm == 'TF' & EcoTypesDF$Biome == 1 & EcoTypesDF$EFG == 3),] # Permanent marshland
Location8 <- EcoTypesDF[which(EcoTypesDF$Realm == 'T' & EcoTypesDF$Biome == 5 & EcoTypesDF$EFG == 3),] # Sclerophyll hot deserts and semi-deserts

# Reformat variables to aid visualization clarity
EcoTypesDF$Realm_2 <- factor(EcoTypesDF$Realm_2,
                           levels = c('Marine', 'Subterranean', 'Freshwater', 'Terrestrial', 'Coastal', 'Wetland'))

# Create plot
ggplot(EcoTypesDF, aes(x=Dmean, y=Rmean)) +
  annotate('rect', xmin = 2.235, xmax = 2.35, ymin = 0.003, ymax = 0.11, fill = 'grey', alpha = 0.5) +
  annotate('rect', xmin = 2.1, xmax = 2.235, ymin = 0.0001, ymax = 0.003, fill = 'grey', alpha = 0.2) +
  geom_hline(aes(yintercept = 0.003), linetype = "dashed", col = 'black', linewidth = 1.5) +
  geom_vline(aes(xintercept = 2.235), linetype = "dashed", col = 'black', linewidth = 1.5) +
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.15, yend = 0.05), 
               data = Location1, color = "#00204DFF", linewidth = 2, alpha = 1, linetype = 'dashed') + # Seamounts
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.245, yend = 0.00021), 
               data = Location2, color = "lightblue", linewidth = 2, alpha = 0.7, linetype = 'dashed') + # Subglacial Lakes
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.31, yend = 0.1), 
               data = Location3, color = "#5F9258FF", linewidth = 2, alpha = 0.7, linetype = 'dashed') + # Oceanic Temperate Forests
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.33, yend = Rmean), 
               data = Location4, color = "#5F9258FF", linewidth = 2, alpha = 0.7, linetype = 'dashed') + # Polar Tundra
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.34, yend = 0.0005), 
               data = Location5, color = "#CBBA69FF", linewidth = 2, alpha = 0.7, linetype = 'dashed') + # Coastal River Deltas
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.14, yend = 0.0003), 
               data = Location6, color = "lightblue", linewidth = 2, alpha = 0.9, linetype = 'dashed') + # Permanent Salt Lakes
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.34, yend = 0.00018), 
               data = Location7, color = "#FDE725FF", linewidth = 2, alpha = 0.9, linetype = 'dashed') + # Permanent Marsh
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.2, yend = 0.0005), 
               data = Location8, color = "#5F9258FF", linewidth = 2, alpha = 0.7, linetype = 'dashed') + # Sclerophyll hot deserts and semi-deserts
  geom_point(aes(color = Realm_2), size = 9) +
  xlab(NULL) + # Fractal Dimension
  ylab(NULL) + # Rugosity
  scale_y_continuous(trans = 'log10', breaks = c(0.0002, 0.001, 0.01, 0.1), labels = function(i) { format(i, scientific = F) }, expand = c(0,0)) +
  scale_x_continuous(limits = c(2.10,2.35), expand = c(0,0)) +
  scale_color_manual(values = c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF"),
                     guide = guide_legend(title = 'Realm', 
                                          reverse = F, label = T,
                                          na.value = "white")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_blank(),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = 0)) +
  theme(axis.line = element_line(color = 'black'))

### Land Use Plot -------------------------------

# Identify the human mediated landscapes
Urban <- LandUseDF[LandUseDF$Cat_Name == 'Urban',] 
Cropland1 <- LandUseDF[LandUseDF$Cat_Name == 'Cropland',]
Cropland2 <- LandUseDF[LandUseDF$Cat_Name == 'Mosaic Vegetation/Cropland',]

# Create plot
ggplot(LandUseDF, aes(x = Dmean, y=Rmean)) +
  annotate('rect', xmin = 2.255, xmax = 2.3, ymin = 0.01, ymax = 0.061, fill = 'grey', alpha = 0.5) +
  annotate('rect', xmin = 2.208, xmax = 2.255, ymin = 0.0018, ymax = 0.01, fill = 'grey', alpha = 0.2) +
  geom_hline(aes(yintercept = 0.01), linetype = "dashed", col = 'black', linewidth = 1.5) +
  geom_vline(aes(xintercept = 2.255), linetype = "dashed", col = 'black', linewidth = 1.5) +
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.225, yend = Rmean), 
               data = Urban, color = "#000004FF", linewidth = 2, alpha = 1, linetype = 'dashed') + # Urban
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.225, yend = Rmean), 
               data = Cropland1, color = "#0D0B2AFF", linewidth = 2, alpha = 1, linetype = 'dashed') + # Cropland
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.225, yend = 0.008), 
               data = Cropland2, color = "#281259FF", linewidth = 2, alpha = 1, linetype = 'dashed') + # Mosaic Vegetation/Cropland
  geom_point(aes(color = Cat_Name), size = 9) +
  xlab(NULL) + # Fractal Dimension
  ylab(NULL) + # Rugosity
  scale_y_continuous(trans = 'log10', breaks = c(0.003, 0.01, 0.03), labels = function(i) { format(i, scientific = F) }, expand = c(0,0)) +
  scale_x_continuous(limits = c(2.208,2.3), expand = c(0,0)) +
  scale_color_manual(values = magma(n = length(unique(LandUseDF$Cat_Name))),
                     limits = c("Urban","Cropland","Mosaic Vegetation/Cropland","Mosaic Vegetation",
                                "Grassland","Shrubland","Lichens & Mosses","Sparse Vegetation",         
                                "Mixed Tree Cover","Tree Cover (Needleleaved)","Tree Cover (Broadleaved)",
                                'Wetlands',"Bare substrate","Permenant Snow & Ice"),
                     guide = guide_legend(title = 'Land Cover Type', 
                                          reverse = F, label = T,
                                          na.value = "white")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_blank(),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = 0)) +
  theme(axis.line = element_line(color = 'black'))


#####################################
# STEP 6: Analyse ecosystem complexity patterns 
######################################

# Implement brms regression to explore the relationship between rugosity and fractal dimension across ecosystem types

# Determine most appropriate model format by running linear and non-linear form
EcoTypesMod <- brm(log10(Rmean) ~ Dmean , data = EcoTypesDF, family = 'gaussian', chains = 20, iter = 5000, warmup = 1000)
EcoTypesMod2 <- brm(log10(Rmean) ~ Dmean + I(Dmean^2), data = EcoTypesDF, family = 'gaussian', chains = 20, iter = 5000, warmup = 1000)
model_weights(EcoTypesMod, EcoTypesMod2, weights = 'loo') # none-linear is the better fit

# Extract estimate of model fit
bayes_R2(EcoTypesMod2)


######################################
# STEP 7: Extract examples of complex environments
######################################

# Define working resolution
L <- res(DRast)[1]
# extract raster projection
Moll <- crs(DRast)

# Isolate the fractal dimension and rugosity regimes of selected global features

## 1. Amazon Basin
# Calculate mollweide coordinates from lat longs obtained from Google Earth
FeatureCoords1 <- matrix(c(-62.2159, -2.4653), ncol = 2)
FeatureCoords1 <- spTransform(SpatialPoints(FeatureCoords1, CRS('+proj=longlat')), CRS(Moll))@coords
# Crop the complexity rasters around these selected points
Feature1D <- dem_crop(DRast, x0 = FeatureCoords1[1], y0 = FeatureCoords1[2], L = L*500, plot = FALSE)
Feature1R <- dem_crop(RRast, x0 = FeatureCoords1[1], y0 = FeatureCoords1[2], L = L*500, plot = FALSE)
# Estimate mean complexity values
mean(values(Feature1D), na.rm = TRUE)
mean(values(Feature1R), na.rm = TRUE)

## 2. Mariana Trench
# Calculate mollweide coordinates from lat longs obtained from Google Earth
FeatureCoords2 <- matrix(c(142.1996, 11.3493), ncol = 2)
FeatureCoords2 <- spTransform(SpatialPoints(FeatureCoords2, CRS('+proj=longlat')), CRS(Moll))@coords
# Crop the complexity rasters around these selected points
Feature2D <- dem_crop(DRast, x0 = FeatureCoords2[1], y0 = FeatureCoords2[2], L = L*500, plot = FALSE)
Feature2R <- dem_crop(RRast, x0 = FeatureCoords2[1], y0 = FeatureCoords2[2], L = L*500, plot = FALSE)
# Estimate mean complexity values
mean(values(Feature2D), na.rm = TRUE)
mean(values(Feature2R), na.rm = TRUE)

## 3. Sahara Desert
# Calculate mollweide coordinates from lat longs obtained from Google Earth
FeatureCoords3 <- matrix(c(25.6628, 23.4162), ncol = 2)
FeatureCoords3 <- spTransform(SpatialPoints(FeatureCoords3, CRS('+proj=longlat')), CRS(Moll))@coords
# Crop the complexity rasters around these selected points
Feature3D <- dem_crop(DRast, x0 = FeatureCoords3[1], y0 = FeatureCoords3[2], L = L*500, plot = FALSE)
Feature3R <- dem_crop(RRast, x0 = FeatureCoords3[1], y0 = FeatureCoords3[2], L = L*500, plot = FALSE)
# Estimate mean complexity values
mean(values(Feature3D), na.rm = TRUE)
mean(values(Feature3R), na.rm = TRUE)

# Convert data into data-frames (required for plotting rasters using ggplot)
Feature1D_df <- st_as_sf(terra::as.points(Feature1D))
Feature1R_df <- st_as_sf(terra::as.points(Feature1R))
Feature2D_df <- st_as_sf(terra::as.points(Feature2D))
Feature2R_df <- st_as_sf(terra::as.points(Feature2R))
Feature3D_df <- st_as_sf(terra::as.points(Feature3D))
Feature3R_df <- st_as_sf(terra::as.points(Feature3R))
colnames(Feature1D_df)[1] <- colnames(Feature1R_df)[1] <- colnames(Feature2D_df)[1] <- colnames(Feature2R_df)[1] <-
  colnames(Feature3D_df)[1] <- colnames(Feature3R_df)[1] <- "value"

# identify the maximum and minimum fractal dimension and rugosity estimates across these selected features (to keep plot colour scales consistent)
maxD <- ceiling(max(c(Feature1D_df$value, Feature2D_df$value, Feature3D_df$value),na.rm = T)*10)/10 # this little trick ensures the value is rounded up (at one decimal place)
maxR <- log10(ceiling(max(c(Feature1R_df$value, Feature2R_df$value, Feature3R_df$value),na.rm = T)*10)/10)
minD <- floor(min(c(Feature1D_df$value, Feature2D_df$value, Feature3D_df$value),na.rm = T)*10)/10
minR <- -12

# 1. Amazon Basin
# Fractal Dimension
ggplot(data = Feature1D_df) +
  geom_sf(aes(col = value)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_colour_gradientn(colours = fish(50, direction = -1, option = 'Variola_louti'),
                       limits = c(minD,maxD),
                       guide = guide_colorbar(ticks = F, title = 'D',
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_sf(expand = F) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), axis.text.y = element_text(size = 30, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.9,0.8),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# Rugosity
ggplot(data = Feature1R_df) +
  geom_sf(aes(col = log10(value))) +
  xlab(NULL) +
  ylab(NULL) +
  scale_colour_gradientn(colours = fish(50, direction = 1, option = 'Ostracion_whitleyi'),
                       guide = guide_colorbar(ticks = F, title = expression('log'[10]*'(R)'),
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_sf(expand = F) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), axis.text.y = element_text(size = 30, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.9,0.8),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# 2. Mariana Trench
# Fractal Dimension
ggplot(data = Feature2D_df) +
  geom_sf(aes(col = value)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_colour_gradientn(colours = fish(50, direction = -1, option = 'Variola_louti'),
                       limits = c(minD,maxD),
                       guide = guide_colorbar(ticks = F, title = 'D',
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_sf(expand = F) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), axis.text.y = element_text(size = 30, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.9,0.8),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# Rugosity
ggplot(data = Feature2R_df) +
  geom_sf(aes(col = log10(value))) +
  xlab(NULL) +
  ylab(NULL) +
  scale_colour_gradientn(colours = fish(100, direction = 1, option = 'Ostracion_whitleyi'),
                       guide = guide_colorbar(ticks = F, title = expression('log'[10]*'(R)'),
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_sf(expand = F) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), axis.text.y = element_text(size = 30, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.9,0.8),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# 3. Sahara Desert
# Fractal Dimension
ggplot(data = Feature3D_df) +
  geom_sf(aes(col = value)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_colour_gradientn(colours = fish(50, direction = -1, option = 'Variola_louti'),
                       limits = c(minD,maxD),
                       guide = guide_colorbar(ticks = F, title = 'D',
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_sf(expand = F) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), axis.text.y = element_text(size = 30, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.9,0.8),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# Rugosity
ggplot(data = Feature3R_df) +
  geom_sf(aes(col = log10(value))) +
  xlab(NULL) +
  ylab(NULL) +
  scale_colour_gradientn(colours = fish(100, direction = 1, option = 'Ostracion_whitleyi'),
                       guide = guide_colorbar(ticks = F, title = expression('log'[10]*'(R)'),
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_sf(expand = F) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), axis.text.y = element_text(size = 30, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.9,0.8),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# ----------------------------------------------------- End of Code ---------------------------------