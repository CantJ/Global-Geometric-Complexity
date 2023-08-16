# This script is for producing schematic plots to showcase how different ecosystem types and land use scenarios align with measures of complexity
# The ecosystem classification used follows the typology outlined in Keith et al. (2022) Nature.
# The Land Use classifications have been derived from the Land Cover 2020 data product from the Copernicus Climate Change service

# Date last modified: August 2023
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

###############################
# STEP 1: Data Download
###############################

# Define file pathways
EcoTypes <- 'FileDirectory1/' # Location of Ecosystem typology maps. These maps can be downloaded from the following respositiory: 
# Keith, D. et al. (2020). Indicative distribution maps for Ecosystem Functional Groups - Level 3 of IUCN Global Ecosystem Typology (Version 2.1.1) [Data set]. Zenodo. DOI: 10.5281/zenodo.3546513
# and unzipped using untar() [see below].
ComplexRast <- 'FileDirectory2/' # Location of Complexity rasters
LURast <- 'FileDirectory3/' # Location of Land Use template

# Extract typology maps from downloaded typology .tar file
# untar(paste0(EcoTypes, "all-maps-raster-geotiff.tar.bz2"), exdir = EcoTypes) # only needed once to unzip downloaded maps.

###############################
# STEP 2: Define Extraction Function
###############################

# Complexity extraction function
ComplexExtract <- function(ii, Ecosystem = Ecosystem){
  # progress read out
  print(ii)
  
  # Are ecosystem classifications being used?
  if(Ecosystem == TRUE) {
    # Open ecosystem raster (and re-project to match complexity raster)
    RastSelect <- project(rast(fileNames[ii]), DRast, method = 'near', align = TRUE, mask = TRUE) 
    # Remove minor occurrences of selected ecosystem
    RastSelect[RastSelect == 2] <- NA
    # Ensure rasters align
    tmpD <- resample(DRast, RastSelect)
    tmpR <- resample(RRast, RastSelect)
  } else {
    # If Land use classifications are being used instead
    # Duplicate Land Use raster
    RastSelect <- LandUse
    # Isolate selected Land Use scenario
    RastSelect[RastSelect != LU_df$Cat1[ii]] <- NA
    # duplicate complexity rasters
    tmpD <- DRast2
    tmpR <- RRast2
  }
  
  # Isolate corresponding complexity values
  tmpD[is.na(RastSelect)] <- NA
  tmpR[is.na(RastSelect)] <- NA
  # Extract values
  Dmax <- max(values(tmpD),na.rm = TRUE)
  Dmin <- min(values(tmpD),na.rm = TRUE)
  Dmean <- mean(values(tmpD),na.rm = TRUE)
  Dsigma <- sd(values(tmpD),na.rm = TRUE)
  Rmax <- max(values(tmpR),na.rm = TRUE)
  Rmin <- min(values(tmpR),na.rm = TRUE)
  Rmean <- mean(values(tmpR),na.rm = TRUE)
  Rsigma <- sd(values(tmpR),na.rm = TRUE)
  # Return outputs
  return(c(Dmax = Dmax, Dmin = Dmin, Dmean = Dmean, Dsigma = Dsigma,
           Rmax = Rmax, Rmin = Rmin, Rmean = Rmean, Rsigma = Rsigma))
}

###############################
# STEP 3: Identify Classification Categories
###############################

### Ecosystem types -----------------------------------

# List avalible ecosystem type maps
fileNames <- list.files(EcoTypes, pattern = '.tif$')
fileNames <- gsub('_',' ', fileNames) # remove underscores
fileNames <- gsub('.tif','', fileNames) # and file extensions
# Extract necessary file details
EcoTypesDF <- data.frame(Realm = gsub("[^a-zA-Z]", "", word(fileNames, 1)),
                         Biome = word(fileNames, 1),
                         EFG = word(fileNames, 3,-1),
                         # Place holder for complexity details
                         Dmean =  rep(NA, length.out = length(fileNames)),
                         Dmin =  rep(NA, length.out = length(fileNames)),
                         Dmax = rep(NA, length.out = length(fileNames)),
                         Dsigma = rep(NA, length.out = length(fileNames)),
                         Rmean =  rep(NA, length.out = length(fileNames)),
                         Rmin =  rep(NA, length.out = length(fileNames)),
                         Rmax = rep(NA, length.out = length(fileNames)),
                         Rsigma = rep(NA, length.out = length(fileNames)))

# regenerate file list for extraction process
fileNames <- list.files(EcoTypes, pattern = '.tif$', full.names = TRUE)

### Land Use types -----------------------------------

# Open land use raster
LandUse <- rast(paste0(LURast, 'LandMask_Mollweide_1122m.tif'))
# List unique scenarios
LU_cats <- unique((values(LandUse)))
LU_cats <- LU_cats[!is.nan(LU_cats)]
# Refine and name land use scenario list
LU_df <- data.frame(Cat1 = LU_cats,
                    Cat_Name = rep(NA, length(LU_cats)))
LU_df[LU_df$Cat1 %in% c(10,11,12,20), 'Cat_Name'] <- 'Cropland'
LU_df[LU_df$Cat1 %in% c(30,40), 'Cat_Name'] <- 'Mosaic Vegetation/Cropland'
LU_df[LU_df$Cat1 %in% c(50:62), 'Cat_Name'] <- 'Tree Cover (Broadleaved)'
LU_df[LU_df$Cat1 %in% c(70:82), 'Cat_Name'] <- 'Tree Cover (Needleleaved)'
LU_df[LU_df$Cat1 %in% c(90), 'Cat_Name'] <- 'Mixed Tree Cover'
LU_df[LU_df$Cat1 %in% c(100:110), 'Cat_Name'] <- 'Mosaic Vegetation'
LU_df[LU_df$Cat1 %in% c(120:122), 'Cat_Name'] <- 'Shrubland'
LU_df[LU_df$Cat1 %in% c(130), 'Cat_Name'] <- 'Grassland'
LU_df[LU_df$Cat1 %in% c(140), 'Cat_Name'] <- 'Lichens & Mosses'
LU_df[LU_df$Cat1 %in% c(150:153), 'Cat_Name'] <- 'Sparse Vegetation'
LU_df[LU_df$Cat1 %in% c(160:180), 'Cat_Name'] <- 'Wetlands'
LU_df[LU_df$Cat1 %in% c(190), 'Cat_Name'] <- 'Urban'
LU_df[LU_df$Cat1 %in% c(200:202), 'Cat_Name'] <- 'Bare substrate'
LU_df[LU_df$Cat1 %in% c(220), 'Cat_Name'] <- 'Permenant Snow & Ice'

# Place holder for complexity details
LU_df$Dmean =  rep(NA, length(LU_cats))
LU_df$Dmin =  rep(NA, length(LU_cats))
LU_df$Dmax = rep(NA, length(LU_cats))
LU_df$Dsigma = rep(NA, length(LU_cats))
LU_df$Rmean =  rep(NA, length(LU_cats))
LU_df$Rmin =  rep(NA, length(LU_cats))
LU_df$Rmax = rep(NA, length(LU_cats))
LU_df$Rsigma = rep(NA, length(LU_cats))

###############################
# STEP 4: Extract Complexities
###############################

# Load Complexity rasters
DRast <- rast(paste0(ComplexRast, 'GlobalFractalDimension.tif'))
RRast <- rast(paste0(ComplexRast, 'GlobalRugosity.tif'))

# Working through each ecosystem type, overlay ecosystem map onto complexity map and extract corresponding complexity details 
EcoComplex <- lapply(1:length(fileNames), ComplexExtract, Ecosystem = TRUE) # takes a while to compute
# Repeat for each Land Use category
# save processing time by ensure complexity and land use raster are aligned ahead of loop.
DRast2 <- resample(DRast, LandUse) 
RRast2 <- resample(RRast, LandUse) 
# run land use extraction
LUComplex <- lapply(1:length(LU_cats), ComplexExtract, Ecosystem = FALSE)

# Assign complexity values to corresponding ecosystem type.
EcoTypesDF$Dmax <- sapply(EcoComplex,'[[', 1)
EcoTypesDF$Dmin <- sapply(EcoComplex,'[[', 2)
EcoTypesDF$Dmean <- sapply(EcoComplex,'[[', 3)
EcoTypesDF$Dsigma <- sapply(EcoComplex,'[[', 4)
EcoTypesDF$Rmax <- sapply(EcoComplex,'[[', 5)
EcoTypesDF$Rmin <- sapply(EcoComplex,'[[', 6)
EcoTypesDF$Rmean <- sapply(EcoComplex,'[[', 7)
EcoTypesDF$Rsigma <- sapply(EcoComplex,'[[', 8)

# Add an additional ecosystem clustering variable (to group some of the realms)
EcoTypesDF$Realm_2 <- NA
EcoTypesDF[EcoTypesDF$Realm %in% c('F'),]$Realm_2 <- 'Freshwater'
EcoTypesDF[EcoTypesDF$Realm %in% c('MFT','MT'),]$Realm_2 <- 'Coastal'
EcoTypesDF[EcoTypesDF$Realm %in% c('S','SF','SM'),]$Realm_2 <- 'Subterranean'
EcoTypesDF[EcoTypesDF$Realm %in% c('T'),]$Realm_2 <- 'Terrestrial'
EcoTypesDF[EcoTypesDF$Realm %in% c('TF'),]$Realm_2 <- 'Wetland'
EcoTypesDF[EcoTypesDF$Realm %in% c('M','FM'),]$Realm_2 <- 'Marine'

# Assign complexity values to corresponding Land Use type.
LU_df$Dmax <- sapply(LUComplex,'[[', 1)
LU_df$Dmin <- sapply(LUComplex,'[[', 2)
LU_df$Dmean <- sapply(LUComplex,'[[', 3)
LU_df$Dsigma <- sapply(LUComplex,'[[', 4)
LU_df$Rmax <- sapply(LUComplex,'[[', 5)
LU_df$Rmin <- sapply(LUComplex,'[[', 6)
LU_df$Rmean <- sapply(LUComplex,'[[', 7)
LU_df$Rsigma <- sapply(LUComplex,'[[', 8)

# Condense together values across repeated categories
LandUseDF <- LU_df %>%
  group_by(Cat_Name) %>%
  summarize(meanD = mean(Dmean, na.rm = TRUE),
            meanR = mean(Rmean, na.rm = TRUE))

# Save files (data checkpoint)
write.csv(EcoTypesDF, paste0(ComplexRast, 'EcosystemComplexity.csv'), row.names = FALSE)
write.csv(LandUseDF, paste0(ComplexRast, 'LandUseComplexity.csv'), row.names = FALSE)

###############################
# STEP 5: Visualize Complexities
###############################

# Load files from checkpoint (if needed)
EcoTypesDF <- read.csv(paste0(ComplexRast, 'EcosystemComplexity.csv'))
LandUseDF <- read.csv(paste0(ComplexRast, 'LandUseComplexity.csv'))

# Clean data ahead of visualization
EcoTypesDF[sapply(EcoTypesDF, is.infinite)] <- NA # remove infinite's
EcoTypesDF <- EcoTypesDF[complete.cases(EcoTypesDF),] # remove categories with no data

### Ecosystem Type Plot -------------------------------

# Select representative environments
Location1 <- EcoTypesDF[EcoTypesDF$Dmean == min(EcoTypesDF$Dmean),] # Seamounts
Location2 <- EcoTypesDF[EcoTypesDF$Rmean == min(EcoTypesDF$Rmean),] # Subglacial Lakes
Location3 <- EcoTypesDF[EcoTypesDF$Rmean == max(EcoTypesDF$Rmean),] # Oceanic temp rainforests
Location4 <- EcoTypesDF[EcoTypesDF$EFG == 'Polar alpine rock',] 
Location5 <- EcoTypesDF[EcoTypesDF$EFG == 'Temp alpine grasslands',] 
Location6 <- EcoTypesDF[EcoTypesDF$EFG == 'Rocky shores',]
Location7 <- EcoTypesDF[EcoTypesDF$EFG == 'Permanent marshes',]
Location8 <- EcoTypesDF[EcoTypesDF$EFG == 'Abyssal plains',]
Location9 <- EcoTypesDF[EcoTypesDF$EFG == 'Perm salt lakes',]

# Create plot
ggplot(EcoTypesDF, aes(x=Dmean, y=Rmean)) +
  geom_segment(aes(x = Dmean, y = Rmean, xend = Dmean, yend = 0.08), 
               data = Location1, color = "#2A788EFF", linewidth = 1.5, alpha = 1, linetype = 'dashed') + # Seamounts
  geom_segment(aes(x = Dmean, y = Rmean, xend = Dmean, yend = 0.0005), 
               data = Location2, color = "#414487FF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Subglacial Lakes
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.31, yend = 0.1), 
               data = Location3, color = "#BCAF6FFF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Oceanic Temperate Forests
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.36, yend = Rmean), 
               data = Location4, color = "#BCAF6FFF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Polar Alpine Rock
  geom_segment(aes(x = Dmean, y = Rmean, xend = Dmean, yend = 0.13), 
               data = Location5, color = "#BCAF6FFF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Temperate Alpine Grasslands
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.32, yend = 0.02), 
               data = Location6, color = "#440154FF", linewidth = 1.5, alpha = 0.9, linetype = 'dashed') + # Rocky shorelines
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.30, yend = 0.0003), 
               data = Location7, color = "#FDE725FF", linewidth = 1.5, alpha = 0.9, linetype = 'dashed') + # Permenant Marshland
  geom_segment(aes(x = Dmean, y = Rmean, xend = Dmean, yend = 0.0005), 
               data = Location8, color = "#2A788EFF", linewidth = 1.5, alpha = 0.9, linetype = 'dashed') + # Abyssal Plains
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.25, yend = Rmean), 
               data = Location9, color = "#414487FF", linewidth = 1.5, alpha = 0.9, linetype = 'dashed') + # Permenant Salt Lakes
  geom_point(aes(color = Realm_2), size = 7) +
  xlab('\nFractal Dimension') +
  ylab('Rugosity\n') +
  scale_y_continuous(trans = 'log10') +
  scale_color_manual(values = viridis(n = length(unique(EcoTypesDF$Realm_2)), option = 'cividis'),
                     guide = guide_legend(title = 'Realm', 
                                          reverse = F, label = T,
                                          na.value = "white")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     panel.border = element_blank(),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = 0)) +
  theme(axis.line = element_line(color = 'black'))

### Land Use Plot -------------------------------

# Create plot
ggplot(LandUseDF, aes(x=meanD, y=meanR)) +
  geom_segment(aes(x = 2.265062, y = 0.002939217, xend = 2.279, yend = 0.0048), 
               data = Location3, color = "#440154FF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Urban landscapes
  geom_segment(aes(x = 2.260332, y = 0.004467514, xend = 2.28, yend = 0.005), 
               data = Location3, color = "#481D6FFF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Cropland
  geom_point(aes(color = Cat_Name), size = 7) +
  xlab('\nFractal Dimension') +
  ylab('Rugosity\n') +
  scale_y_continuous(trans = 'log10') +
  scale_color_manual(values = viridis(n = length(unique(LandUseDF$Cat_Name)), option = 'cividis'),
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
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     panel.border = element_blank(),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = 0)) +
  theme(axis.line = element_line(color = 'black'))

# ----------------------------------------------------- End of Code ---------------------------------