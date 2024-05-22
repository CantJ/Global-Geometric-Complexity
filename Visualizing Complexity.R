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

# is this the first time running the below script? i.e. does the ecosystem typology zipped folder needed unpacking
FirstRun <- FALSE

# Set random number seed
set.seed(458967)

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
plot(DRast, col = magma(100, direction = -1), 
     buffer = FALSE,
     plg = list(x = 'bottom', at = c(2,3), digits = 1, tic = 'none', size = c(1,2.5)),
     box = FALSE, 
     axes = FALSE)

# Rugosity
plot(log10(RRast), col = cividis(100, direction = 1), 
     buffer = FALSE, 
     plg = list(x = 'bottom', at = c(-12,0), tic = 'none', size = c(1,2.5)),
     box = FALSE, 
     axes = FALSE)

# Height Range
plot(log10(HRast), col = cividis(100, direction = -1), 
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
# This process is faster if done in batch within QGIS. Before the reprojected rasters are then saved into the 'Ecotypes' folder ahead of the processing below.

###############################
# STEP 3: Identify Classification Categories (only needed if first time running script)
###############################

### Ecosystem types -----------------------------------

if(FirstRun == TRUE) {
  # List available ecosystem type maps
  fileNames <- list.files(EcoTypes, pattern = '.tif$')
  fileNames <- str_extract(fileNames, '[^.]*.[^.]*') # remove unnecessary file name details
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

  # regenerate file list for extraction process
  fileNames <- list.files(EcoTypes, pattern = '.tif$', full.names = TRUE)

  ### Land Use types -----------------------------------

  # Open land use raster
  LandUse <- rast(paste0(FilePath, 'LandCover2022_Mollweide_1870m.tif'))
  # List unique scenarios
  LU_cats <- unique((values(LandUse)))
  LU_cats <- LU_cats[!is.na(LU_cats)]
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
    ext(RastSelect) <- ext(tmpD)
  
    # Isolate corresponding complexity values
    tmpD[is.na(RastSelect)] <- NA
    tmpR[is.na(RastSelect)] <- NA
  
    # Extract values
    EcoTypesDF$Dmax[ii] <- max(values(tmpD),na.rm = TRUE)
    EcoTypesDF$Dmin[ii] <- min(values(tmpD),na.rm = TRUE)
    EcoTypesDF$Dmean[ii] <- mean(values(tmpD),na.rm = TRUE)
    EcoTypesDF$Dsigma[ii] <- sd(values(tmpD),na.rm = TRUE)
    EcoTypesDF$Rmax[ii] <- max(values(tmpR),na.rm = TRUE)
    EcoTypesDF$Rmin[ii] <- min(values(tmpR),na.rm = TRUE)
    EcoTypesDF$Rmean[ii] <- mean(values(tmpR),na.rm = TRUE)
    EcoTypesDF$Rsigma[ii] <- sd(values(tmpR),na.rm = TRUE)
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
  EcoTypesDF <- EcoTypesDF[complete.cases(EcoTypesDF),] # remove categories with no data
  
  ### 2. Land Use types

  for(ii in 1:dim(LU_df)[1]) { # working through each land use category
    # progress read out
    print(ii)

    # Duplicate Land use raster
    RastSelect <- LandUse
    # Isolate selected Land Use scenario
    RastSelect[values(RastSelect) != LU_df$Cat1[ii]] <- NA 
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
  }
  
  # Condense together values across repeated categories
  LandUseDF <- LU_df %>%
    group_by(Cat_Name) %>%
    summarize(Dmean = mean(Dmean, na.rm = TRUE),
              Rmean = mean(Rmean, na.rm = TRUE))
  
  # Save files (data checkpoint)
  write.csv(EcoTypesDF, paste0(FilePath, 'EcosystemComplexity.csv'), row.names = FALSE)
  write.csv(LandUseDF, paste0(FilePath, 'LandUseComplexity.csv'), row.names = FALSE)
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
Location4 <- EcoTypesDF[which(EcoTypesDF$Realm == 'T' & EcoTypesDF$Biome == 6 & EcoTypesDF$EFG == 2),] # Polar Alpine Cliffs 
Location5 <- EcoTypesDF[which(EcoTypesDF$Realm == 'MT' & EcoTypesDF$Biome == 1 & EcoTypesDF$EFG == 1),] # Rocky Shores
Location6 <- EcoTypesDF[which(EcoTypesDF$Realm == 'TF' & EcoTypesDF$Biome == 1 & EcoTypesDF$EFG == 3),] # Permanent marshes
Location7 <- EcoTypesDF[which(EcoTypesDF$Realm == 'M' & EcoTypesDF$Biome == 3 & EcoTypesDF$EFG == 3),] # Abyssal Plains
Location8 <- EcoTypesDF[which(EcoTypesDF$Realm == 'F' & EcoTypesDF$Biome == 2 & EcoTypesDF$EFG == 6),] # Permanent Salt Lakes

# Create plot
ggplot(EcoTypesDF, aes(x=Dmean, y=Rmean)) +
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.15, yend = 0.05), 
               data = Location1, color = "#2A788EFF", linewidth = 1.5, alpha = 1, linetype = 'dashed') + # Seamounts
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.245, yend = 0.00021), 
               data = Location2, color = "#414487FF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Subglacial Lakes
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.31, yend = 0.1), 
               data = Location3, color = "#BCAF6FFF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Oceanic Temperate Forests
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.33, yend = Rmean), 
               data = Location4, color = "#BCAF6FFF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Polar Alpine Rock
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.31, yend = 0.0095), 
               data = Location5, color = "#440154FF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Rocky shorelines
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.35, yend = 0.00015), 
               data = Location6, color = "#FDE725FF", linewidth = 1.5, alpha = 0.9, linetype = 'dashed') + # Permanent Marshes
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.15, yend = 0.001), 
               data = Location7, color = "#2A788EFF", linewidth = 1.5, alpha = 0.9, linetype = 'dashed') + # Abyssal Plains
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.15, yend = 0.0003), 
               data = Location8, color = "#414487FF", linewidth = 1.5, alpha = 0.9, linetype = 'dashed') + # Permanent Salt Lakes
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

# Identify the human mediated landscapes
Urban <- LandUseDF[LandUseDF$Cat_Name == 'Urban',] 
Cropland1 <- LandUseDF[LandUseDF$Cat_Name == 'Cropland',]
Cropland2 <- LandUseDF[LandUseDF$Cat_Name == 'Mosaic Vegetation/Cropland',]

# Create plot
ggplot(LandUseDF, aes(x = Dmean, y=Rmean)) +
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.250, yend = Rmean), 
               data = Urban, color = "#440154FF", linewidth = 1.5, alpha = 1, linetype = 'dashed') + # Urban
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.250, yend = Rmean), 
               data = Cropland1, color = "#481D6FFF", linewidth = 1.5, alpha = 1, linetype = 'dashed') + # Cropland
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.255, yend = 0.008), 
               data = Cropland2, color = "#481D6FFF", linewidth = 1.5, alpha = 1, linetype = 'dashed') + # Mosaic Vegetation/Cropland
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


######################################
# STEP 6: Analyse ecosystem complexity patterns 
######################################

# Implement brms regression to explore the relationship between rugosity and fractal dimension across ecosystem types

# Determine most appropriate model format by running linear and non-linear form
EcoTypesMod <- brm(log10(Rmean) ~ Dmean , data = EcoTypesDF, family = 'gaussian', chains = 20, iter = 5000, warmup = 1000)
EcoTypesMod2 <- brm(log10(Rmean) ~ Dmean + I(Dmean^2), data = EcoTypesDF, family = 'gaussian', chains = 20, iter = 5000, warmup = 1000)
model_weights(EcoTypesMod, EcoTypesMod2, weights = 'loo') # none-linear is the better fit

# Extract regression predicted lines from selected model
EcoTypesDF$Rpredict <- 10^(predict(EcoTypesMod2)[,1]) # mean fit
# confidence bounds
EcoTypesDF$Rlower <- 10^(predict(EcoTypesMod2)[,3]) # 2.5% 
EcoTypesDF$Rupper <- 10^(predict(EcoTypesMod2)[,4]) # 97.5% 

# Also predict how this model appears on the Land Use plot
LandUseDF$Rpredict <- 10^(predict(EcoTypesMod2, newdata = LandUseDF)[,1]) # mean fit only

# add regression lines to the ecosystem typology plot
ggplot(EcoTypesDF, aes(x=Dmean, y=Rmean)) +
  geom_ribbon(aes(ymin = Rlower, ymax = Rupper), fill = 'gray', alpha = 0.2) +
  geom_line(aes(y = Rpredict, x = Dmean), col = 'black', linetype = 'solid', linewidth = 2, alpha = 0.8) +
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.15, yend = 0.05), 
               data = Location1, color = "#2A788EFF", linewidth = 1.5, alpha = 1, linetype = 'dashed') + # Seamounts
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.245, yend = 0.00021), 
               data = Location2, color = "#414487FF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Subglacial Lakes
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.31, yend = 0.1), 
               data = Location3, color = "#BCAF6FFF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Oceanic Temperate Forests
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.33, yend = Rmean), 
               data = Location4, color = "#BCAF6FFF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Polar Alpine Rock
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.31, yend = 0.0095), 
               data = Location5, color = "#440154FF", linewidth = 1.5, alpha = 0.7, linetype = 'dashed') + # Rocky shorelines
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.35, yend = 0.00015), 
               data = Location6, color = "#FDE725FF", linewidth = 1.5, alpha = 0.9, linetype = 'dashed') + # Permanent Marshes
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.15, yend = 0.001), 
               data = Location7, color = "#2A788EFF", linewidth = 1.5, alpha = 0.9, linetype = 'dashed') + # Abyssal Plains
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.15, yend = 0.0003), 
               data = Location8, color = "#414487FF", linewidth = 1.5, alpha = 0.9, linetype = 'dashed') + # Permanent Salt Lakes
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

# Extract estimate of model fit
bayes_R2(EcoTypesMod2)

# Recreate Land Use plot with fit line
ggplot(LandUseDF, aes(x = Dmean, y=Rmean)) +
  geom_line(aes(x = Dmean, y = Rpredict), col = 'black', linetype = 'solid', linewidth = 1.5, alpha = 1) +
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.250, yend = Rmean), 
               data = Urban, color = "#440154FF", linewidth = 1.5, alpha = 1, linetype = 'dashed') + # Urban
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.250, yend = Rmean), 
               data = Cropland1, color = "#481D6FFF", linewidth = 1.5, alpha = 1, linetype = 'dashed') + # Cropland
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.255, yend = 0.008), 
               data = Cropland2, color = "#481D6FFF", linewidth = 1.5, alpha = 1, linetype = 'dashed') + # Mosaic Vegetation/Cropland
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

## 4. Great Barrier Reef
# Calculate mollweide coordinates from lat longs obtained from Google Earth
FeatureCoords4 <- matrix(c(147.6992, -18.2871), ncol = 2)
FeatureCoords4 <- spTransform(SpatialPoints(FeatureCoords4, CRS('+proj=longlat')), CRS(Moll))@coords
# Crop the complexity rasters around these selected points
Feature4D <- dem_crop(DRast, x0 = FeatureCoords4[1], y0 = FeatureCoords4[2], L = L*500, plot = FALSE)
Feature4R <- dem_crop(RRast, x0 = FeatureCoords4[1], y0 = FeatureCoords4[2], L = L*500, plot = FALSE)
# Estimate mean complexity values
mean(values(Feature4D), na.rm = TRUE)
mean(values(Feature4R), na.rm = TRUE)

## Plot Global Features

# Convert data into data-frames (required for plotting rasters using ggplot)
Feature1D_df <- as.data.frame(as(raster(Feature1D), "SpatialPixelsDataFrame"))
Feature1R_df <- as.data.frame(as(raster(Feature1R), "SpatialPixelsDataFrame"))
Feature2D_df <- as.data.frame(as(raster(Feature2D), "SpatialPixelsDataFrame"))
Feature2R_df <- as.data.frame(as(raster(Feature2R), "SpatialPixelsDataFrame"))
Feature3D_df <- as.data.frame(as(raster(Feature3D), "SpatialPixelsDataFrame"))
Feature3R_df <- as.data.frame(as(raster(Feature3R), "SpatialPixelsDataFrame"))
Feature4D_df <- as.data.frame(as(raster(Feature4D), "SpatialPixelsDataFrame"))
Feature4R_df <- as.data.frame(as(raster(Feature4R), "SpatialPixelsDataFrame"))
colnames(Feature1D_df) <- colnames(Feature1R_df) <- colnames(Feature2D_df) <- colnames(Feature2R_df) <-
  colnames(Feature3D_df) <- colnames(Feature3R_df) <- colnames(Feature4D_df) <- colnames(Feature4R_df) <- c("value", "x", "y")

# identify the maximum and minimum fractal dimension and rugosity estimates across these selected features (to keep plot colour scales consistent)
maxD <- ceiling(max(c(Feature1D_df$value, Feature2D_df$value, Feature3D_df$value, Feature4D_df$value),na.rm = T)*10)/10 # this little trick ensures the value is rounded up (at one decimal place)
maxR <- log10(ceiling(max(c(Feature1R_df$value, Feature2R_df$value, Feature3R_df$value, Feature4R_df$value),na.rm = T)*10)/10)
minD <- floor(min(c(Feature1D_df$value, Feature2D_df$value, Feature3D_df$value, Feature4D_df$value),na.rm = T)*10)/10
minR <- -12

# 1. Amazon Basin
# Fractal Dimension
ggplot(aes(x = x, y = y), data = Feature1D_df) +
  geom_raster(aes(fill = value)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradientn(colours = viridis(n = 100, option = 'magma', direction = -1),
                       limits = c(minD,maxD),
                       guide = guide_colorbar(ticks = F, title = 'D',
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_equal() +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.8,0.6),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# Rugosity
ggplot(aes(x = x, y = y), data = Feature1R_df) +
  geom_raster(aes(fill = log10(value))) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradientn(colours = viridis(n = 100, option = 'magma', direction = -1),
                       limits = c(minR, maxR),
                       guide = guide_colorbar(ticks = F, title = expression('log'[10]*'(R)'),
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_equal() +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.8,0.6),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# 2. Mariana Trench
# Fractal Dimension
ggplot(aes(x = x, y = y), data = Feature2D_df) +
  geom_raster(aes(fill = value)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradientn(colours = viridis(n = 100, option = 'magma', direction = -1),
                       limits = c(minD,maxD),
                       guide = guide_colorbar(ticks = F, title = 'D',
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_equal() +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.8,0.1),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# Rugosity
ggplot(aes(x = x, y = y), data = Feature2R_df) +
  geom_raster(aes(fill = log10(value))) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradientn(colours = viridis(n = 100, option = 'magma', direction = -1),
                       limits = c(minR, maxR),
                       guide = guide_colorbar(ticks = F, title = expression('log'[10]*'(R)'),
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_equal() +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.8,0.1),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# 3. Sahara Desert
# Fractal Dimension
ggplot(aes(x = x, y = y), data = Feature3D_df) +
  geom_raster(aes(fill = value)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradientn(colours = viridis(n = 100, option = 'magma', direction = -1),
                       limits = c(minD,maxD),
                       guide = guide_colorbar(ticks = F, title = 'D',
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_equal() +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.8,0.6),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# Rugosity
ggplot(aes(x = x, y = y), data = Feature3R_df) +
  geom_raster(aes(fill = log10(value))) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradientn(colours = viridis(n = 100, option = 'magma', direction = -1),
                       limits = c(minR, maxR),
                       guide = guide_colorbar(ticks = F, title = expression('log'[10]*'(R)'),
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_equal() +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.8,0.6),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# 4 Great Barrier Reef
# Fractal Dimension
ggplot(aes(x = x, y = y), data = Feature4D_df) +
  geom_raster(aes(fill = value)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradientn(colours = viridis(n = 100, option = 'magma', direction = -1),
                       limits = c(minD,maxD),
                       guide = guide_colorbar(ticks = F, title = 'D',
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_equal() +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.8,0.6),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# Rugosity
ggplot(aes(x = x, y = y), data = Feature4R_df) +
  geom_raster(aes(fill = log10(value))) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradientn(colours = viridis(n = 100, option = 'magma', direction = -1),
                       limits = c(minR, maxR),
                       guide = guide_colorbar(ticks = F, title = expression('log'[10]*'(R)'),
                                              reverse = F, label = T,
                                              na.value = "white")) +
  coord_equal() +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(0, 1, 0.2, 0.2, "cm"),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.justification=c(0,0), legend.position=c(0.8,0.6),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = "black"))

# ----------------------------------------------------- End of Code ---------------------------------