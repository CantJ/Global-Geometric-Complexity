# This script is for statistically testing the similarities and differences between the structural complexity regimes
# associated with differing ecosystem types and land use scenarios

# Date last modified: June 2026
# Primary Author: James Cant
# -----------------------------------------------------------------------------------------

# clear working directory
rm(list = ls())

# load required files
library(terra)
library(raster)
library(stringr)
library(dplyr)
library(sf)
library(vegan)

# Set random number seed
set.seed(458967)
# prevent rounding of small values
options(digits = 22)
# First time running script?
FirstRun <- TRUE
# Set up file save directory
FilePath <- 'E:/Complexity/Visualising Outputs/'

###############################
# STEP 1: Identify Classification Categories (only needed if first time running script)
###############################

# Define file pathways for the desired Land use and ecosystem typology files
EcoTypes <- 'E:/Complexity/Ecotypes/'
ComplexRast <- 'E:/Complexity/Complexity Analyses/Final Data Products/'


if(FirstRun == TRUE) {
  ### Ecosystem types ----------------------------------- 
  # List available ecosystem type maps
  fileNames <- list.files(EcoTypes, pattern = '.tif$', full.names = TRUE)
  EFGs <- str_extract(fileNames, '[^.]*.[^.]*') # remove unnecessary file name details
  EFGs <- sub("_", " ", EFGs)
  EFGs <- sub("_", " ", EFGs) # remove second underscore
  EFGs <- word(EFGs, 3)
  # Extract necessary file metadata
  EcoMetadata <- data.frame(Realm = gsub("[^a-zA-Z]", "", EFGs),
                            Biome = gsub(".*?([0-9]+).*", "\\1", EFGs),
                            EFG = sapply(strsplit(EFGs, '\\D+'),'[',3))
  # match up EFG codes to corresponding ecosystem names
  EFGcodes <- read.csv(paste0(FilePath, 'EFGcodes_names.csv'))
  EcoMetadata <- merge(EcoMetadata, EFGcodes, by = c('Realm', 'Biome', 'EFG'))
  # Add appropriate realm name classification
  EcoMetadata$Realm_Name <- NA
  EcoMetadata[EcoMetadata$Realm %in% c('F'),]$Realm_Name <- 'Freshwater'
  EcoMetadata[EcoMetadata$Realm %in% c('MFT','MT'),]$Realm_Name <- 'Coastal'
  EcoMetadata[EcoMetadata$Realm %in% c('S','SF','SM'),]$Realm_Name <- 'Subterranean'
  EcoMetadata[EcoMetadata$Realm %in% c('T'),]$Realm_Name <- 'Terrestrial'
  EcoMetadata[EcoMetadata$Realm %in% c('TF'),]$Realm_Name <- 'Wetland'
  EcoMetadata[EcoMetadata$Realm %in% c('M','FM'),]$Realm_Name <- 'Marine'
  
  # Output storage for distributions of complexity measures associated with each ecosystem type
  EcoTypesD <- EcoTypesR <- EcoTypesH <- list()
  
  ### Land Use types -----------------------------------
  # Open land use raster
  LandUse <- rast(paste0(ComplexRast, 'LandCover2022_Mollweide_1870m.tif'))
  # List unique scenario codes
  LU_cats <- unique((values(LandUse)))
  LU_cats <- LU_cats[!is.na(LU_cats)]
  # List unique scenario names
  LandUseMetadata <- c('Cropland', 'Mosaic Vegetation/Cropland', 'Tree Cover (Broadleaved)', 'Tree Cover (Needleleaved)',
                       'Mixed Tree Cover', 'Mosaic Vegetation', 'Shrubland', 'Grassland', 'Lichens & Mosses', 'Sparse Vegetation',
                       'Wetlands', 'Urban', 'Bare substrate', 'Permenant Snow & Ice')
  
  # Output storage for density plots
  LandUseD <- LandUseR <- LandUseH <- list()
  
  ###############################
  # STEP 4: Extract Complexities (again only needed if first time running script)
  ###############################
  
  ### 1. Ecosystem Types
  # Working through each ecosystem type, overlay ecosystem map onto complexity map and extract corresponding complexity details 
  
  for(ii in 1:dim(EcoMetadata)[1]) { # working through each ecosystem typology file
    # progress read out
    print(ii)
    
    # Open selected ecosystem raster
    RastSelect <- rast(fileNames[ii])
    # Remove minor occurrences of selected ecosystem
    RastSelect[RastSelect != 1] <- NA
    # re-load required complexity rasters
    tmpD <- rast(paste0(ComplexRast, 'GlobalFractalDimension.tif'))
    tmpR <- rast(paste0(ComplexRast, 'GlobalRugosity.tif'))
    tmpH <- rast(paste0(ComplexRast, 'GlobalHeightRange.tif'))
    # Ensure rasters align
    RastSelect <- resample(RastSelect, tmpR)
    
    # Isolate corresponding complexity values
    tmpD[is.na(RastSelect)] <- NA
    tmpR[is.na(RastSelect)] <- NA
    tmpH[is.na(RastSelect)] <- NA
    tmpR <- log10(tmpR) # Apply data transformation to rugosity and height range values
    tmpH <- log10(tmpH)
    tmpR[is.infinite(tmpR)] <- NA
    tmpH[is.infinite(tmpH)] <- NA
    
    # Extract distribution of complexity values
    try(EcoTypesD[[ii]] <- na.omit(values(tmpD)))
    try(EcoTypesR[[ii]] <- na.omit(values(tmpR)))
    try(EcoTypesH[[ii]] <- na.omit(values(tmpH)))
  }
  
  ### Worked up to here -----------------------------

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
    try(LUDensityD[[ii]] <- density(na.omit(values(tmpD))))
    tmpR <- log10(tmpR) # Apply data transformation
    tmpR[is.infinite(tmpR)] <- NA
    try(LUDensityR[[ii]] <- density(na.omit(values(tmpR))))
  }
  
  # Save files (data checkpoint)
  # Ecosystem types
  write.csv(EcoTypesDF, paste0(FilePath, 'EcosystemComplexity.csv'), row.names = FALSE)
  saveRDS(EcoTypesDensityD, paste0(FilePath, 'EcoTypesDensityD.RData'))
  saveRDS(EcoTypesDensityR, paste0(FilePath, 'EcoTypesDensityR.RData'))
  # Land Use types
  write.csv(LU_df, paste0(FilePath, 'LandUseComplexity.csv'), row.names = FALSE)
  saveRDS(LUDensityD, paste0(FilePath, 'LandUseDensityD.RData'))
  saveRDS(LUDensityR, paste0(FilePath, 'LandUseDensityR.RData'))
}