# This script is for statistically testing the similarities and differences between the structural complexity regimes
# associated with differing ecosystem types and land use scenarios

# Date last modified: June 2026
# Primary Author: James Cant
# -----------------------------------------------------------------------------------------

# Set random number seed
set.seed(458967)

# Load required complexity rasters
DRast <- rast(paste0(FilePath, 'GlobalFractalDimension.tif'))
RRast <- rast(paste0(FilePath, 'GlobalRugosity.tif'))
HRast <- rast(paste0(FilePath, 'GlobalHeightRange.tif'))
# Apply data transformation to rugosity and height range values
RRast <- log10(RRast)
HRast <- log10(HRast)
# Remove error values
RRast[is.infinite(RRast)] <- NA
HRast[is.infinite(HRast)] <- NA

###############################
# STEP 1: Identify Classification Categories (only needed if first time running script)
###############################

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
  LandUse <- rast(paste0(FilePath, 'LandCover2022_Mollweide_1870m.tif'))
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
    try(RastSelect[RastSelect != 1] <- NA)
    # Ensure ecosystem raster aligns with complexity rasters
    try(RastSelect <- resample(RastSelect, RRast))
    # Vectorise ecosystem distribution
    try(RastSelect <- as.polygons(RastSelect))
    # Extract distribution of complexity values
    try(EcoTypesD[[ii]] <- na.omit(extract(DRast, RastSelect)[,2]))
    try(EcoTypesR[[ii]] <- na.omit(extract(RRast, RastSelect)[,2]))
    try(EcoTypesH[[ii]] <- na.omit(extract(HRast, RastSelect)[,2]))
  }

  # match up EFG codes to corresponding ecosystem names
  write.csv(EcoMetadata, paste0(FilePath, 'EFGcodes.csv'), row.names = F)
  EFGcodes <- read.csv(paste0(FilePath, 'EFGcodes_names.csv'))
  EcoMetadata <- merge(EcoMetadata, EFGcodes, by = c('Realm', 'Biome', 'EFG', 'Realm_Name'))
  
  ### 2. Land Use types
  
  for(ii in 1:length(LandUseMetadata)) { # working through each land use category
    # progress read out
    print(ii)
    
    # Identify scenario codes associated with the selected category
    if(LandUseMetadata[ii] == 'Cropland'){ CatCode <- c(10,11,12,20) }
    if(LandUseMetadata[ii] == 'Mosaic Vegetation/Cropland'){ CatCode <-  c(30,40) }
    if(LandUseMetadata[ii] == 'Tree Cover (Broadleaved)'){ CatCode <-  c(50:62) }
    if(LandUseMetadata[ii] == 'Tree Cover (Needleleaved)'){ CatCode <-  c(70:82) }
    if(LandUseMetadata[ii] == 'Mixed Tree Cover'){ CatCode <-  c(90) }
    if(LandUseMetadata[ii] == 'Mosaic Vegetation'){ CatCode <-  c(100:110) }
    if(LandUseMetadata[ii] == 'Shrubland'){ CatCode <-  c(120:122) }
    if(LandUseMetadata[ii] == 'Grassland'){ CatCode <-  c(130) }
    if(LandUseMetadata[ii] == 'Lichens & Mosses'){ CatCode <-  c(140) }
    if(LandUseMetadata[ii] == 'Sparse Vegetation'){ CatCode <-  c(150:153) }
    if(LandUseMetadata[ii] == 'Wetlands'){ CatCode <-  c(160:180) }
    if(LandUseMetadata[ii] == 'Urban'){ CatCode <-  c(190) }
    if(LandUseMetadata[ii] == 'Bare substrate'){ CatCode <-  c(200:202) }
    if(LandUseMetadata[ii] == 'Permenant Snow & Ice'){ CatCode <-  c(220) }
    
    # Duplicate Land use raster
    RastSelect <- LandUse
    # Isolate selected Land Use scenario
    RastSelect[!(values(RastSelect) %in% CatCode)] <- NA 
    # Ensure land use and complexity rasters align
    ext(RastSelect) <- ext(tmpD)
    # Vectorise lselected land use distribution
    try(RastSelect <- as.polygons(RastSelect))
    # Extract distribution of complexity values
    try(LandUseD[[ii]] <- na.omit(extract(DRast, RastSelect)[,2]))
    try(LandUseR[[ii]] <- na.omit(extract(RRast, RastSelect)[,2]))
    try(LandUseH[[ii]] <- na.omit(extract(HRast, RastSelect)[,2]))
  }
  
  # Save files (data checkpoint)
  # Ecosystem types
  write.csv(EcoMetadata, paste0(FilePath, 'EcosystemComplexity.csv'), row.names = FALSE)
  saveRDS(EcoTypesD, paste0(FilePath, 'EcoTypesD.RData'))
  saveRDS(EcoTypesR, paste0(FilePath, 'EcoTypesR.RData'))
  saveRDS(EcoTypesH, paste0(FilePath, 'EcoTypesH.RData'))
  # Land Use types
  saveRDS(LandUseH, paste0(FilePath, 'LandUseH.RData'))
  saveRDS(LandUseD, paste0(FilePath, 'LandUseD.RData'))
  saveRDS(LandUseR, paste0(FilePath, 'LandUseR.RData'))
}
