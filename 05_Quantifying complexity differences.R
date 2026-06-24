# This script is for statistically testing and visualising global patterns in geometric complexity 
# and how different ecosystem types and land use scenarios align with these measures
# The ecosystem classification used follows the typology outlined in Keith et al. (2022) Nature.
# The Land Use classifications have been derived from the Land Cover 2020 data product from the Copernicus Climate Change service

# Date last modified: June 2026
# Primary Author: James Cant
# -----------------------------------------------------------------------------------------


# Set random number seed
set.seed(458967)

# Load required complexity rasters
DRast <- rast(paste0(FilePath, 'GlobalFractalDimension.tif'))
RRast <- rast(paste0(FilePath, 'GlobalRugosity.tif'))
HRast <- rast(paste0(FilePath, 'GlobalHeightRange.tif'))

###############################
# STEP 1: Plot global patterns
###############################

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
# STEP 2: Extract complexities for ecosystem and land Use types
###############################

# Extract typology maps from downloaded typology .tar file
if(FirstRun == TRUE) { untar(paste0(FilePath, "all-maps-raster-geotiff.tar.bz2"), exdir = EcoTypes) }
# only needed once to unzip downloaded maps.
# Note: once these ecosystem typology maps have been extracted then need to be reprojected to match the mollweide projection, resolution, and extent of the complexity rasters.
# This process is faster if done in batch within QGIS, first using the warp function to change the desired resolution and CRS of the IUCN maps before clipping these reprojected maps to the same extent as the complexity maps.
# The reprojected rasters are then saved into the 'Ecotypes' folder ahead of the processing below.

if(FirstRun == TRUE) {
  # Apply data transformation to rugosity and height range values
  RRast <- log10(RRast)
  HRast <- log10(HRast)
  # Remove error values
  RRast[is.infinite(RRast)] <- NA
  HRast[is.infinite(HRast)] <- NA
  
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
  # STEP 3: Extract Complexities
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
  
  # Reformat complexity data
  # Ecosystem typologies
  EcoTypeDat <- as.data.frame(do.call(rbind, lapply(1:dim(EcoMetadata)[1], function(x){
    # Isolate complexity measures for selected land use category
    vecD <- as.numeric(EcoTypesD[[x]])
    vecR <- as.numeric(EcoTypesR[[x]])
    vecH <- as.numeric(EcoTypesH[[x]])
    # Determine maximum vector length
    max_length <- max(length(vecD), length(vecR), length(vecH))
    # Set length of each vector equal to max length
    length(vecD) <- max_length                      
    length(vecR) <- max_length 
    length(vecH) <- max_length
    # Assign Land use category
    Index <- rep(EcoMetadata$Ecosystem_Name[x], length.out = max_length)
    # cbind the vectors together
    return(do.call(cbind, list(Index, vecD, vecR, vecH)))
  })))
  names(EcoTypeDat) <- c('EcosystemType', 'D', 'R', 'H')
  
  # LandUse Data
  LandUseDat <- as.data.frame(do.call(rbind, lapply(1:length(LandUseMetadata), function(x){
    # Isolate complexity measures for selected land use category
    vecD <- as.numeric(LandUseD[[x]])
    vecR <- as.numeric(LandUseR[[x]])
    vecH <- as.numeric(LandUseH[[x]])
    # Determine maximum vector length
    max_length <- max(length(vecD), length(vecR), length(vecH))
    # Set length of each vector equal to max length
    length(vecD) <- max_length                      
    length(vecR) <- max_length 
    length(vecH) <- max_length
    # Assign Land use category
    Index <- rep(LandUseMetadata[x], length.out = max_length)
    # cbind the vectors together
    return(do.call(cbind, list(Index, vecD, vecR, vecH)))
  })))
  names(LandUseDat) <- c('LandUse', 'D', 'R', 'H')
  
  #EcoTypeDat <- readRDS(paste0(FilePath, 'EcoTypeData.RData'))
  #LandUseDat <- readRDS(paste0(FilePath, 'LandUseData.RData'))
  
  # Ensure correct variable formatting
  EcoTypeDat$EcosystemType <- as.factor(EcoTypeDat$EcosystemType)
  EcoTypeDat$D <- as.numeric(EcoTypeDat$D)
  EcoTypeDat$R <- as.numeric(EcoTypeDat$R)
  EcoTypeDat$H <- as.numeric(EcoTypeDat$H)
  LandUseDat$LandUse <- as.factor(LandUseDat$LandUse)
  LandUseDat$D <- as.numeric(LandUseDat$D)
  LandUseDat$R <- as.numeric(LandUseDat$R)
  LandUseDat$H <- as.numeric(LandUseDat$H)

  # Drop complete entries 
  EcoTypeDat <- EcoTypeDat[complete.cases(EcoTypeDat[,c('D', 'R')]),]
  LandUseDat <- LandUseDat[complete.cases(LandUseDat[,c('D', 'R')]),]
  
  # Save files (data checkpoint)
  write_parquet(EcoTypeDat, paste0(FilePath, 'EcoTypeData.parquet'))
  write_parquet(LandUseDat, paste0(FilePath, 'LandUseData.parquet'))
} else {
  # Load files
  EcoTypeDat <- open_dataset(paste0(FilePath, 'EcoTypeData.parquet')) 
  LandUseDat <- open_dataset(paste0(FilePath, 'LandUseData.parquet'))
}


###############################
# STEP 4: Quantify and visualize complexity patterns
#############################

# Determine how structural complexity predicts ecosystem and land use types
## Ecosystem types
# Isolate data samples to not overwhelm computational system
EcoDatSample <- EcoTypeDat %>%
  group_by(EcosystemType) %>% 
  map_batches(~ as_record_batch(sample_frac(as.data.frame(.), 0.001))) %>%   
  collect()                                 # Pull only this tiny slice into RAM

# Removed unused factor levels
EcoDatSample$EcosystemType <- droplevels(EcoDatSample$EcosystemType)

# Fit multinomial logistic model
EcoMod <- multinom(EcosystemType ~ D + R, data = EcoDatSample)
# Explore effect of dropping terms
EcoModD <- multinom(EcosystemType ~ D, data = EcoDatSample)
EcoModR <- multinom(EcosystemType ~ R, data = EcoDatSample)
# Repeat model without structural metrics to generate baseline null model.
null_EcoMod <- multinom(EcosystemType ~ 1, data = EcoDatSample)
# Use Likelihood Ratio Test to determine if structural metrics explain a significant amount of variance compared to null model
lmtest::lrtest(null_EcoMod, EcoModD, EcoModR, EcoMod)
AIC(EcoModD, EcoModR, EcoMod)

## Land Use types
# Isolate data samples to not overwhelm computational system
LUDatSample <- LandUseDat %>%
  group_by(LandUse) %>% 
  map_batches(~ as_record_batch(sample_frac(as.data.frame(.), 0.001))) %>%   
  collect()                                 # Pull only this tiny slice into RAM

# Removed unused factor levels
LUDatSample$LandUse <- droplevels(LUDatSample$LandUse)

# Fit multinomial logistic model
LUMod <- multinom(LandUse ~ D + R, data = LUDatSample)
# Explore effect of dropping terms
LUModD <- multinom(LandUse ~ D, data = LUDatSample)
LUModR <- multinom(LandUse ~ R, data = LUDatSample)
# Repeat model without structural metrics to generate baseline null model.
null_LUMod <- multinom(LandUse ~ 1, data = LUDatSample)
# Use Likelihood Ratio Test to determine if structural metrics explain a significant amount of variance compared to null model
lmtest::lrtest(null_LUMod, LUModD, LUModR, LUMod)
AIC(LUModD, LUModR, LUMod)

# Extract marginal effects slope coefficients (probability) to visualise effects of structural metrics on Ecosystem and Land use types
LU_b <- avg_slopes(LUMod, variables = c('D', 'R'))
Eco_b <- avg_slopes(EcoMod, variables = c('D', 'R'))


### -------------------------------------- WORKS TO HERE 


# Extract post-hoc summary statistics for plotting.
# Ecosystem Types
df %>% group_by(plant_var) %>%  summarise(n = n(), mean = mean(canopy_vol), sd = sd(canopy_vol))
# Land Use
LU_lda <- data.frame(LandUse = LandUseDat[, "LandUse"], lda = predict(LandUsePostHoc)$x)
plotLU <- LU_lda %>% group_by(LandUse) %>% summarise(n = n(), meanLD1 = mean(lda.LD1), meanLD2 = mean(lda.LD2))

# Plot centroid positions of Ecosystem and Land use types on their respective LDA planes
## Ecosystem types


## Land Use
ggplot(aes(x = meanLD1, y = meanLD2), data = plotLU) + 
  geom_point(aes(colour = LandUse), size = 6) +
  scale_color_manual(values = magma(n = length(unique(plotLU$LandUse))),
                     limits = c("Urban","Cropland","Mosaic Vegetation/Cropland","Mosaic Vegetation",
                                "Grassland","Shrubland","Lichens & Mosses","Sparse Vegetation",         
                                "Mixed Tree Cover","Tree Cover (Needleleaved)","Tree Cover (Broadleaved)",
                                'Wetlands',"Bare substrate","Permenant Snow & Ice"),
                     guide = guide_legend(title = 'Land Cover Type', 
                                          reverse = F, label = T)) +
  scale_x_continuous(labels = function(x) { format(x, digits = 2) }) +
  scale_y_continuous(labels = function(x) { format(x, digits = 2) }) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"),
                     panel.border = element_blank(),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = 0),
                     axis.line = element_line(color = 'black'),
                     plot.margin = margin(15, 10, 10, 10))
  

# Extract densities of complexity distributions
