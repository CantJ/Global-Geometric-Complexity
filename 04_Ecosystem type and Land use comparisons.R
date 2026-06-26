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
# STEP 4: Quantify complexity patterns
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


###############################
# STEP 5: Visualize complexity patterns
###############################

## Isolate mean ecosystem complexities to not overwhelm computational system
EcoDatSample <- EcoTypeDat %>%
  group_by(EcosystemType) %>% 
  summarise(Dmean = mean(D),
            Rmean = mean(R)) %>%   
  collect()
# Add realm variable back into data
EFGcodes <- read.csv(paste0(FilePath, 'EFGcodes_names.csv'))
EcoDatSample <- merge(EcoDatSample, EFGcodes[,c('Ecosystem_Name', 'Realm_Name')], by.x = 'EcosystemType', by.y = 'Ecosystem_Name')
names(EcoDatSample) <- c("EcosystemType", "Dmean", "Rmean", 'Realm')
# Reformat variables to aid visualization clarity
EcoDatSample$Realm <- factor(EcoDatSample$Realm,
                        levels = c('Marine', 'Subterranean', 'Freshwater', 'Terrestrial', 'Coastal', 'Wetland'))

### Mean Ecosystem complexity plot -------------------------------
# Select representative environments
Location1 <- EcoDatSample[which(EcoDatSample$EcosystemType == 'Seamounts, ridges & plateaus'),]
Location2 <- EcoDatSample[which(EcoDatSample$EcosystemType == 'Permenant marshes'),]
Location3 <- EcoDatSample[which(EcoDatSample$EcosystemType == 'Oceanic cool temperate rainforests'),]
Location4 <- EcoDatSample[which(EcoDatSample$EcosystemType == 'Polar tundra'),]
Location5 <- EcoDatSample[which(EcoDatSample$EcosystemType == 'Coastal river deltas'),] 
Location6 <- EcoDatSample[which(EcoDatSample$EcosystemType == 'Permenant salt & soda lakes'),]
Location7 <- EcoDatSample[which(EcoDatSample$EcosystemType == 'Sclerophyll hot deserts & semi-deserts'),]

# Create mean complexity plot
ggplot(EcoDatSample, aes(x=Dmean, y=Rmean)) +
  annotate('rect', xmin = 2.225, xmax = 2.35, ymin = -2.85, ymax = -0.5, fill = 'grey', alpha = 0.5) +
  annotate('rect', xmin = 2.1, xmax = 2.225, ymin = -5.2, ymax = -2.85, fill = 'grey', alpha = 0.2) +
  geom_hline(aes(yintercept = -2.85), linetype = "dashed", col = 'gray50', linewidth = 1) +
  geom_vline(aes(xintercept = 2.225), linetype = "dashed", col = 'gray50', linewidth = 1) +
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.13, yend = -1.8), 
               data = Location1, color = "#00204DFF", linewidth = 1.5, alpha = 0.7, linetype = 'solid') + # Seamounts
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.33, yend = -3.7), 
               data = Location2, color = "#FDE725FF", linewidth = 1.5, alpha = 0.7, linetype = 'solid') + # Subglacial Lakes
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.28, yend = Rmean), 
               data = Location3, color = "#5F9258FF", linewidth = 1.5, alpha = 0.7, linetype = 'solid') + # Oceanic Temperate Forests
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.28, yend = Rmean), 
               data = Location4, color = "#5F9258FF", linewidth = 1.5, alpha = 0.7, linetype = 'solid') + # Polar Tundra
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.35, yend = -4.5), 
               data = Location5, color = "#CBBA69FF", linewidth = 1.5, alpha = 0.7, linetype = 'solid') + # Coastal River Deltas
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.14, yend = -4.9), 
               data = Location6, color = "lightblue", linewidth = 1.5, alpha = 0.9, linetype = 'solid') + # Permanent Salt Lakes
  geom_segment(aes(x = Dmean, y = Rmean, xend = 2.16, yend = -4.4), 
               data = Location7, color = "#5F9258FF", linewidth = 1.5, alpha = 0.7, linetype = 'solid') + # Sclerophyll hot deserts and semi-deserts
  geom_point(aes(color = Realm), size = 8) +
  xlab(NULL) + # Fractal Dimension
  ylab(NULL) + # Rugosity
  scale_y_continuous(limits = c(-5.2,-0.5), breaks = c(-4.9,-2.85,-0.8), labels = function(i) { format(exp(i), scientific = F, digits = 3) }, expand = c(0,0)) +
  scale_x_continuous(limits = c(2.10,2.35), labels = function(x) { format(x, digits = 3) }, expand = c(0,0)) +
  scale_color_manual(values = c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF"),
                     guide = guide_legend(title = 'Realm', 
                                          reverse = F, label = T)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_blank(),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(colour = 0)) +
  theme(axis.line = element_line(color = 'black'))

### Complexity distributions across ecosystems and biomes plots -------------------------------

# Compute density distributions across ecosystem types
# Rugosity
EcoRugosity <- EcoTypeDat %>%
  dplyr::select(EcosystemType) %>%
  distinct() %>%
  collect()
# Add Biome and Realm Data
EcoRugosity <- merge(setDT(EcoRugosity), setDT(EFGcodes[,c('Ecosystem_Name', 'Realm_Name', 'Biome_Name')]), by.x = 'EcosystemType', by.y = 'Ecosystem_Name')
# Reformat variables to aid visualization clarity
names(EcoRugosity) <- c('EcosystemType', 'Realm', 'Biome')
EcoRugosity$Realm <- factor(EcoRugosity$Realm,
                             levels = rev(c('Marine', 'Subterranean', 'Freshwater', 'Terrestrial', 'Coastal', 'Wetland')))
EcoRugosity$Biome <- factor(EcoRugosity$Biome, levels = rev(c('Continental shelf','Pelagic & deep sea','Transitional inlets & bays','Anthropogenic marine', # Marine
                                                              'Subterranean caves & pools', # Subterranean
                                                              'Rivers & streams','Lakes','Artificial wetland', # Freshwater
                                                              'Tropical-subtropical forest','Temperate forests & woodland','Shrublands & shrubby woodland','Savanna & grassland','Desert & semi-desert','Polar/alpine','Urban & Farmland', # Terrestrial
                                                              'Brackish tidal systems','Natural intertidal shoreline','Supralittoral coastal systems','Anthropogenic shorelines', # Coastal
                                                              'Wetland'))) # Wetlands
EcoRugosity$EcosystemType <- factor(EcoRugosity$EcosystemType, levels = rev(c(
                                                              # Marine
                                                              'Seagrass meadows','Rhodolith/Maerl beds','Kelp forests','Photic coral reefs','Shellfish beds & reefs','Upwelling zones','Epipelagic ocean waters','Mesopelagic ocean waters',
                                                              'Bathypelagic ocean waters','Abyssopelagic ocean waters','Sea ice','Continental & island slopes','Submarine canyons',
                                                              'Abyssal plains','Seamounts, ridges & plateaus','Hadal trenches & troughs','Chemosynthetic-based-ecosystems','Deepwater coastal inlets','Permenantly open riverine estuaries and bays',
                                                              'Intermittently closed and open lakes & lagoons','Marine aquafarms',
                                                              # Subterranean
                                                              'Aerobic caves','Underground streams & pools','Groundwater ecosystems','Anchialine caves','Anchialine pools',
                                                              # Freshwater
                                                              'Permenant upland streams','Permenant lowland rivers','Freeze-thaw rivers & streams','Seasonal upland streams','Seasonal lowland rivers','Eposodic arid rivers','Large lowland rivers','Large permenant freshwater lakes','Small permenant freshwater lakes','Seasonal freshwater lakes',
                                                              'Freeze-thaw freshwater lakes','Ephemeral freshwater lakes','Permenant salt & soda lakes','Ephemeral salt lakes',
                                                              'Artesian springs & oases','Geothermal pools & wetlands','Large reservoirs','Constructed lacustrine wetlands','Rice paddies','Canals, ditches & drains',
                                                              # Terrestrial
                                                              'Tropical/Subtropical lowland rainforests','Tropical/Subtropical dry forests & thickets','Tropical/Subtropical montane forests','Tropical heath forests',
                                                              'Boreal and temperate high montane forests & woodlands','Deciduous temperate forests','Oceanic cool temperate rainforests','Warm temperate laurophyll forests','Temperate pyric humid forests','Temperate pyric sclerophyll forests & woodlands',
                                                              'Seasonally dry tropical shrublands','Seasonally dry temperate heath & shrublands','Cool temperate heathlands','Young rocky pavements, lava flows & screes','Trophic savannas','Pyric tussock savannas',
                                                              'Hummock savannas','Temperate savannas','Temperate subhumid grasslands','Semi-desert steppe',
                                                              'Succulent/Thorny deserts & semi-deserts','Sclerophyll hot deserts & semi-deserts','Cool deserts & semi-deserts','Hyper-arid deserts','Ice sheets, glaciers & perennial snowfields','Polar/alpine cliffs, screes, outcrops & lava flows','Polar tundra','Temperate alpine grasslands & shrublands',
                                                              'Tropical alpine grasslands & herbfields','Annual croplands','Sown pastures & fields','Plantations','Urban & industrial ecosystems','Derived semi-natural pastures & old fields',
                                                              # Coastal
                                                              'Coastal river deltas','Intertidal forests & shrublands','Coastal saltmarshes & reedbeds','Rocky shorelines','Muddy shorelines','Sandy shorelines','Boulder & cobble shores',
                                                              'Coastal shrublands & grasslands','Large seabird & pinniped colonies','Artificial shorelines',
                                                              # Wetlands
                                                              'Tropical flooded forests & peat forests','Permenant marshes',
                                                              'Seasonal floodplain marshes','Episodic arid floodplains','Boreal, temperate & montane peat bogs','Boreal & temperate fens')))
# Duplicate output for Fractal Dimension also
EcoFD <- EcoRugosity
# set additional output storage
EcoRugosity$xx <- EcoRugosity$yy <- EcoFD$xx <- EcoFD$yy <- I(list())

# compute density distributions for complexity estimates associated with each ecosystem typology.
for(ii in seq_along(EcoRugosity$EcosystemType)){
  print(ii) # progress read out
  # Isolate Rugosity values for each ecosystem type 
  TmpR <- EcoTypeDat %>% 
    filter(EcosystemType == EcoRugosity$EcosystemType[ii]) %>%
    dplyr::select(R) %>%
    # Pull the raw values into memory (ensuring data fits in RAM)
    collect() %>% 
    reframe(# standard density returns an object; we wrap it in a data frame
      as.data.frame(density(R)[c("x", "y")])) %>%
    rename(value = x, density = y)
  # Repeat for fractal dimension 
  TmpD <- EcoTypeDat %>% 
    filter(EcosystemType == EcoRugosity$EcosystemType[ii]) %>%
    dplyr::select(D) %>%
    # Pull the raw values into memory (ensuring data fits in RAM)
    collect() %>% 
    reframe(# standard density returns an object; we wrap it in a data frame
      as.data.frame(density(D)[c("x", "y")])) %>%
    rename(value = x, density = y)
  # Save desired outputs
  # Rugosity
  EcoRugosity$xx[[ii]] <- TmpR$value
  EcoRugosity$yy[[ii]] <- TmpR$density
  # Fractal dimension
  EcoFD$xx[[ii]] <- TmpD$value
  EcoFD$yy[[ii]] <- TmpD$density
} # close loop

# Reshape to a long format
EcoRugosity <- data.frame(Realm = rep(EcoRugosity$Realm, each = 512),
                   Biome = rep(EcoRugosity$Biome, each = 512),
                   EcosystemType = rep(EcoRugosity$EcosystemType, each = 512),
                   xx = unlist(EcoRugosity$xx),
                   yy = unlist(EcoRugosity$yy))
EcoFD <- data.frame(Realm = rep(EcoFD$Realm, each = 512),
                          Biome = rep(EcoFD$Biome, each = 512),
                          EcosystemType = rep(EcoFD$EcosystemType, each = 512),
                          xx = unlist(EcoFD$xx),
                          yy = unlist(EcoFD$yy))

# Generate density plots (Realm level)
# Rugosity
ggplot(EcoRugosity) +
  geom_ridgeline(aes(x = xx, y = Realm, height = yy, fill = Realm, colour = Realm, scale = 1)) +
  scale_x_continuous(breaks = c(-12, -6, 0), labels = function(i) { format(exp(i), scientific = F, digits = 1) }) +
  scale_fill_manual(values = rev(c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF"))) +
  scale_color_manual(values = rev(c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF"))) +
  xlab(NULL) +
  ylab(NULL) +
  theme_ridges() +
  theme(axis.text.x = element_text(size = 20, colour = "black"), 
        axis.text.y = element_text(size = 20, colour = "black"),
        panel.border = element_rect(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 0.9),
        legend.position = 'none',
        plot.margin = margin(5,15,5,5),
        axis.line = element_line(color = 'black'))

# Fractal Dimension
ggplot(EcoFD) +
  geom_ridgeline(aes(x = xx, y = Realm, height = yy, fill = Realm, colour = Realm, scale = 0.2)) +
  scale_x_continuous(labels = function(i) { format(i, digits = 3) }) +
  scale_fill_manual(values = rev(c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF"))) +
  scale_color_manual(values = rev(c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF"))) +
  xlab(NULL) +
  ylab(NULL) +
  theme_ridges() +
  theme(axis.text.x = element_text(size = 20, colour = "black"), 
        axis.text.y = element_text(size = 20, colour = "black"),
        panel.border = element_rect(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 0.9),
        legend.position = 'none',
        plot.margin = margin(5,15,5,5),
        axis.line = element_line(color = 'black'))

# Generate plots (Biome level)
# Rugosity
ggplot(EcoRugosity) +
  geom_ridgeline(aes(x = xx, y = Biome, height = yy, fill = Realm, colour = Realm, scale = 1)) +
  scale_x_continuous(labels = function(i) { format(exp(i), scientific = F, digits = 1) }) +
  scale_fill_manual(values = c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF")) +
  scale_color_manual(values = c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF")) +
  xlab(NULL) +
  ylab(NULL) +
  theme_ridges() +
  theme(axis.text.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        panel.border = element_rect(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 0.9),
        legend.position = 'none',
        plot.margin = margin(5,15,5,5),
        axis.line = element_line(color = 'black'))

# Fractal dimension
ggplot(EcoFD) +
  geom_ridgeline(aes(x = xx, y = Biome, height = yy, fill = Realm, colour = Realm, scale = 0.15)) +
  scale_x_continuous(labels = function(i) { format(i, digits = 3) }) +
  scale_fill_manual(values = c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF")) +
  scale_color_manual(values = c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF")) +
  xlab(NULL) +
  ylab(NULL) +
  theme_ridges() +
  theme(axis.text.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        panel.border = element_rect(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 0.9),
        legend.position = 'none',
        plot.margin = margin(5,15,5,5),
        axis.line = element_line(color = 'black'))

# Generate plots (Ecosystem level)
# Rugosity
ggplot(EcoRugosity) +
  geom_ridgeline(aes(x = xx, y = EcosystemType, height = yy, fill = Realm, colour = Realm, scale = 1)) +
  scale_x_continuous(labels = function(i) { format(exp(i), scientific = F, digits = 1) }) +
  scale_fill_manual(values = c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF")) +
  scale_color_manual(values = c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF")) +
  xlab(NULL) +
  ylab(NULL) +
  theme_ridges() +
  theme(axis.text.x = element_text(size = 20, colour = "black"), 
        axis.text.y = element_text(size = 9, colour = "black"), 
        panel.border = element_rect(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 0.9),
        legend.position = 'none',
        plot.margin = margin(5,15,5,5),
        axis.line = element_line(color = 'black'))

# Fractal dimension
ggplot(EcoFD) +
  geom_ridgeline(aes(x = xx, y = EcosystemType, height = yy, fill = Realm, colour = Realm, scale = 0.15)) +
  scale_x_continuous(labels = function(i) { format(i, digits = 3) }) +
  scale_fill_manual(values = c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF")) +
  scale_color_manual(values = c("#00204DFF", "#31446BFF", 'lightblue', "#5F9258FF", "#CBBA69FF", "#FFEA46FF")) +
  xlab(NULL) +
  ylab(NULL) +
  theme_ridges() +
  theme(axis.text.x = element_text(size = 20, colour = "black"), 
        axis.text.y = element_text(size = 9, colour = "black"), 
        panel.border = element_rect(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 0.9),
        legend.position = 'none',
        plot.margin = margin(5,15,5,5),
        axis.line = element_line(color = 'black'))

### Complexity distributions across Land Use types -------------------------------
# Compute density distributions across Land use types
LandUseTypes <- LandUseDat %>%
  dplyr::select(LandUse) %>%
  distinct() %>%
  collect()
# set output storage
LURugosity <- data.frame(LandUse = character(),
                         xx = I(list()),
                         yy = I(list()))
LUFD <- data.frame(LandUse = character(),
                   xx = I(list()),
                   yy = I(list()))

# compute density distributions for complexity estimates associated with each Land Use type.
for(ii in seq_along(LandUseTypes$LandUse)){
  print(ii) # progress read out
  # Isolate Rugosity values for each category
  TmpR <- LandUseDat %>% 
    filter(LandUse == LandUseTypes$LandUse[ii]) %>%
    dplyr::select(R) %>%
    # Pull the raw values into memory (ensuring data fits in RAM)
    collect() %>% 
    reframe(# standard density returns an object; we wrap it in a data frame
      as.data.frame(density(R)[c("x", "y")])) %>%
    rename(value = x, density = y)
  # Repeat for fractal dimension 
  TmpD <- LandUseDat %>% 
    filter(LandUse == LandUseTypes$LandUse[ii]) %>%
    dplyr::select(D) %>%
    # Pull the raw values into memory (ensuring data fits in RAM)
    collect() %>% 
    reframe(# standard density returns an object; we wrap it in a data frame
      as.data.frame(density(D)[c("x", "y")])) %>%
    rename(value = x, density = y)
  # Save desired outputs
  # Rugosity
  LURugosity[ii, 'LandUse'] <- as.character(LandUseTypes$LandUse[ii])
  LURugosity[[ii, 'xx']] <- TmpR$value
  LURugosity[[ii, 'yy']] <- TmpR$density
  # Fractal dimension
  LUFD[ii, 'LandUse'] <- as.character(LandUseTypes$LandUse[ii])
  LUFD[[ii, 'xx']] <- TmpD$value
  LUFD[[ii, 'yy']] <- TmpD$density
} # close loop

# Reshape to a long format
LURugosity <- data.frame(LandUse = rep(LURugosity$LandUse, each = 512),
                          xx = unlist(LURugosity$xx),
                          yy = unlist(LURugosity$yy))
LUFD <- data.frame(LandUse = rep(LUFD$LandUse, each = 512),
                    xx = unlist(LUFD$xx),
                    yy = unlist(LUFD$yy))

# Ensure appropriate variable formats
LURugosity$LandUse <- factor(LURugosity$LandUse, levels = rev(c("Urban","Cropland","Mosaic Vegetation/Cropland","Mosaic Vegetation",
                                                            "Grassland","Shrubland","Lichens & Mosses","Sparse Vegetation",         
                                                            "Mixed Tree Cover","Tree Cover (Needleleaved)","Tree Cover (Broadleaved)",
                                                            "Wetlands","Bare substrate","Permenant Snow & Ice")))
LUFD$LandUse <- factor(LUFD$LandUse, levels = rev(c("Urban","Cropland","Mosaic Vegetation/Cropland","Mosaic Vegetation",
                                                "Grassland","Shrubland","Lichens & Mosses","Sparse Vegetation",         
                                                "Mixed Tree Cover","Tree Cover (Needleleaved)","Tree Cover (Broadleaved)",
                                                "Wetlands","Bare substrate","Permenant Snow & Ice")))

# Generate plots
# Rugosity
ggplot(LURugosity) +
  geom_ridgeline(aes(x = xx, y = LandUse, height = yy, fill = LandUse, colour = LandUse, scale = 4)) +
  scale_x_continuous(labels = function(i) { format(exp(i), scientific = F, digits = 1) }) +
  scale_fill_manual(values = magma(n = length(LandUseTypes$LandUse)),
                    limits = c("Urban","Cropland","Mosaic Vegetation/Cropland","Mosaic Vegetation",
                               "Grassland","Shrubland","Lichens & Mosses","Sparse Vegetation",         
                               "Mixed Tree Cover","Tree Cover (Needleleaved)","Tree Cover (Broadleaved)",
                               'Wetlands',"Bare substrate","Permenant Snow & Ice")) +
  scale_color_manual(values = magma(n = length(LandUseTypes$LandUse)),
                     limits = c("Urban","Cropland","Mosaic Vegetation/Cropland","Mosaic Vegetation",
                                "Grassland","Shrubland","Lichens & Mosses","Sparse Vegetation",         
                                "Mixed Tree Cover","Tree Cover (Needleleaved)","Tree Cover (Broadleaved)",
                                'Wetlands',"Bare substrate","Permenant Snow & Ice")) +
  xlab(NULL) +
  ylab(NULL) +
  theme_ridges() +
  theme(axis.text.x = element_text(size = 20, colour = "black"), 
        axis.text.y = element_text(size = 15, colour = "black"), 
        panel.border = element_rect(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 0.9),
        legend.position = 'none',
        plot.margin = margin(5,15,5,5),
        axis.line = element_line(color = 'black'))

# Fractal Dimension
ggplot(LUFD) +
  geom_ridgeline(aes(x = xx, y = LandUse, height = yy, fill = LandUse, colour = LandUse, scale = 0.4)) +
  scale_fill_manual(values = magma(n = length(LandUseTypes$LandUse)),
                    limits = c("Urban","Cropland","Mosaic Vegetation/Cropland","Mosaic Vegetation",
                               "Grassland","Shrubland","Lichens & Mosses","Sparse Vegetation",         
                               "Mixed Tree Cover","Tree Cover (Needleleaved)","Tree Cover (Broadleaved)",
                               'Wetlands',"Bare substrate","Permenant Snow & Ice")) +
  scale_color_manual(values = magma(n = length(LandUseTypes$LandUse)),
                     limits = c("Urban","Cropland","Mosaic Vegetation/Cropland","Mosaic Vegetation",
                                "Grassland","Shrubland","Lichens & Mosses","Sparse Vegetation",         
                                "Mixed Tree Cover","Tree Cover (Needleleaved)","Tree Cover (Broadleaved)",
                                'Wetlands',"Bare substrate","Permenant Snow & Ice")) +
  xlab(NULL) +
  ylab(NULL) +
  theme_ridges() +
  theme(axis.text.x = element_text(size = 20, colour = "black"), 
        axis.text.y = element_text(size = 15, colour = "black"), 
        panel.border = element_rect(linewidth = 1.5),
        axis.ticks = element_line(linewidth = 0.9),
        legend.position = 'none',
        plot.margin = margin(5,15,5,5),
        axis.line = element_line(color = 'black'))


######################################
# STEP 6: Extract examples of complex environments
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
