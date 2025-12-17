# This script evaluates how differing measures of complexity and geodiversity correspond to one another
# and explores the mechanisms through which structural complexity mediates biodiversity.

# Date last modified: Oct 2025
# Primary Author: James Cant
# -----------------------------------------------------------------------------------------

# Clear working directory
rm(list = ls())

# Load required packages
library(sf) 
library(terra) 
library(rcompanion) 
library(ggplot2) 
library(dplyr) 
library(tidyr) 
library(data.table)
library(mFD) 
library(ggdist) 
library(ggeffects)
library(pbapply) 
library(spsUtil) 
library(brms) 
library(tidybayes) 
library(ggh4x) 
library(forcats) 
library(data.table)
library(reshape2)
library(fishualize)

# Set random number seed
set.seed(56570)

# Modify R digit printing to prevent the rounding of very close values
options(digits = 22)

# Define file pathways
ComplexPath <- '/FILE _DIRECTORY 1/' # Containing complexity raster files
IslandPath <- '/FILE _DIRECTORY 2/' # containing downloaded island shapefiles
DataPath <- '/FILE _DIRECTORY 3/' # Bird diversity data

# Setup time saving index
FirstRun <- FALSE

# Define function for calculating Island isolation following the method presented in Weigelt & Kreft (2013) Ecography 36: 417-429.
calcIsolation <- function(island, b) {
  
  # Identify the landmasses (Mainland and large islands) situated within distance buffers  from selected island
  LM_A <- pbsapply(seq_along(b), function(x) {
    # determine buffer radius
    buffer <- st_buffer(island, dist = b[x]) 
    # Isolate landmass regions within buffer area
    Intersects <- vect(c(quiet(st_intersection(st_buffer(Mainland,0), buffer)$geometry), quiet(st_intersection(st_buffer(BigIslands,0), buffer)$geometry)))
    Intersects <- terra::intersect(Intersects, vect(buffer))
    # ensure the selected island is not duplicated within overlapping landmasses
    Intersects <- erase(Intersects, vect(island))
    # calculate landmass area covered by the overlap (Km2) within selected buffer
    A <- sum(st_area(st_as_sf(Intersects))/1e+06)
    return(A) })
  
  # Sum all land mass areas covered across the buffer
  return(sum(LM_A))
}


######################################
# STEP 1: Import & Clean Raw Data
######################################

## ----- Island Biodiversity and Complexity Data

if(FirstRun == TRUE) { # this component can take a period to run so is only necessary if this is the first time running the script
  
  ### Island Complexity  --------------------
  # Load the Complexity and Geodiversity Rasters
  DRast <- rast(paste0(ComplexPath, 'GlobalFractalDimension.tif'))
  RRast <- rast(paste0(ComplexPath, 'GlobalRugosity.tif'))
  HRast <- rast(paste0(ComplexPath, 'GlobalHeightRange.tif'))
  
  # Ensure resolutions and extents match across the geometric complexity rasters
  ext(RRast) == ext(DRast)
  ext(HRast) == ext(DRast)
  
  # Remove pixel from the complexity rasters corresponding with non-fractal locations (i.e. D = 2)
  RRast <- mask(RRast, DRast, maskvalues = 2) 
  HRast <- mask(HRast, DRast, maskvalues = 2)
  DRast <- mask(DRast, DRast, maskvalues = 2) 
  
  # Transform complexity values as needed
  RRast <- app(RRast, log10)
  HRast <- app(HRast, log10)
  
  # Load raster of Human Population size
  PopRast <- rast(paste0(ComplexPath, 'PopCount_Mollweide_1870m.tif'))
  
  # Load global island shape files (confirming their validity simultaneously)
  BigIslands <- st_read(paste0(IslandPath, 'GlobalIslands_BigIslands.shp'))
  SmallIslands <- st_read(paste0(IslandPath, 'GlobalIslands_SmallIslands.shp'))
  VSmallIslands <- st_read(paste0(IslandPath, 'GlobalIslands_VerySmallIslands.shp'))
  GlobalIslands <- do.call(rbind, list(BigIslands,SmallIslands,VSmallIslands))
  # remove NA's from name variables
  GlobalIslands <- GlobalIslands %>% mutate_at(c('Name_USGSO','NAME_wcmcI'), ~replace_na(.,''))
  
  # Load mainland shapefiles
  Mainland <- st_read(paste0(IslandPath, 'GlobalIslands_MainLand.shp'))
  
  # Load Raw Islands data
  IslandData <- read.csv(paste0(DataPath, 'Raw Island Data_Combined.csv'), encoding = "UTF-8", stringsAsFactors = F)
  # Remove entries with no corresponding shape file names and those that represent coral atolls
  IslandData <- IslandData[!(IslandData$NameUSGS == "" & IslandData$NameWCMC == ""),]
  IslandData <- IslandData[is.na(IslandData$Atoll),]
  
  # Add in blank variables for populating
  IslandData$Hsd <- IslandData$H <- IslandData$Rsd <-
    IslandData$R <- IslandData$Dsd <- IslandData$D <- IslandData$Isolation <- IslandData$HPI <- IslandData$AreaKM2 <- 
    IslandData$MeanLat <-  NA
  
  # Work through each island to extract complexity and size details using corresponding shapefile.
  for (ii in 1:dim(IslandData)[1]){
    # progress print out
    cat(ii, IslandData$AuthorIslandName[ii], '\n')
    
    # Isolate corresponding island shapefile
    islandSelect <- GlobalIslands %>% filter(Name_USGSO == IslandData$NameUSGS[ii] & 
                                               NAME_wcmcI == IslandData$NameWCMC[ii])
    
    # Only continue if data found
    if(dim(islandSelect)[1] == 1) { # This place holder is required because there are some island name clashes that need to be resolved manually.
      # Extract Island Size
      IslandData$AreaKM2[ii] <- islandSelect$IslandArea
      
      # Extract central island latitude
      IslandCoords <- suppressWarnings(st_centroid(islandSelect)$geometry) %>% st_transform(crs = st_crs('+proj=longlat +datum=WGS84 +no_defs +type=crs'))
      IslandData$MeanLat[ii] <- IslandCoords[[1]][2]
      
      # Estimate Island Isolation. Isolation is calculated as the proportional land mass area within incremental distance buffers of 100, 1000, 3000km
      IslandData$Isolation[ii] <- as.numeric(calcIsolation(island = islandSelect, b = c(1e+05, 1e+06, 2e+06)))
      
      # Isolate corresponding complexity and Geodiversity values
      IslandMask <- vect(islandSelect)
      # Human Population Index
      IslandP <- mask(crop(PopRast, IslandMask), IslandMask)
      IslandData$HPI[ii] <- sum(values(IslandP), na.rm = TRUE)
      # Fractal Dimension
      IslandD <- mask(crop(DRast, IslandMask), IslandMask)
      IslandData$D[ii] <- mean(values(IslandD), na.rm = TRUE)
      IslandData$Dsd[ii] <- sd(values(IslandD), na.rm = TRUE)
      # Rugosity
      IslandR <- mask(crop(RRast, IslandMask), IslandMask)
      IslandData$R[ii] <- mean(values(IslandR), na.rm = TRUE)
      IslandData$Rsd[ii] <- sd(values(IslandR), na.rm = TRUE)
      # Height Range
      IslandH <- mask(crop(HRast, IslandMask), IslandMask)
      IslandData$H[ii] <- mean(values(IslandH), na.rm = TRUE)
      IslandData$Hsd[ii] <- sd(values(IslandH), na.rm = TRUE)
    }
    
  } # end of loop
  
  # Manually add the details for the Islands skipped due to naming clashes.
  Index <-  as.numeric(rownames(IslandData[is.na(IslandData$MeanLat) & is.na(IslandData$AreaKM2),])) # Row index of clashing islands or islands with name codes that did not open correctly.
  ID <- c(282526,285393,274063,324152,278543,278542,277887,277165,273844,273975,274257,273982,274121,278530,278561,120251,274005,274103,273781,273770,NA,277045,273839,280058,280057,340374,206954,206929,199388,274684,279685,274669,279686,279689,279688,279692,280673,280682,280707,280708,280799,
          280836,273907,274386,274466,280394,280366,280360,280329,280265,274271,280322,280318,280300,280293,280269,280266,280249,280273,280271,280294,280298,280279,280286,280281,280283,280282,280239,280308,280299,280291,280296,280290,280292,280311,280261,280253,280264,280244,280238,209770,209766,
          209763,209726,209725,280226,209761,280225,280227,209742,280235,280230,280232,280231,209735,209704,280228,122955) # manually identified corresponding unique island IDs
  
  # Repeat extraction loop
  for (ii in 1:length(Index)){  
    # progress print out
    cat(ii, IslandData[rownames(IslandData) == Index[ii],]$AuthorIslandName, '\n')
    
    # Isolate corresponding island shapefile
    islandSelect <- GlobalIslands %>% filter(ALL_Uniq == ID[ii])
    
    # In case it is still not possible to identify the corresponding island
    if(dim(islandSelect)[1] == 1) {
      # Extract Island Size
      IslandData[rownames(IslandData) == Index[ii],]$AreaKM2 <- islandSelect$IslandArea
      
      # Extract central island latitude
      IslandCoords <- suppressWarnings(st_centroid(islandSelect)$geometry) %>% st_transform(crs = st_crs('+proj=longlat +datum=WGS84 +no_defs +type=crs'))
      IslandData[rownames(IslandData) == Index[ii],]$MeanLat <- IslandCoords[[1]][2]
      
      # Estimate Island Isolation. 
      IslandData[rownames(IslandData) == Index[ii],]$Isolation <- as.numeric(calcIsolation(island = islandSelect, b = c(1e+05, 1e+06, 2e+06)))
      
      # Isolate corresponding complexity and Geodiversity values
      IslandMask <- vect(islandSelect)
      # Human Population Index
      IslandP <- mask(crop(PopRast, IslandMask), IslandMask)
      IslandData[rownames(IslandData) == Index[ii],]$HPI <- sum(values(IslandP), na.rm = TRUE)
      # Fractal Dimension
      IslandD <- mask(crop(DRast, IslandMask), IslandMask)
      IslandData[rownames(IslandData) == Index[ii],]$D <- mean(values(IslandD), na.rm = TRUE)
      IslandData[rownames(IslandData) == Index[ii],]$Dsd <- sd(values(IslandD), na.rm = TRUE)
      # Rugosity
      IslandR <- mask(crop(RRast, IslandMask), IslandMask)
      IslandData[rownames(IslandData) == Index[ii],]$R <- mean(values(IslandR), na.rm = TRUE)
      IslandData[rownames(IslandData) == Index[ii],]$Rsd <- sd(values(IslandR), na.rm = TRUE)
      # Height Range
      IslandH <- mask(crop(HRast, IslandMask), IslandMask)
      IslandData[rownames(IslandData) == Index[ii],]$H <- mean(values(IslandH), na.rm = TRUE)
      IslandData[rownames(IslandData) == Index[ii],]$Hsd <- sd(values(IslandH), na.rm = TRUE)
    }
  }
  
  # Clean up any NaNs
  IslandData <- IslandData %>% mutate_all(~ifelse(is.nan(.), NA, .))
  
  # Remove entries for which no corresponding island data can be found
  IslandData <- IslandData[!(is.na(IslandData$Isolation)),]
  
  ### Island Functional diversity --------------------
  
  # Load trait matrix
  traitDF <- read.csv(paste0(DataPath, 'SpeciesTraits.csv'))
  speciesNames <- traitDF$Species.2 
  
  # Isolate traits of interest
  traitDF <- traitDF[,c('Beak.Length.culmen','Beak.Length.nares','Beak.Width','Beak.Depth',
                        'Tarsus.Length','Wing.Length','Hand.Wing.Index','Tail.Length','Mass',
                        'Trophic.Niche','Primary.Lifestyle')]
  traitDF$Trophic.Niche <- as.factor(traitDF$Trophic.Niche)
  traitDF$Primary.Lifestyle <- as.factor(traitDF$Primary.Lifestyle)
  rownames(traitDF) <- speciesNames
  # remove, incomplete entries as functional indices are sensitive to missing trait values
  traitDF <- traitDF[complete.cases(traitDF),]
  speciesNames <- rownames(traitDF) # update species name vector
  
  # Define matrix outlining trait types
  typesDF <- data.frame(trait_name = c('Beak.Length.culmen','Beak.Length.nares','Beak.Width','Beak.Depth',
                                       'Tarsus.Length','Wing.Length','Hand.Wing.Index','Tail.Length','Mass',
                                       'Trophic.Niche','Primary.Lifestyle'),
                        trait_type = c('Q','Q','Q','Q','Q','Q','Q','Q','Q','N','N'))
  
  # Calculate functional trait-based distances between each species.
  funDist <- funct.dist(sp_tr = traitDF, tr_cat = typesDF, # requires time to compute.
                        metric        = "gower",
                        scale_euclid  = "scale_center",
                        ordinal_var   = "classic",
                        weight_type   = "equal",
                        stop_if_NA    = TRUE)
  
  # Use the functional distances between species to estimate a measure of the functional space characterised by the species list. 
  # This function also determines which multidimensional space offers the best characterisation of the functional diversity present.
  fspaces <- quality.fspaces( # again requires time to compute.
    sp_dist             = funDist,
    maxdim_pcoa         = 10,
    deviation_weighting = "absolute",
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # compare multidimensional space 
  fspaces$"quality_fspaces"
  # with the lowest value the D multidimensional approach offers the best representation of the functional diversity plane
  # extract axis coordinates for the computed functional space
  axisCoords <- fspaces$"details_fspaces"$"sp_pc_coord"
  
  # Add new data variable
  IslandData$FRich <- NA
  
  # Number of different studies (for which functional diversity can be computed)
  studyList <- unique(IslandData[IslandData$Source == 'Matthews et al', c('Archipelago')])
  
  # Working through each island read in the corresponding assemblage data and estimate its functional richness
  for(ii in 1:length(studyList)){
    # progress
    cat(ii, studyList[ii], '\n')
    
    # Isolate file containing corresponding assemblage data
    fileSelect <- grep(studyList[ii], list.files(DataPath, full.names = TRUE), value = TRUE)
    
    # Open file and identify corresponding assemblage data (presented as occurrence)
    dat <- read.csv(fileSelect, stringsAsFactors = F, check.names = FALSE)
    dat <- dat[!dat$species %in% c('Area (ha)', 'sp.r', ''),]
    # extract island names ahead of transposing the data
    IslandNames <- names(dat[-1]) 
    dat <- as.matrix(transpose(setDT(dat), make.names = 'species'))
    rownames(dat) <- IslandNames
    
    # Remove Island assemblages with less species than the number of desired pcoa axes.  
    dat <- dat[rowSums(dat) > 6,]
    # Remove Islands for which there is no corresponding complexity data
    dat <- dat[rownames(dat) %in% IslandData[IslandData$Archipelago == studyList[ii], 'AuthorIslandName'],]
    # reorder dat rows to match island ordering in main dataframe
    x <- c(IslandData[IslandData$Archipelago == studyList[ii] &
                        IslandData$AuthorIslandName %in% rownames(dat), 'AuthorIslandName']) # create a vector with letters in the desired order
    # apply reordering
    dat <- as.data.frame(dat)
    dat <- as.matrix(dat %>% slice(match(x, rownames(dat))))
    
    # Estimate functional richness
    funValues <- alpha.fd.multidim(
      sp_faxes_coord   = axisCoords[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")],
      asb_sp_w         = dat,
      ind_vect         = "fric",
      scaling          = TRUE,
      check_input      = TRUE,
      details_returned = TRUE)
    
    # Insert functional richness values  
    IslandData[IslandData$Archipelago == studyList[ii] &
                 IslandData$AuthorIslandName %in% rownames(dat), 'FRich'] <- funValues$functional_diversity_indices$fric
  } 
  
  # Save data file
  write.csv(IslandData, paste0(ComplexPath, 'IslandComplexityData.csv'), row.names = F)
  # clear memory space
  rm(funDist, fspaces, traitDF, typesDF)
  gc()
  
} else { # load already processed data
  IslandData <- read.csv(paste0(ComplexPath, 'IslandComplexityData.csv'))
}

# Estimate absolute latitude
IslandData$absLat <- abs(IslandData$MeanLat)
hist(IslandData$absLat)

# Check if data transformations required
# Species Richness
hist(IslandData$Extant.natives)
hist(log10(IslandData$Extant.natives)) # log transformation appropriate
IslandData$SpRich <- log10(IslandData$Extant.natives)

# Functional Richness
hist(IslandData$FRich)
hist(log10(IslandData$FRich)) 
transformTukey(IslandData$FRich, plotit = FALSE)
hist(IslandData$FRich^0.275)
IslandData$FRich2 <- IslandData$FRich^0.275

# Island human population index
hist(IslandData$HPI)
hist(log10(IslandData$HPI+1)) # log transformation appropriate 
IslandData$HPI <- log10(IslandData$HPI+1)

# Island Isolation
hist(IslandData$Isolation)
hist(log10(IslandData$Isolation))
transformTukey(IslandData$Isolation, plotit = FALSE)
hist(IslandData$Isolation^0.5)
IslandData$Isolation <- IslandData$Isolation^0.5

# Island Size
hist(IslandData$AreaKM2)
hist(log10(IslandData$AreaKM2))
IslandData$AreaKM2 <- log10(IslandData$AreaKM2)

# Fractal Dimension
hist(IslandData$D)

# Rugosity
hist(IslandData$R)

# Height Range
hist(IslandData$H) # No transformations needed for complexity variables.

# Convert Archipelago into a factor variable
IslandData$Archipelago <- as.factor(IslandData$Archipelago)

######################################
# STEP 2: Run Analyses
######################################  

## -----------------------------------------------------------
##### Biodiversity and complexity
# Implement species-area/complexity models.
#### Species Richness
# baseline model
SpMod <- brm(formula = SpRich ~ poly(AreaKM2, 2) + (AreaKM2|Archipelago) + absLat + Isolation + HPI, # accounting for island isolation, latitude, and associated human population
             data = IslandData,                                                                     # and the random effect of Archipelago.
             family = 'gaussian',
             iter = 3000, chains = 4, seed = 1200)
SpModR <- brm(formula = SpRich ~ poly(R, 2) + (R|Archipelago) + AreaKM2 + absLat + Isolation + HPI, # accounting for the confounding effects of island area.
              data = IslandData,
              family = 'gaussian',
              iter = 3000, chains = 4, seed = 15999)
SpModD <- brm(formula = SpRich ~ poly(D, 2) + (D|Archipelago) + AreaKM2 + absLat + Isolation + HPI,
              data = IslandData,
              family = 'gaussian',
              iter = 3000, chains = 4, seed = 8765)
SpModH <- brm(formula = SpRich ~ poly(H, 2) + (H|Archipelago) + AreaKM2 + absLat + Isolation + HPI, 
              data = IslandData,
              family = 'gaussian',
              iter = 3000, chains = 4, seed = 1256)
SpModG <- brm(formula = SpRich ~ poly(G, 2) + (G|Archipelago) + AreaKM2 + absLat + Isolation + HPI, 
              data = IslandData[!(is.na(IslandData$G)),],
              family = 'gaussian',
              iter = 3000, chains = 4, seed = 1256)

# Model fits
bayes_R2(SpMod)
bayes_R2(SpModH)
bayes_R2(SpModR)
bayes_R2(SpModD)
bayes_R2(SpModG)

#### Functional Richness
# baseline model
FMod <- brm(formula = FRich2 ~ poly(AreaKM2, 2) + (AreaKM2|Archipelago) + absLat + Isolation + HPI, # accounting for island isolation, latitude, and associated human population
            data = IslandData,
            family = 'gaussian',
            iter = 3000, chains = 4, seed = 4333)
FModR <- brm(formula = FRich2 ~ poly(R, 2) + (R|Archipelago) + AreaKM2 + absLat + Isolation + HPI, 
             data = IslandData,
             family = 'gaussian',
             iter = 3000, chains = 4, seed = 2167)
FModD <- brm(formula = FRich2 ~ poly(D, 2) + (D|Archipelago) + AreaKM2 + absLat + Isolation + HPI,
             data = IslandData,
             family = 'gaussian',
             iter = 3000, chains = 4, seed = 3412)
FModH <- brm(formula = FRich2 ~ poly(H, 2) + (H|Archipelago) + AreaKM2 + absLat + Isolation + HPI, 
             data = IslandData,
             family = 'gaussian',
             iter = 3000, chains = 4, seed = 4987)
FModG <- brm(formula = FRich2 ~ poly(G, 2) + (G|Archipelago) + AreaKM2 + absLat + Isolation + HPI, 
             data = IslandData[!(is.na(IslandData$G)),],
             family = 'gaussian',
             iter = 3000, chains = 4, seed = 4987)
# Model fits
bayes_R2(FMod)
bayes_R2(FModH)
bayes_R2(FModR)
bayes_R2(FModD)
bayes_R2(FModG)

# Extract posterior coefficients of the first polynomial slope coefficient
# Species Richness
SpDat <- data.frame(A = spread_draws(SpMod, b_polyAreaKM221)[4],
                    R = spread_draws(SpModR, b_polyR21)[4],
                    H = spread_draws(SpModH, b_polyH21)[4],
                    D = spread_draws(SpModD, b_polyD21)[4])
names(SpDat) <- c('A', 'R', 'H', 'D')
SpDat <- melt(SpDat, variable.name = 'Metric', value.name = 'Slope')
SpDat$Metric <- factor(SpDat$Metric, levels = c('D', 'H', 'R', 'A'))
# Functional Richness
FRDat <- data.frame(A = spread_draws(FMod, b_polyAreaKM221)[4],
                    R = spread_draws(FModR, b_polyR21)[4],
                    H = spread_draws(FModH, b_polyH21)[4],
                    D = spread_draws(FModD, b_polyD21)[4])
names(FRDat) <- c('A', 'R', 'H', 'D')
FRDat <- melt(FRDat, variable.name = 'Metric', value.name = 'Slope')
FRDat$Metric <- factor(FRDat$Metric, levels = c('D', 'H', 'R', 'A'))

# Plot posterior coefficients for the first polynomial
# Species Richness
ggplot(SpDat, aes(y = Metric, x = Slope, fill = Metric)) +
  stat_halfeye(alpha = 0.8) +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', linewidth = 1) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_manual(values = c('#73D8FEFF','#1D373AFF','#3E938BFF','#FDC718FF')) + #fish(n = 5, option = 'Acanthurus_olivaceus', direction = -1)
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))
# Functional Richness
ggplot(FRDat, aes(y = Metric, x = Slope, fill = Metric)) +
  stat_halfeye(alpha = 0.8) +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', linewidth = 1) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_manual(values = c('#73D8FEFF','#1D373AFF','#3E938BFF','#FDC718FF')) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))

# Visualise affects of complexity on species and functional richness
SpA <- conditional_effects(SpMod)[[1]][,c('effect1__', 'estimate__', 'lower__', 'upper__')]; names(SpA) <- c('Area', 'mean', 'lower', 'upper')
SpR <- conditional_effects(SpModR)[[1]][,c('effect1__', 'estimate__', 'lower__', 'upper__')]; names(SpR) <- c('R', 'mean', 'lower', 'upper')
SpD <- conditional_effects(SpModD)[[1]][,c('effect1__', 'estimate__', 'lower__', 'upper__')]; names(SpD) <- c('D', 'mean', 'lower', 'upper')
SpH <- conditional_effects(SpModH)[[1]][,c('effect1__', 'estimate__', 'lower__', 'upper__')]; names(SpH) <- c('H', 'mean', 'lower', 'upper')
SpG <- conditional_effects(SpModG)[[1]][,c('effect1__', 'estimate__', 'lower__', 'upper__')]; names(SpG) <- c('G', 'mean', 'lower', 'upper')
FA <- conditional_effects(FMod)[[1]][,c('effect1__', 'estimate__', 'lower__', 'upper__')]; names(FA) <- c('Area', 'mean', 'lower', 'upper')
FR <- conditional_effects(FModR)[[1]][,c('effect1__', 'estimate__', 'lower__', 'upper__')]; names(FR) <- c('R', 'mean', 'lower', 'upper')
FD <- conditional_effects(FModD)[[1]][,c('effect1__', 'estimate__', 'lower__', 'upper__')]; names(FD) <- c('D', 'mean', 'lower', 'upper')
FH <- conditional_effects(FModH)[[1]][,c('effect1__', 'estimate__', 'lower__', 'upper__')]; names(FH) <- c('H', 'mean', 'lower', 'upper')
FG <- conditional_effects(FModG)[[1]][,c('effect1__', 'estimate__', 'lower__', 'upper__')]; names(FG) <- c('G', 'mean', 'lower', 'upper')

#### Species Richness versus 
# Area
ggplot() +
  geom_ribbon(aes(x = Area, ymin = lower, ymax = upper), data = SpA, alpha = 0.4, fill = '#FDC718FF') +
  geom_point(aes(x = AreaKM2, y = SpRich), data = IslandData, size = 5, col = '#FDC718FF') + 
  geom_line(aes(x = Area, y = mean), data = SpA, linetype = 'solid', linewidth = 4, col = '#191819FF') +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6), labels = function(i) format(exp(i), digits = 2)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2), labels = function(i) format(exp(i), digits = 2)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))
# Rugosity
ggplot() +
  geom_ribbon(aes(x = R, ymin = lower, ymax = upper), data = SpR, alpha = 0.4, fill = '#3E938BFF') +
  geom_point(aes(x = R, y = SpRich), data = IslandData, size = 5, col = '#3E938BFF') + 
  geom_line(aes(x = R, y = mean), data = SpR, linetype = 'solid', linewidth = 4, col = '#191819FF') +
  scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0), labels = function(i) format(exp(i), digits = 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2), labels = function(i) format(exp(i), digits = 2)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))
# Height Range
ggplot() +
  geom_ribbon(aes(x = H, ymin = lower, ymax = upper), data = SpH, alpha = 0.4, fill = '#1D373AFF') +
  geom_point(aes(x = H, y = SpRich), data = IslandData, size = 5, col = '#1D373AFF') + 
  geom_line(aes(x = H, y = mean), data = SpH, linetype = 'solid', linewidth = 4, col = '#191819FF') +
  scale_x_continuous(breaks = c(-1.8, -1.2, -0.6, -0, 0.6), labels = function(i) format(exp(i), digits = 2)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2), labels = function(i) format(exp(i), digits = 2)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))
# Geodiversity
ggplot() +
  geom_ribbon(aes(x = G, ymin = lower, ymax = upper), data = SpG, alpha = 0.4, fill = '#1942CDFF') +
  geom_point(aes(x = G, y = SpRich), data = IslandData, size = 5, col = '#1942CDFF') + 
  geom_line(aes(x = G, y = mean), data = SpG, linetype = 'solid', linewidth = 4, col = '#191819FF') +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2), labels = function(i) format(exp(i), digits = 2)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))
# Fractal Dimension
ggplot() +
  geom_ribbon(aes(x = D, ymin = lower, ymax = upper), data = SpD, alpha = 0.4, fill = '#73D8FEFF') +
  geom_point(aes(x = D, y = SpRich), data = IslandData, size = 5, col = '#73D8FEFF') + 
  geom_line(aes(x = D, y = mean), data = SpD, linetype = 'solid', linewidth = 4, col = '#191819FF') +
  scale_x_continuous(breaks = c(2, 2.25, 2.5), labels = function(i) format(i, digits = 3)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2), labels = function(i) format(exp(i), digits = 2)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))

# Functional Richness versus
# Area
ggplot() +
  geom_ribbon(aes(x = Area, ymin = lower, ymax = upper), data = FA, alpha = 0.7, fill = '#FDC718FF') +
  geom_point(aes(x = AreaKM2, y = FRich2), data = IslandData, size = 5, col = '#FDC718FF') + 
  geom_line(aes(x = Area, y = mean), data = FA, linetype = 'solid', linewidth = 4, col = '#490D07FF') +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6), labels = function(i) format(exp(i), digits = 2)) +
  scale_y_continuous(breaks = c(0,0.1, 0.2, 0.3, 0.4), labels = function(i) format(i^(1/0.275), digits = 2)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))
# Rugosity
ggplot() +
  geom_ribbon(aes(x = R, ymin = lower, ymax = upper), data = FR, alpha = 0.7, fill = '#3E938BFF') +
  geom_point(aes(x = R, y = FRich2), data = IslandData, size = 5, col = '#3E938BFF') + 
  geom_line(aes(x = R, y = mean), data = FR, linetype = 'solid', linewidth = 4, col = '#490D07FF') +
  scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0), labels = function(i) format(exp(i), digits = 1)) +
  scale_y_continuous(breaks = c(0,0.1, 0.2, 0.3, 0.4), labels = function(i) format(i^(1/0.275), digits = 2)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))
# Height Range
ggplot() +
  geom_ribbon(aes(x = H, ymin = lower, ymax = upper), data = FH, alpha = 0.7, fill = '#1D373AFF') +
  geom_point(aes(x = H, y = FRich2), data = IslandData, size = 5, col = '#1D373AFF') + 
  geom_line(aes(x = H, y = mean), data = FH, linetype = 'solid', linewidth = 4, col = '#490D07FF') +
  scale_x_continuous(breaks = c(-1.8, -1.2, -0.6, -0, 0.6), labels = function(i) format(exp(i), digits = 2)) +
  scale_y_continuous(breaks = c(0,0.1, 0.2, 0.3, 0.4), labels = function(i) format(i^(1/0.275), digits = 2)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))
# Geodiversity
ggplot() +
  geom_ribbon(aes(x = G, ymin = lower, ymax = upper), data = FG, alpha = 0.4, fill = '#1942CDFF') +
  geom_point(aes(x = G, y = FRich2), data = IslandData, size = 5, col = '#1942CDFF') + 
  geom_line(aes(x = G, y = mean), data = FG, linetype = 'solid', linewidth = 4, col = '#191819FF') +
  scale_y_continuous(breaks = c(0,0.1, 0.2, 0.3, 0.4), labels = function(i) format(i^(1/0.275), digits = 2)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))
# Fractal Dimension
ggplot() +
  geom_ribbon(aes(x = D, ymin = lower, ymax = upper), data = FD, alpha = 0.7, fill = '#73D8FEFF') +
  geom_point(aes(x = D, y = FRich2), data = IslandData, size = 5, col = '#73D8FEFF') + 
  geom_line(aes(x = D, y = mean), data = FD, linetype = 'solid', linewidth = 4, col = '#490D07FF') +
  scale_x_continuous(breaks = c(2, 2.25, 2.5), labels = function(i) format(i, digits = 3)) +
  scale_y_continuous(breaks = c(0,0.1, 0.2, 0.3, 0.4), labels = function(i) format(i^(1/0.275), digits = 2)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_rect(linewidth = 1.5),
                     plot.margin = margin(10, 15, 10, 10),
                     legend.position = 'none',
                     axis.line = element_line(color = 'black'))

##### ------------------------------------------------- End of Code ----------------------------------