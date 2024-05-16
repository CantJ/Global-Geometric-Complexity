# This script is a compilation of analyses exploring the relationship between variables of structural 
# complexity across terrestrial and marine habitats.
#
# Author: James Cant
# Last Modified: May 2024
# ------------------------------------------------------------------------------------

# Clear working directory
rm(list = ls())

# load required files
library(terra)
library(arrow)
library(ggplot2)
library(hexbin)
library(ggdist)
library(viridis)
library(BayesFactor)
library(brms)
library(dplyr)
library(cmdstanr)
library(sf)
library(spsUtil)
library(gmodels)
library(ggeffects)
library(pbapply)

# Define monte carlo brms regression function
MCbrm <- function(RHbf, DHbf, DRbf, dat, n) {
  
  # Extract random data sample
  TmpData <- dat %>% select(Lon, Lat, D, R2, H2) %>% collect() %>% sample_n(size = n, replace = FALSE) 
  
  # Run Bayesian regression model for each pairwise combination
  mod1 <- quiet(brm(RHbf, family = 'gaussian', data = TmpData, chains = 1, iter = 2000, warmup = 1000, backend = 'cmdstanr'))
  mod2 <- quiet(brm(DHbf, family = 'gaussian', data = TmpData, chains = 1, iter = 2000, warmup = 1000, backend = 'cmdstanr'))
  mod3 <- quiet(brm(DRbf, family = 'gaussian', data = TmpData, chains = 1, iter = 2000, warmup = 1000, backend = 'cmdstanr'))
  
  # Extract model fits
  mod1R <- bayes_R2(mod1)[1]
  mod2R <- bayes_R2(mod2)[1]
  mod3R <- bayes_R2(mod3)[1]
  
  # Predict complexity values using assessed relationship
  pred1 <- ggpredict(mod1, terms = list(H2 = seq(-7,2, 0.01)))$predicted # covers the full range of H2 
  pred2 <- ggpredict(mod2, terms = list(D = seq(2,3,0.001)))$predicted # covers the full range of D
  pred3 <- ggpredict(mod3, terms = list(D = seq(2,3,0.001)))$predicted
  
  # Collate and return outputs
  return(list(RH = pred1, DH = pred2, DR = pred3, RH_R = mod1R, DH_R = mod2R, DR_R = mod3R))
}

######################################
# STEP 1: Data Cleaning
######################################

# Set random number seed
set.seed(45560)

# Modify R digit printing to prevent the rounding of very close values
options(digits = 22)



# Define number of Monte Carlo simulations
mcSim <- 1000

# Define sample size for plotting (This is to speed up plot processing times)
Nsize <- 2e+07 # 10% sample
# Define sample size for mcmc Bayesian regression analysis
N <- 1000

# Define file pathways
FilePath <- '/FILE_DIRECTORY 1/'
DEMPath <- '/FILE_DIRECTORY 2/'

# Define spatial resolution and CRS
L0 <- res(rast(paste0(DEMPath, 'GlobalDEM_Mollweide_187m.tif')))[1]
Moll <- crs(rast(paste0(DEMPath, 'GlobalDEM_Mollweide_187m.tif')))

# Is this the first time this script has been run? TRUE or FALSE?
FirstRun <- FALSE

# If this script is being run for the first time there are a few steps of data cleaning required. However, these can take a while to complete so if 
# these have been run in the past this step can be skipped with the already processed data loaded directly instead.
if(FirstRun == TRUE) { 
  
  # Load in the Complexity Data (needs initial processing using functions that aren't compatible with arrow)
  RawData <- read_parquet(paste0(FilePath, 'GlobalComplexity.parquet'), as_data_frame = T)
  # Omit entries for which R = 0, H = 0 and D = 2. These locations are not truly fractal so will be excluded from further analysis.
  RawData %>% filter(D==2) %>% dim()
  RawData <- RawData %>% mutate(R = replace(R, D==2, NA)) # drop corresponding rugosity entries
  RawData <- RawData %>% mutate(H = replace(H, D==2, NA)) # drop corresponding height range entries
  RawData <- RawData %>% mutate(D = replace(D, D==2, NA)) # finally drop fractal dimension entries
  # Remove all pixels with missing complexity data
  # This will also remove pixels outside mollweide projection boundaries for which complexity estimates were never estimated.
  RawData <- RawData %>% filter(complete.cases(H,R,D))
  
  # Confirm data ranges of complexity variables
  # Also a good way to check there are no NA's
  RawData %>% 
    summarise(Hmin = min(H), # shouldn't be below or equal to 0
              Rmin = min(R), # shouldn't be below or equal to 1
              Dmin = min(D), # shouldn't be below or equal to 2
              Hmax = max(H),
              Rmax = max(R),
              Dmax = max(D)) # shouldn't be higher than 3
  
  # Height Range and Rugosity need to be converted into their standardized form to ensure they are on a consistent plane to fractal dimension.
  RawData <- RawData %>%
    mutate(H2 = H/(sqrt(2)*L0),
           R2 = (R^2)-1) 
  
  # Check distributions
  # Height range
  RawData %>% with(hist(H2)) # Height Range requires transformation
  RawData %>% with(hist(log10(H2)))
  # Rugosity
  RawData %>% with(hist(R2)) # Rugosity requires transformation
  RawData %>% with(hist(log10(R2)))
  # Fractal Dimension
  RawData %>% with(hist(D)) # slight skew
  
  # Apply necessary transformations
  RawData <- RawData %>%
    mutate(H2 = log10(H2),
           R2 = log10(R2))
  
  # Return data to memory saving parquet format
  write_parquet(RawData, paste0(FilePath, 'ProcessedData.parquet')) # overwrite output to memory efficient object
  # Reload updated parquet file
  RawData <- read_parquet(paste0(FilePath,'ProcessedData.parquet'), as_data_frame = F)
  # clean memory space
  gc()
  
} else {
  
  # Load the already processed data frame in a memory efficient format.
  RawData <- read_parquet(paste0(FilePath,'ProcessedData.parquet'), as_data_frame = F)
  
}


######################################
# STEP 2: Test for Spatial Autocorrelation
######################################

# Estimating spatial autocorrelation takes a short while so can be skipped if already done.
if(FirstRun == TRUE) {
  
  # Load complexity rasters
  DRast <- rast(paste0(FilePath,'GlobalFractalDimension.tif'))
  RRast <- rast(paste0(FilePath,'GlobalRugosity.tif'))
  HRast <- rast(paste0(FilePath,'GlobalHeightRange.tif'))
  
  # Define nearest neighbor matrix (queens format - diagonals included)
  w <- matrix(c(1,1,1,1,0,1,1,1,1), nrow = 3)
  # test for spatial autocorrelation using Moran's I (-1 < I > 1)
  # Estimate Moran's I (can take a while to complete)
  D_I <- autocor(DRast, w, method = 'moran', global = TRUE)
  R_I <- autocor(RRast, w, method = 'moran', global = TRUE)  
  H_I <- autocor(HRast, w, method = 'moran', global = TRUE) 
  
}

######################################
# STEP 3: Quantify geometric relationships
######################################

# 1. Compare complexities across Terrestrial and Marine ecoregions.
# Extract desired variables
plotdata <- RawData %>% dplyr::select(Realm, D, R2, H2) %>% collect() %>% sample_n(size = Nsize, replace = FALSE)

# generate plots

# Fractal Dimension
ggplot(plotdata, aes(x = factor(Realm), y = D, fill = factor(Realm))) +
  stat_eye(adjust = 0.5, show.legend = FALSE, .width = 0, point_colour = NA) +
  scale_x_discrete(labels = c('Marine', 'Terrestrial')) +
  scale_fill_manual(values = c('#00008B','#BCAF6FFF'),
                    guide = NULL) +
  xlab('\nRealm') +
  ylab('Fractal Dimension\n') +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), 
                     axis.text.y = element_text(size = 15, colour = "black"),
                     panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))

# Rugosity
ggplot(plotdata, aes(x = factor(Realm), y = R2, fill = factor(Realm))) +
  stat_eye(adjust = 0.5, show.legend = FALSE, .width = 0, point_colour = NA) +
  scale_x_discrete(labels = c('Marine', 'Terrestrial')) +
  scale_y_continuous(labels = function(i) 10^i) +
  scale_fill_manual(values = c('#00008B','#BCAF6FFF'),
                    guide = NULL) +
  xlab('\nRealm') +
  ylab('Rugosity\n') +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), 
                     axis.text.y = element_text(size = 15, colour = "black"),
                     panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))

# Height Range
ggplot(plotdata, aes(x = factor(Realm), y = H2, fill = factor(Realm))) +
  stat_eye(adjust = 0.5, show.legend = FALSE, .width = 0, point_colour = NA) +
  scale_x_discrete(labels = c('Marine', 'Terrestrial')) +
  scale_y_continuous(labels = function(i) format(10^i, digits = 1, scientific = TRUE)) +
  scale_fill_manual(values = c('#00008B','#BCAF6FFF'),
                    guide = NULL) +
  xlab('\nRealm') +
  ylab('Height Range\n') +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 15, colour = "black"), 
                     axis.text.y = element_text(size = 15, colour = "black"),
                     panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))

# Determine statistical significance of differences between complexity distributions
# Collate required data
dat <- rbind(data.frame(D = RawData %>% filter(Realm == 1) %>% pull(D, as_vector = TRUE), R = RawData %>% filter(Realm == 1) %>% pull(R2, as_vector = TRUE), H = RawData %>% filter(Realm == 1) %>% pull(H2, as_vector = TRUE), Realm = 'T'),
             data.frame(D = RawData %>% filter(Realm == 0) %>% pull(D, as_vector = TRUE), R = RawData %>% filter(Realm == 0) %>% pull(R2, as_vector = TRUE), H = RawData %>% filter(Realm == 0) %>% pull(H2, as_vector = TRUE), Realm = 'M'))
# ensure data is in correct format
dat$Realm <- as.factor(dat$Realm)

# Run analyses
# Fractal Dimension
# Bayesian anova with resampling
FDMod <- ttestBF(formula = D~Realm, data = dat)
FDMod

# Rugosity
# Bayesian anova with resampling
RMod <- ttestBF(formula = R~Realm, data = dat)
RMod

# Height Range
# Bayesian anova with resampling
HRMod <- ttestBF(formula = H~Realm, data = dat)
HRMod

# clear memory
rm(dat)
gc()

# 2. Quantify the pairwise relationships between each complexity variable

# Firstly, visualise the relationships.
# Rugosity and Height range
(RH_plot <- ggplot(plotdata, aes(x=H2, y=R2)) +
    geom_hex(bins = 250) +
    xlab('\nHeight Range') +
    ylab('Rugosity\n') +
    scale_y_continuous(labels = function(i) 10^i) +
    scale_x_continuous(labels = function(i) format(10^i, digits = 1, scientific = TRUE)) +
    scale_fill_gradientn(colours = viridis(n = 50, option = 'cividis'),
                         guide = guide_colorbar(ticks = F, title = 'Density',
                                                reverse = F, label = T,
                                                na.value = "white")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title = element_text(size = 15, colour = 'black'),
                       axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                       legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black")))

# Height range and Fractal Dimension
(DH_plot <- ggplot(plotdata, aes(x=D, y=H2)) +
    geom_hex(bins = 250) +
    xlab('\nFractal Dimension') +
    ylab('Height Range\n') +
    scale_y_continuous(labels = function(i) format(10^i, digits = 1, scientific = TRUE)) +
    scale_fill_gradientn(colours = viridis(n = 50, option = 'cividis'),
                         guide = guide_colorbar(ticks = F, title = 'Density',
                                                reverse = F, label = T,
                                                na.value = "white")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title = element_text(size = 15, colour = 'black'),
                       axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"), 
                       legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black")))

# Rugosity and Fractal Dimension
(DR_plot <- ggplot(plotdata, aes(x=D, y=R2)) +
    geom_hex(bins = 250) +
    xlab('\nFractal Dimension') +
    ylab('Rugosity\n') +
    scale_y_continuous(labels = function(i) 10^i) +
    scale_fill_gradientn(colours = viridis(n = 50, option = 'cividis'),
                         guide = guide_colorbar(ticks = F, title = 'Density',
                                                reverse = F, label = T,
                                                na.value = "white")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title = element_text(size = 15, colour = 'black'),
                       axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                       legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black")))

# Next determine the most appropriate model structure for each of the pairwise analytical comparisons.
# Extract random data sample
TmpData <- RawData %>% select(Lon, Lat, D, R2, H2) %>% collect() %>% sample_n(size = N, replace = FALSE) 

# Compare linear and non-linear formats for each pairwise complexity variable combination (using a quick brms run-through)
# Rugosity versus Height Range
RH_brm_mod <- add_criterion(brm(R2 ~ H2 + Lat + Lon, family = 'gaussian', # GPS details included as fixed effects variables
                                data = TmpData, chains = 1, iter = 1000, warmup = 500, backend = 'cmdstanr'), 'loo')
RH_brm_mod2 <- add_criterion(brm(R2 ~ H2 + I(H2^2) + s(Lat, Lon), family = 'gaussian', 
                                 data = TmpData, chains = 1, iter = 1000, warmup = 500, backend = 'cmdstanr'), 'loo')
# Height Range versus Fractal Dimension
DH_brm_mod <- add_criterion(brm(H2 ~ D + Lat + Lon, family = 'gaussian', 
                                data = TmpData, chains = 1, iter = 1000, warmup = 500, backend = 'cmdstanr'), 'loo')
DH_brm_mod2 <- add_criterion(brm(H2 ~ D + I(D^2) + s(Lat, Lon), family = 'gaussian', 
                                 data = TmpData, chains = 1, iter = 1000, warmup = 500, backend = 'cmdstanr'), 'loo')
# Rugosity versus Fractal Dimension
DR_brm_mod <- add_criterion(brm(R2 ~ D + + Lat + Lon, family = 'gaussian', 
                                data = TmpData, chains = 1, iter = 1000, warmup = 500, backend = 'cmdstanr'), 'loo')
DR_brm_mod2 <- add_criterion(brm(R2 ~ D + I(D^2) + s(Lat, Lon), family = 'gaussian', 
                                 data = TmpData, chains = 1, iter = 1000, warmup = 500, backend = 'cmdstanr'), 'loo')

# Access model fit
loo_compare(RH_brm_mod, RH_brm_mod2) 
loo_compare(DH_brm_mod, DH_brm_mod2)
loo_compare(DR_brm_mod, DR_brm_mod2) # In all cases the none-linear model is afforded more predictive power.

# Store most appropriate model structure 
RH_mod <- bf(R2 ~ H2 + I(H2^2) + s(Lat, Lon))
DH_mod <- bf(H2 ~ D + I(D^2) + s(Lat, Lon))
DR_mod <- bf(R2 ~ D + I(D^2) + s(Lat, Lon))

#--- Run pairwise analyses
# Evaluate relationship using Monte Carlo Resampling
modList <- pblapply(1:mcSim, function(x) { MCbrm(RHbf = RH_mod, DHbf = DH_mod, DRbf = DR_mod, dat = RawData, n = N) })
saveRDS(modList, paste0(DEMPath, "modList.rds")) # save as a checkpoint

# Extract R squared statistics
RH_R <- quiet(ci(unlist(lapply(modList,'[[', 'RH_R'))))
DH_R <- quiet(ci(unlist(lapply(modList,'[[', 'DH_R'))))
DR_R <- quiet(ci(unlist(lapply(modList,'[[', 'DR_R'))))

# Extract predicted relationships for plotting
RH <- quiet(as.data.frame(t(apply(sapply(modList, '[[', 'RH'), 1, ci)))); RH$H2 <- seq(-7,2,0.01); names(RH) <- c('R2','Lower','Upper','SE','H2') # Rugosity and Height Range
DH <- quiet(as.data.frame(t(apply(sapply(modList, '[[', 'DH'), 1, ci)))); DH$D <- seq(2,3,0.001); names(DH) <- c('H2','Lower','Upper','SE','D') # Fractal Dimension and Height Range
DR <- quiet(as.data.frame(t(apply(sapply(modList, '[[', 'DR'), 1, ci)))); DR$D <- seq(2,3,0.001); names(DR) <- c('R2','Lower','Upper','SE','D') # Fractal Dimension and Rugosity

# Add details to plots
# Rugosity and Height Range
RH_plot + 
  geom_ribbon(aes(x = H2, ymin = Lower, ymax = Upper), fill = 'red', alpha = 1, data = RH, inherit.aes = FALSE) +
  geom_line(aes(x = H2, y = R2), col = 'red', linetype = 'solid', linewidth = 0.5, data = RH, inherit.aes = FALSE) # adding the mean line ensures a clean line across the plot.

# Fractal Dimension and Height range
DH_plot +
  geom_ribbon(aes(x = D, ymin = Lower, ymax = Upper), fill = 'red', alpha = 1, data = DH, inherit.aes = FALSE) +
  geom_line(aes(x = D, y = H2), col = 'red', linetype = 'solid', linewidth = 0.5, data = DH, inherit.aes = FALSE)

# Fractal Dimension and Height Range
DR_plot +
  geom_ribbon(aes(x = D, ymin = Lower, ymax = Upper), fill = 'red', alpha = 1, data = DR, inherit.aes = FALSE) +
  geom_line(aes(x = D, y = R2), col = 'red', linetype = 'solid', linewidth = 0.5, data = DR, inherit.aes = FALSE)

# ----------------------------------------------------- End of Code ---------------------------------