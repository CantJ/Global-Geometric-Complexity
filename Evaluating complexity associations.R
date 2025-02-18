# This script will explore the correlation between measures of Geometric complexity and Geodiversity, Human population density, and climate variability

# Primary Author: James Cant
# Date last modified: July 2024
# ---------------------------------------------------------------

# Clear working directory
rm(list = ls())

# Load required packages
library(terra)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(pbapply)
library(ggeffects)
library(gmodels)
library(spsUtil)
library(brms)
library(ggdist)
library(mgcv)

# Set random number seed
set.seed(457034)

# Modify R digit printing to prevent the rounding of very close values
options(digits = 22)

# Define file pathways
RastPath <- '~/Documents/Complexity Analyses/Final Data Products/'

# Define sample size/and iteration number for resampling analyses
N <- 1000
iter <- 1000

# Define Monte Carlo resampling Generalized Additive Model function
mc_GAM <- function(dat, n) {
  
  # Extract random data sample
  TmpData <- dat %>% select(Lon, Lat, D, R2, H2, G) %>% collect() %>% sample_n(size = n, replace = FALSE) 
  
  # Run GAM regression comparing the relationship between each measure of complexity and geodiversity
  mod1 <- bam(R2 ~ s(G) + s(Lon,Lat), family = gaussian(), data = TmpData)
  mod2 <- bam(D ~ s(G) + s(Lon,Lat), family = gaussian(), data = TmpData)
  mod3 <- bam(H2 ~ s(G) + s(Lon,Lat), family = gaussian(), data = TmpData)
  
  # Extract model fits
  mod1R <- summary(mod1)$r.sq
  mod2R <- summary(mod2)$r.sq
  mod3R <- summary(mod3)$r.sq
  
  # Predict complexity values using assessed relationship
  pred1 <- ggpredict(mod1, terms = c('G[1:21]'))$predicted
  pred2 <- ggpredict(mod2, terms = c('G[1:21]'))$predicted
  pred3 <- ggpredict(mod3, terms = c('G[1:21]'))$predicted
  
  # and return outputs
  return(list(R = pred1, D = pred2, H = pred3, R_R = mod1R, D_R = mod2R, H_R = mod3R))
}

# Define monte carlo brms regression function
MCbrm <- function(Hbf, Rbf, Dbf, dat, n, measure) {
  
  # Extract random data sample
  TmpData <- dat %>% select(Lon, Lat, D, R2, H2, PD2, CV2) %>% collect() %>% sample_n(size = n, replace = FALSE) 
  
  if(measure == 'PD'){
    # Run Bayesian regression model for each pairwise combination
    mod1 <- quiet(brm(Hbf, family = hurdle_gamma(), data = TmpData, chains = 1, iter = 2000, warmup = 1000, backend = 'cmdstanr'))
    mod2 <- quiet(brm(Rbf, family = hurdle_gamma(), data = TmpData, chains = 1, iter = 2000, warmup = 1000, backend = 'cmdstanr'))
    mod3 <- quiet(brm(Dbf, family = hurdle_gamma(), data = TmpData, chains = 1, iter = 2000, warmup = 1000, backend = 'cmdstanr'))
  } else {
    mod1 <- quiet(brm(Hbf, family = gaussian(), data = TmpData, chains = 1, iter = 2000, warmup = 1000, backend = 'cmdstanr'))
    mod2 <- quiet(brm(Rbf, family = gaussian(), data = TmpData, chains = 1, iter = 2000, warmup = 1000, backend = 'cmdstanr'))
    mod3 <- quiet(brm(Dbf, family = gaussian(), data = TmpData, chains = 1, iter = 2000, warmup = 1000, backend = 'cmdstanr'))
  }
  
  # Extract model fits
  mod1R <- bayes_R2(mod1)[1]
  mod2R <- bayes_R2(mod2)[1]
  mod3R <- bayes_R2(mod3)[1]
  
  if(measure == 'PD'){
    # Extract hurdle coefficients (probability of entry being zero)
    H_hu <- c(fixef(mod1)[2], fixef(mod1)[4])
    R_hu <- c(fixef(mod2)[2], fixef(mod2)[4])
    D_hu <- c(fixef(mod3)[2], fixef(mod3)[4])
  }
  
  # Predict population densities (when not zero) using assessed relationship
  pred1 <- ggpredict(mod1, terms = list(H2 = seq(-7,1.2, 0.01)))$predicted # covers the full range of H2 
  pred2 <- ggpredict(mod2, terms = list(R2 = seq(-13.7,0.8, 0.01)))$predicted # covers the full range of R2
  pred3 <- ggpredict(mod3, terms = list(D = seq(2,3,0.01)))$predicted # covers full range of D
  
  if(measure == 'PD'){
    # Collate and return outputs
    return(list(H = pred1, R = pred2, D = pred3, H_hu = H_hu, R_hu = R_hu, D_hu = D_hu, H_R = mod1R, R_R = mod2R, D_R = mod3R))
  } else {
    return(list(H = pred1, R = pred2, D = pred3, H_R = mod1R, R_R = mod2R, D_R = mod3R))
  }
}

######################################
# STEP 1: Import & Clean Raw Data
######################################

# Load in the Complexity and Geodiversity Rasters
DRast <- rast(paste0(RastPath, 'GlobalFractalDimension.tif'))
RRast <- rast(paste0(RastPath, 'GlobalRugosity.tif'))
HRast <- rast(paste0(RastPath, 'GlobalHeightRange.tif'))
GRast <- rast(paste0(RastPath, 'Geodiversity_Mollweide_1870m.tif'))  
# Load raster of Human Population size
PopRast <- rast(paste0(RastPath, 'PopCount_Mollweide_1870m.tif'))
# and climate variability
ClimRast <- rast(paste0(RastPath, 'ClimVar_Mollweide_1870m.tif'))

# Ensure resolutions and extents match across the geodiversity and geometric complexity rasters
ext(GRast) <- ext(DRast)
ext(PopRast) <- ext(DRast)
ext(ClimRast) == ext(DRast)
res(GRast) == res(DRast)
res(PopRast) == res(DRast)
res(ClimRast) == res(DRast)
ext(RRast) == ext(DRast)
ext(HRast) == ext(DRast)

# Remove values from the geometric complexity rasters for which there are no corresponding geodiversity values (Focus on terrestrial locations)
DRast <- mask(DRast, GRast, maskvalues = NA)
RRast <- mask(RRast, GRast, maskvalues = NA)
HRast <- mask(HRast, GRast, maskvalues = NA)

# # Estimate spatial autocorrelation of terrestrial complexity measures.
# Define nearest neighbor matrix (queens format - diagonals included)
w <- matrix(c(1,1,1,1,0,1,1,1,1), nrow = 3)

#test for spatial autocorrelation using Moran's I (-1 < I > 1)
# Estimate Moran's I (can take a while to complete)
D_I <- terra::autocor(DRast, w, method = 'moran', global = TRUE)
R_I <- terra::autocor(RRast, w, method = 'moran', global = TRUE)  
H_I <- terra::autocor(HRast, w, method = 'moran', global = TRUE) 


# collate raster variables together in a single data frame
# Extract GPS coordinates
RawData <- crds(DRast, df = TRUE, na.rm = FALSE)
# rename columns 
names(RawData) <- c('Lon', 'Lat')
# add in complexity estimates
RawData$D <- as.numeric(values(DRast))
RawData$R <- as.numeric(values(RRast))
RawData$H <- as.numeric(values(HRast))
RawData$G <- as.numeric(values(GRast))
RawData$PD <- as.numeric(values(PopRast))
RawData$CV <- as.numeric(values(ClimRast))
# remove incomplete entries
RawData <- RawData[complete.cases(RawData),]

# Clean data
RawData %>% 
  summarise(Hmin = min(H), # shouldn't be below 0
            Rmin = min(R), # shouldn't be below 0
            Dmin = min(D), # shouldn't be below 2
            Gmin = min(G), # shouldn't be below 0
            PDmin = min(PD), # shouldn't be below 0
            CVmin = min(CV), # shouldn't be below 0
            Hmax = max(H),
            Rmax = max(R),
            Dmax = max(D), # shouldn't be higher than 3
            Gmax = max(G),
            PDmax = max(PD),
            CVmax = max(CV))

# Remove non-fractal entries
RawData %>% filter(D==2) %>% dim()
RawData <- RawData %>% mutate(R = replace(R, D==2, NA)) # drop corresponding rugosity entries
RawData <- RawData %>% mutate(H = replace(H, D==2, NA)) # drop corresponding height range entries
RawData <- RawData %>% mutate(D = replace(D, D==2, NA)) # finally drop fractal dimension entries
# again tidy away incomplete entries
RawData <- RawData[complete.cases(RawData),]

# Height range
RawData %>% with(hist(H)) # Height Range requires transformation
RawData %>% with(hist(log10(H)))
# Rugosity
RawData %>% with(hist(R)) # Rugosity requires transformation
RawData %>% with(hist(log10(R)))
# Fractal Dimension
RawData %>% with(hist(D)) 
# Geodiversity
RawData %>% with(hist(G))
# Population density
RawData %>% with(hist(PD)) # Population density requires transformation
# this needs to be done carefully to maintain the position of zero entries.
RawData %>% filter(PD>0) %>% with(hist(PD^0.05))
# Climate variability
RawData %>% with(hist(CV))
RawData %>% with(hist(log10(CV+1)))

# Apply necessary transformations
RawData <- RawData %>%
  mutate(H2 = log10(H),
         R2 = log10(R),
         PD2 = PD,
         CV2 = log10(CV+1))

#transform non-zero entries in the population density variable
c = 0.05
RawData$PD2[RawData$PD2 > 0] <- RawData$PD2[RawData$PD2 > 0]^c
# confirm desired effect
RawData %>% 
  summarise(PD2min = min(PD2, na.rm = T), # should still be zero
            PD2max = max(PD2, na.rm = T))

######################################
# STEP 2: Compare Geodiversity and Geometric complexity
######################################

# Plot relationships
# Define sample size for plotting (This is to speed up plot processing times)
Nsize <- 2e+07 
tmpData <- RawData %>% select(D, R2, H2, G) %>% collect() %>% sample_n(size = Nsize, replace = FALSE) # isolate plotting sample
# Fractal Dimension
(GDplot <- ggplot(tmpData, aes(x = G, y = D, fill = factor(G))) +
    stat_halfeye(adjust = 0.5, show.legend = FALSE, .width = 0, point_colour = NA) +
    xlab(NULL) +
    ylab(NULL) +
    coord_cartesian(ylim = c(2, 2.8)) + 
    scale_y_continuous(breaks = c(2.0,2.2,2.4,2.6,2.8), labels = c(2.0,2.2,2.4,2.6,2.8)) +
    scale_fill_manual(values = magma(n = 20)) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title = element_text(size = 15, colour = 'black'),
                       axis.text.x = element_text(size = 30, colour = "black"), 
                       axis.text.y = element_text(size = 30, colour = "black"),
                       panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')))

# Rugosity
(GRplot <- ggplot(tmpData, aes(x = G, y = R2, fill = factor(G))) +
    stat_halfeye(adjust = 0.5, show.legend = FALSE, .width = 0, point_colour = NA) +
    xlab(NULL) +
    ylab(NULL) +
    scale_y_continuous(breaks = c(log10(1), log10(1e-05), log10(1e-10)), labels = function(i) 10^i) +
    coord_cartesian(ylim = c(-11, 2)) + 
    scale_fill_manual(values = magma(n = 20)) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title = element_text(size = 15, colour = 'black'),
                       axis.text.x = element_text(size = 30, colour = "black"), 
                       axis.text.y = element_text(size = 30, colour = "black"),
                       panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')))

# Height Range
(GHplot <- ggplot(tmpData, aes(x = G, y = H2, fill = factor(G))) +
    stat_halfeye(adjust = 0.5, show.legend = FALSE, .width = 0, point_colour = NA) +
    xlab(NULL) +
    ylab(NULL) +
    scale_y_continuous(breaks = c(log10(10), log10(1), log10(1e-02), log10(1e-4)), labels = function(i) 10^i) +
    coord_cartesian(ylim = c(-4.1, 1.6)) + 
    scale_fill_manual(values = magma(n = 20)) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title = element_text(size = 15, colour = 'black'),
                       axis.text.x = element_text(size = 30, colour = "black"), 
                       axis.text.y = element_text(size = 30, colour = "black"),
                       panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')))

# Quantify relationships

# Run Generalized Additive Models on full data set using resampling to evaluate the relationships between each complexity variable and geodiversity
GAMlist <- pblapply(1:iter, function(i) { mc_GAM(dat = RawData, n = N) })

# Extract R squared statistics
R_R <- quiet(ci(unlist(lapply(GAMlist,'[[', 'R_R'))))
D_R <- quiet(ci(unlist(lapply(GAMlist,'[[', 'D_R'))))
H_R <- quiet(ci(unlist(lapply(GAMlist,'[[', 'H_R'))))

# Extract predicted relationships for plotting
RG <- quiet(as.data.frame(t(apply(sapply(GAMlist, '[[', 'R'), 1, ci)))); RG$G <- as.numeric(1:21); names(RG) <- c('R2','Lower','Upper','SE','G') # Rugosity
DG <- quiet(as.data.frame(t(apply(sapply(GAMlist, '[[', 'D'), 1, ci)))); DG$G <- as.numeric(1:21); names(DG) <- c('D','Lower','Upper','SE','G') # Fractal Dimension
HG <- quiet(as.data.frame(t(apply(sapply(GAMlist, '[[', 'H'), 1, ci)))); HG$G <- as.numeric(1:21); names(HG) <- c('H2','Lower','Upper','SE','G') # Height Range

# Add details to plots
# Rugosity
GRplot + 
  geom_ribbon(aes(x = G, ymin = Lower, ymax = Upper), fill = 'Black', alpha = 0.9, data = RG, inherit.aes = FALSE)
# Fractal Dimension
GDplot + 
  geom_ribbon(aes(x = G, ymin = Lower, ymax = Upper), fill = 'Black', alpha = 0.9, data = DG, inherit.aes = FALSE) 
# Height Range
GHplot + 
  geom_ribbon(aes(x = G, ymin = Lower, ymax = Upper), fill = 'Black', alpha = 0.9, data = HG, inherit.aes = FALSE) 


######################################
# STEP 3: Population density and Geometric complexity
######################################

# Model the relationship between population density and each of the geometric complexity variables
# Define model structures for implementing hurdle gamma regression
Hbf <- bf(PD2 ~ H2 + s(Lat,Lon), # pattern across non-zero entries
          hu ~ H2) # probability of zero entries
Rbf <- bf(PD2 ~ R2 + s(Lat,Lon),
          hu ~ R2)
Dbf <- bf(PD2 ~ D + s(Lat,Lon),
          hu ~ D)

# Run Bayesian approach with resampling used to evaluate the relationships between each complexity variable and population density
PDlist <- pblapply(1:iter, function(i) { MCbrm(Hbf = Hbf, Rbf = Rbf, Dbf = Dbf, dat = RawData, n = N, measure = 'PD') })
saveRDS(PDlist, file = paste0(RastPath, "PDlist.rds")) # checkpoint

# Extract R squared statistics
H_R <- quiet(ci(unlist(lapply(PDlist,'[[', 'H_R'))))
R_R <- quiet(ci(unlist(lapply(PDlist,'[[', 'R_R'))))
D_R <- quiet(ci(unlist(lapply(PDlist,'[[', 'D_R'))))

#Extract probability of zero coefficients
coefR <- quiet(as.data.frame(t(apply(sapply(PDlist, '[[', 'R_hu'), 1, ci))))
coefH <- quiet(as.data.frame(t(apply(sapply(PDlist, '[[', 'H_hu'), 1, ci))))
coefD <- quiet(as.data.frame(t(apply(sapply(PDlist, '[[', 'D_hu'), 1, ci))))
# Compute mean and variance bounds in each modeled relationship (transform to probability of >0 )
Rvec <- seq(-13.7,0.8, 0.01); Dvec <- seq(2,3, 0.01); Hvec <- seq(-7,1.2, 0.01)
R_hu <- data.frame(R = Rvec, Mean = 1-(exp(coefR[1,1]+coefR[2,1]*Rvec)), Lower = 1-(exp(coefR[1,2]+coefR[2,2]*Rvec)), Upper = 1-(exp(coefR[1,3]+coefR[2,3]*Rvec)))
D_hu <- data.frame(D = Dvec, Mean = 1-(exp(coefD[1,1]+coefD[2,1]*Dvec)), Lower = 1-(exp(coefD[1,2]+coefD[2,2]*Dvec)), Upper = 1-(exp(coefD[1,3]+coefD[2,3]*Dvec)))
H_hu <- data.frame(H = Hvec, Mean = 1-(exp(coefH[1,1]+coefH[2,1]*Hvec)), Lower = 1-(exp(coefH[1,2]+coefH[2,2]*Hvec)), Upper = 1-(exp(coefH[1,3]+coefH[2,3]*Hvec)))

# Extract predicted relationships for plotting
RPD <- quiet(as.data.frame(t(apply(sapply(PDlist, '[[', 'R'), 1, ci)))); RPD$R <- Rvec; names(RPD) <- c('PD','Lower','Upper','SE','R') # Rugosity
DPD <- quiet(as.data.frame(t(apply(sapply(PDlist, '[[', 'D'), 1, ci)))); DPD$D <- Dvec; names(DPD) <- c('PD','Lower','Upper','SE','D') # Fractal Dimension
HPD <- quiet(as.data.frame(t(apply(sapply(PDlist, '[[', 'H'), 1, ci)))); HPD$H <- Hvec; names(HPD) <- c('PD','Lower','Upper','SE','H') # Height Range

# Plot modeled relationships
# Rugosity
# determine axis scaling factor
dR <- max(RPD$Upper)/max(R_hu$Mean)
# generate plot
ggplot(data = RPD) +
  geom_ribbon(aes(x = R, ymin = Lower, ymax = Upper), fill = '#0000FF', alpha = 0.2) +
  geom_line(aes(x = R, y = PD), col = '#000099', linewidth = 1.1) +
  geom_ribbon(aes(x = R, ymin = Lower*dR, ymax = Upper*dR), fill = '#F90707', alpha = 0.2, data = R_hu) +
  geom_line(aes(x = R, y = Mean*dR), col = '#A80000', linewidth = 1.1, data = R_hu) +
  xlab(NULL) +
  scale_x_continuous(expand = c(0,0), breaks = c(log10(1), log10(1e-05), log10(1e-10)), labels = function(i) 10^i) +
  scale_y_continuous(name = NULL,
                     labels = function(i) {format(i, digits = 2, scientific = FALSE)},
                     sec.axis = sec_axis(~./dR,
                                         name = NULL,
                                         labels = function(i) round(i, 3))) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))

# Fractal Dimension
dD <- max(DPD$Upper)/max(D_hu$Mean)
# generate plot
ggplot(data = DPD) +
  geom_ribbon(aes(x = D, ymin = Lower*dD, ymax = Upper*dD), fill = '#F90707', alpha = 0.2, data = D_hu) +
  geom_line(aes(x = D, y = Mean*dD), col = '#A80000', linewidth = 1.1, data = D_hu) +
  geom_ribbon(aes(x = D, ymin = Lower, ymax = Upper), fill = '#0000FF', alpha = 0.2) +
  geom_line(aes(x = D, y = PD), col = '#000099', linewidth = 1.1) +
  xlab(NULL) +
  scale_y_continuous(name = NULL,
                     labels = function(i) {format(i, digits = 2, scientific = FALSE)},
                     sec.axis = sec_axis(~./dD,
                                         name = NULL,
                                         labels = function(i) round(i, 3))) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))

# Height Range
dH <- max(HPD$Upper)/max(H_hu$Mean)
# generate plot
ggplot(data = HPD) +
  geom_ribbon(aes(x = H, ymin = Lower, ymax = Upper), fill = '#0000FF', alpha = 0.2) +
  geom_line(aes(x = H, y = PD), col = '#000099', linewidth = 1.1) +
  geom_ribbon(aes(x = H, ymin = Lower*dH, ymax = Upper*dH), fill = '#F90707', alpha = 0.2, data = H_hu) +
  geom_line(aes(x = H, y = Mean*dH), col = '#A80000', linewidth = 1.1, data = H_hu) +
  xlab(NULL) +
  scale_x_continuous(expand = c(0,0), breaks = c(log10(10), log10(1), log10(1e-02), log10(1e-4), log10(1e-6)), labels = function(i) 10^i) +
  scale_y_continuous(name = NULL,
                     labels = function(i) {format(i, digits = 2, scientific = FALSE)},
                     sec.axis = sec_axis(~./dH,
                                         name = NULL,
                                         labels = function(i) round(i, 3))) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))

######################################
# STEP 4: Compare Climate variability and Geometric complexity
######################################

# Model the relationship between climate variability and each of the geometric complexity variables
Hcvbf <- bf(CV2 ~ H2 + s(Lat,Lon))
Rcvbf <- bf(CV2 ~ R2 + s(Lat,Lon))
Dcvbf <- bf(CV2 ~ D + s(Lat,Lon))

# Run Bayesian approach with resampling used to evaluate the relationships between each complexity variable and population density
CVmods <- pblapply(1:iter, function(i) { MCbrm(Hbf = Hcvbf, Rbf = Rcvbf, Dbf = Dcvbf, dat = RawData, n = N, measure = 'CV') })
saveRDS(CVmods, file = paste0(FilePath, "CVmods.rds")) # checkpoint

# Extract R squared statistics
(cvR_R <- quiet(ci(unlist(lapply(CVmods,'[[', 'R_R')))))
(cvD_R <- quiet(ci(unlist(lapply(CVmods,'[[', 'D_R')))))
(cvH_R <- quiet(ci(unlist(lapply(CVmods,'[[', 'H_R')))))

# Extract predicted relationships for plotting
Rvec <- seq(-13.7,0.8, 0.01); Dvec <- seq(2,3, 0.01); Hvec <- seq(-7,1.2, 0.01)
RCV <- quiet(as.data.frame(t(apply(sapply(CVmods, '[[', 'R'), 1, ci)))); RCV$R <- Rvec; names(RCV) <- c('CV','Lower','Upper','SE','R') # Rugosity
DCV <- quiet(as.data.frame(t(apply(sapply(CVmods, '[[', 'D'), 1, ci)))); DCV$D <- Dvec; names(DCV) <- c('CV','Lower','Upper','SE','D') # Fractal Dimension
HCV <- quiet(as.data.frame(t(apply(sapply(CVmods, '[[', 'H'), 1, ci)))); HCV$H <- Hvec; names(HCV) <- c('CV','Lower','Upper','SE','H') # Height Range

# Plot relationships
# Rugosity
(CVRplot <- ggplot(RCV) +
    geom_ribbon(aes(x = R, ymin = Lower, ymax = Upper), fill = '#fc8961', alpha = 0.5) +
    geom_line(aes(x = R, y = CV), col = '#000004', linewidth = 1.1, linetype = 'dashed') + 
    xlab('\nRugosity') +
    ylab('Climate Variability\n') +
    scale_x_continuous(expand = c(0,0), labels = function(i) format((10^i), digits = 2)) +
    scale_y_continuous(limits = c(log10(0.65+1),log10(1.13+1)), labels = function(i) format((10^i)-1, digits = 2)) +
    scale_fill_gradientn(colours = viridis(n = 50, option = 'cividis'),
                         guide = guide_colorbar(ticks = F, title = 'Density',
                                                reverse = F, label = T,
                                                na.value = "white")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title = element_text(size = 15, colour = 'black'),
                       axis.text.x = element_text(size = 15, colour = "black"), 
                       axis.text.y = element_text(size = 15, colour = "black"),
                       panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')))

# Height Range
(CVHplot <- ggplot(HCV) +
    geom_ribbon(aes(x = H, ymin = Lower, ymax = Upper), fill = '#fc8961', alpha = 0.7) +
    geom_line(aes(x = H, y = CV), col = '#000004', linewidth = 1.1, linetype = 'dashed') + 
    xlab('\nHeight Range') +
    ylab('Climate Variability\n') +
    scale_x_continuous(expand = c(0,0), breaks = c(log10(10), log10(1), log10(1e-02), log10(1e-4), log10(1e-6)), labels = function(i) 10^i) +
    scale_y_continuous(limits = c(log10(0.65+1),log10(1.13+1)), labels = function(i) format((10^i)-1, digits = 2)) +
    scale_fill_gradientn(colours = viridis(n = 50, option = 'cividis'),
                         guide = guide_colorbar(ticks = F, title = 'Density',
                                                reverse = F, label = T,
                                                na.value = "white")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title = element_text(size = 15, colour = 'black'),
                       axis.text.x = element_text(size = 15, colour = "black"), 
                       axis.text.y = element_text(size = 15, colour = "black"),
                       panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')))

# Fractal Dimension
(CVDplot <- ggplot(DCV) +
    geom_ribbon(aes(x = D, ymin = Lower, ymax = Upper), fill = '#fc8961', alpha = 0.5) +
    geom_line(aes(x = D, y = CV), col = '#000004', linewidth = 1.1, linetype = 'dashed') + 
    xlab('\nFractal Dimension') +
    ylab('Climate Variability\n') +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(limits = c(log10(0.6+1),log10(1.13+1)), labels = function(i) format((10^i)-1, digits = 2)) +
    scale_fill_gradientn(colours = viridis(n = 50, option = 'cividis'),
                         guide = guide_colorbar(ticks = F, title = 'Density',
                                                reverse = F, label = T,
                                                na.value = "white")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title = element_text(size = 15, colour = 'black'),
                       axis.text.x = element_text(size = 15, colour = "black"), 
                       axis.text.y = element_text(size = 15, colour = "black"),
                       panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')))

# ----------------------------------------------------- End of Code ---------------------------------