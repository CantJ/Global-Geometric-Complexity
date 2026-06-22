# This script will explore the correlation between measures of Geometric complexity and Geodiversity, Human population density, and climate variability

# Primary Author: James Cant
# Date last modified: July 2024
# ---------------------------------------------------------------

# Set random number seed
set.seed(457034)

# Define sample size/and iteration number for resampling analyses
N <- 100000
iter <- 100

# Define Monte Carlo brms regression function
MCbrm <- function(dat, n, measure) {
  
  # Extract random data sample
  TmpData <- dat %>% dplyr::select(Lon, Lat, D, R2, H2, G, PD2) %>% collect() %>% sample_n(size = n, replace = FALSE) 
  if(measure == 'PD'){
  # Separate zero and non-zero entries
    PDdat <- TmpData %>% filter(PD2>0)
    TmpData$Zero <- NA
    TmpData$Zero[TmpData$PD2 > 0] <- 0
    TmpData$Zero[TmpData$PD2 == 0] <- 1
  }
  
  # Run Bayesian regression model for each pairwise combination
  if(measure == 'PD'){
    mod1 <- quiet(qbrms(formula = PD2 ~ H2 + Lat + Lon, family = gaussian(), data = PDdat))
    mod2 <- quiet(qbrms(formula = PD2 ~ R2 + Lat + Lon, family = gaussian(), data = PDdat))
    mod3 <- quiet(qbrms(formula = PD2 ~ D + Lat + Lon, family = gaussian(), data = PDdat))
  }
  if(measure == 'GD'){
    mod1 <- quiet(qbrms(formula = H2 ~ G + Lon + Lat, family = gaussian(), data = TmpData))
    mod2 <- quiet(qbrms(formula = R2 ~ G + Lon + Lat, family = gaussian(), data = TmpData))
    mod3 <- quiet(qbrms(formula = D ~ G + Lon + Lat, family = gaussian(), data = TmpData))
  } 
  
  # Extract model fits (R2)
  mod1R <- bayes_R2(mod1, verbose = F)[1]
  mod2R <- bayes_R2(mod2, verbose = F)[1]
  mod3R <- bayes_R2(mod3, verbose = F)[1]

  # Extract hurdle coefficients as required (probability of entry being zero)
  if(measure == 'PD'){
    zero1 <- quiet(qbrms(formula = Zero ~ H2 + Lat + Lon, family = binomial(), data = TmpData))
    zero2 <- quiet(qbrms(formula = Zero ~ R2 + Lat + Lon, family = binomial(), data = TmpData))
    zero3 <- quiet(qbrms(formula = Zero ~ D + Lat + Lon, family = binomial(), data = TmpData))
    # Extract coefficients
    H_hu <- c(coef(zero1)[1], coef(zero1)[2])
    R_hu <- c(coef(zero2)[1], coef(zero2)[2])
    D_hu <- c(coef(zero3)[1], coef(zero3)[2])
  }
  
  # Predict mean relationship
  if(measure == 'GD'){
    pred1 <- conditional_effects(mod1)$G 
    pred2 <- conditional_effects(mod2)$G
    pred3 <- conditional_effects(mod3)$G
  }
  if(measure == 'PD'){
    pred1 <- conditional_effects(mod1)$H2 
    pred2 <- conditional_effects(mod2)$R2
    pred3 <- conditional_effects(mod3)$D
  }

  # Collate and return outputs
  if(measure == 'PD'){
    return(list(H = pred1, R = pred2, D = pred3, H_hu = H_hu, R_hu = R_hu, D_hu = D_hu, H_R = mod1R, R_R = mod2R, D_R = mod3R))
  } 
  if(measure == 'PD') {
    return(list(H = pred1, R = pred2, D = pred3, H_R = mod1R, R_R = mod2R, D_R = mod3R))
  }
}

# Function for computing vector confidence intervals.
CI95 <- function(vector) {
  # Standard deviation of sample
  s <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vmean <- mean(vector)
  # Error according to t distribution
  error <- qt(0.975, df = n - 1) * s / sqrt(n)
  # summary stats as a vector
  result <- c('Mean' = vmean, '2_5%' = vmean - error, '97_5%' = vmean + error)
  return(result)
}

######################################
# STEP 1: Import & Clean Raw Data
######################################

# Load in the Complexity and Geodiversity Rasters
DRast <- rast(paste0(FilePath, 'GlobalFractalDimension.tif'))
RRast <- rast(paste0(FilePath, 'GlobalRugosity.tif'))
HRast <- rast(paste0(FilePath, 'GlobalHeightRange.tif'))
GRast <- rast(paste0(FilePath, 'Geodiversity_Mollweide_1870m.tif'))  
# Load raster of Human Population size
PopRast <- rast(paste0(FilePath, 'PopCount_Mollweide_1870m.tif'))
# and climate variability
ClimRast <- rast(paste0(FilePath, 'ClimVar_Mollweide_1870m.tif'))

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
tmpData <- RawData %>% dplyr::select(D, R2, H2, G) %>% collect() %>% sample_n(size = Nsize, replace = FALSE) # isolate plotting sample
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
# Determine the most appropriate model structure for the relationships between each complexity measure and geodiversity 
if(FirstRun == TRUE){
  # Extract random data sample
  TmpData <- RawData %>% dplyr::select(Lon, Lat, D, R2, H2, G) %>% collect() %>% sample_n(size = N, replace = FALSE)
  
  # Compare linear and non-linear formats for each complexity variable GEODIVERSITY combination
  # Rugosity
  RG_brm_mod <- qbrms(formula = R2 ~ G + Lat + Lon, family = gaussian(), # GPS details included as fixed effects variables to accommodate spatial autocorrelation
                      data = TmpData)
  RG_brm_mod2 <- qbrms(formula = R2 ~ G + I(G^2) + Lat + Lon, family = gaussian(),
                       data = TmpData)
  # Height Range
  HG_brm_mod <- qbrms(formula = H2 ~ G + Lat + Lon, family = gaussian(),
                      data = TmpData)
  HG_brm_mod2 <- qbrms(formula = H2 ~ G + I(G^2) + Lat + Lon, family = gaussian(),
                       data = TmpData)
  # Fractal Dimension
  DG_brm_mod <- qbrms(formula = D ~ G + Lat + Lon, family = gaussian(),
                      data = TmpData)
  DG_brm_mod2 <- qbrms(formula = D ~ G + I(G^2) + Lat + Lon, family = gaussian(),
                       data = TmpData)
  
  # Access model fit
  loo_compare(RG_brm_mod, RG_brm_mod2) 
  loo_compare(HG_brm_mod, HG_brm_mod2) 
  loo_compare(DG_brm_mod, DG_brm_mod2) # in all cases there is no significant difference between fits, linear will be used as it is a simpler model
}

# Run Bayesian models to evaluate the relationships between each complexity variable and geodiversity
GeoList <- lapply(1:iter, function(i) { print(i)
  MCbrm(dat = RawData, n = N, measure = 'GD') })
saveRDS(GeoList, paste0(FilePath, "GeoList.rds")) # save as a checkpoint

# Extract R squared statistics
(R_R <- CI95(unlist(lapply(GeoList,'[[', 'R_R'))))
(D_R <- CI95(unlist(lapply(GeoList,'[[', 'D_R'))))
(H_R <- CI95(unlist(lapply(GeoList,'[[', 'H_R'))))

# Extract predicted relationships for plotting
RG <- do.call(rbind, lapply(GeoList, '[[', 'R')) %>% arrange(G); names(RG) <- c('G','R2','Lower','Upper') 
DG <- do.call(rbind, lapply(GeoList, '[[', 'D')) %>% arrange(G); names(DG) <- c('G','D','Lower','Upper') 
HG <- do.call(rbind, lapply(GeoList, '[[', 'H')) %>% arrange(G); names(HG) <- c('G','H2','Lower','Upper')

# Add details to plots
# Rugosity
# determine smoothed relationships
RG$LowerSmooth <- predict(loess(Lower ~ G, data = RG), RG$G)
RG$UpperSmooth <- predict(loess(Upper ~ G, data = RG), RG$G)
RG$MeanSmooth <- predict(loess(R2 ~ G, data = RG), RG$G)
GRplot + 
  geom_ribbon(aes(x = G, ymin = LowerSmooth, ymax = UpperSmooth), fill = 'black', alpha = 0.2, data = RG, inherit.aes = FALSE) +
  geom_line(aes(x = G, y = MeanSmooth), col = 'black', linewidth = 2, data = RG, inherit.aes = FALSE)
# Fractal Dimension
# determine smoothed relationships
DG$LowerSmooth <- predict(loess(Lower ~ G, data = DG), DG$G)
DG$UpperSmooth <- predict(loess(Upper ~ G, data = DG), DG$G)
DG$MeanSmooth <- predict(loess(D ~ G, data = DG), DG$G)
GDplot + 
  geom_ribbon(aes(x = G, ymin = LowerSmooth, ymax = UpperSmooth), fill = 'black', alpha = 0.2, data = DG, inherit.aes = FALSE) +
  geom_line(aes(x = G, y = MeanSmooth), col = 'black', linewidth = 2, data = DG, inherit.aes = FALSE)
# Height Range
# determine smoothed relationships
HG$LowerSmooth <- predict(loess(Lower ~ G, data = HG), HG$G)
HG$UpperSmooth <- predict(loess(Upper ~ G, data = HG), HG$G)
HG$MeanSmooth <- predict(loess(H2 ~ G, data = HG), HG$G)
GHplot + 
  geom_ribbon(aes(x = G, ymin = LowerSmooth, ymax = UpperSmooth), fill = 'black', alpha = 0.2, data = HG, inherit.aes = FALSE) +
  geom_line(aes(x = G, y = MeanSmooth), col = 'black', linewidth = 2, data = HG, inherit.aes = FALSE)

######################################
# STEP 3: Population density and Geometric complexity
######################################

# Model the relationship between population density and each of the geometric complexity variables
# A modified hurdle modelling approach will be used to deal with the large numbers of zeros in the data.
# Firstly however, what is the most appropriate model fit for the non-zero data.
if(FirstRun == TRUE){
  # Extract random data sample whilst dropping zero entries
  TmpData <- RawData %>% dplyr::select(Lon, Lat, D, R2, H2, PD2) %>% filter(PD2>0) %>% collect() %>% sample_n(size = N, replace = FALSE)
  
  # Compare linear and non-linear formats for each pairwise complexity variable-population density combination.
  # Rugosity
  R_mod <- qbrms(formula = PD2 ~ R2 + Lat + Lon, family = gaussian(),
                      data = TmpData)
  R_mod2 <- qbrms(formula = PD2 ~ R2 + I(R2^2) + Lat + Lon, family = gaussian(),
                       data = TmpData)
  # Height Range
  H_mod <- qbrms(formula = PD2 ~ H2 + Lat + Lon, family = gaussian(),
                      data = TmpData)
  H_mod2 <- qbrms(formula = PD2 ~ H2 + I(D^2) + Lat + Lon, family = gaussian(),
                       data = TmpData)
  # Fractal Dimension
  D_mod <- qbrms(formula = PD2 ~ D + Lat + Lon, family = gaussian(),
                      data = TmpData)
  D_mod2 <- qbrms(formula = PD2 ~ D + I(D^2) + Lat + Lon, family = gaussian(),
                       data = TmpData)
  
  # Access model fit
  loo_compare(R_mod, R_mod2) 
  loo_compare(H_mod, H_mod2) 
  loo_compare(D_mod, D_mod2) # in all cases there is no significant difference between fits, linear will be used as it is a simpler model
}

# Run Bayesian approach with resampling used to evaluate the relationships between each complexity variable and population density
PDlist <- lapply(1:iter, function(i) { print(i)
  MCbrm(dat = RawData, n = N, measure = 'PD') })
saveRDS(PDlist, file = paste0(FilePath, "PDlist.rds")) # checkpoint

# Extract R squared statistics
(PDR_R <- CI95(unlist(lapply(PDlist,'[[', 'R_R'))))
(PDD_R <- CI95(unlist(lapply(PDlist,'[[', 'D_R'))))
(PDH_R <- CI95(unlist(lapply(PDlist,'[[', 'H_R'))))

# Extract predicted relationships for plotting
RPD <- do.call(rbind, lapply(PDlist, '[[', 'R')) %>% arrange(R2); names(RPD) <- c('R','PD','Lower','Upper') 
DPD <- do.call(rbind, lapply(PDlist, '[[', 'D')) %>% arrange(D); names(DPD) <- c('D','PD','Lower','Upper') 
HPD <- do.call(rbind, lapply(PDlist, '[[', 'H')) %>% arrange(H2); names(HPD) <- c('H','PD','Lower','Upper')

#Extract probability of zero coefficients
coefR <- data.frame(Intercept = CI95(sapply(PDlist,'[[', 'R_hu')[1,]), b = CI95(sapply(PDlist,'[[', 'R_hu')[2,]))
coefH <- data.frame(Intercept = CI95(sapply(PDlist,'[[', 'H_hu')[1,]), b = CI95(sapply(PDlist,'[[', 'H_hu')[2,]))
coefD <- data.frame(Intercept = CI95(sapply(PDlist,'[[', 'D_hu')[1,]), b = CI95(sapply(PDlist,'[[', 'D_hu')[2,]))
# Compute mean and variance bounds in each modeled relationship (transform to probability of >0 )
Rvec <- seq(min(RPD$R),max(RPD$R), 0.01); Dvec <- seq(min(DPD$D),max(DPD$D), 0.01); Hvec <- seq(min(HPD$H),max(HPD$H), 0.01)
R_hu <- data.frame(R = Rvec, Mean = 1-(exp(coefR[1,1] + coefR[1,2]*Rvec)), Lower = 1-(exp(coefR[2,1] + coefR[2,2]*Rvec)), Upper = 1-(exp(coefR[3,1] + coefR[3,2]*Rvec)))
D_hu <- data.frame(D = Dvec, Mean = 1-(exp(coefD[1,1] + coefD[1,2]*Dvec)), Lower = 1-(exp(coefD[2,1] + coefD[2,2]*Dvec)), Upper = 1-(exp(coefD[3,1] + coefD[3,2]*Dvec)))
H_hu <- data.frame(H = Hvec, Mean = 1-(exp(coefH[1,1] + coefH[1,2]*Hvec)), Lower = 1-(exp(coefH[2,1] + coefH[2,2]*Hvec)), Upper = 1-(exp(coefH[3,1] + coefH[3,2]*Hvec)))

# Plot modeled relationships
# Rugosity
# determine smoothed relationships
RPD$LowerSmooth <- predict(loess(Lower ~ R, data = RPD), RPD$R)
RPD$UpperSmooth <- predict(loess(Upper ~ R, data = RPD), RPD$R)
RPD$MeanSmooth <- predict(loess(PD ~ R, data = RPD), RPD$R)
# determine axis scaling factor
dR <- max(RPD$UpperSmooth)/max(R_hu$Mean)
# generate plot
ggplot(data = RPD) +
  geom_ribbon(aes(x = R, ymin = LowerSmooth, ymax = UpperSmooth), fill = '#0000FF', alpha = 0.2) +
  geom_line(aes(x = R, y = MeanSmooth), col = '#000099', linewidth = 2) +
  geom_ribbon(aes(x = R, ymin = Lower*dR, ymax = Upper*dR), fill = '#F90707', alpha = 0.2, data = R_hu) +
  geom_line(aes(x = R, y = Mean*dR), col = '#A80000', linewidth = 1.1, data = R_hu) +
  xlab(NULL) +
  scale_x_continuous(expand = c(0,0), breaks = c(log10(1), log10(1e-05), log10(1e-10)), labels = function(i) 10^i) +
  scale_y_continuous(name = NULL,
                     labels = function(i) {format(i, digits = 4, scientific = FALSE)},
                     sec.axis = sec_axis(~./dR,
                                         name = NULL,
                                         labels = function(i) round(i, 3))) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     plot.margin = margin(t = 12, l = 5, b = 5, r = 5),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))

# Fractal Dimension
# determine smoothed relationships
DPD$LowerSmooth <- predict(loess(Lower ~ D, data = DPD), DPD$D)
DPD$UpperSmooth <- predict(loess(Upper ~ D, data = DPD), DPD$D)
DPD$MeanSmooth <- predict(loess(PD ~ D, data = DPD), DPD$D)
# determine axis scaling factor
dD <- max(DPD$UpperSmooth)/max(D_hu$Mean)
# generate plot
ggplot(data = DPD) +
  geom_ribbon(aes(x = D, ymin = Lower*dD, ymax = Upper*dD), fill = '#F90707', alpha = 0.2, data = D_hu) +
  geom_line(aes(x = D, y = Mean*dD), col = '#A80000', linewidth = 1.1, data = D_hu) +
  geom_ribbon(aes(x = D, ymin = LowerSmooth, ymax = UpperSmooth), fill = '#0000FF', alpha = 0.2) +
  geom_line(aes(x = D, y = MeanSmooth), col = '#000099', linewidth = 2) +
  xlab(NULL) +
  scale_y_continuous(name = NULL,
                     labels = function(i) {format(i, digits = 2, scientific = FALSE)},
                     sec.axis = sec_axis(~./dD,
                                         name = NULL,
                                         labels = function(i) round(i, 3))) +
  scale_x_continuous(expand = c(0,0), labels = function(i) {format(i, digits = 3, scientific = FALSE)},
                     breaks = c(2.1, 2.45, 2.8)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = margin(t = 12, l = 5, b = 5, r = 5),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))

# Height Range
# determine smoothed relationships
HPD$LowerSmooth <- predict(loess(Lower ~ H, data = HPD), HPD$H)
HPD$UpperSmooth <- predict(loess(Upper ~ H, data = HPD), HPD$H)
HPD$MeanSmooth <- predict(loess(PD ~ H, data = HPD), HPD$H)
# determine axis scaling factor
dH <- max(HPD$UpperSmooth)/max(H_hu$Mean)
# generate plot
ggplot(data = HPD) +
  geom_ribbon(aes(x = H, ymin = LowerSmooth, ymax = UpperSmooth), fill = '#0000FF', alpha = 0.2) +
  geom_line(aes(x = H, y = MeanSmooth), col = '#000099', linewidth = 2) +
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
                     plot.margin = margin(t = 12, l = 5, b = 5, r = 5),
                     axis.title = element_text(size = 15, colour = 'black'),
                     axis.text.x = element_text(size = 30, colour = "black"), 
                     axis.text.y = element_text(size = 30, colour = "black"),
                     panel.border = element_blank()) +
  theme(axis.line = element_line(color = 'black'))

# ----------------------------------------------------- End of Code ---------------------------------