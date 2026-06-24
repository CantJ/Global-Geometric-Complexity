# Script to centralize file and package details across an analysis exploring how global patterns in geometric landscape complexity
# mediate the configuration of natural and human communities.

# Primary Author: James Cant
# Date last modified: June 2026
# ------------------------------------------------------------------------------

# clear working directory
rm(list = ls())

# Modify R digit printing to prevent the rounding of very small values
options(digits = 22)

# Is this the first time this code is being implemented? I.e. is DEM tiling required
FirstRun <- FALSE

############################################
# STEP 1: Load required packages
############################################

# install INLA backend
inlaSetup <- FALSE # is installation required?
if(inlaSetup == TRUE) {
  install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}

# Define required packages
pkges <- c('tidyr', 'dplyr', 'reshape2', 'terra', 'sf', 'sp', 'raster', 'data.table', 'rcompanion', 'arrow', 'readr',
           'habtools', 'parallel', 'doParallel', 'mFD', 'forcats', 'ggh4x', 'xpectr', 'pbapply', 'gmodels', 'spsUtil', 
           'stringr', 'ggplot2', 'ggridges', 'hexbin', 'fishualize', 'viridis', 'brms', 'tidybayes', 'ggdist', 'gtools', 
           'ggeffects', 'bayestestR', 'BayesFactor', 'mgcv', 'qbrms', 'nnet', 'lmtest', 'marginaleffects')

# Install all packages required across the scripts (unless already installed on the system)
install.packages(setdiff(pkges, rownames(installed.packages())))
# Load installed packages
lapply(pkges, require, character.only = TRUE)

# Are cmdstanr and cmdstan installed and set up?
cmdstanSetup <- FALSE # is set up required?
if(cmdstanSetup == TRUE){
  install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
  library(cmdstanr)
  install_cmdstan()
  set_cmdstan_path(cmdstan_path())
} else {
  library(cmdstanr)
  set_cmdstan_path(cmdstan_path())
}

############################################
# STEP 2: Set up file paths
############################################

# The following file path is the main directory used by all scripts for key outputs
FilePath <- '/MAIN_FILE_DIRECTORY/'
# A few scripts require additional storage directories to keep interim temporary files seperate from main outputs or for storing 
# large quantities of raw data downloaded from external sources.
TilePath <- '/Temporary_File_Directory_1/' # directory for storing temporary DEM tiles during complexity indices computation
fileSave <- '/Temporary_File_Directory_2/' # directory for saving interim .csv files during complexity indices computation
EcoTypes <- '/Raw_File_Directory_1/' # initial storage for ecosystem typology maps downloaded from the following repository: 
                                     # Keith, D. et al. (2020). Indicative distribution maps for Ecosystem Functional Groups - Level 3 of IUCN Global Ecosystem Typology (Version 2.1.1) [Data set]. Zenodo. DOI: 10.5281/zenodo.3546513
ClimatePath <- '/Raw_File_Directory_2/' # Raw climate data (original data sources can be found in the manuscript associated with these analyses scripts)
tmpSave <- '/Temporary_File_Directory_3/' # directory for storing interim climate files during computation of climate variability
IslandPath <- '/Raw_File_Directory_3/' # island shape files downloaded from:
                                       # Sayre, R., 2023, Global Islands: U.S. Geological Survey data release, https://doi.org/10.5066/P91ZCSGM.
BirdData <- '/Raw_File_Directory_4/' # Bird diversity data (original data sources can be found in the manuscript associated with these analyses scripts)

# The data extraction and analysis scripts can now be implemented as desired.
# Each individual script requires user implementation to ensure required files exist (typically achieved by previous scripts in the sequence have been completed 
# or through downloading external data files as directed in the associated script).

# -------------------------------------------------------- End of Script -------