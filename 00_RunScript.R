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

# Define required packages
pkges <- c('tidyr', 'dplyr', 'reshape2', 'terra', 'sf', 'sp', 'raster', 'data.table', 'rcompanion', 'arrow', 'readr',
           'habtools', 'parallel', 'doParallel', 'mFD', 'forcats', 'ggh4x', 'xpectr', 'pbapply', 'gmodels', 'spsUtil', 
           'stringr', 'ggplot2', 'ggridges', 'hexbin', 'fishualize', 'viridis', 'brms', 'tidybayes', 'ggdist', 'gtools', 
           'ggeffects', 'vegan', 'bayestestR', 'BayesFactor', 'mgcv')
qbrms 'INLA'

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

# The data extraction and analysis script can now be implemented as desired.
# Each individual script requires user implementation to ensure appropriate file paths are set up and 
# required files exist (achieved by previous scripts in the sequence have been completed).

# -------------------------------------------------------- End of Script -------