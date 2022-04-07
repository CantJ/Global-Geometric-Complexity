# This script is for downloading multi-year (2015-2019) satellite data obtained by the Copernicus Land Cover Service.
# Further details on the product can be found at: https://land.copernicus.eu/global/index.html

# Primary Author: James Cant
# Date Last Modified: April 2022
# ------------------------------------------------------

# Clear the workspace
rm(list=ls(all=TRUE))

# source required packages
library(curl)

###################################################
# STEP 1: Specify download URL, file destination and file requirements
###################################################

destfile <- c("/Volumes/Pocillopora/Geodiversity Data/Raw/Land Cover/2015",
              "/Volumes/Pocillopora/Geodiversity Data/Raw/Land Cover/2016",
              "/Volumes/Pocillopora/Geodiversity Data/Raw/Land Cover/2017",
              "/Volumes/Pocillopora/Geodiversity Data/Raw/Land Cover/2018",
              "/Volumes/Pocillopora/Geodiversity Data/Raw/Land Cover/2019")

URL_date_code <- c(3939038,
                   3518026,
                   3518036,
                   3518038,
                   3939050)

URL_mode <- c("base", "conso",  "conso", "conso", "nrt")

Year <- c(2015,2016,2017,2018,2019)

file_select <- c("Bare-CoverFraction",
                 "BuiltUp-CoverFraction",
                 "Crops-CoverFraction",
                 "Grass-CoverFraction",
                 "MossLichen-CoverFraction",
                 "PermanentWater-CoverFraction",
                 "SeasonalWater-CoverFraction",
                 "Shrub-CoverFraction",
                 "Snow-CoverFraction",
                 "Tree-CoverFraction",
                 "Forest-Type")

###################################################
# STEP 2: Define download function
###################################################

download_func <- function(ii){
  # progress
  print(ii)
  
  # work through each land cover type
  for (xx in 1:length(file_select)) {
    # additional progress read out
    print(file_select[xx])
    
    # stitch together desired URL
    url_use <- paste0("https://zenodo.org/record/", URL_date_code[ii], "/files/PROBAV_LC100_global_v3.0.1_", Year[ii],
                      "-", URL_mode[ii], "_", file_select[xx], "-layer_EPSG-4326.tif?download=1")
    
    # define name for saved file
    destfile_use <- paste0(destfile[ii], "/", file_select[xx], ".tif")
    
    # download desired file
    curl_download(url = url_use, destfile = destfile_use, quiet = FALSE, mode = "wb")
  }
}

###################################################
# STEP 3: Run download sequence
###################################################

lapply(1:length(Year), download_func)

# ------------------------------------------------------ End of Code -------------------------------------------------