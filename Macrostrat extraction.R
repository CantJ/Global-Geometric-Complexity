# This script is for extracting and downloading lithographic (rock) data from the Macrostrat geological database (https://macrostrat.org/)

# Primary Author: James Cant
# Date last modified: April 2022.
# ------------------------------------------------------

# Clear the workspace
rm(list=ls(all=TRUE))

# source required packages
library(RCurl)
library(curl)
library(geojsonsf)
library(sf)
library(httr)

###################################################
# STEP 1: Specify download URL and file destination
###################################################

# Data from Macrostrat is downloaded through an Application Programming Interface (API).
# First determine the different bedrock classifications present within the database.
strat_types <- read.csv("https://macrostrat.org/api/v2/defs/strat_names?all&format=csv")
num_strat <- dim(strat_types)[1] # how many classification is this?

# Define the permanent sections of URL for downloading data for each bedrock class - the Macrostrat API allows these different files to be called simultaneously. 
# However the file size is too large to call in a single URL so it is nessecary to divide the calls.
URL_main <- "https://macrostrat.org/api/v2/geologic_units/map?strat_name_id="
URL_format <- "&format=geojson_bare"
  
###################################################
# STEP 2: Define download function
###################################################

download_func <- function(ii){
  # progress read out
  print(ii)
  
  # define full URL
  URL_full <- paste0(URL_main, strat_types$strat_name_id[ii], URL_format)
  
  # obtain shapefile from URL link
  sf_map <- suppressWarnings(geojson_sf(URL_full))
  
  # return output
  return(sf_map)
}

###################################################
# STEP 3: Run download 
###################################################

# run download
Strat_data <- sapply(1:100, download_func, simplify = TRUE)

###################################################
# STEP 4: Save the extracted file
###################################################

# Convert file into a single shapefile


# Specify destination
fileDest <- "/Volumes/Pocillopora/Geodiversity Data/Raw/Macrostrat/"

# write to file
st_write(sf_map, paste0(fileDest, "Macrostrat.shp"))

# ------------------------------------------------------ End of Code -------------------------------------------------
