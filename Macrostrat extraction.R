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

# There is a lot of strat classifications so the download below needs to be broken down into sections.
index_start <- seq(1, num_strat, by = 500)
index_end <- c(seq(500, num_strat, by = 500), num_strat)

# Now define the permanent sections of URL for downloading data for each bedrock class - the Macrostrat API allows these different files to be called simultaneously. 
# However the file size is too large to call in a single URL so it is nessecary to divide the calls.
URL_main <- "https://macrostrat.org/api/v2/geologic_units/map?strat_name_id="
URL_format <- "&format=geojson_bare"
  
###################################################
# STEP 2: Run download
###################################################

# The first download section will be run manually - with all later downloads added on top.
# define full URL
URL_full <- paste0(URL_main, paste0(strat_types[index_start[1]:index_end[1],]$strat_name_id, collapse = ','), URL_format)

# obtain shapefile from URL link
sf_map <- geojson_sf(URL_full)

# a little place holder
jar <- 2

# The remainder of the download will be done as a loop
for (ii in jar:length(index_start)) {
  # for dealing with errors 
  jar <- ii+1
  
  # progress read out
  print(ii)
  
  # define full URL
  URL_full <- paste0(URL_main, paste0(strat_types[index_start[ii]:index_end[ii],]$strat_name_id, collapse = ','), URL_format)
  
  # Obtain associated shapefile
  sf_temp <- suppressWarnings(geojson_sf(URL_full))
  
  # and add to existing data
  sf_map <- rbind(sf_map, sf_temp)
}

###################################################
# STEP 3: Save the extracted file
###################################################

# Specify destination
fileDest <- "/Volumes/Pocillopora/Geodiversity Data/Raw/Macrostrat/"

# write to file
st_write(sf_map, paste0(fileDest, "Macrostrat.shp"))

# ------------------------------------------------------ End of Code -------------------------------------------------
