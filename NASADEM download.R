# This script is for downloading high resolution digital elevation models from the NASADEM program 
# Further details on the product can be found at: https://earthdata.nasa.gov/esds/competitive-programs/measures/nasadem

# Primary Author: James Cant
# Date Last Modified: April 2022
# ------------------------------------------------------

# Clear the workspace
rm(list=ls(all=TRUE))

# source required packages
library(XML)
library(RCurl)
library(curl)
library(R.utils)
library(httr)

# Store username and password details registered with LP DAAC system
myusername <- "james.cant91@gmail.com"
mypassword <- "Sharky@74"

###################################################
# STEP 1: Specify download URL, file destination and file requirements
###################################################

# Specify URL location of files for download
URL <- "https://e4ftl01.cr.usgs.gov/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n00e006.zip/n00e006.hgt"
# This URL connects through the USGS Earth resources observation and Science Center LP DAAC Data Pool (https://lpdaac.usgs.gov/about/)
# There are many files listed in this directory - this script is only interested in the .nc files.
all_files <- GET(URL, config = list(authenticate(myusername, mypassword)), content_type(".zip"))
#convert response into HTML text string
all_files <- content(all_files, as = "text")
# extract names of the desired .nc files
HGT_files <- getHTMLLinks(all_files, xpQuery = "//a/@href['.zip'=substring(.,string-length(.) - 2)]") # subsets those that are .nc files

# Specify destination file
destfile <- "/Volumes/Pocillopora/Geodiversity Data/Raw/NASADEM/"

# clear memory space
rm(all_files)
curl_download("https://e4ftl01.cr.usgs.gov/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n00e006.zip", destfile = '/Volumes/Pocillopora/Geodiversity Data/TempFile/NASADEM_HGT_n00e006.zip/n00e006.hgt')
###################################################
# STEP 2: Define download function
###################################################

download_func <- function(ii){
  # progress print out
  print(ii)
  
  # define file name
  url_use <- paste0(URL, nc_files[ii])
    
  # download desired file
  downloadFile(url = URL, path = destfile, username = myusername, password = mypassword)
}

###################################################
# STEP 3: Run download sequence
###################################################

lapply(1:length(nc_files), download_func)

# ------------------------------------------------------ End of Code -------------------------------------------------
