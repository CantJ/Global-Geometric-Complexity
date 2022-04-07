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

###################################################
# STEP 1: Specify download URL, file destination and file requirements
###################################################

# Specify URL location of files for download
URL <- "https://e4ftl01.cr.usgs.gov/MEASURES/NASADEM_NC.001/2000.02.11/"
# This URL connects through the USGS Earth resources observation and Science Center LP DAAC Data Pool

# There are many files listed in this directory - this script is only interested in the .nc files.
all_files <- getURL(URL,verbose=TRUE,ftp.use.epsv=TRUE, dirlistonly = TRUE) # identifies all available files
nc_files <- getHTMLLinks(all_files, xpQuery = "//a/@href['.nc'=substring(.,string-length(.) - 2)]") # subsets those that are .nc files

# Specify destination file
destfile <- "/Volumes/Pocillopora/Geodiversity Data/Raw/NASADEM/"

# clear memory space
rm(all_files)

###################################################
# STEP 2: Define download function
###################################################

download_func <- function(ii){
  # progress print out
  print(ii)
  
  # define file name
  url_use <- paste0(URL, nc_files[ii])

  # define name for saved file
  destfile_use <- paste0(destfile, nc_files[ii])
    
  # download desired file
  curl_download(url = url_use, destfile = destfile_use, quiet = FALSE, mode = "wb")
}

###################################################
# STEP 3: Run download sequence
###################################################

lapply(1:length(nc_files), download_func)

# ------------------------------------------------------ End of Code -------------------------------------------------
