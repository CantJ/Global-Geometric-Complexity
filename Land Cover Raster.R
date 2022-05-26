# This script is for processing global landcover data into a single raster file.
# These files will then be used for determining a measure of habitat complexity across space.
# the data used in this script has been downloaded via the USGS EarthData portal from the MODIS/Terra+Aqua Land Cover survey

# Date last modified: May 2022
# Primary Author: James Cant
# ------------------------------------------------------------------------------------------------------------

#Clear the workspace
rm(list=ls(all=TRUE))

# load required packages
library(raster)
library(terra)
library(gdalUtilities)
library(sp)
library(sf)
library(gdalUtils)
library(rgdal)

# Define file directories
MODIS_dir <- "/Volumes/Pocillopora/Geodiversity_Data/Raw/Land Cover/MODIS/"

# Define temporary file directories
TmpDir <- '/Users/james/Documents/R/R_TempFiles/' # (to be cleared at the end of each session)

# Save directory
fileDest <- '/Volumes/Pocillopora/Geodiversity_Data/Processed/'

#################################################
# STEP 1: LOAD AND REFORMAT DATA
#################################################

# Re-set working directory
setwd(MODIS_dir)

# Determine the files available and which years each file corresponds with
HDFfiles <- list.files(pattern = '.hdf')
fileYear <- substr(HDFfiles, 10, 13)

# Load global raster template (consistent with the template used across the other geodiversity variables)
worldRaster <- raster('/Volumes/Pocillopora/Geodiversity_Data/worldRaster.grd')

# Each .hdf file contains a series of subdatasets
gdalinfo(HDFfiles[1])


test <- h5file(HDFfiles[1])


# work through each file to translate the available data into a raster file.
# Tells me what subdatasets are within my hdf4 MODIS files and makes them into a list

sds <- get_subdatasets(HDFfiles[1])
sds

[1] "HDF4_EOS:EOS_GRID:MOD17A3H.A2000001.h21v09.006.2015141183401.hdf:MOD_Grid_MOD17A3H:Npp_500m"   
[2] "HDF4_EOS:EOS_GRID:MOD17A3H.A2000001.h21v09.006.2015141183401.hdf:MOD_Grid_MOD17A3H:Npp_QC_500m"

# I'm only interested in the first subdataset and I can use gdal_translate to convert it to a .tif

gdal_translate(src_dataset = HDFfiles[1], dst_dataset = paste0(fileDest,"NPP2000.tif"), sds = TRUE, of = 'GTiff')

# Load and plot the new .tif

rast <- raster("NPP2000.tif")
plot(rast)

# If you have lots of files then you can make a loop to do all this for you

files <- dir(pattern = ".hdf")
files

[1] "MOD17A3H.A2000001.h21v09.006.2015141183401.hdf" "MOD17A3H.A2001001.h21v09.006.2015148124025.hdf"
[3] "MOD17A3H.A2002001.h21v09.006.2015153182349.hdf" "MOD17A3H.A2003001.h21v09.006.2015166203852.hdf"
[5] "MOD17A3H.A2004001.h21v09.006.2015099031743.hdf" "MOD17A3H.A2005001.h21v09.006.2015113012334.hdf"
[7] "MOD17A3H.A2006001.h21v09.006.2015125163852.hdf" "MOD17A3H.A2007001.h21v09.006.2015169164508.hdf"
[9] "MOD17A3H.A2008001.h21v09.006.2015186104744.hdf" "MOD17A3H.A2009001.h21v09.006.2015198113503.hdf"
[11] "MOD17A3H.A2010001.h21v09.006.2015216071137.hdf" "MOD17A3H.A2011001.h21v09.006.2015230092603.hdf"
[13] "MOD17A3H.A2012001.h21v09.006.2015254070417.hdf" "MOD17A3H.A2013001.h21v09.006.2015272075433.hdf"
[15] "MOD17A3H.A2014001.h21v09.006.2015295062210.hdf"

filename <- substr(files,11,14)
filename <- paste0("NPP", filename, ".tif")
filename

[1] "NPP2000.tif" "NPP2001.tif" "NPP2002.tif" "NPP2003.tif" "NPP2004.tif" "NPP2005.tif" "NPP2006.tif" "NPP2007.tif" "NPP2008.tif"
[10] "NPP2009.tif" "NPP2010.tif" "NPP2011.tif" "NPP2012.tif" "NPP2013.tif" "NPP2014.tif"

i <- 1

for (i in 1:15){
  sds <- get_subdatasets(files[i])
  gdal_translate(sds[1], dst_dataset = filename[i])
}