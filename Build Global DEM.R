# This script is for constructing a global-scale digital elevation model from satellite imagery products downloaded from GEBCO (https://www.gebco.net/).
# The specific map product selected was GEBCO's 2023 Global terrain model for ocean and land (including ice surface elevation).

# Date last modified: March 2024
# Primary Author: James Cant
# ------------------------------------------------------------------------------------------------------------

# Clear workspace
rm(list=ls(all=TRUE))

# Define directory pathways
fileLoad <- '/File_Directory_1/' # file directory containing downloaded GEBCO products.
fileSave <- '/File_Directory_2/'

# load necessary packages
library(terra)


#################################################
# STEP 1: Generate Global Map
#################################################

# Locate raw GEBCO products
mapList <- list.files(fileLoad, pattern = '.tif', full.names = TRUE)
# these files are a series of tiles outlining elevation patterns of the earths surface.
# They are currently projected in WGS84 format and to a resolution of 15 arc-seconds.

# Mosaic together the map tiles and open in raster format.
RastList <- lapply(mapList, rast)
GlobalRast <- terra::merge(sprc(RastList)) # this is a large raster and can take a period of time to compute.

# Save raster file
writeRaster(GlobalRast, paste0(fileSave, 'GlobalDEM_WGS84.tif'))

# The size of this raster file necessitates that reprojection is carried out using QGIS (due to the removal of rGDAL related packages from CRAN)
# Using QGIS this raster file is to be reprojected to a Mollweide projection (EPSG: 9001, ESR1: 54009) with a resolution of 187m using bilinear interpolation.
# The python code associated with this reprojection is as below:

# gdalwarp -s_srs EPSG:4326 -t_srs ESRI:54009 -tr 187.0 187.0 -r bilinear -of GTiff FILE_DIRECTORY_PATH_TO/GlobalDEM_WGS84.tif "FILE_DIRECTORY_PATH_FOR_SAVING/GlobalDEM_Mollweide_187m.tif"

# --------------------------------------------------------- End of code ---------------------------------------------- #
