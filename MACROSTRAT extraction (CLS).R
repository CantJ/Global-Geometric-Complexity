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
num_strat <- dim(strat_types)[1] # how many classifications is this?

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

# run download, there is a lot to be done here, trycatch ensures that any download errors don't end the process.
Strat_data <- sapply(1:num_strat,function(x) tryCatch(download_func(x), error=function(e) NULL), simplify = TRUE)

###################################################
# STEP 4: Save the extracted file/Open the extracted shape file if restarting from the check point
###################################################

# Specify destination
fileDest <- "D:/MACROSTRAT/"

# write to file
#save(Strat_data, file = paste0(fileDest, "Macrostrat.RData"))
# this works as a little check point as the file is so large!

# load file
load(paste0(fileDest, "Macrostrat.RData")) # this will automatically open the file as 'Strat_data'

# define file dimension
num_strat <- length(Strat_data)

###################################################
# STEP 5: Reshape the extracted file into a shapefile
###################################################

# how many rows are being dealt with across the list
num_entries <- sum(unlist(Map(nrow, Strat_data)))
# for strata types for which there is data 21 different variables have been extracted
num_cols <- 21
# extract variable names
col_names <- colnames(Strat_data[[19]])

# define an empty storage output
data_output <- as.data.frame(matrix(NA, nrow = num_entries, ncol = num_cols))
# rename column IDs
colnames(data_output) <- col_names

# define indexing count
index <- 1

# run the reshape
for (ii in 1:num_strat) {

  # print progress
  print(ii)
  
  # select next shape file within the list
  sfExtract <- Strat_data[[ii]]
  
  # determine dimensions
  sfSize <- dim(sfExtract)[1]
  
  # store extracted shapefile (if there is anything to store)
  if(sfSize == 0 || is.null(sfSize)){
    
    # if there is not data - inform user
    cat("There is no data here \n")
    
  } else {
    
    # define row index(s)
    index_seq <- seq(index, ((index-1)+sfSize), by = 1)
    
    # format data for inserting
    class(sfExtract) <- "data.frame"
    
    # insert the extracted file into the blank frame
    # being sure to use the correct indexing
    for(x in 1:dim(sfExtract)[1]){
      data_output[index_seq[x],] <- sfExtract[x,]
    }

    # finally advance indexing for next run through
    index <- index+sfSize
    
  }
 # end of extraction loop 
}

###################################################
# STEP 6: Save the reformatted shapefile
###################################################

# write to file
st_write(data_output, paste0(fileDest, "Macrostrat.shp"), driver = "ESRI Shapefile")

# ------------------------------------------------------ End of Code -------------------------------------------------