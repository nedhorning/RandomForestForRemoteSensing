##########################################################################################
#Updated 11 September 2018
#
# This script is used to edit the output classified image from the SegmentRF_Classification.R 
# script. When running the SegmentRF_Classification script you need to make sure you set the 
# variable "outClassMapCSV" to a file. This script uses a point Shapefile with the each
# point representing a segment that will have its class label changed. The attribute table for
# the point file must have an integer attribute ("newClassNum") to hold the new class value. 
# The output is a new classified image with the new class labels. 

#Set the variables below in the "SET VARIABLES HERE" section of the script. 

#This script was written by Ned Horning [horning@amnh.org] and Wayne Walker [wwalker@whrc.org]

#Support for writing this script came from Google.org.

#This script is free software; you can redistribute it and/or modify it under the terms of the GNU General 
# Public License as published by the Free Software Foundation either version 2 of the License, or ( at 
# your option ) any later version.            
# 
##########################################################################################
#Load libraries
require(raster)
require(rgdal)
#############################  SET VARIABLES HERE  #######################################
# Set working directory
setwd("/media/ned/Data1/AMNH/WHRC_CarbonProject/GoogleProject/Tutorials/MappingForestCover/Data")

# Name and location of the segment raster image 
segImage <- "MeanShift15_10_50Rasterize.tif"

# Segment raster nodata value.
nd <- 1

# Name and location of the edited classified map
outImage <- "classImageEdited.tif"

# Input class mapping CSV file
inClassMapCSV <- "classMapping.csv"

# Data set name for the vector point file containing locations and class assignments for segments to be modified. 
# This is often a file name or directory. This and "layer" are defined by the ORG drivers. 
# Look at http://www.gdal.org/ogr/ogr_formats.html for more info
editPointsDsn <- "editSegs.shp"

# Enter EITHER the name (case sensitive and in quotes) or the column number of the 
# field containing class (forest or non-forest) number
newClassNum <- "id"
###########################################################################################
## Start processing
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")

# Read the vector file 
cat("Reading the vector file\n")
editPointsLayer <- strsplit(tail(unlist(strsplit(editPointsDsn, "/")), n=1), "\\.")[[1]] [1]
vec <- readOGR(editPointsDsn, editPointsLayer)
newClass <- slot(vec, "data")

# Load the segment raster image
segImg <- raster(segImage)

#open the CSV file
classMap <- read.csv(inClassMapCSV)

# Extract segment IDs under the points to edit
cat("Extracting segment IDs under the points to modify\n")
changeSegIDs <- cbind(as.integer(as.character(newClass[,newClassNum])), extract(segImg, vec))

# Change class mapping based on point data
for (i in 1:nrow(changeSegIDs)) {
  classMap[which(classMap[,1]==changeSegIDs[i, 2]), 2] <- changeSegIDs[i,1]
}
 
# Write the output forest/nonforest raster map.
# Reload the raster package.
bs <- blockSize(segImg)

# Create the output raster and begin writing to it.
img.out <- raster(segImg)
img.out <- writeStart(img.out, outImage, overwrite=TRUE, datatype='INT1U')

# Loop over blocks of the output raster from eCognition cland write the new classified value.
# This looping method will allow for the input of larger rasters without memory problems.
for (i in 1:bs$n) {
  cat("processing block", i, "of", bs$n, "\r")
  # require(raster)
	img <- getValues(segImg, row=bs$row[i], nrows=bs$nrows[i])
	# Set the no data value to NA so it doesn't get converted to a predicted value
	is.na(img) <- img == nd
	# Convert the segment ID to the predicted (numeric) class so that a nodata value can be set.
  img.match <- as.numeric(classMap$pred[match(img, classMap$segAtr.segnumber)])
  # Set the no data value to the default value for the output image
  img.match[is.na(img.match) == TRUE] <- NAvalue(img.out)
  writeValues(img.out, img.match, bs$row[i])
}

# Finish saving and close the connection to the image.
img.out <- writeStop(img.out)
  
# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("Processing time", format(timeDiff), "\n")