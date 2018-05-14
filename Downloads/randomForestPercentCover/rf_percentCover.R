#############################################################################
# This script reads training data from the CSV file created using the "percentCoverResample.R 
# script.  The script then uses the X and Y coordinates from the training data file to select
# the pixel values (predictor values) for each sample point in the input image. The predictor 
# values and the percent cover data from the training data file (response variable) are 
# combined and used as input to the random forests model. After the model is created percent 
# cover predictions are made on the input image to create an output image with percent cover 
# values ranging from 0 to 1. 
# 
# Set the variables below in the "SET VARIABLES HERE" section of the script. 
#
# This script was written by Ned Horning [horning@amnh.org]
# Support for writing and maintaining this script comes from The John D. and 
# Catherine T. MacArthur Foundation.
#
# This script is free software; you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software Foundation
# either version 2 of the License, or ( at your option ) any later version.                                   *
#
#############################################################################
#Load libraries
require(maptools)
require(sp)
require(randomForest)
require(raster)
require(rgdal)
#
#############################   SET VARIABLES HERE  ###################################
#
# The CSV file containing X, Y, and percent cover point data created by the percentCoverResample.R script.
pointData <- '/media/nedhorning/684EE5FF4EE5C642/AMNH/R_Project/TestData/outSamples2.csv'
# Name and path for the input satellite image 
inImage <-'/media/nedhorning/684EE5FF4EE5C642/AMNH/R_Project/TestData/alaska1_dem_slo_asp_v2_subset641.img'
# Name and path of the output GeoTiff image
outImage <- '/media/nedhorning/684EE5FF4EE5C642/AMNH/R_Project/TestData/out_percent_dual_ArcGIS_v4.tif'
# No data value for satellite image
nd <- 0
######################################################################################
#
# Start processing
print("Set variables and start processing")
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")

pointTable <- read.csv(pointData, header=TRUE)
xy <- SpatialPoints(pointTable[,1:2])
response <- as.numeric(pointTable[,3])

# Load the moderate resolution image
satImage <- stack(inImage)
for (b in 1:nlayers(satImage)) { NAvalue(satImage@layers[[b]]) <- nd }

# Get pixel DNs from the input image for each sample point
print("Getting the pixel values under each point")
trainvals <- cbind(response, extract(satImage, xy)) 

# Remove NA values from trainvals
trainvals_no_na <- na.omit(trainvals)

# Run Random Forest
print("Starting to calculate random forest object")
randfor <- randomForest(response ~. , data=trainvals_no_na)

# Start predictions
print("Starting predictions")
predict(satImage, randfor, filename=outImage, progress='text', format='GTiff', datatype='FLT4S', type='response', overwrite=TRUE)
#
# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("Processing time", format(timeDiff), "\n")
