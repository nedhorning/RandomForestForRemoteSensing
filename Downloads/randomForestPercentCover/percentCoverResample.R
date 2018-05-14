#############################################################################
# The script gets the pixel values in a high resolution classified image that correspond to 
# individual randomly selected moderate resolution pixels and then calculates the percent of 
# the classified image pixels that represent your cover type of interest. In other words, if 
# your high resolution image has a pixel size of 1m and your moderate resolution image has a 
# pixel size of 30m the sampling process would take a block of 900 of the 1m resolution pixels 
# that correspond to a single 30m pixel and calculate the percentage of the 1m pixels that 
# are forest. For example, if there were 600 forest pixels and 300 non-forest pixels the value 
# given for the output pixel would be 0.67 since 67% of the block of 1m pixels were forest. If 
# there are clouds or other no-data values in the high resolution image the following logic will 
# apply. If the total no-data values for a block of high resolution pixels is greater than or equal 
# to a user defined threshold (we will use 10% i.e., 90 or more pixels in our example above) then 
# it will not be included in the training data set since there is too much missing data to provide 
# a reliable cover percentage. If the cloud cover is less then 10% the no-data pixels are removed 
# from the total number of pixels when calculating the percent forest cover. 
#
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
require(raster)

#############################   SET VARIABLES HERE  ###################################
# Number of samples to be selected
numSamps <- 2000
# Name and path for the input classified image 
inClassImage <-'/media/684EE5FF4EE5C642/AMNH/R_Project/TestData/spot_641_v2_color_aea_reg.tif'
# Name and path for the input image that will be used for predictions 
inPredImage <-'/media/684EE5FF4EE5C642/AMNH/R_Project/TestData/alaska1_dem_slo_asp_v2_subset641.img'
# No data value for the prediction image
ndPred <- 0
# Class numbers that will be mapped using the following scheme:
# 0 = no data such as background, clouds and shadow
# 1 = class for which percent cover is being calculated 
# 2 = all other land cover classes
fromVals <- c(0,1, 2, 3, 4, 5, 6)
toVals <-   c(0,1, 1, 2, 1, 2, 0)
# Threshold for no-data processing
noDataPct <- 0.1
# Output file path and name
outFile <- '/media/684EE5FF4EE5C642/AMNH/R_Project/TestData/outSamples2.csv'
#############################################################################
#
# Start processing
print("Set variables and start processing")
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")

# Load the classified image
classImage <- raster(inClassImage)

# Load the first band of the image that will be used for predicitons. 
predImage <- raster(inPredImage, band=1)
NAvalue(predImage) <- ndPred

# Calculate the resolution of the two images
predImageRes <- res(predImage)[1]
classImageRes <- res(classImage)[1]
# Calculate the distance from the center of a raster cell to the edge - assuming a square pixel
halfCell <- predImageRes/2

# Check to make sure fromVals and toVals are the same length
if (length(fromVals) != length(toVals)) {
   stop("fromVals and toVals must have the same number of values \n\n", call.=FALSE)
}

# Select the class values that will be used to create the percent cover map
percentCoverValues <- fromVals[toVals == 1]

# Calculate the common extent to ensure sample cover both images
commonExt <- intersect(extent(predImage), extent(classImage))

# Select random samples from the prediciton image
cat("\nSelecting", numSamps, "samples from the prediction image\n")
sampleCells <- sampleRandom(crop(predImage,commonExt), size=numSamps, xy=TRUE, ext=commonExt, na.rm=TRUE)
lenSampCells <- nrow(sampleCells)
# Create the matrix that will hold the output data and label the columns
outputMatrix <- matrix(nrow=lenSampCells, ncol=3)
colnames(outputMatrix) <- c("X", "Y", "pct_cover")

fromTo <- data.frame(from=fromVals, to=as.integer(toVals==1))

# Get the pixels values from the classified image that fit inside the selected course-resolution pixel
cat("Calculating percent cover values for the output matrix\n\n")
for (i in 1:lenSampCells) {
   # Get the x and y coordinate from the center of the predition image pixel indicated by the sample point
   centerCoords <- sampleCells[i,1:2]
   # Insert x and y coordinate into the output matrix
   outputMatrix[i,1:2] <- c(centerCoords)
   # Calculate the extent of the selected prediction image pixel by adding/subtracting half the resolution of the pixel
   ext <- extent(centerCoords[1] - halfCell, centerCoords[1] + halfCell, centerCoords[2] - halfCell, centerCoords[2] + halfCell)
   # Get the cell numbers of all the classified image pixels that fall inside the extent of the selected prediction image pixel
   classCellNumbers <- cellsFromExtent(classImage, ext)
   # Get the class number of all of the selected classified image pixels
   classCellValues <- extract(classImage, classCellNumbers)
   # Convert classCellValues no-data pixels (0 in toVals) to NA
   classCellValues[classCellValues %in% fromVals[toVals==0]] <- NA
   # If the number of no-data pixels from the classified image is less than 'noDataPct' then calculate the percent cover
   # otherwise output NA
   if (sum(is.na(classCellValues))/length(classCellValues) < noDataPct) {
      outputMatrix[i,3] <- sum(!is.na(match(classCellValues, fromVals[toVals==1])))/length(classCellNumbers)
   } else {
      outputMatrix[i,3] <- NA
   }
   if (i %% 25 == 0)    cat("Processing sample", i, "of", lenSampCells, "\r")
}

# Write out the non-NA values in the output matrix to a CSV file
write.csv(outputMatrix[which(!is.na(outputMatrix[,3])),], outFile, row.names=FALSE)
#
# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("Processing time", format(timeDiff), "\n")
