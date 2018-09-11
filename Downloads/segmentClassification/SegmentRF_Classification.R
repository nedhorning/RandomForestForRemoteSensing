##########################################################################################
# Updated 11 September 2018
#
# This script is used to classify a segmented image using a random forests classifier. Before running this 
# script it is necessary to generate a comma separated Values (CSV) file with the segment number in the 
# first column and the predictor variables for each segment in the remaining columns. The R script 
# “CreateSegemntFeatures” can be used to create the CSV file. 
#
# In addition to the CSV file you will also need to have a raster image with the segments. Each segment in 
# the image should have a unique pixel value. You will also need a vector file that contains the locations 
# of the training (response variable) data such as a class number. Point, line and polygon vector files can 
# be used although at this time it is not possible to use multiple files or files with a mix of points, lines, 
# and polygons. The more pixels covered by the vector features the slow the scrip will run so in general points are 
# fastest and polygons are slowest. You will need to specify the column name or number of the column containing 
# the response variable (i.e., a unique number for each class that is to be mapped). It is important that the 
# projection of the image and vector data are the same or no matches will be found.
#
# The output classified image is written out to a user defined file (name and location is defined by the 
# "outImage variable) using the GeoTIFF format. A variable importance plot is displayed to provide information 
# about the influence of each variable. 
#
# There is also an option to assess the quality of the training data. The metric for thist 
# is the “margin”. The margin of a training segment is the proportion of votes for the correct 
# class minus maximum proportion of votes for the other classes for that segment. Positive margin 
# values represent correct classification, and vice versa. The margin data are written to a 
# point ESRI Shapefile so they can be overlaid on the image and segment polygons to assess which 
# segments need to be removed and relabeled in the training data and it can help determine which 
# classes needs additional training segments. If this output is not needed you can enter two 
# double or single-quotes (“” or '') for the variable outPointsFile.
# 
#Set the variables below in the "SET VARIABLES HERE" section of the script. 

#This script was written by Ned Horning [horning@amnh.org] and Wayne Walker [wwalker@whrc.org]

#Support for writing and maintaining this script comes from the Gordon and Betty Moore Foundation, 
# Google.org, and the David and Lucile Packard Foundation.

#This script is free software; you can redistribute it and/or modify it under the terms of the GNU General 
# Public License as published by the Free Software Foundation either version 2 of the License, or ( at 
# your option ) any later version.            
# 
##########################################################################################
#Load libraries
require(randomForest)
require(raster)
require(rgdal)
require(sp)
#############################  SET VARIABLES HERE  #######################################
# Set working directory
setwd("/media/ned/Data1/AMNH/WHRC_CarbonProject/GoogleProject/Tutorials/MappingForestCover/Data")

#Name and location of segment feature data CSV file output from the CreateSegmentFeatures_v1.R script
segCsv <- "segmentFeatures_SubRas3.csv"

# Name and location of the segment raster image 
segImage <- "MeanShift15_10_50Rasterize.tif"

# Segment raster nodata value.
nd <- 0

# Name and location of the classified image
outImage <- 'classImage.tif'

# Name and location of the output Shapefile point file that will be created. If this output 
# is not needed you can enter two double or single-quotes (“” or '')
# Note that if this file exists the write will fail with the message "Creation of output file failed" Enter 
outMarginFile <- 'margin.shp'

# Data set name for the vector file containing training data. This is often a file name or directory. 
# This and "layer" are defined by the ORG drivers. Look at http://www.gdal.org/ogr/ogr_formats.html for more info
trainingDsn <- 'clippedTrainingData.shp'

# Enter EITHER the name (case sensitive and in quotes) or the column number of the 
# field containing class (forest or non-forest) number
classNum <- "class_int"

# Output CSV file with class mapping information. If this output is not needed you can enter two double or single-quotes (“” or '')
outClassMapCSV <- 'classMapping.csv'

# Output classification without applying threshold (enter TRUE or FALSE)
classImage <- TRUE

# Output probability image (enter TRUE or FALSE)
probImage <- TRUE

# Output classification and set pixel under threshold to 0 (enter TRUE or FALSE)
threshImage <- TRUE

# Enter threshold probability in percent (values must be between 0 and 100) only used if threshImage=TRUE
probThreshold <- 75
###########################################################################################
## Start processing
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")

# Read the vector file 
cat("Reading the vector file\n")
trainingLayer <- strsplit(tail(unlist(strsplit(trainingDsn, "/")), n=1), "\\.")[[1]] [1]
vec <- readOGR(trainingDsn, trainingLayer)
pts <- slot(vec, "data")

# Load the segment raster 
segImg <- raster(segImage)

# Extract segment IDs under the point, line or polygon features
cat("Extracting segment IDs under the vector features\n")
exSegNums <- extract(segImg, vec, cellnumbers=TRUE)
if (is.matrix(exSegNums)) {
  exSegNums <- as.list(data.frame(t(exSegNums)))
}

# Remove NULL values from the list
exSegNums <- exSegNums[!sapply(exSegNums, is.null)]

# Select unique segment IDs under each vector features and associate the response 
# variable ("classNum") to the segment ID
trainSegs <- matrix(ncol=3, nrow=0)
for (i in 1:length(exSegNums)) {
  lineResponse <- pts[i,classNum]
  if (is.matrix(exSegNums[[i]])) { 
    segNums <- exSegNums[[i]][which(duplicated(exSegNums[[i]][,2]) == FALSE),]
  } else {
    segNums <- exSegNums[[i]]
  }
  
  if (is.matrix(segNums)) {
    trainSegs <- rbind(trainSegs, cbind(lineResponse, segNums))
  }
  else {
    trainSegs <- rbind(trainSegs, cbind(lineResponse, segNums[1], segNums[2]))
  }  
}

# Remove row names and add column names 
rownames(trainSegs) <- NULL
colnames(trainSegs) <- c("response", "cellNums", "segNums")

# Remove NA values from the list of unique segment IDs
trainSegs_no_na <- as.data.frame(na.omit(trainSegs))

# Read the CSV file with the segment feature information
segAtr <- read.csv(segCsv, header=TRUE)

# Remove NAs from the feature table
segAtr <- na.omit(segAtr)

#Create training data set by matching unique training segment IDs with segment feature information 
train <- segAtr[match(trainSegs_no_na$segNums, segAtr$segnumber),]

# Remove NAs from the training data frame
train_no_na <- as.data.frame(na.omit(train))

# Create response variable data frame
response_no_na <- trainSegs_no_na[match(train_no_na$segnumber, trainSegs_no_na$segNums), c(1,3)]

# Run Random Forest classification algorithm
cat("Starting to calculate random forest object \n)")
randfor <- randomForest(as.factor(response_no_na$response) ~. , data=train_no_na[,-1], proximity=TRUE)

# Write the output forest/nonforest raster map.
# Reload the raster package.
bs <- blockSize(segImg)

# Calculate how many bands the output image should have
numBands <- classImage + probImage + threshImage

# Create the output raster and begin writing to it.
img.out <- brick(segImg, values=FALSE, nl=numBands)
img.out <- writeStart(img.out, outImage, overwrite=TRUE, datatype='INT1U')

# Prediction
predValues <- predict(randfor, segAtr, type='response')
predValuesDF <- data.frame(segAtr$segnumber, predValues)

if (probImage || threshImage) { 
  predProbs <- predict(randfor, segAtr, type='prob')
  maxProb <- round(apply(predProbs, 1, max) * 100)
  maxProbDF <- data.frame(segAtr$segnumber, maxProb)
}

# Loop over blocks of the output raster from eCognition and write the new classified value.
# This looping method will allow for the input of larger rasters without memory problems.
for (i in 1:bs$n) {
  cat("processing block", i, "of", bs$n, "\r")
  img <- getValues(segImg, row=bs$row[i], nrows=bs$nrows[i])
  outMatrix <- matrix(nrow=length(img), ncol=0)
  # Set the no data value to NA so it doesn't get converted to a predicted value
  is.na(img) <- img == nd
  if (classImage) {
    # Convert the segment ID to the predicted (numeric) class so that a nodata value can be set.
    outMatrix <- cbind(outMatrix, predValuesDF$predValues[match(img, predValuesDF$segAtr.segnumber)])
  }
  if (probImage) { 
    outMatrix <- cbind(outMatrix, maxProbDF$maxProb[match(img, maxProbDF$segAtr.segnumber)])
  }
  if (threshImage) {
    threshValues <- as.numeric(predValuesDF$predValues[match(img, predValuesDF$segAtr.segnumber)])
    threshValues[which(maxProbDF$maxProb[match(img, maxProbDF$segAtr.segnumber)] <= probThreshold)] <- 0
    outMatrix <- cbind(outMatrix,threshValues)
  }
  writeValues(img.out, outMatrix, bs$row[i])
}

# Finish saving and close the connection to the image.
img.out <- writeStop(img.out)

# Output class mapping CSV file if a filename was provided for outClassMapCSV
if (outClassMapCSV != "") {
  write.csv(predValuesDF, file=outClassMapCSV, row.names=FALSE)
}

# View variable importance plot.
varImpPlot(randfor)

# Print error rate and confusion matrix for this classification
confMatrix <- randfor$confusion
cat("\n#################################################################################\n")
cat("OOB error rate estimate\n", 1 - (sum(diag(confMatrix)) / sum(confMatrix[,1:ncol(confMatrix)-1])), "%\n\n", sep="")
cat("Confusion matrix\n")
print(randfor$confusion)
cat("\n")

if (outMarginFile != "") {
  # Calculate margin (proportion of votes for correct class minus maximum proportion of votes for other classes)
  marginData <- margin(randfor)
  trainingAccuracy <- cbind(marginData[order(marginData)], trainSegs_no_na[order(marginData),])
  
  # Add column names to attributes table
  colnames(trainingAccuracy) <- c("margin", "classNum", "cellNum", "segID")
  
  # Calculate X and Y coordinates for training data points
  xyCoords <- matrix(ncol=2, nrow=0)
  for (z in 1:nrow(trainingAccuracy)) {
    xyCoords <- rbind(xyCoords, xyFromCell(segImg, trainingAccuracy[z,3]))
  }
  
  # Create and write point Shapefile with margin information to help improve training data
  pointVector <- SpatialPointsDataFrame(xyCoords, trainingAccuracy[, c(1,2,4)], coords.nrs = numeric(0), proj4string = segImg@crs)
  writeOGR(pointVector, outMarginFile, "layer", driver="ESRI Shapefile", check_exists=TRUE)
}

# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("Processing time", format(timeDiff), "\n")