#############################################################################
# Updated 11 September 2018

# This script is used to calculate features from image segment using a zonal statistics approach. The 
# output is a Comma Separated Variable (CSV) file that can be used as input to the SegmentRF_Classification 
# R script to classify and a segmented image. The features that are calculated are:
#    Mean for each image band
#    Standard deviation for each image band
#    Mean, standard deviation, and coefficient of variation for each band after removing 40% of the pixels by clipping 20% of the pixels 
#         and 20% of the pixels with the highest DNs to reduce effects from outlier pixels included in the segment
#    Meadian for each band

# The script requires a multi-band image and a segment image that will typically be created from the multi-band 
# image using software such as eCognition or the open source Orfeo Toolbox (OTB) software. The output from the 
# script is a CSV file with the first column containing segment IDs for each segment in the raster segment image 
# and then several columns containing the feature information (e.g., mean and standard deviation for each image 
# band) for each segment. 

# Set the variables below in the "SET VARIABLES HERE" section of the script. 

# This script was written by Ned Horning [horning@amnh.org] 

# Support for writing and maintaining this script comes from Google.org.

#This script is free software; you can redistribute it and/or modify it under the terms of the GNU General 
# Public License as published by the Free Software Foundation either version 2 of the License, or ( at 
# your option ) any later version.            
#  
#
###########################################################################################
#Load libraries
require(raster)
require(data.table)
require(maptools)
#
print("Set variables and start processing")
#
#
#############################  SET VARIABLES HERE  ########################################
# Set variables here
# Set working directory
setwd("/media/ned/Data1/AMNH/WHRC_CarbonProject/GoogleProject/Tutorials/MappingForestCover/Data")
# Name and path for the image that will be used to calculate feature information image
satImage <- "RadarSubset_2.tif"

# Name and path for the segment polygon vector file. If using a raster image for segments use two single or double quotes
# with no space between them
#segVector <- "MeanShift15_10_50.shp"
segVector <- ''

# The label for the attribute field in the Shapefile that has segment IDs. This is ignored if segment data comes from an image
idAttributeLabel <- "DN"

# Name and path for the segment image. If using a vecotor file for segments use two single or double quotes
# with no space between them
segImage <- "MeanShift15_10_50Rasterize.tif"
#segImage <- ""

# Name and path of the output CSV file
outFeatures <- "segmentFeatures_SubRas3.csv"

# No-data value for the input satellite image
ndSat <- 0
""
# No-data value for the input segment image
ndSeg <- -1
###########################################################################################

## Start processing
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")

# Define functions for calculating features from segments
# Function to calculate mean after clipping upper and lower 20% of values
clippedMeanFun <- function(nums,...) {
  a=quantile(nums, probs=seq(0,1,0.2), na.rm = TRUE)
  mean(nums[which(nums >= a[2] & nums <= a[5])])
}

# Function to calculate standard deviation after clipping upper and lower 20% of values
clippedSDFun <- function(nums,...) {
  a=quantile(nums, probs=seq(0,1,0.2), na.rm = TRUE)
  sd(nums[which(nums >= a[2] & nums <= a[5])])
}

# Function to calculate coeficient of variation after clipping upper and lower 20% of values
clippedCVFun <- function(nums,...) {
  a=quantile(nums, probs=seq(0,1,0.2), na.rm = TRUE)
  as.double(cv(nums[which(nums >= a[2] & nums <= a[5])]))
}

# Function to calculate median since data.table has problems with meadian of integer date
medianFun <- function(nums,...) {
  median(as.double(nums), na.rm = TRUE)
}

# Function to create column labels to match the number of image bands
createColLabels <- function(numBands, baseName) {
  labels <- "segnumber"
  for (i in 1:numBands) {
    labels[i+1] <- paste(baseName, "Band", i, sep='')
  }
  labels
}

# Load images that will be used to create segment features
cat("Reading satellite image\n")
satelliteImage <- stack(satImage)
for (b in 1:nlayers(satelliteImage)) { NAvalue(satelliteImage@layers[[b]]) <- ndSat }

# Determine if per-segment or full image processing will be used
# If segVector has a path and filename process this block of code
if (segVector != '' & segImage == '') { 
  cat("Reading vector segments\n")
  segVectorLayerName <- strsplit(tail(unlist(strsplit(segVector, "/")), n=1), "\\.")[[1]] [1]
  vec <- readOGR(segVector, segVectorLayerName)
  allAtt <- slot(vec, "data")
  numRows <- nrow(allAtt)
  
  zonalMean <- matrix(nrow=0, ncol=nlayers(satelliteImage)+1)  
  zonalSD <- matrix(nrow=0, ncol=nlayers(satelliteImage)+1)
  zonalClippedMean <- matrix(nrow=0, ncol=nlayers(satelliteImage)+1)
  zonalClippedSD <- matrix(nrow=0, ncol=nlayers(satelliteImage)+1)
  zonalMedian <- matrix(nrow=0, ncol=nlayers(satelliteImage)+1)
  zonalCV <- matrix(nrow=0, ncol=nlayers(satelliteImage)+1)
  zonalClippedCV  <- matrix(nrow=0, ncol=nlayers(satelliteImage)+1)
  
  for (x in 1:numRows) {
  #for (x in 1:5) {
    segPoly<- vec[vec[[idAttributeLabel]]==allAtt[x,idAttributeLabel],]
    
    # Calculate zonal statistics
    zonalMean <- rbind(zonalMean, c(x, extract(satelliteImage, segPoly, fun=mean, na.rm = TRUE)))
    zonalSD <- rbind(zonalSD, c(x,extract(satelliteImage, segPoly, fun=sd, na.rm = TRUE)))
    zonalClippedMean <- rbind(zonalClippedMean, c(x,extract(satelliteImage, segPoly, fun=clippedMeanFun, na.rm = TRUE)))
    zonalClippedSD <- rbind(zonalClippedSD, c(x,extract(satelliteImage, segPoly, fun=clippedSDFun, na.rm = TRUE)))
    zonalMedian <- rbind(zonalMedian, c(x,extract(satelliteImage, segPoly, fun=medianFun, na.rm = TRUE)))
    zonalCV <- rbind(zonalSD, c(x,extract(satelliteImage, segPoly, fun=cv, na.rm = TRUE)))
    zonalClippedCV <- rbind(zonalClippedMean, c(x,extract(satelliteImage, segPoly, fun=clippedCVFun, na.rm = TRUE)))

    cat("Processing segment", x, "of",numRows, "\r")
  }
  zonalMean <- as.data.frame(zonalMean)
  zonalSD <- as.data.frame(zonalSD)
  zonalClippedMean <- as.data.frame(zonalClippedMean)
  zonalClippedSD <- as.data.frame(zonalClippedSD)
  zonalMedian <- as.data.frame(zonalMedian)
  zonalCV <- as.data.frame(zonalCV)
  zonalClippedCV <- as.data.frame(zonalClippedCV)
  cat("\nFinished processing polygons\n")
  
  # Relabel data table column headings
  meanLab <- createColLabels(nlayers(satelliteImage), "mean")
  colnames(zonalMean) <- meanLab
  sdLab <- createColLabels(nlayers(satelliteImage), "sd")
  colnames(zonalSD) <- sdLab
  clipMeanLab <- createColLabels(nlayers(satelliteImage), "clpMean")
  colnames(zonalClippedMean) <- clipMeanLab
  clipSDLab <- createColLabels(nlayers(satelliteImage), "clpSD")
  colnames(zonalClippedSD) <- clipSDLab
  medianLab <- createColLabels(nlayers(satelliteImage), "median")
  colnames(zonalMedian) <- medianLab
  cvLab <- createColLabels(nlayers(satelliteImage), "cv")
  colnames(zonalCV) <- cvLab
  clipCVLab <- createColLabels(nlayers(satelliteImage), "clpCV")
  colnames(zonalClippedCV) <- clipCVLab

# If segImage has an image path and filename process this block of code
} else if (segVector == '' & segImage != '') {
  cat("Reading raster segments\n")
  segmentImage <- raster(segImage)
  NAvalue(segmentImage) <- ndSeg
  
  # Prepare segment values for zonal statistics calculations
  cat("Calculating zonal statistics\n")
  vals <- getValues(satelliteImage)
  zones <- round(getValues(segmentImage), digits = 0)
  rasterDT <- data.table(vals, z=zones)
  setkey(rasterDT, z)
  
  # Calculate zonal statistics
  zonalMean <- rasterDT[, lapply(.SD, "mean", na.rm = TRUE), by=z]
  zonalSD <- rasterDT[, lapply(.SD, "sd", na.rm = TRUE), by=z]
  zonalClippedMean <- rasterDT[,lapply(.SD, clippedMeanFun), by =z]
  zonalClippedSD <- rasterDT[,lapply(.SD, clippedSDFun), by =z]
  zonalMedian <- rasterDT[, lapply(.SD, medianFun), by=z]
  zonalCV <- rasterDT[, lapply(.SD, "cv", na.rm = TRUE), by=z]
  zonalClippedCV <- rasterDT[,lapply(.SD, clippedCVFun), by =z]
  
  # Relabel data table column headings
  meanLab <- createColLabels(nlayers(satelliteImage), "mean")
  setnames(zonalMean, 1:(nlayers(satelliteImage)+1), meanLab)
  sdLab <- createColLabels(nlayers(satelliteImage), "sd")
  setnames(zonalSD, 1:(nlayers(satelliteImage)+1), sdLab)
  clipMeanLab <- createColLabels(nlayers(satelliteImage), "clpMean")
  setnames(zonalClippedMean, 1:(nlayers(satelliteImage)+1), clipMeanLab)
  clipSDLab <- createColLabels(nlayers(satelliteImage), "clpSD")
  setnames(zonalClippedSD, 1:(nlayers(satelliteImage)+1), clipSDLab)
  medianLab <- createColLabels(nlayers(satelliteImage), "median")
  setnames(zonalMedian, 1:(nlayers(satelliteImage)+1), medianLab)
  cvLab <- createColLabels(nlayers(satelliteImage), "cv")
  setnames(zonalCV, 1:(nlayers(satelliteImage)+1), cvLab)
  clipCVLab <- createColLabels(nlayers(satelliteImage), "clpCV")
  setnames(zonalClippedCV, 1:(nlayers(satelliteImage)+1), clipCVLab)
 
# If the segImage and or segVector parameters are not set properly display a message and stop processig
} else {
  cat("\n\n*****************************************************************************")
  cat("\nOne and only one of the variables segImage or segVector can have a file path and name\n")
  cat("and the other variable must be two quotes (single or double) with no space between them.\n\n")
  stop("Change the values for segImage or segVector and then restart the script\n\n", call.=FALSE)
}

# Merge results into a single table of training variables
meanSDMerge <- merge(zonalMean, zonalSD)
clippedmeanSDMerge <- merge(zonalClippedMean, zonalClippedSD)
trainValsTemp1 <- merge(meanSDMerge, clippedmeanSDMerge)
trainValsTemp2 <- merge(trainValsTemp1, zonalMedian)
cvMerge <- merge(zonalCV, zonalClippedCV)
trainVals <- merge(trainValsTemp2, cvMerge)

# Write the CSV file
write.csv(trainVals, file=outFeatures, row.names=FALSE)
# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("Processing time", format(timeDiff), "\n")
