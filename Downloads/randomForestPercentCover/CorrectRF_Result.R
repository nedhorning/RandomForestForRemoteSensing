#############################################################################
# When used for regression the random forests algorithm will overestimate low values 
# and underestimate high values. This script adjusts this effect by calculating 
# regression coefficients using the response variable from the training data and the
# corresponding predictor variables from the image output using X/Y coordinates provided
# by the input file. A new adjusted image image is output after applying gain and offset 
# (slope and intercept) values to the original predicted image values. The user can 
# select the type of regression to apply and set the minimum and maximum values for the 
# output image. This is useful when processing percent cover images when the valid range 
# is between 0 and 1. 
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
require(mblm)

#############################   SET VARIABLES HERE  ###################################
#
# Enter a 1 if the response variable training data are in a .dbf file or 2 if the data are in a CSV format
# with a header
fileType <- 2
# The CSV or DBF file containing X, Y, and response variable (i.e., biomass, % cover...) data. 
pointData <- '/home/nedhorning/AMNH/WHRC_CarbonProject/R_Project/TestData/outSamples2.csv'
# Name and path for the input predicted image 
inImage <-'/home/nedhorning/AMNH/WHRC_CarbonProject/R_Project/TestData/out_percent_dual_ArcGIS_v4.tif'
# Name and path of the output adjusted image 
outImage <-'/home/nedhorning/AMNH/WHRC_CarbonProject/R_Project/TestData/out_percent_dual_ArcGIS_v4_adj3.tif'
# No-data value for the input image
nd <- -1
# Enter EITHER the name (case sensitive and in quotes) or the column number of the 
# field containing X coordinates
x_coord <- "X"
# Enter EITHER the name (case sensitive and in quotes) or the column number of the 
# field containing Y coordinates
y_coord <- "Y"
# Enter EITHER the name (case sensitive and in quotes) or the column number of the 
# field containing the response variable values
responseVar <- "pct_cover"
# Minimum valid output value
minValue <- 0
# Maximum valid output value
maxValue <- 1
# Number of points to be randomly sampled. Enter -1 to use all sample points. It 
# may be necessary to use a subset of the sample points to avoid memory problems.
numSamps <- -1
#
# Type of regression to be applied- 1 = Theil-Sen Siegel repeated medians, 2 = linear, 
#                                   3 = linear with intercept of 0
regType = 2
#
# Display regression graphs and compare different models (TRUE or FALSE)?
# Solid line shows regression for Theil-Sen Siegel repeated medians
# Dashed line shows regression for linear regression
# Dotted line shows regression for linear regression with intercept of 0
# Display a graph?
dispGraphs = TRUE
#############################################################################
# Start processing
cat("Set variables and start processing\n")
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")

# Read the data in the dbf or CSV file and create a table of XY values and another 
# with the response variable values
if (fileType==1) {
   xy <- read.dbf(pointData)[,c(x_coord, y_coord)]
   responseVar <- as.numeric(read.dbf(pointData)[,responseVar])
} else if (fileType==2) {
   pointTable <- read.table(pointData, header=TRUE, sep=",")
   xy <- SpatialPoints(pointTable[,c(x_coord, y_coord)])
   responseVar <- as.numeric(pointTable[,responseVar])
}

# Sample xy and responce if numSamps is not -1
#
if (numSamps != -1) {
   if (numSamps > length(xy[,1])) {
      cat("\n\n*****************************************************************************")
      cat("\nThe variable numSamps is greater than the total number of points available to be randomly sampled\n")
      cat("Please change numSamps to be smaller than", length(xy[,1]), "which is the total number of sample points.\n\n")
      stop("Change the value for numSamps and then restart the script\n\n", call.=FALSE)
   }
   sampIdx = sample(seq(1, length(xy[,1])), size=numSamps)
   xy = xy[sampIdx,]
   responseVar = responseVar[sampIdx]
}

# Load the input predected image then flag all no-data values (nd) so they are not processed
predImage <- raster(inImage)
NAvalue(predImage) <- nd 

# Get pixel responseVar values from the image under each sample point and create a table with 
# observed and predicted responseVar values
cat("Getting the pixel values under each point\n")
samples <- cbind(responseVar, extract(predImage, xy)) 

# Remove NA values from trainvals table created above
samples <- na.omit(samples)

# Calculate blocksize for image writing at the end of the script
bs <- blockSize(predImage)

# Get output data type from the input image
dataType <- dataType(predImage)

samplesX <- samples[,2] # Predicted values from the input image
samplesY <- samples[,1] # Observed values from training data
   
cat("Plot graphs\n")
if (regType==1) {
   fit <- mblm(samplesY ~ samplesX, repeated=TRUE)
   slope <- fit$coefficients[2]
   intercept <- fit$coefficients[1]
   if (dispGraphs) {
      plot(samplesX, samplesY, xlab="Predicted", ylab="Observed")
      abline(fit, lty=1)
      abline(lm(samplesY ~ samplesX), lty=2)
      abline(lm(samplesY ~ samplesX+0), lty=3)
   }
} else if (regType==2) {
   fit <- lm(samplesY ~ samplesX)
   slope <- fit$coefficients[2]
   intercept <- fit$coefficients[1]
   if (dispGraphs) {
      plot(samplesX, samplesY, xlab="Predicted", ylab="Observed")
      abline(fit, lty=2)
      abline(lm(samplesY ~ samplesX+0), lty=3)
      abline(mblm(samplesY ~ samplesX), lty=1)
   }
} else if (regType==3) {
   fit <- lm(samplesY ~ samplesX+0)
   slope <- fit$coefficients[1]
   intercept <- 0
   if (dispGraphs) {
      plot(samplesX, samplesY, xlab="Predicted", ylab="Observed")
      abline(fit, lty=3)
      abline(lm(samplesY ~ samplesX), lty=2)
      abline(mblm(samplesY ~ samplesX), lty=1)
   }
} else {
   cat("\n   The variable regType must be an integer between 1 and 3 \n\nEnter Q to quit\n\n")
  browser()
}

print(fit)   

cat("Writing output image\n")
#Create output image and start writing to it.
img.out <- raster(predImage)
img.out <- writeStart(img.out, outImage, overwrite=TRUE, datatype=dataType)

#Loop over blocks of the image and write the adjusted values

for (i in 1:bs$n) {
   cat("Processing block", i, "of", bs$n, "\r")
   img <- getValues(predImage, row=bs$row[i], nrows=bs$nrows[i])
   img.pred <- img * slope + intercept
   # Set the no data value to the default value for the output image
   img.pred[is.na(img.pred) == TRUE] <- nd
   img.pred[img.pred < minValue] <- minValue
   img.pred[img.pred > maxValue] <- maxValue
   writeValues(img.out, img.pred, bs$row[i])
}
cat("\n")

#Finish saving and close connection to image.
img.out <- writeStop(img.out)
#
# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("Processing time", format(timeDiff), "\n")
