#############################################################################
# The script reads an ESRI Shapefile (defined by the "shapefile" variable) with 
# training polygons and then either select all pixels or randomly select a 
# user-determined number of samples (defied by the "numsamps" variable) from 
# each land cover type. A multilayer image that contains spectral, other 
# continuous data or categorical data is also input (defined by the inImage 
# variable). For each randomly selected sample the data values for that pixel 
# are determined and these data are used to run the Random Forest model. 
#
# After building the model the multilayer image is read, and up to three output 
# layers (classImage, probImage, threshImage) can be selected for the output image. 
#      "classImage" classifies all of the pixels.
#
#     "probImage" outputs the class probability of the class that got the most votes 
#      (i.e., the class that was selected for the classImage layer). 
#
#     "threshImage" is the same as "classImage" except all pixels with a class probability 
#      of the class that got the most votes below the "probThreshold" parameter are set to 0. 
#      This is useful to identify pixels with inter-class confusion.
#
# The image is written out (name and location is defined by the "outImage variable) 
# using the GeoTIFF format. A variable importance plot is displayed to provide information 
# about the influence of each variable. An error rate estimate and confusion matrix are also 
# printed to provide information about classification accuracy.
#
# There is an option to assess the quality of the training data. The metric for thist 
# is the “margin”. The margin of a training point is the proportion of votes for the correct 
# class minus maximum proportion of votes for the other classes for that segment. Positive margin 
# values represent correct classification, and vice versa. The margin data are written to a 
# point ESRI Shapefile so they can be overlaid on the image and training polygons to assess which 
# points need to be removed and relabeled in the training data and it can help determine which 
# classes needs additional training segments. If this output is not needed you can enter two 
# double or single-quotes (“” or '') for the variable outPointsFile.
#
# There is also an option to output a feature space plot using two bands of your choice.
# If a feature space plot is not needed then enter "0" for the variables xBand and/or yBand.
# 
# Set the variables below in the "SET VARIABLES HERE" section of the script. 
#
# This script was written by Ned Horning [horning@amnh.org]
# Support for writing and maintaining this script comes from The John D. and 
# Catherine T. MacArthur Foundation and Google.org.
#
# This script is free software; you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software Foundation
# either version 2 of the Licenase, or ( at your option ) any later version.                               *
#
#############################################################################
#Load libraries
require(maptools)
require(sp)
require(randomForest)
require(raster)
#
cat("Set variables and start processing\n")
#
#############################   SET VARIABLES HERE  ###################################
# Name and path for the Shapefile (don't need the .shp extension)
shapefile <- '/home/nedhorning/AMNH/R_Project/TestData/spot_400_train/spot_400_train.shp'
# Approximate number of training samples to be randomly selected for each land cover class
# If numsamps is set to "0" then all pixels in all of the polygons will be used as training samples
numsamps <- 200
# Name of the attribute that holds the integer land cover type identifyer
attName <- 'type_id'
# No-data value for the input image
nd <- 0
# Name and path for the input satellite image 
inImage <-'/home/nedhorning/AMNH/R_Project/TestData/spot_subset_400.tif'
# Name and path of the output GeoTiff image
outImageName <- '/home/nedhorning/AMNH/R_Project/TestData/testProb_v6.tif'
# Name and location of the output Shapefile point file that will be created. If this output 
# is not needed you can enter two double or single-quotes (“” or '')
# Note that if this file exists the write will fail with the message "Creation of output file failed"  
outMarginFile <- '/home/nedhorning/AMNH/R_Project/TestData/junkPoint2.shp'
# Output classification layer without applying threshold (enter TRUE or FALSE)
classImage <- TRUE
# Output probability image layer (enter TRUE or FALSE)
probImage <- TRUE
# Output classification layer and set pixels with probability less than "probThreshold" to 0 (enter TRUE or FALSE)
threshImage <- TRUE
# Enter threshold probability in percent (values must be between 0 and 100) only used if threshImage=TRUE
probThreshold <- 75
# Layer number (band number) for the X and Y axis of the feature space plot. 
# If you do not want to calculate a feature plot enter 0 as the layer number
xBand <- 0
yBand <- 2
#######################################################################################
#
# Start processing
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")

# Read the Shapefile
vec <- readShapePoly(shapefile) 
#
# Load the image then flag all no-data values (nd) so they are not processed
satImage <- stack(inImage)
for (b in 1:nlayers(satImage)) { NAvalue(satImage@layers[[b]]) <- nd }
#
# Create vector of unique land cover attribute values
allAtt <- slot(vec, "data")
tabAtt <-table(allAtt[[attName]])
uniqueAtt <-as.numeric(names(tabAtt))

# Check if there are data in uniqueAtt
if (is.na(uniqueAtt[1])) {
  cat("\n*************************No attributes were found**************************** \n")
  stop("Check the attName variable in the variable settings\n", call.=FALSE)
}
# If all pixels in a polygon are to be used process this block
if (numsamps == 0) {
  # Create input data from a Shapefile using all training data 
  cat("Create training data using all pixels in training polygons\n")
  predictors <- data.frame()
  response <- numeric()
  for (x in 1:length(uniqueAtt)) {
    # Get the metadata for all polygons for a particular class (based on the uniqueAtt variable)
    class_data<- vec[vec[[attName]]==uniqueAtt[x],]
    # Extract and combine predictor and response variables for each polygon within a class
    for (i in 1:dim(class_data)[1]) {
      satValues <- extract(satImage, class_data[i,])
      satValues <- as.data.frame(do.call(rbind,satValues))
      attributeVector <- rep.int(uniqueAtt[x],nrow(satValues))
      
      predictors <- rbind(predictors, satValues)
      response <- c(response, attributeVector)
    }
  }
  trainvals <- cbind(response, predictors)
} else {
  # Create input data from a Shapefile by sampling training data 
  cat("Create training data by sampling", numsamps, "pixels for each class\n")
  for (x in 1:length(uniqueAtt)) {
    # Get the metadata for all polygons for a particular class (based on the uniqueAtt variable)
    class_data<- vec[vec[[attName]]==uniqueAtt[x],]
    # Get the area of each polygon for a particular class
    areas <- sapply(slot(class_data, "polygons"), slot, "area")
    # Calculate the number of samples for each polygon based on the area in proportion to total area for a class
    nsamps <- ceiling(numsamps*(areas/sum(areas)))
    # Use random sampling to select training points (proportial based on area) from each polygon for a given class 
    for (i in 1:dim(class_data)[1]) {
      xy_class <- spsample(class_data[i,], type="random", n=nsamps[i])
      # Add coordinates to create a list of random points for all polygons
      if (i == 1) cpts <- xy_class
      else cpts <- rbind(cpts, xy_class)
    }
    # The number of points might not match numsamps exactly.
    classpts <- cpts
    if (x == 1) {
      xy_allClasses<- classpts
    } else {
      xy_allClasses<- rbind(xy_allClasses, classpts)
    }
  } 
  # Get class number for each sample point for responce variable
  response <- over(xy_allClasses, vec)[[attName]]
  # Get pixel DNs from the image for each sample point
  trainvals <- cbind(response, extract(satImage, xy_allClasses))
}
# Test if feature space plot is needed
if (xBand != 0 & yBand != 0) {
  #Plot feature space and samples
  #
  continue <- "c"
  while (continue == "c") {
    plotImage <- stack(satImage[[xBand]], satImage[[yBand]])
    # Get pixel values from the image under each sample point and create a table with 
    # observed and predicted values
    cat("Getting pixel values to create feature space plot\n\n")
    featurePlotPoints <- sampleRegular(plotImage,100000 )
  
    # Remove NA values from trainvals table created above
    featurePlotPoints <- na.omit(featurePlotPoints)
  
    minBand1 <- min(featurePlotPoints[,1])
    maxBand1 <- max(featurePlotPoints[,1])
    minBand2 <- min(featurePlotPoints[,2])
    maxBand2 <- max(featurePlotPoints[,2])
    rangeBand1 <- maxBand1 - minBand1 + 1
    rangeBand2 <- maxBand2 - minBand2 + 1
  
    xAxisLabel <- paste("Layer", xBand, sep=" ")
    yAxisLabel <- paste("Layer", yBand, sep=" ")
  
    plot(featurePlotPoints[,1], featurePlotPoints[,2], col="lightgrey", xlab=xAxisLabel, ylab=yAxisLabel)
  
    uniqueValues <- unique(trainvals[,1])
    for (v in 1:length(uniqueValues)) {
      points(trainvals[which(trainvals[,1]==uniqueValues[v]), xBand+1], trainvals[which(trainvals[,1]==uniqueValues[v]), yBand+1], col=v, pch=20)
    }
  
    legend(minBand1, maxBand2, col=1:v, pch=20, title="Classes", legend=as.character(uniqueValues))
  
    continue <- readline(prompt="Type n to stop, c to change feature space bands or any other key to continue with randome forests model creation and prediciton: \n\n")
  
    if (substr(continue, 1,1) == "n") {
      stop("Processing stopped at users request \n\n", call.=FALSE)
    }
    if (substr(continue, 1,1) == "c") {
      xBand <- as.numeric(readline(prompt="Enter the band number for the x axis: \n"))
      yBand <- as.numeric(readline(prompt="Enter the band number for the y axis: \n"))
    }
  }
}

# Remove NA values 
trainvals <- na.omit(trainvals)

# Check to make sure Shapefile and input image are in the same projection
if (nrow(trainvals) == 0) {
  cat("\n*************************No training data found**************************** \n")
  stop("It is possible the projection of the Shapefile with training data and input image are different\nCheck projections and run again", call.=FALSE)
}

# Run Random Forest
cat("Calculating random forest object\n")
randfor <- randomForest(as.factor(response) ~., data=trainvals, importance=TRUE, na.action=na.omit)

# Start predictions
cat("Starting predictions\n")

# Calculate how many bands the output image should have
numBands <- classImage + probImage + threshImage

# Calculate the image block size for processing
bs <- blockSize(satImage)

# Create the output raster block
outImage <- brick(satImage, values=FALSE, nl=numBands)
outImage <- writeStart(outImage, filename=outImageName, progress='text', format='GTiff', datatype='INT1U', overwrite=TRUE)

# Loop though each of the image blocks to calculate the output layers selected in the variables section
for (i in 1:bs$n) {
  cat("processing block", i, "of", bs$n, "\r")
  imageBlock <-  getValuesBlock(satImage, row=bs$row[i], nrows=bs$nrows[i])
  predValues <- predict(randfor, imageBlock, type='response')
  classValues <- as.numeric(levels(predValues))[predValues]
  outMatrix <- matrix(nrow=nrow(imageBlock), ncol=0)
  if (classImage) {
    outMatrix <- cbind(outMatrix, classValues)
  }
  if (probImage || threshImage) { 
    predProbs <- as.data.frame(predict(randfor, imageBlock, type='prob'))
    maxProb <- round(apply(predProbs, 1, max) * 100)
    if (probImage) { 
      outMatrix <- cbind(outMatrix, maxProb)
    }
    if (threshImage) {
      threshValues <- classValues
      threshValues[which(maxProb <= probThreshold)] <- 0
      outMatrix <- cbind(outMatrix,threshValues)
    }
  }
  writeValues(outImage, outMatrix, bs$row[i])
}

# Stop writing and close the file
outImage <- writeStop(outImage)

# Plotting variable importance plot
varImpPlot(randfor)

# Print error rate and confusion matrix for this classification
confMatrix <- randfor$confusion
cat("#################################################################################\n")
cat("OOB error rate estimate\n", 1 - (sum(diag(confMatrix)) / sum(confMatrix[,1:ncol(confMatrix)-1])), "%\n\n", sep="")
cat("Confusion matrix\n")
print(randfor$confusion)
cat("\n")

if (outMarginFile != "") {
  # Calculate margin (proportion of votes for correct class minus maximum proportion of votes for other classes)
  marginData <- margin(randfor)
  trainingAccuracy <- cbind(marginData[order(marginData)], trainvals[order(marginData),1])
  
  # Add column names to attributes table
  colnames(trainingAccuracy) <- c("margin", "classNum")  
  # Calculate X and Y coordinates for training data points
  xyCoords <- xy_allClasses@coords
  xyCoords <- xyCoords[order(marginData),]
  
  # Create and write point Shapefile with margin information to help improve training data
  pointVector <- SpatialPointsDataFrame(xyCoords, as.data.frame(trainingAccuracy), coords.nrs = numeric(0), proj4string = satImage@crs)
  writeOGR(pointVector, outMarginFile, "layer", driver="ESRI Shapefile", check_exists=TRUE)
}

# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("\nProcessing time", format(timeDiff), "\n")