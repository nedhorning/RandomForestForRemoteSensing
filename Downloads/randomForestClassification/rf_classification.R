#############################################################################
#
# The script reads a Vector File (for example an ESRI Shapefile or Geopackage) with 
# training polygons and then either selects all pixels contained in the polygons or 
# randomly selects a user-determined number of samples (defined using classNums and
# classSampNums) for each land cover type. A multilayer image that contains spectral, other 
# continuous data or categorical data is also input (defined by the inImage 
# variable). For each randomly selected sample the data values for that pixel 
# are determined and these data are used to run the Random Forest model. 
#
# After building the model the multilayer image is read, and up to three output 
# images (classImage, probImage, threshImage) can be selected. 
#     "classImage" classifies all of the pixels.
#
#     "probImage" outputs the class probability of the class that got the most votes 
#      (i.e., the class that was selected for the classImage layer). 
#
#     "threshImage" is the same as "classImage" except all pixels with a class probability 
#      of the class that got the most votes below the "probThreshold" parameter are set to 0. 
#      This is useful to identify pixels with inter-class confusion.
#
# The images are written out using the GeoTIFF format and the file name is created by appending 
# "_Class" to the input image file name and it is written to the same directory as the input 
# image. A variable importance plot is displayed to provide information 
# about the influence of each variable. An error rate estimate and confusion matrix are also 
# printed to provide information about classification accuracy.
#
# There is an option to assess the quality of the training data. The metric for this
# is the margin. The margin of a training point is the proportion of votes for the correct 
# class minus maximum proportion of votes for the other classes for that segment. Positive margin 
# values represent correct classification, and vice versa. The margin data are written to a 
# point vector file with the same type as the input file (Using Geopackages is strongly recommended). 
# This vector file can be overlaid on the image and training polygons to assess which 
# points need to be removed and relabeled in the training data and it can help determine which 
# classes needs additional training segments.
#
# There is also an option to output a feature space plot using two bands of your choice.
# If a feature space plot is not needed then enter "0" for the variables xBand and/or yBand.
# When a feature space plot is drawn it is possible to define a rectange on the plot to highlight
# pixels in the image that are not well represented in the trianing data. 
# 
# Set the variables below in the "SET VARIABLES HERE" section of the script. 
#
# This script was written by Ned Horning [horning@amnh.org] and modified by Jonas Viehweger
# Support for writing and maintaining this script comes from The John D. and 
# Catherine T. MacArthur Foundation and Google.org.
#
# This script is free software; you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software Foundation
# either version 2 of the License, or ( at your option ) any later version.
#
#############################################################################
#Load libraries
library(snow)
library(maptools)
library(sf)
library(randomForest)
library(raster)
library(rgdal)
library(lwgeom)

cat("Set variables and start processing\n")
#
#############################  SET VARIABLES HERE  ###################################
# Set working directory
setwd("C:/Users/USERNAME/Scriptlocation")
# To get reproducible results, set the seed to a constant number
set.seed(13579)
# Should output files be overwritten if they already exist (enter TRUE or FALSE)
overwrite = F
# Toggle for using multi-core processing
multiprocessing = F
cores = 4
# Name and path for the Shapefile (with extension)
shapefile <- "./Sample_Data/BALATON_FNF.gpkg"
# Class numbers that you want to select training sample from
classNums <- c(1,2,3)
# For each land cover class the approximate number of training samples to be randomly selected 
# If a value is "0" then all pixels in all of the polygons for that class will be used 
classSampNums <- c(100,100,100)
# Name of the attribute that holds the integer land cover type identifyer
attName <- "DN"
# No-data value for the input image
nd <- 0
# Name and path for the input satellite image
inImageName <- "./Sample_Data/BALATON_L8.tif"
# Output point file with margins of sampled points (enter TRUE or FALSE)
marginFile <- TRUE
# Output classification image (enter TRUE or FALSE)
classImage <- TRUE
# Output probability image layer (enter TRUE or FALSE)
probImage <- TRUE
# Output classification layer and set pixels with probability less than "probThreshold" to 0 (enter TRUE or FALSE)
threshImage <- TRUE
# Enter threshold probability in percent (values must be between 0 and 100) only used if threshImage=TRUE
probThreshold <- 75
# Layer number (band number) for the X and Y axis of the feature space plot. 
# If you do not want to calculate a feature plot enter 0 as the layer number
xBand <- 2
yBand <- 3
#######################################################################################

# Start processing
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")

# Initialize Cluster and set Raster options
if(multiprocessing) beginCluster(cores)
rasterOptions(progress = "text", timer = TRUE, overwrite = TRUE, datatype = "INT2U")

# Helper function to initialize all paths for output and check if the files exist
createOutPath <- function(appString){
  out <- tryCatch(
    {
      outFileBaseName <- tools::file_path_sans_ext(inImageName)
      outPath <- paste0(outFileBaseName, appString)
      
      if(file.exists(outPath) && overwrite){
        warning(paste("File",outPath,"already exists \n"))
      } else if(file.exists(outPath)){
        stop(paste("File",outPath,"already exists \n"))
      } else {
        return(outPath)
      }
    },
    warning = function(cond){
      message(cond)
      message("It will be overwritten.")
      return(outPath)
    },
    error = function(cond){
      message(cond)
      if(multiprocessing) endCluster()
      stop("Please move or rename any output files that already exist.
       If the output files should automatically be overwritten you can set the variable 'overwrite' to TRUE.", call. = F)
    }
  )
}

# Initializing necessary outputs
if(classImage) outClassImage <- createOutPath("_Class.tif")
if(probImage) outProbImage <- createOutPath("_Prob.tif")
if(threshImage) outThreshImage <- createOutPath("_Thresh.tif")
if(marginFile) outMarginFile <- createOutPath(paste0("_Margin.",tools::file_ext(shapefile)))

# Read the Shapefile
vec <- st_read(shapefile)

# Load the image then flag all no-data values(nd) so they are not processed
satImage <- brick(inImageName)
NAvalue(satImage) <- nd

# Create vector of unique land cover attribute values
uniqueAtt <-as.numeric(vec[[attName]])

# Check if length of classNums and classSampNums is equal
if (length(classNums) != length(classSampNums)) {
  cat("\n***************length of classNums and classSampNums not equal***************** \n")
  if(multiprocessing) endCluster()
  stop("Check the classNums and classSampNums variable\n", call.=FALSE)
}

# Check if classNums and classSampNums are numeric vectors
if (!is.numeric(classNums) && !is.numeric(classSampNums)) {
  cat("\n***************classNums and classSampNums can only contain numbers***************** \n")
  if(multiprocessing) endCluster()
  stop("Check the classNums and classSampNums variable\n", call.=FALSE)
}

# Check if all classNums exist in uniqueAtt
#### CHECK THIS FUNCTION TO SEE IF classNums ARE IN uniqueAtt  ################
if (sum(classNums %in% uniqueAtt) != length(classNums)) {
  cat("\n*******not all classes in classNums are defined in the vector file******* \n")
  if(multiprocessing) endCluster()
  stop("Check classNums and vector attribute table\n", call.=FALSE)
}

## Scalable and quick implementation of Sampling using sf instead of sp, omitting
## multiple for loops. Not a big difference in speed for Ground Truth with small amounts
## of polygons or sampled points, but significantly faster for a lot of polygons and a lot of sampled points.
## This implementation will not guarantee that every polygon in a class will get sampled,
## since the sampled points are generated over the overall extent of the polygons.
## The number of sampled points will be exactly the amount specified

trainvals <- data.frame()
xyCoords <- data.frame()

cat("Create training data to train model\n")
for(i in 1:length(classNums)){
  if (classSampNums[i] == 0) {
    cat("Create training data using all pixels in polygons with class",classNums[i],"\n")
    onlyAvailable <- vec[which(vec$class %in% classNums[i]),]
    
    if (marginFile){
      extracted <- within(extract(satImage, as(onlyAvailable, "Spatial"), df=TRUE, cellnumbers=TRUE), rm('ID'))
      coords <- xyFromCell(satImage, extracted$cell)
      colnames(coords) <- c("X", "Y")
      xyCoords <- rbind(xyCoords, coords)
      extracted <- within(extracted, rm('cell'))
    } else {
      extracted <- within(extract(satImage, as(onlyAvailable, "Spatial"), df=TRUE), rm('ID'))
    }
    
    response <- rep.int(classNums[i],nrow(extracted))
    extracted <- cbind(response,extracted)
    trainvals <- rbind(trainvals, na.omit(extracted))
  } else {
    cat("Create training data using",classSampNums[i],"pixels in polygons with class",classNums[i],"\n")
    onlyAvailable <- vec[which(vec[[attName]] %in% classNums[i]),]
    coords <- st_coordinates(st_sample(onlyAvailable, classSampNums[i]))
    response <- rep.int(classNums[i], classSampNums[i])
    extracted <- within(extract(satImage, coords, df=TRUE), rm('ID'))
    temp <- na.omit(cbind(response, extracted))
    
    if (marginFile) xyCoords <- rbind(xyCoords, coords)
    
    trainvals <- rbind(trainvals, temp)
  }
}

timeDiff <- Sys.time() - startTime
cat("\nProcessing time for training data", format(timeDiff), "\n\n")

# Test if feature space plot is needed
if (xBand != 0 & yBand != 0) {
  #Plot feature space and samples
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
    
    continue <- readline(prompt="Type n to stop, c to change feature space bands, s to define a rectangle to locate gaps in feature space, or any other key to continue with random forests model creation and prediciton: \n\n")
    
    if (substr(continue, 1,1) == "n") {
      if(multiprocessing) endCluster()
      stop("Processing stopped at users request \n\n", call.=FALSE)
    }
    if (substr(continue, 1,1) == "s") {
      cat("Click two points to define the area on the feature space plot that you want to highlight\n")
      coords <- locator(n=2)
      coords <- unlist(coords)
      xvals <- coords[1:2]
      yvals <- coords[3:4]
      
      # Print out the corner coordinates for the rectangle
      cat("min X =", min(xvals), "\n")
      cat("max X =", max(xvals), "\n")
      cat("min y =", min(yvals), "\n")
      cat("max y =", max(yvals), "\n")
      
      # Draw the rectangle on the feature space plot
      rectangle <- matrix(nrow=5, ncol=2)
      rectangle[1,] <- c(min(xvals), max(yvals))
      rectangle[2,] <- c(max(xvals), max(yvals))
      rectangle[3,] <- c(max(xvals), min(yvals))
      rectangle[4,] <- c(min(xvals), min(yvals))
      rectangle[5,] <- c(min(xvals), max(yvals))
      lines(rectangle[,1], rectangle[,2])
      
      # Get the bands used to calculate the feature space plot
      b1 <- raster(plotImage, layer=1)
      b2 <- raster(plotImage, layer=2)
      
      # Threshold satImage so all values selected in the rectangle on the feature space plot are set to 255
      satImage[(b1 > min(xvals)) & (b1 < max(xvals)) & (b2 > min(yvals)) & (b2 < max(yvals))] <- 255
      
      # Plot the thresholded image with selected pixels displayed as white pixels
      plotRGB(satImage, r=1,g=2,b=3,stretch='hist')
      cat("White pixels in the plotted image were selected in the rectangle drawn on the feature space plot")
      if(multiprocessing) endCluster()
      stop("Add new training data and re-run the script \n\n", call.=FALSE)
    }
    if (substr(continue, 1,1) == "c") {
      xBand <- as.numeric(readline(prompt="Enter the band number for the x axis: \n"))
      yBand <- as.numeric(readline(prompt="Enter the band number for the y axis: \n"))
    }
  }
}

# Check to make sure Shapefile and input image are in the same projection
if (nrow(trainvals) == 0) {
  cat("\n*************************No training data found**************************** \n")
  stop("It is possible the projection of the Shapefile with training data and input image are different\nCheck projections and run again", call.=FALSE)
}

# Run Random Forest
cat("Calculating random forest object\n")
randfor <- randomForest(as.factor(response) ~., data=trainvals, importance=TRUE, na.action=na.omit)

# Start predictions

classPred <- function(x){
  return(x %/% 1000)
}

probPred <- function(x){
  return(x %% 1000)
}

classAndProb <- function(x){
  require(randomForest)
  voteValues <- predict(randfor, x, type='vote', norm.votes=TRUE)
  #classValues <- apply(voteValues, 1, my.max)*1000
  classValues <- as.integer(colnames(voteValues)[max.col(voteValues,ties.method="first")])*1000
  return(round(apply(voteValues, 1, max) * 100)+classValues)
}

threshPred <- function(x,y){
  result <- x
  result[y <= probThreshold] <- NA
  return(result)
}

clusterPred <- function(x, fun=NULL) {
  calc(x, fun)
}

if(multiprocessing){

  cl <- getCluster()
  clusterExport(cl, list('randfor','probThreshold'))

  if (probImage || threshImage || classImage){
    cat("Starting predictions\n")
    pred_calc <- clusterR(satImage, clusterPred, args = list(fun=classAndProb), export = 'classAndProb')
  }
  
  if (classImage || threshImage) {
    cat("Starting class calculation\n")
    class_calc <- calc(pred_calc, fun=classPred, filename = outClassImage)
    # Multiprocessing Implementation:
    # Works, but is often slower
    #class_calc <- clusterR(satImage, clusterPred, args = list(fun=classPred), filename = outClassImage, export = 'classPred')
  }
  
  if (probImage || threshImage){
    cat("Starting probability calculation\n")
    prob_calc <- calc(pred_calc, fun=probPred, filename = outProbImage)
    # Multiprocessing Implementation:
    # Works, but is often slower
    # prob_calc <- clusterR(satImage, clusterPred, args = list(fun=probPred), filename = outProbImage, export = 'probPred')
  }
  
  if(threshImage){
    cat("Starting threshold calculation\n")
    thresh_calc <- overlay(class_calc, prob_calc, fun = threshPred, filename = outThreshImage)
  }
  
  endCluster()

} else {
  if (probImage || threshImage || classImage){
    cat("Starting predictions\n")
    pred_calc <- calc(satImage, fun=classAndProb)
  }
  
  if (classImage || threshImage) {
    cat("Starting class calculation\n")
    class_calc <- calc(pred_calc, fun=classPred, filename = outClassImage)
  }
    
  if (probImage || threshImage){
    cat("Starting probability calculation\n")
    prob_calc <- calc(pred_calc, fun=probPred, filename = outProbImage)
  }
  
  if(threshImage){
    cat("Starting threshold calculation\n")
    thresh_calc <- overlay(class_calc, prob_calc, fun = threshPred, filename = outThreshImage)
  }
}
 
# Print error rate and confusion matrix for this classification
confMatrix <- randfor$confusion
cat("\n#################################################################################\n")
cat("OOB error rate estimate\n", 1 - (sum(diag(confMatrix)) / sum(confMatrix[,1:ncol(confMatrix)-1])), "%\n\n", sep="")
cat("Confusion matrix\n")
print(randfor$confusion)
cat("\n")

if (marginFile) {
  # Calculate margin (proportion of votes for correct class minus maximum proportion of votes for other classes)
  marginData <- margin(randfor)
  trainingAccuracy <- cbind(marginData[order(marginData)], trainvals[order(marginData),1])

  # Add column names to attributes table
  colnames(trainingAccuracy) <- c("margin", "classNum")
  # Order X and Y coordinates
  xyCoords <- xyCoords[order(marginData),]

  # Create and write point Shapefile with margin information to help improve training data
  xyCoords <- cbind(xyCoords,trainingAccuracy)
  mpoints <- st_as_sf(xyCoords, crs = satImage@crs, coords = c("X", "Y"), 
                      agr = "constant", dim = "XY", remove = TRUE, na.fail = FALSE,
                      sf_column_name = NULL)
  
  st_write(mpoints, outMarginFile, update = FALSE, delete_layer = TRUE)
}

# Plotting variable importance plot
varImpPlot(randfor)

# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("\nProcessing time", format(timeDiff), "\n")
