################################################################################
#
# R_call_csv.R
#
# 
#
# Input: fileIn: String with raw .csv file
#
# Output: 
#
# Author:  Brandon Jordon-Thaden
# Date:    15 Nov 2011
# Version: 15 Nov 2011 Draft
#
################################################################################
library(lattice)


# Data file
fileIn <- "combined.csv"

# Output file
#fileOut <- "combined_linked.csv"

# Filter parameters
cutTheta <- 10

thresholdTheta1 <- 0.0
thresholdTheta2 <- 22.5
thresholdTheta3 <- 45.0
thresholdTheta4 <- 67.5
thresholdTheta5 <- 90.0

################################################################################

# What am I doing?
print("Running R_call_csv.R script")

print(paste("  -> Working directory:", getwd(), "..."))

#Open files
tableIn <- read.csv(fileIn, header = TRUE)
print(paste("  -> Reading", fileIn, "..."))

vectorAssayID <- unique(as.character(tableIn$Assay.Id))
#numAssayID <- length(vectorAssayID)
numAssayID <- 100 # For testing purposes
print(paste("  -> Found", numAssayID, " Assay ID's ..."))

tableIn$radius <- sqrt(tableIn$Area^2 + tableIn$Area.2^2)
tableIn$theta <- atan(tableIn$Area/(tableIn$Area.2 + 0.00001)) / pi * 180

#tableIn$New.Call <- "Not Found"

#
# Loop over the assays
#


filePlot <- paste("Plot_plate_data", ".pdf", sep = "")
pdf(file = filePlot)
for (i in 1:numAssayID) {

  vectorNewcalls <- vector()
  # Subset the assay of interest
  tableSubset <- subset(tableIn, Assay.Id == vectorAssayID[i])

  # Determine the original calls for the assay of interest
  vectorCalls <- unique(as.character(tableSubset$Call))
  #numCalls <- nchar(vectorCalls)

  # Determine the number of parents and children
  vectorParents <- vector()
  vectorChildren <- vector()
  numBad <- 0;
  
  for (j in 1:length(vectorCalls)) {
    if (nchar(vectorCalls[j]) == 1) {
      vectorParents <- rbind(vectorParents, vectorCalls[j])} 
    else if (nchar(vectorCalls[j]) == 2) {
      vectorChildren <- rbind(vectorChildren, vectorCalls[j])}
    else numBad <- numBad + 1
  }

  numParents <- length(vectorParents)
  numChildren <- length(vectorChildren)

  print(paste("    -> Found",
              numParents,
              "parents and",
              numChildren,
              "children in Assay",
              vectorAssayID[i]))

  # Find the correlation between Peak 1 or Peak 2 based upon the original call


  flagfail <- FALSE
      stringoutput <- ""
  
  # Compare average of parent 1 calls theta to 0 and 90 deg.
  # Parent 1 goes to peak 1.
  if ( (numParents == 2) && (numChildren == 1) ){
    tableTemp1 <- subset(tableSubset, Call == vectorParents[1])
    tableTemp2 <- subset(tableSubset, Call == vectorParents[2])
  
    if      (mean(tableTemp1$theta) > (90 - cutTheta)) flagPeak1 = TRUE
    else if (mean(tableTemp2$theta) > (90 - cutTheta)) flagPeak1 = FALSE
    else if (mean(tableTemp1$theta) < cutTheta       ) flagPeak1 = FALSE
    else if (mean(tableTemp2$theta) < cutTheta       ) flagPeak1 = TRUE
    else flagfail <- TRUE
    
    stringoutput <- (paste("Found two parents", vectorParents[1],
              "and", vectorParents[2]))}

  else if( (numParents == 1) && (numChildren == 1) ) {
    tableTemp1 <- subset(tableSubset, Call == vectorParents[1])
    if      (mean(tableTemp1$theta) > (90 - cutTheta)) flagPeak1 = TRUE
    else if (mean(tableTemp1$theta) < cutTheta       ) flagPeak1 = FALSE
    else flagfail <- TRUE

    if (vectorParents[1] == substr(vectorChildren[1],1,1)) {
        vectorParents[2] <- substr(vectorChildren[1],2,2) }
    else vectorParents[2] <- substr(vectorChildren[1],1,1)
       
    stringoutput <- paste("Found one parent", vectorParents[1],
                "and infered the other from", vectorChildren[1])}
  
  else flagfail <- TRUE

  if (flagfail) stringoutput <- "Failed to find the parents."

  print(paste("    ->",stringoutput))

  # Make call
  if (!flagPeak1) vectorParents <- rev(vectorParents)
  
  for (j in 1:nrow(tableSubset)) {
    if ((!flagfail) && (tableSubset$theta[j] < (thresholdTheta1 + cutTheta))) {
      stringCall <- paste(vectorParents[2],
                          vectorParents[2],
                          vectorParents[2],
                          vectorParents[2], sep = "")}
    
    else if ((!flagfail) &&
             (tableSubset$theta[j] > (thresholdTheta2 - cutTheta)) &&
             (tableSubset$theta[j] < (thresholdTheta2 + cutTheta)) ) {
      stringCall <- paste(vectorParents[1],
                          vectorParents[2],
                          vectorParents[2],
                          vectorParents[2], sep = "")}

    else if ((!flagfail) &&
             (tableSubset$theta[j] > (thresholdTheta3 - cutTheta)) &&
             (tableSubset$theta[j] < (thresholdTheta3 + cutTheta)) ) {
      stringCall <- paste(vectorParents[1],
                          vectorParents[1],
                          vectorParents[2],
                          vectorParents[2], sep = "")}

    else if ((!flagfail) &&
             (tableSubset$theta[j] > (thresholdTheta4 - cutTheta)) &&
             (tableSubset$theta[j] < (thresholdTheta4 + cutTheta)) ) {
      stringCall <- paste(vectorParents[1],
                          vectorParents[1],
                          vectorParents[1],
                          vectorParents[2], sep = "")}
    
     else if ((!flagfail) &&
              (tableSubset$theta[j] > (thresholdTheta5 - cutTheta)) ) {
      stringCall <- paste(vectorParents[1],
                          vectorParents[1],
                          vectorParents[1],
                          vectorParents[1], sep = "")}

    else stringCall <- "NA"

    vectorNewcalls <- rbind(vectorNewcalls,stringCall)
    #tableSubset$New.Call[j] <- as.character(stringCall)
  }

  tableSubset <- cbind(tableSubset, New.Calls = vectorNewcalls)
  # Make plots

  #  print(histogram(tableSubset$Description,
#                  type = "count",
#                  main = paste(vectorAssayID[i],"Description"),
#                  xlab = "Description",
#                  ylab = "Number"))
  print( xyplot(Area ~ Area.2 | Description,
                data = tableSubset,
                main = paste(vectorAssayID[i],"\n",stringoutput),
                aspect = 1,
                xlab = "Area 2",
                ylab = "Area 1",
                layout = c(2,3)))
  
  print( histogram( ~ theta | Description,
                   data = tableSubset,
                   type = "count",
                   main = paste(vectorAssayID[i],"\n",stringoutput),
                   xlab = "theta (deg)",
                   ylab = "Number",
                   layout = c(2,3),
                   nint = 46
                   ))

    print( histogram( ~ Call | Description,
                   data = tableSubset,
                   type = "count",
                   main = paste(vectorAssayID[i],"\n",stringoutput),
                   xlab = "Original Call",
                   ylab = "Number",
                   layout = c(2,3)
                   ))

      print( histogram( ~ New.Calls | Description,
                   data = tableSubset,
                   type = "count",
                   main = paste(vectorAssayID[i],"\n",stringoutput),
                   xlab = "New Call",
                   ylab = "Number",
                   layout = c(2,3)
                   ))  
}

dev.off()

#
#for (i in 1:numIn) {
#  stringID <- c(strsplit(as.character(tableIn$Sample.Id[i]),"_"))
#
#  vecDNA[i] <- stringID[[1]] [1]
#
#  vecPopulation[i] <- stringID[[1]] [2]
#
#  vecSpecies[i] <- stringID[[1]] [length(stringID[[1]])]
#
#  strtemp <- stringID[[1]] [3]
#  for (j in 4:(length(stringID[[1]])-1)) {
#    strtemp <- paste(strtemp,stringID[[1]][j],sep="_")}
#  vecCoord[i] <- strtemp 
#}

#  if (numParents == 2) {
#    tableTemp1 <- subset(tableSubset, Call == vectorParents[1])
#    tableTemp2 <- subset(tableSubset, Call == vectorParents[2])
#
#    if (mean(tableTemp1$Area) > mean(tableTemp2$Area)) {
#      print(paste("    -> Identified parent",
#                  vectorParents[1],
#                  "with mass spec peak 1"))}
#    else {
#      print(paste("    -> Identified parent",
#                  vectorParents[2],
#                  "with mass spec peak 1"))
#      vectorParents <- rev(vectorParents)}}
#  else print("    -> I can't deal with this")


