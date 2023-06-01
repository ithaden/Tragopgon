library('tcltk')

# Data out file
fileOut <- "combined_called.csv.gz"

#Call
strHigh <- "H"
strLow <- "L"

# Call Limits
vecDiploidcuts <- c(22.5, 67.5)
vecTetraploidcuts <- c(11.25, 33.75, 56.25, 78.75)

# Chose input files
directoryIn <- tk_choose.dir(caption = "Choose working directory")

setwd(directoryIn)
fileIn <- tk_choose.files(default = "combined_clustered.csv.gz", caption = "Choose input file", multi = FALSE)

#Open data file
print(paste("  -> Reading", fileIn, "..."))
tableIn <- read.csv(fileIn, header = TRUE)
print("     -> ... Done")

numSamples <- nrow(tableIn)
vecCall <- vector()
#numSamples <- 10000

for (i in 1: numSamples) {
  if (!(i%%1000)) print(paste(" -> working on",i,"of",numSamples,"..."))
  
  intCluster <- as.numeric(tableIn[i, ]$cluster)
  numClustertheta <- as.numeric(tableIn[i, ]$clustertheta)
  numClusterdtheta <- as.numeric(tableIn[i, ]$clusterdtheta)
  strPloidy <- tableIn[i, ]$Ploidy
  #  strHigh <- tableOut[i, ]$HighCall
  #  strLow <- tableOut[i, ]$LowCall
  
  #  print(paste(intCluster, numClustertheta, (tableIn[i, ]$clusterx)))
  
  if (intCluster >= 1) {
    if (strPloidy == "diploid"){
      if (numClustertheta <= vecDiploidcuts[1]) {
        new_call <- paste(strHigh, strHigh, sep = "")}
      else if ( (numClustertheta > vecDiploidcuts[1]) & (numClustertheta < vecDiploidcuts[2])) {
        new_call <- paste(strLow, strHigh, sep = "")}
      else {
        new_call <- paste(strLow, strLow, sep = "")}
    }
    else {
      if (numClustertheta < vecTetraploidcuts[1]) {
        new_call <- paste(strHigh, strHigh, strHigh, strHigh, sep = "")}
      else if ( (numClustertheta > vecTetraploidcuts[1]) & (numClustertheta < vecTetraploidcuts[2])) {
        new_call <- paste(strLow, strHigh, strHigh, strHigh, sep = "")}
      else if ( (numClustertheta > vecTetraploidcuts[2]) & (numClustertheta < vecTetraploidcuts[3])) {
        new_call <- paste(strLow, strLow, strHigh, strHigh, sep = "")}
      else if ( (numClustertheta > vecTetraploidcuts[3]) & (numClustertheta < vecTetraploidcuts[4])) {
        new_call <- paste(strLow, strLow, strLow, strHigh, sep = "")}
      else {
        new_call <- paste(strLow, strLow, strLow, strLow, sep = "")}
    }
  }
  else if (intCluster == 0) new_call <- "NoAlleles"
  else if (intCluster == -1) new_call <- "BadSpectrum"
  else if (intCluster == -2) new_call <- "BadTheta"
  else new_call <- NA
  
  vecCall[i] <- new_call
}

tableIn$New_Call <- vecCall

# Save data file
print(paste("  -> Saving", fileOut, "..."))
outStream <- gzfile(fileOut)
write.csv(tableIn, outStream)
print("     -> ... Done")
