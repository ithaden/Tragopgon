################################################################################
#
# R_combine_csv.R
#
# Script to combine Ingrid's .csv mass spec data into one
# really big .csv file.
#
# Author:  Brandon Jordon-Thaden
# Date:    13 Nov 2011
# Version: 13 Nov 2011 Draft
#          15 Nov 2011 Added a check on the number of data frame columns
#                      to make sure that the incoming csv file is ok.
#          29 Nov 2011 Minor edits
#
################################################################################

# Output file
fileOut <- "combined.csv"

# Number of datafram columns
numCol <- 19

# Number of lines at the top (other than the header) to strip out
numStrip <- 1

################################################################################

# What am I doing?
print("Running combine_csv.R script")
print(paste("  -> Working directory:", getwd(), "..."))

# List the .csv files in the directory
listFiles <- list.files(pattern = "\\.csv$",ignore.case = TRUE)

# Number of files in the directory
numFiles <- length(listFiles)

print(paste("  -> Found", numFiles, "files ..."))

# Loop over the .csv file in the directory and merge the data frames.
flagFirst <- 0

for (i in 1:numFiles) {
  print(paste("      ", listFiles[i], "..."))
  tableTemp <- read.csv(listFiles[i], header = TRUE,
                        skip = numStrip, check.names = TRUE)

  # Add a column to the data frame with the original csv file name
  tableTemp <- transform(tableTemp,file = listFiles[i])
  
  # In the orginal csv file the last column end is a comma, R interprets
  # that as an additional column and calls it 'X'. Drop it.
  tableTemp$X <- NULL 

  if ( ncol(tableTemp) != numCol)
    {print("    -> Wrong number of columns, not added")}
  else if (flagFirst == 0) {
    tableBig <- tableTemp
    flagFirst <- 1 }
  else tableBig <- rbind(tableBig,tableTemp)               
}

# Close things up
print(paste("  -> Found", nrow(tableBig), "entries ..."))
write.csv(tableBig,file = fileOut,quote=FALSE)
print(paste("  -> Wrote", fileOut, " ..."))
print("Done")
