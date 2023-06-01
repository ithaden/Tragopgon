################################################################################
#
# R_link_csv.R
#
# Given a .csv file and a lookup csv file. Link the two data file through
# the ID_chamber_pot.
#
# Input: fileIn: String with raw .csv file
#        fileLookup: String with lookup .csv file
#        fileOut: String with output file
#
# Output: Write a .csv file (fileOut) with linked information 
#
# Author:  Brandon Jordon-Thaden
# Date:    15 Nov 2011
# Version: 15 Nov 2011 Draft
#
################################################################################

# Data file
fileIn <- "combined.csv"

# Lookup table
fileLookup <- ""

# Output file
fileOut <- "combined_linked.csv"

# Magic numbers
vecDNA <- c("gD","cD")
vecDNAreplace <- c("gDNA","cDNA")

#List of locals
vecLocal <- c("A","B","C","D","E","F","G","H","J","K","L","M","N","O","P")
vecLocalreplace <- c("A","B","C","D","E","F",
                     "G","H","J","K","L","M","N","O","P")

#List (vector) of lines
vecLine <- c("i","ii","iii","iv","v","vi","vii","viii","ix","x",
             "I","II","III","IV","V","VI","VII","VIII","IX","X")

#List (vector) of species
vecSpecies <- c("m","r","prt","prr","db")
vecSpeciesreplace <- c("miscellus","mirus","pratensis","porrifolius","dubius")
################################################################################

# What am I doing?
print("Running R_link_csv.R script")

print(paste("  -> Working directory:", getwd(), "..."))

#Open files
print(paste("  -> Reading", fileIn, "..."))
tableIn <- read.csv(fileIn, header = TRUE)

#tableLookup <- read.csv(fileIn, header = TRUE)
#print(paste("  -> Reading", fileLookup, "..."))

#numLookup <- nrow(tableLookup)
numIn <- nrow(tableIn)
#numIn <- 500 # Used for testing

# Preallocate vectors
vecDNAID <- mat.or.vec(numIn,1)
vecSynNatID <- mat.or.vec(numIn,1)
vecLocalID <- mat.or.vec(numIn,1)
vecLineID <- mat.or.vec(numIn,1)
vecSpeciesID <- mat.or.vec(numIn,1)

# Loop over all of the elements in tableIn parsing the Sample.Id
for (i in 1:numIn) {
  print(i)
  # Break up Sample.Id at the _
  strSampleID <- c(strsplit(as.character(tableIn$Sample.Id[i]),"_"))

  # Number of elements in stringID
#  numSampleID <- length(stringID[[1]])
  numSampleID <- length(strSampleID[[1]])
  
  # Last string in the stringID vector
#  strLaststringID <- strID[[1]][numStrID] 
  
  # Identify the type of DNA
  indexDNA <- pmatch(strSampleID[[1]][1], vecDNA, nomatch = 0)
  
  if (indexDNA != 0) vecDNAID[i] <- vecDNAreplace[indexDNA]
  else { print(paste("    -> Unknown DNA type on line", i+1, "... Skipping"))
       vecDNAID[i] <- NA }

  # Indentify the natural species through local identifier
  indexLocal <- pmatch(strSampleID[[1]][2], vecLocal, nomatch = 0)
  
  if (indexLocal != 0) {
    vecLocalID[i] <- vecLocalreplace[indexLocal]
    vecSynNatID[i] <- "Natural"

    indexSpecies <- pmatch(strSampleID[[1]][numSampleID],
                           vecSpecies, nomatch = 0)
    if (indexSpecies !=0) vecSpeciesID[i] <- vecSpeciesreplace[indexSpecies]
    else {
      print(paste("    -> Unknown species on line", i+1, "... Skipping"))
      vecSpecies[i] <- NA
    }

    vecLineID[i] <- NA
  }

  # Indentify the synthetic species through local identifier
  indexLine <- pmatch(strSampleID[[1]][2], vecLine, nomatch = 0)
  if (indexLine != 0) {
    vecLineID[i] <- vecLine[indexLine]
    vecSynNatID[i] <- "Synthetic"
  }
  else vecLineID[i] <- NA
}

tableIn$DNA <- vecDNAID
tableIn$Species <- vecSpeciesID
tableIn$SynNat <- vecSynNatID
tableIn$Local <- vecLocalID
tableIn$Line <- vecLineID

write.csv(tableIn,file = fileOut,quote=FALSE)
print(paste("  -> Wrote", fileOut, " ..."))
print("Done")
  
#  if      (stringID[[1]][1] == "gD") vecDNA[i] <- "gDNA"
#  else if (stringID[[1]][1] == "cD") vecDNA[i] <- "cDNA"
#  else {
#    print(paste("    -> Can't indentify DNA type for",Sample.Id[i]))
#    vecDNA[i] <- NA
#  }

  #
#  if      (stringID[[1]][2] == "A"){
#    vecType[i] <- "Natural"
#    vecLocal[i] <- "A"
#    vecLineID[i] <- NA
#    vecSpecies[i] <- stringID[[1]][numStringID]
#  }
#  else if (stringID[[1]][2] == "B"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "B"
#    vecLineID[i] <- NA
#    vecSpecies[i] <- stringID[[1]][numStringID]
#  }
#  else if (stringID[[1]][2] == "C"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "C"
#    vecLineID[i] <- NA
#    vecSpecies[i] <- stringID[[1]][numStringID]
#  }
#  else if (stringID[[1]][2] == "D"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "D"
#    vecLineID[i] <- NA
#    vecSpecies[i] <- stringID[[1]][numStringID]
#  }
#  else if (stringID[[1]][2] == "E"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "E"
#    vecLineID[i] <- NA
#    vecSpecies[i] <- stringID[[1]][numStringID]
#  }
#  else if (stringID[[1]][2] == "F"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "F"
#    vecLineID[i] <- NA
#    vecSpecies[i] <- stringID[[1]][numStringID]
#  }
#  else if (stringID[[1]][2] == "G"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "G"
#    vecLineID[i] <- NA
#    vecSpecies[i] <- stringID[[1]][numStringID]
#  }
#  else if (stringID[[1]][2] == "H"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "H"
#    vecLineID[i] <- NA
#    vecSpecies[i] <- stringID[[1]][numStringID]
#  }
#  else if (stringID[[1]][2] == "J"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "J"
#    vecLineID[i] <- NA
#    vecSpecies[i] <- stringID[[1]][numStringID]
#  }
#  else if (stringID[[1]][2] == "K"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "K"
#    vecLineID[i] <- NA
#    vecSpecies[i] <- stringID[[1]][numStringID]
#  }
#  else if (stringID[[1]][2] == "L"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "L"
#  }
#  else if (stringID[[1]][2] == "M"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "N"
#  }
#  else if (stringID[[1]][2] == "O"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "O"
#  }
#  else if (stringID[[1]][2] == "P"){
#    vecType[i] <- "Natural"
#    vecPopulation[i] <- "P"
#  } 
#  else if (stringID[[1]][2] == "i"){
#    vecType[i] <- "Synthetic"
#    vecPopulation[i] <- NA
#    vecLineID[i] <- "i"
#  }
#  else if (stringID[[1]][2] == "ii"){
#    vecType[i] <- "Synthetic"
#    vecPopulation[i] <- NA
#    vecLineId[i] <- "ii"
#  }
#  else if (stringID[[1]][2] == "iii"){
#    vecType[i] <- "Synthetic"
#    vecPopulation[i] <- NA
#        vecLineId[i] <- "ii"
#  }
#  else if (stringID[[1]][2] == "i"){
#    vecType[i] <- "Synthetic"
#    vecPopulation[i] <- NA
#        vecLineId[i] <- "ii"
#  }
#  else if (stringID[[1]][2] == "ii"){
#    vecType[i] <- "Synthetic"
#    vecPopulation[i] <- NA
#        vecLineId[i] <- "ii"
#  }
#  else if (stringID[[1]][2] == "iii"){
#    vecType[i] <- "Synthetic"
#    vecPopulation[i] <- NA
#        vecLineId[i] <- "ii"
#  }
#
#
#  
#  else print("ARG")
#
#  
#  
#          
#  vecPopulation[i] <- stringID[[1]][2]
#
#  vecSpecies[i] <- stringID[[1]] [length(stringID[[1]])]
#
#                              
#  print(stringID[[1]][1])
#}
