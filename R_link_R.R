#///////////////////////////////////////////////////////////////////////////////
# Information #####
# R_link_R.R 
#
# Given a .csv file and a lookup csv file. Link the two data file through
# the ID_chamber_pot.
# 
# Input: fileIn: String with raw .csv file
#
# Output: Write a compressed data file (fileOut) with linked information
#
# TODO: 
# 
# Author: jordobr
# Date: Jan 8, 2012
# 
#Edit: 20 May 2013 Edit for usability
#      19 Feb 2014 First commit
#      21 Feb 2014 Change output file to .csv.gz
#      12 March 2014 Edit to fit Rstudio
#
#///////////////////////////////////////////////////////////////////////////////

##### Load libraries ###########################################################
library(tcltk)
library(foreach)
library(doParallel)

##### User input for file io ###################################################
# Suggested input file
stringIn <- 'combined.csv'

# Suggested lookup file for high and low calls
stringLookup <- 'TragopogonChamberData_fallspring_051211.csv'

# Suggested output file
stringOut <- 'combined_linked.csv.gz'

# Log file
fileSink <- 'log_link.txt'

#### User input parameters for all samples #####################################
# List of DNA types
tableDNA <- matrix( c(
				'gD', 'gDNA', 
				'cD', 'cDNA'
        ), nrow = 2, ncol = 2, byrow = TRUE)

# Lookup table for species and ploidy
tablePloidy <- matrix(c(
				'miscellus', 	 'tetraploid',
				'mirus', 			 'tetraploid',
				'pratensis', 	 'diploid',
				'porrifolius', 'diploid',
				'dubius', 		 'diploid'
        ), nrow = 5, ncol = 2, byrow = TRUE)

#### User input parameters for natural samples #################################
strNatural <- 'Natural'

#List of locals with replacements 
tableLocal <- matrix(c(
				'A', 'Spokane-5th&Walnut',
				'B', 'Spokane-5th-Sher',
				'C', 'Spokane-Canon',
				'D', 'Spangle',
				'E', 'Rosalia',
				'F', 'Tekoa',
				'G', 'Oakesdale2nd  Perl',
				'H', 'Garfield',
				'J', 'Palouse',
				'K', 'Pullman1-Finches',
				'L', 'Pullman2-Elm',
				'M', 'Pullman3-N.Grand',
				'N', 'Pullman-(long-lig)PostOffice',
				'O', 'Moscow',
				'P', 'Troy'
        ), nrow = 15, ncol = 2, byrow = TRUE)

# List of natural species abbreviations
tableSpecies <- matrix(c(
				'm', 	  'miscellus', 
				'r', 		'mirus', 	
				'prt',	'pratensis', 	
				'prr', 	'porrifolius',	
				'db',  	'dubius'
        ), nrow = 5, ncol = 2, byrow = TRUE)

#### User input parameters for synthetic samples ###############################
strSynthetic <- 'Synthetic'

# List (vector) of lines
tableLine <- matrix(c(
				'i', 		'mirus',
				'ii', 	'mirus',
				'iii', 	'mirus',
				'iv', 	'mirus',
				'v', 		'mirus',
				'vi', 	'mirus',
				'vii',	'mirus',
				'viii',	'mirus',
				'ix', 	'mirus',
				'x',		'mirus',
				'I', 		'miscellus',
				'II', 	'miscellus',
				'III', 	'miscellus',
				'IV', 	'miscellus',
				'V', 		'miscellus',
				'VI', 	'miscellus',
				'VII', 	'miscellus',
				'VIII', 'miscellus',
				'IX',		'miscellus',
				'X',		'miscellus'
        ), nrow = 20, ncol = 2, byrow = TRUE)

# List of sythentic locales
vecLineLocal <- c('S3','S1','S2')

#### User input parameters for control samples #################################
# List of controls and specials and negative controls in the order of
# SynNat Line Locale extra Species PlateID and Polodiy
tableParentControl <- matrix(c(
				'Parent Control', 'PM1', 'Troy', 			NA, 'dubius and porrifolius',	'2887_7a_and_Troy22Ba', 'diploid', 
				'Parent Control', 'PM2', 'Pullman2', 	NA, 'dubius and porrifolius',	'2879_4b_and_2878_3a', 	'diploid', 
				'Parent Control', 'PM3', 'Moscow', 		NA, 'dubius and pratensis', 	'2886_5b_and_2884_4a', 	'diploid', 
				'Parent Control', 'PM4', 'Garfield', 	NA, 'dubius and pratensis', 	'2895_1a_and_2893_4a', 	'diploid', 
				'Parent Control', 'PM5', 'Oakesdale', NA, 'dubius and pratensis', 	'2874_7a_and_2873_1a' , 'diploid', 
				'Parent Control', 'PM6', 'Spangle', 	NA, 'dubius and pratensis', 	'2860_4a_and_2861_2a',	'diploid'
        ), nrow = 6, ncol = 7, byrow = TRUE)

tableSynControl <- matrix(c(
				'Synthetic Control', 'S1B', NA, NA, 'mirus', '121_1_7',	  'tetraploid',
				'Synthetic Control', 'S1C', NA, NA, 'mirus', '77_3_5',		'tetraploid',
				'Synthetic Control', 'S1D', NA, NA, 'mirus', '116_15_2',	'tetraploid',
				'Synthetic Control', 'S1E', NA, NA, 'mirus', '99_3_9', 		'tetraploid',
				'Synthetic Control', 'S1F', NA, NA, 'mirus', '84_2_5', 		'tetraploid'
		), nrow = 5, ncol = 7, byrow = TRUE)

tableNegativeControl <- matrix(c(
				'Negative Control', 'n',			NA, NA, NA, NA, NA,
				'Negative Control', 'water',	NA, NA, NA, NA, NA
		), nrow = 2, ncol = 7, byrow = TRUE)
strNegativeControlDNA <- NA

# Special plants from other samples
tableSpecials <- matrix(c(
				'Special Natural', NA,'Spokane', 			NA, 'miscellus', '2730_3_1', 'tetraploid',
				'Special Natural', NA, 'Spokane', 			NA, 'miscellus', '2731_9_2', 'tetraploid',
				'Special Natural', NA, 'Post Falls', 			NA, 'miscellus', '2736_1_1', 'tetraploid',
				'Special Natural', NA, 'Post Falls', 			NA, 'miscellus', '2736_3_1', 'tetraploid',
				'Special Natural', NA, 'Coeur d Alene', 	NA, 'miscellus', '2738_1_2', 'tetraploid',
				'Special Natural', NA, '??', 						NA, 'miscellus', '2738_4_2', 'tetraploid'
		), nrow = 6, ncol = 7, byrow = TRUE)

#### User pick the input files #################################################
directoryIn <- tk_choose.dir(caption = 'Choose working directory')
setwd(directoryIn)

fileIn <- tk_choose.files(default = stringIn, 
		caption = 'Choose input file', multi = FALSE)
fileOut <- tk_choose.files(default = stringOut, 
                           caption = 'Choose output file', multi = FALSE)
fileLookup <- tk_choose.files(
		default = stringLookup, 
		caption = 'Choose chamber data lookup file', multi = FALSE)
directoryLog <- tk_choose.dir(caption = 'Choose logging directory')

#### Start log file ############################################################
sink(file.path(directoryLog,fileSink), append = FALSE, split = FALSE)

#### What am I doing? ##########################################################
print('Running R_link_R.R script ...')
print(paste('  -> Working directory ', getwd(), '...'))
print(paste('  -> Logging in file ', fileSink, '...'))

#### Open files for input ######################################################
# Open input file and get the number of rows
print(paste('  -> Reading', fileIn, '...'))
tableIn <- read.csv(fileIn, header = TRUE)
print('Done')

# Open lookup table
print(paste('  -> Reading', fileLookup, '...'))
tableLookup <- read.csv(fileLookup, header = TRUE)
print('Done')

#### Start processing input file ###############################################
# Loop over all of the rows in tableIn deconstructing the Sample.Id
# then using the elements of Sample.Id to generate ID_chamber as well as
# other parameters

# Number of rows in tableIn
numIn <- nrow(tableIn)
#numIn <- 2000 # Used for testing comment out for the real thing
print(paste('     -> Found', numIn, 'entries ...'))

#Populate new vectors to dataframe
tableOut <- tableIn
tableOut[c('DNA', 'Species', 'SynNat', 'Local', 'Line', 'PlantID_from_plate', 
           'Ploidy', 'Coord', 'Block')] <- 'Not_found'

# cluster <- makeCluster(4)
# registerDoParallel(cluster)
# start <-Sys.time()
#### Loop over all of the elements in tableIn parsing the Sample.Id ############
for (i in 1:numIn) {
#foreach (i=1:numIn) %dopar% {
	if (i%%10000 == 0) print(paste('     -> Working on entry', i, 'out of',numIn, '...'))
	# Break up Sample.Id at the _
	strSampleID <- c( strsplit( as.character( tableIn$Sample.Id[i] ), '_') )
	
	# Number of elements in stringID
	numSampleID <- length( strSampleID[[1]] )
	
	# Look at the first element in the Sample.ID
	# if it is 'n' or 'water' then it is a Negative Control
	indexNegativeControl <- pmatch(strSampleID[[1]][1], tableNegativeControl[ ,2], nomatch = 0)
	if (indexNegativeControl) {
		tableOut$DNA[i] <- strNegativeControlDNA
		tableOut$SynNat[i] <- tableNegativeControl[indexNegativeControl, 1]
		tableOut$Line[i]<- tableNegativeControl[indexNegativeControl, 2]
		tableOut$Local[i] <-  tableNegativeControl[indexNegativeControl, 3]
		tableOut$Species[i] <- tableNegativeControl[indexNegativeControl, 5]
		tableOut$PlantID_from_plate[i] <- tableNegativeControl[indexNegativeControl, 6]
		tableOut$Ploidy[i] <- tableNegativeControl[indexNegativeControl, 7]
		tableOut$Coord[i] <- NA
		tableOut$Block[i] <- NA
		next
	}
	
	# if it is a g or c type of DNA then continue
	indexDNA <- pmatch(strSampleID[[1]][1], tableDNA[ ,1], nomatch = 0)
	if (indexDNA) {
		tableOut$DNA[i] <- tableDNA[ indexDNA,2 ]

		# Look at the second element in Sample.ID
		# if it has a local then it is Natural
		# if it is a Roman numeral it is Synthetic
		# if it is a Parental control ...
		# if it is a Synthentic control ...
		indexLocal <- pmatch( strSampleID[[1]][2], tableLocal[ ,1], nomatch = 0 )
		indexLine  <- pmatch( strSampleID[[1]][2], tableLine[ ,1], nomatch = 0)
		indexParentControl <- pmatch( strSampleID[[1]][2], tableParentControl[ ,2], nomatch = 0)
		indexSynControl <- pmatch( strSampleID[[1]][2], tableSynControl[ ,2], nomatch = 0)
		
		# To find the special plants, 
		# in this case the special plants have a 4 digit _ 1 digit _ 1 digit
		# we strip out the gD_ and the last part to compare to the list of specials
		indexSpecial <- pmatch( substr(tableIn$Sample.Id[i],4,11), tableSpecials[ ,6], nomatch  = 0)
		
		if (indexLocal) {
			tableOut$SynNat[i] <- strNatural
			tableOut$Local[i] <- tableLocal[indexLocal, 2]
			tableOut$Line[i] <- NA
			
			# Indentify the species code from the last element
			# and replace the species name (Natural only)
			indexSpecies <- pmatch(strSampleID[[1]][numSampleID],
					tableSpecies[,1], nomatch = 0)
			if (indexSpecies) tableOut$Species[i] <- tableSpecies[indexSpecies,2]
			else {
				print(paste('    -> Unknown species on line', i+1,
								'with', tableIn$Sample.Id[i], 
								'from', tableIn$file[i]))
			}
			indexPloidy <- pmatch(tableOut$Species[i], tablePloidy[,1])
			tableOut$Ploidy[i] <- tablePloidy[indexPloidy,2]
			
			# Reconstruct the plant ID from the Sample ID by dropping the 
			# first 2 elements of SampleID and pasting with _
			strPlantID <- strSampleID[[1]][3]
			if (numSampleID > 3) {
				for ( j in 4:(numSampleID) ) {
					strPlantID <- paste(strPlantID, strSampleID[[1]][j],sep='_')
				}
			}
			tableOut$PlantID_from_plate[i] <- strPlantID
		}
		else if (indexLine) {
			tableOut$SynNat[i]  <- strSynthetic
			tableOut$Local[i] <- NA
			tableOut$Line[i]    <- tableLine[indexLine,1]
			tableOut$Species[i] <- tableLine[indexLine,2]
			indexPloidy <- pmatch(tableOut$Species[i], tablePloidy[,1])
			tableOut$Ploidy[i] <- tablePloidy[indexPloidy,2]
			
			# The number of plant id digits will tell me the synthetic local
			# If the plant is synthetic the number of terms to the plant is
			# (number of segments in the SampleID) - 2
			# 2 because the first segment is the DNA type and the second is the
			# Line ID
			tableOut$Local[i] <- vecLineLocal[numSampleID - 2]
			
			# Reconstruct the plant ID from the Sample ID by dropping the 
			# first 2 elements of SampleID and pasting with _
			strPlantID <- strSampleID[[1]][3]
			if (numSampleID > 3) {
				for ( j in 4:(numSampleID) ) {
					strPlantID <- paste(strPlantID, strSampleID[[1]][j],sep='_')
				}
			}
			tableOut$PlantID_from_plate[i] <- strPlantID
		}
		else if (indexParentControl) {
			tableOut$SynNat[i] <- tableParentControl[indexParentControl,1]
			tableOut$Local[i] <- tableParentControl[indexParentControl,3]
			tableOut$Line[i] <- tableParentControl[indexParentControl,2]
			tableOut$Species[i] <- tableParentControl[indexParentControl,5]
			tableOut$PlantID_from_plate[i] <- tableParentControl[indexParentControl,6]
			tableOut$Ploidy[i] <- tableParentControl[indexParentControl,7]
		}
		else if (indexSynControl) {
			tableOut$SynNat[i] <- tableSynControl[indexSynControl,1]
			tableOut$Local[i] <- tableSynControl[indexSynControl,3]
			tableOut$Line[i] <- tableSynControl[indexSynControl,2]
			tableOut$Species[i] <- tableSynControl[indexSynControl,5]
			tableOut$PlantID_from_plate[i] <- tableSynControl[indexSynControl,6]
			tableOut$Ploidy[i] <- tableSynControl[indexSynControl,7]
		}
		else if (indexSpecial) {
			tableOut$SynNat[i] <- tableSpecials[indexSpecial,1]
			tableOut$Local[i] <- tableSpecials[indexSpecial,3]
			tableOut$Line[i] <- tableSpecials[indexSpecial,2]
			tableOut$Species[i] <- tableSpecials[indexSpecial,5]
			tableOut$PlantID_from_plate[i] <- tableSpecials[indexSpecial,6]
			tableOut$Ploidy[i] <- tableSpecials[indexSpecial,7]
		}
		else {
			print(paste('    -> Unknown plant type on line', i+1, as.character(tableIn$Sample.Id[i])))
		}
	}
	else {
		print(paste('    -> Unknown DNA type on line', i+1))
	}
	
	indexCoord <- pmatch(tableOut$PlantID_from_plate[i], 
			tableLookup$ID.chamber.pot, nomatch = 0)
	
	if (indexCoord) {
		tableOut$Coord[i] <- as.character(tableLookup$coord[indexCoord])
		tableOut$Block[i] <- tableLookup$BLOCK[indexCoord]
	}
	else {
		tableOut$Coord[i] <- NA
		tableOut$Block[i] <- NA
	}
}
#stopCluster(cluster)
#### Factor the new data frame columns #########################################
tableOut$DNA <- as.factor(tableOut$DNA)
tableOut$Species <- as.factor(tableOut$Species)
tableOut$SynNat <- as.factor(tableOut$SynNat)
tableOut$Local <- as.factor(tableOut$Local)
tableOut$Line <- as.factor(tableOut$Line)
tableOut$PlantID_from_plate <- as.factor(tableOut$PlantID_from_plate)
tableOut$Ploidy <- as.factor(tableOut$Ploidy)
tableOut$Coord <- as.factor(tableOut$Coord)
tableOut$Block <- as.factor(tableOut$Block)

print(Sys.time()-start)
#///////////////////////////////////////////////////////////////////////////////
#
# Save files #####
#
#///////////////////////////////////////////////////////////////////////////////
outStream <- gzfile(fileOut)
write.csv(tableOut, outStream, row.names=FALSE)
print(paste('  -> Wrote', fileOut, ' ...'))
print('Done')
sink()

