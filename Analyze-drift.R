# TODO: Add comment
# 
# Author: jordobr
###############################################################################
# Working directory
directoryIn <- "data"

# Data file
fileIn <-  'calc_drift.csv'

# Datafile columns to drop after subsetting
listNamesDrop <- c('Sample.Id', 'X', 'X.1', 'Sample.Description', 
		'Well.Position', 'Description', 'Entry.Operator', 'Calibration', 
		'Mass.Shift', 'Rasters', 'Fractional.UEP', 'Fractional.Pausing', 
		'Resolution', 'Resolution.2', 'Area', 'Area.2', 'Area..Variation', 
		'Area.2.Variation', 'DNA', 'file')

# Datafile columns to prepend DNA to after subsetting
listNamesRename <- c('Call', 
		'radius', 'theta', 
		'dradius', 'dtheta',
		'cluster', 'clusterx', 'clustery', 'clusterdx', 'clusterdy', 
		'New_Call', 'Inheritance', 'Dosage')

###############################################################################
# Go to workind directory
setwd(directoryIn)

dataIn <- read.csv(fileIn, header = TRUE)

dataGDNA <- subset(dataIn, DNA == "gDNA")
dataGDNA <- dataGDNA[,!(names(dataGDNA) %in% listNamesDrop)]
for (i in 1:length(listNamesRename)) {
	names(dataGDNA)[names(dataGDNA)==listNamesRename[i]] <- 
			paste(listNamesRename,"GDNA",sep = '_')[i]
}

dataCDNA <- subset(dataIn, DNA == "cDNA")
dataCDNA$Inheritance <- as.character(dataCDNA$Inheritance)
dataCDNA$Inheritance[dataCDNA$Description == "N.No-Alleles "] <- "Lost"
dataCDNA$Inheritance[dataCDNA$Description == "N.No_Alleles "] <- "Lost"
dataCDNA$Inheritance <- as.factor(dataCDNA$Inheritance)

dataCDNA <- dataCDNA[,!(names(dataCDNA) %in% listNamesDrop)]
for (i in 1:length(listNamesRename)) {
	names(dataCDNA)[names(dataCDNA)==listNamesRename[i]] <- 
			paste(listNamesRename,"CDNA",sep = '_')[i]
}

dataOut <- merge(dataGDNA, dataCDNA)
dataOut$Drift <- paste(dataOut$Inheritance_GDNA, 'to', dataOut$Inheritance_CDNA)

tableDriftbyContig <- table(dataOut$Assay.Id, dataOut$Drift)
tableDriftbyPlant <- table(dataOut$PlantID_from_plate, dataOut$Drift)

listPlants <- rownames(tableDriftbyPlant)
numPlants <- length(listPlants)
dataPlants <- data.frame(PlantID=character(numPlants), 
		Species = character(numPlants),
		SynNat = character(numPlants),
		Local = character(numPlants),
		Line = character(numPlants),
		stringsAsFactors=FALSE)

for (i in 1:nrow(tableDriftbyPlant)) {
	intRow <- which(dataOut$PlantID_from_plate == listPlants[i])[1]
	
	dataPlants$PlantID[i] <- toString(dataOut$PlantID_from_plate[intRow])
	dataPlants$Species[i] <- toString(dataOut$Species[intRow])
	dataPlants$SynNat[i] <- toString(dataOut$SynNat[intRow])
	dataPlants$Local[i] <- toString(dataOut$Local[intRow])
	dataPlants$Line[i] <- toString(dataOut$Line[intRow])
}

write.csv(tableDriftbyPlant, file = 'DriftbyPlant.csv', quote = FALSE, row.names = TRUE)
write.csv(dataPlants, file = 'PlantInfo.csv', quote = FALSE, row.names = TRUE)

tableGDNAInherit <- table(dataGDNA$PlantID_from_plate, dataGDNA$Inheritance)
tableGDNADosage <- table(dataGDNA$PlantID_from_plate, dataGDNA$Dosage)

tableCDNAInherit <- table(dataCDNA$PlantID_from_plate, dataCDNA$Inheritance)
tableCDNADosage <- table(dataCDNA$PlantID_from_plate, dataCDNA$Dosage)

speciesTable <- table(dataGDNA$Assay.Id, dataGDNA$Species)

rm(bigTable)  
bigTable <- cbind(speciesTable, tableGDNAInherit)
bigTable <- cbind(tableGDNAInherit, tableGDNADosage)
bigTable <- cbind(bigTable,tableCDNAInherit)
bigTable <- cbind(bigTable,tableCDNADosage )

write.csv(bigTable, file = 'bigTable.csv', quote = FALSE, row.names = TRUE)
write.csv(speciesTable, file = 'temp2.csv', quote = FALSE, row.names = TRUE)

print(ftable(contigTable))
print(speciesTable)

temp <- table(tableIn$Assay.Id,tableIn$Local,tableIn$Species)
sink(file = 'out.txt')
print(ftable(temp))
sink()

print(temp)
write.csv(tableTemp, file = 'temp.csv', quote=FALSE, row.names = TRUE)

summary(dataIn$New_Call)





