#///////////////////////////////////////////////////////////////////////////////
##### Load libraries ###########################################################
library('tcltk')

##### Data entry #############################################################
# DNA
listDNA <- c("cDNA", "gDNA")

# Experiment triange parent -> children
tableExperiment <- matrix( c(
  "dubius","porrifolius","mirus", 
  "dubius", "pratensis", "miscellus"), nrow = 2, ncol = 3, byrow = TRUE)

# Filter on proportion of parents at 0 and 90 deg
cutParent <- 0.8

# Data in file
fileIn <- "combined_called.csv.gz"

# Data out file
fileOut <- "combined_inheritance.csv.gz"

# Log file
fileLog <- "log_inheritance.txt"

##### Open data files ##########################################################
# Chose input files
directoryIn <- tk_choose.dir(caption = "Choose working directory")
setwd(directoryIn)

fileIn <- tk_choose.files(default = fileIn, 
                          caption = "Choose input file", multi = FALSE)
#fileLookup <- tk_choose.files(default = "combined_parents.csv.gz", 
#                              caption = "Choose input file", multi = FALSE)

# Start log file
sink(fileLog)

#Open data file
print(paste("  -> Reading", fileIn, "..."))
tableIn <- read.csv(fileIn, header = TRUE)
print("     -> ... Done")

#Open parents file
#print(paste("  -> Reading", fileIn, "..."))
#tableLookup <- read.csv(fileLookup, header = TRUE)
#print("     -> ... Done")

##### Processing ###############################################################
# Prune incoming data set (X, X.1 ... ) is cruft from previous saves
tableIn <- subset(tableIn, select=-c(X, X.1, X.2, X.3))

# Generate list of contigs
vecAssayID <- unique(tableIn$Assay.Id)

##### Create the call frequency table ##########################################
# Subset out > 0 clusters for the parents table
tableTemp <- subset(tableIn, ((cluster > 0)))
tableCallFreq <- table( tableTemp$Assay.Id, 
                        tableTemp$DNA, 
                        tableTemp$Species, 
                        tableTemp$New_Call)
ftable(tableCallFreq)

# Create empty dataframe 
tableContig <- data.frame()

##### Loop over Assay ID's to calculate inheritance ############################
for (i in 1: length(vecAssayID)) {
#for (i in 1:10) {
  for (j in 1: length(listDNA)) {
    for (k in 1: nrow(tableExperiment) ) {
      strAssayID <- vecAssayID[i]
      
      strDNA <- listDNA[j]
      
      strParent1 <- tableExperiment[k,1] 
      strParent2 <- tableExperiment[k,2] 
      strChild <- tableExperiment[k,3]
      
      indexAssayID <- which(levels(tableIn$Assay.Id)==strAssayID)
      indexDNA <- which(levels(tableIn$DNA)==strDNA)
      indexParent1 <- which(levels(tableIn$Species)==strParent1)
      indexParent2 <- which(levels(tableIn$Species)==strParent2)
      indexChild <- which(levels(tableIn$Species)==strChild)
      
      tableChild <- subset(tableIn, ((Assay.Id==vecAssayID[i]) 
                                     & (DNA==listDNA[j]) 
                                     & (Species==strChild)), 
                           select=c(Assay.Id, Sample.Id, New_Call))
      if (nrow(tableChild) == 0) next
      
      tableChild$Dosage <- NA
      tableChild$Dosage_Description <- NA
      
      tableTemp <- as.data.frame(subset(
        prop.table(tableCallFreq[indexAssayID, indexDNA, ,],1), 
        select=c(HH,LL)))
      
      ##### Calculate parent 1 and parent 2 calls with proportions #############
      numParent1Prop <- as.numeric(apply(tableTemp[indexParent1, ], 1, max))
      strParent1Call <- toString(colnames(tableTemp)[
        as.numeric(apply(tableTemp[indexParent1, ], 1, which.max))])
      
      numParent2Prop<- as.numeric(apply(tableTemp[indexParent2, ], 1, max))
      strParent2Call <- toString(colnames(tableTemp)[
        as.numeric(apply(tableTemp[indexParent2, ], 1, which.max))])
      
      
      print(paste("     -> Working on ", vecAssayID[i], listDNA[j], "and", strChild, "..."))
      print(paste("     -> Found", strParent1, strParent1Call, numParent1Prop, 
                  strParent2, strParent2Call, numParent2Prop))
      
      ##### Set good parent flags ##############################################
      flagParent1 <- ((numParent1Prop > cutParent) & !is.na(numParent1Prop))
      flagParent2 <- ((numParent2Prop > cutParent) & !is.na(numParent2Prop))
      flagParentBoth <- !(toString(strParent1Call)==toString(strParent2Call))
      
      
      ###### Caluclate the inheritance #########################################
      if ((flagParent1) & (flagParent2) & (flagParentBoth)) {
        charParent1 <- substr(strParent1Call,1,1)
        for (l in 1: nrow(tableChild)) {
          if (length(unlist(strsplit(as.character(tableChild$New_Call[l]),"")))== 4) { 
            tableChild$Dosage[l] <- length(grep(charParent1, unlist(strsplit(as.character(tableChild$New_Call[l]),""))))
            if (tableChild$Dosage[l] == 2) tableChild$Dosage_Description[l] <- "Equal"
            else if (tableChild$Dosage[l] < 2) tableChild$Dosage_Description[l] <- paste(strParent1, "-like", sep = "")
            else if (tableChild$Dosage[l] > 2) tableChild$Dosage_Description[l] <- paste(strParent2, "-like", sep = "")    
          }
          else {
            tableChild$Dosage_Description[l] <- as.character(tableChild$New_Call[l])
          }
        }
      }
      else if (!(flagParentBoth) & !is.na(flagParentBoth)) {
        tableChild$Dosage_Description <- paste("Same_contig_for", strParent1, "and", strParent2, sep="_")
      }
      else if (!(flagParent1) & (flagParent2)) {
        tableChild$Dosage_Description <- paste("Prop_bad_for", strParent1, sep="_")
        print(paste("Prop_bad_for", strParent1, sep="_"))
      }
      else if ((flagParent1) & !(flagParent2)) {
        tableChild$Dosage_Description <- paste("Prop_bad_for", strParent2, sep="_")
        print( paste("Prop_bad_for", strParent2, sep="_"))
      }
      else {
        tableChild$Dosage_Description <- paste("Prop_bad_for", strParent1, "and", strParent2, sep="_")
        print(paste("Prop_bad_for", strParent1, "and", strParent2, sep="_"))
      }
      
      tableChild$Parent1_Call <- strParent1Call
      tableChild$Parent2_Call <- strParent2Call
      tableChild$Parent1_Prop <- numParent1Prop
      tableChild$Parent2_Prop <- numParent2Prop
   
      tableContig <- rbind(tableContig, tableChild)
    }      
  } 
}
tableContig$Dosage_Description <- as.factor(tableContig$Dosage_Description)
tableContig$Dosage <- as.numeric(tableContig$Dosage)

tableOut <- merge(tableIn, tableContig, all=TRUE, by=c("Assay.Id", "Sample.Id", "New_Call"))

tableGDNA <- subset(tableOut, DNA == "gDNA", 
                   select=c(Assay.Id, PlantID_from_plate, Dosage_Description))
names(tableGDNA)[names(tableGDNA)=="Dosage_Description"] <- "Dosage_GDNA"

tableCDNA <- subset(tableOut, DNA == "cDNA", 
                   select=c(Assay.Id, Sample.Id, PlantID_from_plate, Description, Dosage_Description))

tableCDNA$Dosage_Description <- as.character(tableCDNA$Dosage_Description)
tableCDNA$Dosage_Description[tableCDNA$Description == "N.No-Alleles "] <- "Lost"
tableCDNA$Dosage_Description[tableCDNA$Description == "N.No_Alleles "] <- "Lost"
tableCDNA$Dosage_Description <- as.factor(tableCDNA$Dosage_Description)

tableCDNA <- merge(tableCDNA, tableGDNA, by=c("Assay.Id", "PlantID_from_plate"))

tableCDNA$Drift <- as.factor(paste(tableCDNA$Dosage_GDNA, 'to', tableCDNA$Dosage_Description))

tableCDNA <- subset(tableCDNA, select=c('Assay.Id', 'Sample.Id', 'Drift'))

tableOut <- merge(tableOut, tableCDNA, all=TRUE, by=c("Assay.Id", "Sample.Id"))

#///////////////////////////////////////////////////////////////////////////////
#
# Save files #####
#
#///////////////////////////////////////////////////////////////////////////////
# Save data file
print(paste("  -> Saving", fileOut, "..."))
outStream <- gzfile(fileOut)
write.csv(tableOut, outStream)
print("     -> ... Done")

sink()

