library('tcltk')

# Lookup table for low mass and hight mass with the appropriate 
# experiment triangle

listDNA <- c("cDNA", "gDNA")
vecDNA <- c(1, 2)

# Parents to examine
listSpecies <- c("dubius", "porrifolius", "pratensis")
vecSpecies <-c(1, 4, 5)

# Output file
fileOut <- "combined_parents.csv.gz"

# Chose input files
directoryIn <- tk_choose.dir(caption = "Choose working directory")
setwd(directoryIn)
fileIn <- tk_choose.files(default = "combined_called.csv.gz", 
                          caption = "Choose input file", multi = FALSE)

#Open data file
print(paste("  -> Reading", fileIn, "..."))
tableIn <- read.csv(fileIn, header = TRUE)
print("     -> ... Done")

# Drop row numbers that pollute dataset
tableIn <- subset(tableIn, select=-c(X, X.1, X.2, X.3))

# Subset out > 0 clusters
tableData <- subset(tableIn, ((cluster > 0)))
freqData <- table( tableData$Assay.Id, tableData$DNA, tableData$Species, tableData$New_Call)
ftable(freqData)

tableOut <- data.frame()
Assay.Id <- levels(tableData$Assay.Id)

for (i in 1: length(listDNA)){
  for (k in 1: length(listSpecies)) {
    strDNA <- listDNA[i]
    strSpecies <- listSpecies[k]
    intDNA <- vecDNA[i]
    intSpecies <- vecSpecies[k]
    
    tableTemp <- as.data.frame(subset(prop.table(freqData[,intDNA, intSpecies,],1), select=c(HH,LL)))
    tableContig <- as.data.frame(Assay.Id)
    tableContig$Contig_Prop <- as.numeric(apply(tableTemp, 1, max))
    tableContig$Contig_Call <- as.factor(colnames(tableTemp)[as.numeric(apply(tableTemp, 1, which.max))])
    tableContig$Species <- strSpecies
    tableContig$DNA <- strDNA
    #    names(tableContig)[names(tableContig) == 'Call'] <- paste(strDNA, strSpecies, "call", sep="_")
    #    names(tableContig)[names(tableContig) == 'Prop'] <- paste(strDNA, strSpecies, "prop", sep="_")
    tableOut <- rbind(tableOut,tableContig)
  }
}

# Save data file
print(paste("  -> Saving", fileOut, "..."))
outStream <- gzfile(fileOut)
write.csv(tableOut, outStream)
print("     -> ... Done")
