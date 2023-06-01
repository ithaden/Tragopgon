###############################################################################
# R_cluster.R
#
# Short summary
# 
# Input: Add inputs
#
# Output: Add outputs
#
# TODO: Add comment
# 
# Author: jordobr
# Date: Nov 7, 2012
###############################################################################
library('lattice') 
library('fpc')
library('sm')
library('mixtools')
library('tcltk')

################################################################################
#
# User needs to add information here
#
################################################################################
# Output file
fileOut <- "combined_clustered.csv"

# Lookup table for low mass and hight mass with the appropriate 
# experiment triangle
listSpeciesLookupmir <- c("dubius","porrifolius","mirus")
listSpeciesLookupmisc <- c("dubius", "pratensis", "miscellus")

#Plot output file
filePlots <- "R_diagnostics.pdf"

# Log file
fileSink <- "log_cluster.txt"

# Experment parameters
listDNA <- c("gDNA", "cDNA")
listSynNat <- c("Natural","Synthetic")

# List of species
listSpecies <- c("dubius", "mirus", "miscellus", "porrifolius", "pratensis" )
listSpeciesploidy <- c(2,2,4,4,2,2) 
listColors <- c("red", "blue", "green", "magenta", "purple")

# Bad Allelles Description
strNoAlleles <- c('N.No-Alleles','N.No_Alleles')

# Bad Spectrum Description
strBadSpectrum <- "I.Bad Spectrum"

# Filter on dtheta, the is related to A1 A2 variance
cutdtheta <- 15.0

# Window on 0 deg and 90 deg axis. Automatically call these 
# as parents
cutParent <- 5.0

# Pre-call clusters
vecHomo <- c(0,90)

# Call Limits
vecDipoloidcuts <- c(22.5, 67.5)
vecTetraploid <- c(11.25, 33.75, 56.25, 78.75)
			
################################################################################
#
# Parameters
#
################################################################################
# Really small number
eps <- 0.0000001

# Number of DNA types
numDNA <- length(listDNA)

# Number of species
numSpecies <- length(listSpecies)

################################################################################
#
# Preliminaries
#
################################################################################

# Chose input files
directoryIn <- tk_choose.dir(caption = "Choose working directory")

setwd(directoryIn)
fileIn <- tk_choose.files(default = ""combined_highlow.csv.gz"", 
		caption = "Choose input file", multi = FALSE)
fileLookupmir <- tk_choose.files(
		default = "Mirus_loss_and_de_SNPs.csv", 
		caption = "Choose Mirus SNPs lookup file", multi = FALSE)
fileLookupmisc <- tk_choose.files(
		default = "Misc_loss_and_de_snps.csv", 
		caption = "Choose Miscellus SNPs lookup file", multi = FALSE)
directoryLogs <- tk_choose.dir(caption = "Choose logging directory")
#directoryLog <- 'logs'
directoryPlots <- tk_choose.dir(caption = "Choose directory for plots")
#directoryPlots <- 'plots'

# Set plotting file
#filePlots <- file.path(directoryPlots, filePlots)

# Open log file
#sink(file.path(directoryLogs,fileSink), append = FALSE, split = FALSE)


# What am I doing?
print("Running R_cluster_R.R script ...")
print(paste("  -> Working directory", getwd(), "..."))
print(paste("  -> Logging in file", fileSink, "..."))

################################################################################
#
# Open files
#
################################################################################
#Open data file

print(paste("  -> Reading", fileIn, "..."))
tableIn <- read.csv(fileIn, header = TRUE)
print("     -> ... Done")

fileMir <- fileLookupmir
print(paste("  -> Reading", fileMir, "..."))
tableLookupmir <- read.csv(fileMir, header = TRUE)
print("     -> ... Done")

fileMisc <- fileLookupmisc
print(paste("  -> Reading", fileMisc, "..."))
tableLookupmisc <- read.csv(fileMisc, header = TRUE)
print("     -> ... Done")

################################################################################
#
# Start analysis
#
################################################################################
# With Area 2 on the x-axis and Area 1 on the y axis, transform to 
# polar coodinates (r and theta) where theta is in degrees.
# Calculate the uncertainty in radius and angle theta

print("  -> Transforming to Area 1 and Area 2 polar coordinates ...")
tableIn$radius <- sqrt(tableIn$Area^2 + tableIn$Area.2^2)
tableIn$theta <- atan(tableIn$Area / (tableIn$Area.2 + eps)) * 180 / pi
tableIn$dradius <- sqrt( (tableIn$Area^2 * tableIn$Area..Variation + tableIn$Area.2^2 * tableIn$Area.2.Variation) / (tableIn$Area^2 + tableIn$Area.2^2 + eps))
tableIn$dtheta <- sqrt(tableIn$Area.2^2 * tableIn$Area..Variation + tableIn$Area^2 * tableIn$Area.2.Variation) / (tableIn$Area^2 + tableIn$Area.2^2 + eps) * 180 / pi

# Generate list of AssayID's from the data file
listAssayID <- unique(as.character(tableIn$Assay.Id))
numAssayID <- length(listAssayID)
print(paste("     -> Found", numAssayID, "Assay ID's ..."))

# Used for testing
#numAssayID <- 10
tableOut <- data.frame()

#print(paste("     -> Saving plots to", filePlots, "..."))
#pdf(filePlots)

for (i in 1:numAssayID) {
	strAssayID <- listAssayID[i]
	print("----------------------------------------------------------")
	print(paste("-> Working on Assay ID", strAssayID, 
					"number", i, "out of", numAssayID, "..."))
	
	# Make a table with the Assay ID of interest
	tableSubset <- subset(tableIn, (( Assay.Id == strAssayID ) 
                                  & ( SynNat %in% listSynNat)))
  	
	for (j in 1:length(listDNA)) {
	  tablePlot <- data.frame()
    
		strDNA <- listDNA[j]
		print(paste("  -> Working on", strDNA, "..."))
		
		# Make a table with DNA of interest
		tableSubsetDNA <- subset(tableSubset, DNA == strDNA)
		
		# Make a list of the species in the subset
		listSubsetSpecies <- unique(as.character(tableSubsetDNA$Species))
		listSubsetSpecies <- intersect(listSubsetSpecies, listSpecies)
		numSubsetSpecies <- length(listSubsetSpecies)
		print(paste("    -> Found", numSubsetSpecies, "species ..."))
		
		for (k in 1:numSubsetSpecies) {
			# Zero species
			tableBadSpectrum <- data.frame()
			tableBadTheta <- data.frame()
			tableNoAllele <- data.frame()
			tableLowmass <- data.frame()
			tableHighmass <- data.frame()
			
			strSpecies <- listSubsetSpecies[k]
			print(paste("    -> Working on species", strSpecies, "..."))
			
			tableSubsetSpecies <- subset(tableSubsetDNA, Species == strSpecies)
			
			# Search the order form for high and low mass for calls based upon
			# the experiment triangle
			indexAssayIDmir <- match(strAssayID, 
					tableLookupmir$SNP_ID, nomatch = 0)
			indexAssayIDmisc <- match(strAssayID, 
					tableLookupmisc$SNP_ID, nomatch = 0)
			
			if ((indexAssayIDmisc == 0) & 
					(strSpecies %in% listSpeciesLookupmir)) {
				strMasslow  <- tableLookupmir$EXT1_CALL[indexAssayIDmir]
				strMasshigh <- tableLookupmir$EXT2_CALL[indexAssayIDmir]
			}
			else if ((indexAssayIDmir == 0) & 
					(strSpecies %in% listSpeciesLookupmisc)) {
				strMasslow  <- tableLookupmisc$EXT1_CALL[indexAssayIDmisc]
				strMasshigh <- tableLookupmisc$EXT2_CALL[indexAssayIDmisc]
			}
			else if ((indexAssayIDmisc != 0) & indexAssayIDmir != 0) {
				print(paste("   -> Assay", strAssayID,
								"is in both mirus and miscellus list ..."))
				
				strMasslowmir  <- tableLookupmir$EXT1_CALL[indexAssayIDmir]
				strMasshighmir <- tableLookupmir$EXT2_CALL[indexAssayIDmir]
				
				strMasslowmisc  <- tableLookupmisc$EXT1_CALL[indexAssayIDmisc]
				strMasshighmisc <- tableLookupmisc$EXT2_CALL[indexAssayIDmisc]
				
				if ((strMasslowmir == strMasslowmisc) & 
						(strMasshighmir == strMasshighmisc)) {
					print("    -> No worries, masses are the same ...")
					strMasslow <- strMasslowmir
					strMasshigh <- strMasshighmir }
				else {
					print(paste("   -> Error Assay", strAssayID, 
									"species", strSpecies, 
									"needs further review 
											setting masses to L and H..."))
					strMasslow  <- "L"
					strMasshigh <- "H"	
					
				}
			}
			
			else {
				print(paste("   -> Error did not find SNP for", strAssayID,
								"and", strDNA,
								"setting masses to L and H..."))
				
				strMasslow  <- "L"
				strMasshigh <- "H"
			}
			
			tableSubsetSpecies$HighCall <- strMasshigh
			tableSubsetSpecies$LowCall <- strMasslow
			print(paste("      -> Using",strMasslow, "for the low mass and", 
							strMasshigh, "for the high mass ..." ))
			
			
			print(paste("      -> Found", 
							length(unique(as.character(tableSubsetSpecies$PlantID_from_plate))),
							"plant ID's ..."))
			
			# Remove the samples with bad spectrum
			tableBadSpectrum <- subset(tableSubsetSpecies, 
					Description == strBadSpectrum)
			numBadSpectrum <- nrow(tableBadSpectrum)
			print(paste("      -> Found", numBadSpectrum, 
							"samples with a bad spectrum" ))
			tableSubsetSpecies <- subset(tableSubsetSpecies, 
					Description != strBadSpectrum)
			if (numBadSpectrum > 0) {
				tableBadSpectrum$cluster <- -2
				tableBadSpectrum$clusterx <- NA
				tableBadSpectrum$clustery <- NA
				tableBadSpectrum$clusterdx <- NA
				tableBadSpectrum$clusterdy <- NA
			}
			
			# Remove the samples with large dtheta
			tableBadTheta <- subset(tableSubsetSpecies, 
					dtheta > cutdtheta)
			numBadTheta <- nrow(tableBadTheta)
			print(paste("      -> Found", numBadTheta, 
							"samples with a bad theta" ))
			tableSubsetSpecies <- subset(tableSubsetSpecies, 
					dtheta < cutdtheta)
			if (numBadTheta > 0) {
				tableBadTheta$cluster <- -1
				tableBadTheta$clusterx <- NA
				tableBadTheta$clustery <- NA
				tableBadTheta$clusterdx <- NA
				tableBadTheta$clusterdy <- NA
			}
			
			# Remove the N.No−Alleles and N.No_Alleles
			tableNoAllele <- subset(tableSubsetSpecies, Description %in% strNoAlleles )
			numNoAllele <- nrow(tableNoAllele)
#			print(summary(tableSubsetSpecies$Description))
			print(paste("      -> Found", numNoAllele, 
							"samples with No Alleles" ))
			tableSubsetSpecies <- subset(tableSubsetSpecies, !Description %in% strNoAlleles)
			if (numNoAllele > 0) {
				tableNoAllele$cluster <- 0
				tableNoAllele$clusterx <- NA
				tableNoAllele$clustery <- NA
				tableNoAllele$clusterdx <- NA
				tableNoAllele$clusterdy <- NA
			}
			# How many samples are left?
			print(paste("      -> Found", nrow(tableSubsetSpecies), 
							"remaining samples" ))
			
			# Remove the samples at 0 and 90 deg calling them parents
			tableHighmass <- subset(tableSubsetSpecies, theta < cutParent)
			numHighmass <- nrow(tableHighmass)
			print(paste("      -> Found", numHighmass, 
							"samples at 0 deg" ))
			if (numHighmass > 0) {
				tableHighmass$cluster <- 1
				tableHighmass$clusterx <- 0
				tableHighmass$clustery <- mean(tableHighmass$radius)
				tableHighmass$clusterdx <- 1
				tableHighmass$clusterdy <- var(tableHighmass$radius)
			}
			
			tableLowmass <- subset(tableSubsetSpecies, (theta > (90 - cutParent)))
			numLowmass <- nrow(tableLowmass)
			print(paste("      -> Found", numLowmass, 
							"samples at 90 deg" ))
			if (numLowmass > 0){
				tableLowmass$cluster <- 2
				tableLowmass$clusterx <- 90
				tableLowmass$clustery <- mean(tableLowmass$radius)
				tableLowmass$clusterdx <- 1
				tableLowmass$clusterdy <- var(tableLowmass$radius)
			}
			
			tableSubsetSpecies <- subset(tableSubsetSpecies, 
					((theta > cutParent) & (theta < (90 - cutParent))))
			
			#Start clustering to look for natural groupings for 
			# the off axis
			if (nrow(tableSubsetSpecies) > 1) {
				listFit <- cbind(tableSubsetSpecies$theta,tableSubsetSpecies$radius)
				fit <- Mclust(listFit)
				#fitA <- Mclust(cbind(tableSubsetSpecies$Area,tableSubsetSpecies$Area.2))
				#print(summary(fit))
				print(paste("      -> Found", fit$G, "clusters" ))
				
				# Add cluster info to tableSubset dataframe
				tableSubsetSpecies$cluster <- fit$classification + 2
				tableSubsetSpecies$clusterx <- fit$parameter$mean[1, fit$classification]
				tableSubsetSpecies$clustery <- fit$parameter$mean[2, fit$classification]
				tableSubsetSpecies$clusterdx <- fit$parameter$variance$sigma[1, 1, fit$classification]
				tableSubsetSpecies$clusterdy <- fit$parameter$variance$sigma[2, 2, fit$classification]
			}
			else if (nrow(tableSubsetSpecies) == 1) {
				tableSubsetSpecies$cluster <- 3
				tableSubsetSpecies$clusterx <- tableSubsetSpecies$theta
				tableSubsetSpecies$clustery <- tableSubsetSpecies$radius
				tableSubsetSpecies$clusterdx <- tableSubsetSpecies$dtheta
				tableSubsetSpecies$clusterdy <- tableSubsetSpecies$dradius
			}
      tablePlot <- rbind(tablePlot, tableBadSpectrum, tableBadTheta, tableNoAllele,
                  tableLowmass, tableHighmass, tableSubsetSpecies)

			# add new data to tableOut
			tableOut <- rbind(tableOut, tableBadSpectrum, tableBadTheta, tableNoAllele,
					tableLowmass, tableHighmass, tableSubsetSpecies)
		}
pdf(file.path(directoryPlots, paste('cluster ',strAssayID,' ', strDNA,'.pdf', sep = '')))  
print(
  xyplot(radius~theta | Species, 
         group = cluster, 
         data = tablePlot,
         main = paste(strAssayID,strDNA),
         xlab = "theta (deg)",
         ylab = "radius (arb)",
         xlim = c(0,90),
         scales = list(tick.number = 10),
         auto.key=TRUE
         )
)
dev.off()
	}
}
#dev.off()

################################################################################
#
# Save files
#
################################################################################
# Save data file
fileOut <- file.path(directoryIn, fileOut)
print(paste("  -> Saving", fileOut, "..."))
write.csv(tableOut, fileOut)
print("     -> ... Done")

sink()

################################################################################
#
# Old Cruft
#
################################################################################
#print(xyplot(radius ~ theta | Call * Description,
#				groups = cluster,
#				data = tableSubsetSpecies,
#				type = "p",
#				auto.key = TRUE,
#				main = paste(strAssayID,strDNA,strSpecies),
#				xlab = "Theta (deg)",
#				ylab = "r (arb)",
#				xlim = c(0,90),
#				scales = list(tick.number = 10),
#				panel = function(...) {
#					panel.grid(h=0, v= 8, col="blue")
#					panel.xyplot(...)}))
#
#print(xyplot(radius ~ theta | SynNat,
#				groups = cluster,
#				data = tableSubsetSpecies,
#				type = "p",
#				auto.key = TRUE,
#				main = paste(strAssayID,strDNA,strSpecies),
#				xlab = "Theta (deg)",
#				ylab = "r (arb)",
#				xlim = c(-4,94),
#				scales = list(tick.number = 10),
#				panel = function(...) {
#					panel.grid(h=0, v= 8, col="blue")
#					panel.xyplot(...)}))
#
#for (i in 1:numAssayID) {
#	print(paste("     -> Working on Assay ID", listAssayID[i], 
#					"number", i, "out of", numAssayID, "..."))
#	
#	for (j in 1:numSpecies){
#		tableSubset <- subset(tableIn, Assay.Id == listAssayID[i] &
#						Species == listSpecies[j])
#		
#		numFound <- nrow(tableSubset)
#		print(paste("          -> Number of", listSpecies[j], 
#						"found:", numFound, "..."))
#		
#		if (numFound != 0) {
#			# Reduce data frame to x,y elements for cluster analysis
#			x  <- data.frame(tableSubset$theta, tableSubset$radius)
#			
#			# Run the clustering routine
#			tableResults <- dbscan(x, 3, 5, method = c("raw"), showplot = 0)
#			
#			print(paste("               -> Found", max(tableResults$cluster), 
#							"clusters ..."))
#			
#			tableSubset$cluster <- tableResults$cluster
#			
#			print(histogram(~dtheta | Description,
#							data = tableSubset,
#							type = "count",
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "d theta",
#							ylab = "Frequency"))
#			print(histogram(~Resolution | Description,
#							data = tableSubset,
#							type = "count",
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "Resolution 1 (m/dm)",
#							ylab = "Frequency"))
#			print(histogram(~Resolution.2 | Description,
#							data = tableSubset,
#							type = "count",
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "Resolution 2 (m/dm)",
#							ylab = "Frequency"))
#			print(histogram(~sqrt(Resolution^2 + Resolution.2^2) | Description,
#							data = tableSubset,
#							type = "count",
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "sqrt(Resolution 1^2 + Resolution 2^2) (m/dm)",
#							ylab = "Frequency"))
#			print(histogram(~dradius | Description,
#							data = tableSubset,
#							type = "count",
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "d radius",
#							ylab = "Frequency"))
#			print(histogram(~dtheta * dradius * radius | Description,
#							data = tableSubset,
#							type = "count",
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "dtheta * dradius * radius",
#							ylab = "Frequency"))
#			print(histogram(~Area | Description,
#							data = tableSubset,
#							type = "count",
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "Area",
#							ylab = "Frequency"))
#			print(histogram(~Area.2 | Description,
#							data = tableSubset,
#							type = "count",
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "Area 2",
#							ylab = "Frequency"))
#			print(histogram(~(Area + Area.2) | Description,
#							data = tableSubset,
#							type = "count",
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "Area 1 + Area 2",
#							ylab = "Frequency"))
#			
#			print(xyplot(Area..Variation/(Area+eps) ~ 
#									Area.2.Variation/(Area.2+eps),
#							groups = Description,
#							data = tableSubset,
#							type = "p",
#							auto.key = TRUE,
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "dA2/A2",
#							ylab = "dA1/A1",
#							aspect = "iso"))
#			
#			print(xyplot(Area ~ Area.2 | Description,
#							groups = cluster,
#							data = tableSubset,
#							type = "p",
#							auto.key = TRUE,
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "A2",
#							ylab = "A1",
#							aspect = "iso"))
#			
#			print(xyplot(radius ~ theta | Description,
#							groups = cluster,
#							data = tableSubset,
#							type = "p",
#							auto.key = TRUE,
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "Theta (deg)",
#							ylab = "r (arb)",
#							xlim = c(-4,94)))
#			
#			print(xyplot(Resolution ~ Resolution.2 | Description,
#							groups = cluster,
#							data = tableSubset,
#							type = "p",
#							auto.key = TRUE,
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							xlab = "Resolution 2 (arb)",
#							ylab = "Resolution 1 (arb)",
#							xlim = c(-4,94)))
#			
#			print(densityplot(~theta,
#							data = tableSubset,
#							auto.key = TRUE,
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							from = -20,
#							to = 110,
#							xlab = "Theta (deg)",
#							xlim = c(-4,94)))
#			
#			print(densityplot(~theta,
#							group = Description,
#							data = tableSubset,
#							auto.key = TRUE,
#							main = listAssayID[i],
#							sub = listSpecies[j],
#							from = -20,
#							to = 110,
#							xlab = "Theta (deg)",
#							xlim = c(-4,94)))
#			
#			out <- density(tableSubset$theta,bw = 2,
#					from = -20,
#					to = 110)
#			x <- out$x
#			y <- out$y
#			
#			fit <- nls(y ~
#							(f1/s1)*exp(-(x-u1)^2/(2*s1^2)) +
#							(f2/s2)*exp(-(x-u2)^2/(2*s2^2)) +
#							(f3/s3)*exp(-(x-u3)^2/(2*s3^2)) +
#							(f4/s4)*exp(-(x-u4)^2/(2*s4^2)) +
#							(f5/s5)*exp(-(x-u5)^2/(2*s5^2)),
#					start = list(f1 = 1, s1 = 5, u1 = 0,
#							f2 = 1, s2 = 5, u2 = 27,
#							f3 = 1, s3 = 5, u3 = 45,
#							f4 = 1, s4 = 5, u4 = 72,
#							f5 = 1, s5 = 5, u5 = 90), trace=TRUE,
#					nls.control(maxiter = 50, tol = 1e-03, minFactor = 1/4096,
#							printEval = TRUE, warnOnly = TRUE))
#			print(fit)
#		}
#	}
#}

#tableSubset$Species <- factor(tableSubset$Species, levels = listspecies)
#
#data <- cbind(tableSubset$theta,tableSubset$radius)
#
#fit <- Mclust(data)
#
#print(mclust2Dplot(data = data, what = "classification", identify = TRUE,
#				parameters = fit$parameters, z = fit$z))
#
#print(xyplot(radius ~ theta | Ploidy * DNA, 
#				data = tableSubset
#		))

## Lookup table
#fileLookupmir  <- "Mirus_loss_and_de_SNPs.csv"
#fileLookupmisc <- "Misc_loss_and_de_snps.csv"
#
## Filter parameters
#cutTheta <- 10
#
#thresholdTheta1 <- 0.0
#thresholdTheta2 <- 22.5
#thresholdTheta3 <- 45.0
#thresholdTheta4 <- 67.5
#thresholdTheta5 <- 90.0
##Plot colors and points
##Now make some pretty plots
#listspecies <- c("dubius", "pratensis", "porrifolius", "miscellus", "mirus")
#listcolors <- c("blue", "purple", "green", "red", "yellow")
#listpoints  <- c(1, 1, 1, 2, 2)
## Open lookup tables
#tableLookupmir <- read.csv(fileLookupmir, header = TRUE)
#print(paste("  -> Reading", fileLookupmir, "..."))
#tableLookupmisc <- read.csv(fileLookupmisc, header = TRUE)
#print(paste("  -> Reading", fileLookupmisc, "..."))

#vecNumberSpecies <- mat.or.vec(numAssayID,1)
#
##Plot uncertainty in theta versus theta for diagnostic purposes
#print(xyplot(Area ~ Area.2 | DNA * SynNat,
#				data = tableSubset,
#				group = Species,
#				main = paste(vectorAssayID[i]),
#				aspect = 1,
#				xlab = "Area 2 (arb)",
#				ylab = "Area 1 (arb)",
#				par.settings = list(superpose.symbol = list(
#								pch = listpoints, col = listcolors)),
#				auto.key=TRUE
#		))
#
#print(xyplot(radius ~ theta | DNA,
#				data = tableSubset,
#				group = Species,
#				main = paste(vectorAssayID[i]),
#				aspect = 1,
#				xlab = "sqrt(Area1^2 + Area2^2)",
#				ylab = "Atan(Area1 / Area2)",
#				par.settings = list(superpose.symbol = list(
#								pch = listpoints, col = listcolors)),
#				auto.key=TRUE
#		))
#
#print(histogram(~ radius | DNA * Species,
#				data = tableSubset,
#				group = Species,
#				main = paste(vectorAssayID[i]),
#				aspect = "fill",
#				xlab = "sqrt(Area1^2 + Area2^2)",
#				ylab = "Frequency"
#		))
#
#print(histogram(~ radius_D | DNA * Species,
#				data = tableSubset,
#				group = Species,
#				main = paste(vectorAssayID[i]),
#				aspect = "fill",
#				xlab = "Delta sqrt(Area1^2 + Area2^2)",
#				ylab = "Frequency"
#		))
#
#print(xyplot(radius_D ~ radius | DNA,
#				data = tableSubset,
#				group = Species,
#				main = paste(vectorAssayID[i]),
#				xlab = "sqrt(Area1^2 + Area2^2)",
#				ylab = "Delta sqrt(Area1^2 + Area2^2)",
#				par.settings = list(superpose.symbol = list(
#								pch = listpoints, col = listcolors)),
#				auto.key=TRUE
#		))
#
#print(histogram(~ theta | DNA * Species,
#				data = tableSubset,
#				group = Species,
#				main = paste(vectorAssayID[i]),
#				aspect = "fill",
#				xlim = c(0, 90),
#				xlab = "Atan(Area1 / Area2) (deg)",
#				ylab = "Frequency"
#		))
#
#print(histogram(~ theta_D | DNA * Species,
#				data = tableSubset,
#				group = Species,
#				main = paste(vectorAssayID[i]),
#				aspect = "fill",
#				xlim = c(0, 90),
#				xlab = "Delta Atan(Area1 / Area2) (deg)",
#				ylab = "Frequency"
#		))
#
#print(xyplot(theta_D ~ theta | DNA,
#				data = tableSubset,
#				group = Species,
#				main = paste(vectorAssayID[i]),
#				xlim = c(0, 90),
#				ylim = c(0,90),
#				xlab = "Atan(Area1 / Area2) (deg)",
#				ylab = "Delta Atan(Area1 / Area2) (deg)",
#				par.settings = list(superpose.symbol = list(
#								pch = listpoints, col = listcolors)),
#				auto.key=TRUE
#		))
#	tableSubset <- subset(tableIn, Assay.Id == vectorAssayID[i] &
#					(Species == "dubius" | Species == "pratensis" | Species == "porrifolius" | Species == "miscellus" | Species == "mirus") &
#					(DNA == "cDNA" | DNA == "gDNA")  
#	)
#	# Refactor Species column for pretty plots
#	tableSubset$Species <- factor(tableSubset$Species, levels = listspecies)
#listPairs <- matrix(c("mirus",     "porrifolius", "dubius", 
#				"miscellus", "dubius",      "pratensis"),
#		nrow = 2, ncol = 3, byrow = TRUE) 