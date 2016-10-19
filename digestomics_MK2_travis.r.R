###############################################################################
# loadMascot () read and order a mascot file
# 
# Author: Travis
###############################################################################
loadMascot <- function(sMascotFilename = "")
{
	## Read Mascot file
	if (sMascotFilename == "")
		dfMascot <- read.csv( file.choose(), head = TRUE, stringsAsFactors = FALSE)
	else
		dfMascot <- read.csv( sMascotFilename, head = TRUE, stringsAsFactors = FALSE)	
	
	dfMascot <- dfMascot[order(dfMascot$compound),]	## order by compound  (Polypeptide) 128s for 0.95 GB file
	
	return(dfMascot)
}

###############################################################################
# preparePeptideList() takes a data.frame object filled with data read from
#						a mascot file. generates a list of unique peptides
#						and records all associated information
# 
# Author: Travis
###############################################################################
##	system.time(lPep <- preparePeptideList("D:\\Project Source Files\\mascot_map.csv"))
preparePeptideList <- function(dfMascot)
{
	## generate a list of all unique polypeptides
	uCompound <- unique(dfMascot$compound)
	
	# create list of correct size first!!
	#lPeptides <- list() 1.75 mascot_map2.csv
	lPeptides <- vector(mode = "list", length = length(uCompound)) # 1.70
	names(lPeptides) <- uCompound
	
	#set up list of unique polypeptides
	for (i in 1:length(uCompound))
	{
		vPepIndices <- which(dfMascot$compound == uCompound[i])
		iLen <- length(vPepIndices)
		
	## Find all incidences of a peptide in the file and add the values			
		lPeptides[[i]]$count			<- iLen
		lPeptides[[i]]$sample[1:iLen]	<- dfMascot$sample[vPepIndices]
		lPeptides[[i]]$score[1:iLen]	<- dfMascot$score[vPepIndices]
		lPeptides[[i]]$mScore[1:iLen]	<- dfMascot$mScore[vPepIndices]
		lPeptides[[i]]$mRec[1:iLen]		<- dfMascot$mRec[vPepIndices]
		lPeptides[[i]]$mod[1:iLen]		<- dfMascot$mod[vPepIndices]
		lPeptides[[i]]$file[1:iLen]		<- dfMascot$file[vPepIndices]
	}
	return(lPeptides)
}


loadPf <- function(sPfFileName = "")
{
	if (sPfFileName == "")
		dfPf <- read.csv( file.choose(), head = TRUE, stringsAsFactors = FALSE)
	else
		dfPf <- read.csv( sPfFileName, head = TRUE, stringsAsFactors = FALSE)	
	
	## Order pf map data
	dfPf <- dfPf[order(dfPf$start),]
	dfPf <- dfPf[order(dfPf$chrom),]
	dfPf <- dfPf[order(suppressWarnings(as.numeric(dfPf$chrom))),]
	
	dfPf$Index <- 1:length(dfPf$Index)
	
	return(dfPf)
}

loadSnippet <- function(sSnippetName = "")
{
	if (sSnippetName == "")
		dfSnippet <- read.csv( file.choose(), head = TRUE, stringsAsFactors = FALSE)
	else
		dfSnippet <- read.csv( sSnippetName, head = TRUE, stringsAsFactors = FALSE)	
	
	## Order pf map data
	dfSnippet <- dfSnippet[order(dfSnippet$start),]
	dfSnippet <- dfSnippet[order(dfSnippet$chrom),]
	dfSnippet <- dfSnippet[order(suppressWarnings(as.numeric(dfSnippet$chrom))),]
	
	dfSnippet$Index <- 1:length(dfSnippet$Index)
	
	return(dfSnippet)
}

###############################################################################
#	prepareHumanMap()
# 
# Author: Travis
###############################################################################
##	system.time(dfHuman <- prepareHumanMap("D:\\Project Source Files\\chrHs.csv"))
loadHumanFile <- function(sHumanFilename = "")
{
	## Read Human Proteome to Genome Map file
	if (sHumanFilename == "")
		dfHuman <- read.csv( file.choose(), head = TRUE, stringsAsFactors = FALSE)
	else
		dfHuman <- read.csv( sHumanFilename, head = TRUE, stringsAsFactors = FALSE)
	
	## rename "chromosome" to "chrom" for consistency and brevity
	names(dfHuman)[names(dfHuman) == 'chromosome'] <- "chrom"
	
	## Fix chromosome name problems
	tChr <- dfHuman$chrom
	
	## make sure all x and y chromosomes are in same case
	tChr[which(tChr == "Y")] <- "y"
	tChr[which(tChr == "X")] <- "x"
	tChr [grep("CHR", tChr )] <- "UKN"
	tChr [grep("GL", tChr )] <- "UKN"
	tChr [grep("KI", tChr )] <- "UKN"
	tChr [tChr == ""] <- "UKN"
	dfHuman$chrom <- tChr	

	## Order human map data
	dfHuman <- dfHuman[order(dfHuman$start),]
	dfHuman <- dfHuman[order(dfHuman$chrom),]
	dfHuman <- dfHuman[order(suppressWarnings(as.numeric(dfHuman$chrom))),]
	## after initial ordering assign index to simplify future ordering operations
	## Recycle "X" it was the index before ordering"
	names(dfHuman)[names(dfHuman) == 'X'] <- 'Index'	## give it an informative name
	dfHuman$Index <- 1:length(dfHuman$Index)	## assign an ascending, incremental index starting at 1
	
	return(dfHuman)
}

#prepareSpeciesMap <- function(dfUnknown, sName)
#{
#	lSpecies <- list()
#	lSpecies$name <- sName
#	dfUnknown <- dfUnknown[order(dfUnknown$start),]
#	dfUnknown <- dfUnknown[order(dfUnknown$chrom),]
#	dfUnknown <- dfUnknown[order(suppressWarnings(as.numeric(dfUnknown$chrom))),]	
#	
#	lSpecies$protSeq <- paste(dfUnknown$seq, collapse = '', sep = '')
#	vPepLengths <- nchar(dfUnknown$seq)
#	
#	vIndex <- 1
#	vStartIndices <- vector(mode = "double", length = length(vPepLengths))
#	map <- vector(mode = "integer", length = nchar(lSpecies$protSeq))
#	
#	for (i in 1:length(vPepLengths))
#	{
#		vStartIndices[i] <- vIndex
#		vIndex <- vIndex + vPepLengths[i]
#	}
#	
#	for (i in 1:length(vPepLengths))
#	{	
#		endIndex <- vStartIndices[i] + vPepLengths[i] - 1
#		map[vStartIndices[i]:endIndex] <- i
#	}
#	rm(vStartIndices)
#	lSpecies$map <- map
#	lSpecies$genes <- dfUnknown$GeneName
#	lSpecies$chromes <- dfUnknown$chrom
#	lSpecies$startLoc <- dfUnknown$start
#	lSpecies$endLoc <- dfUnknown$end
#	lSpecies$lCoincidenceMaps <- list()	
#	return(lSpecies)
#}

prepareReverseSeqSpeciesMap <- function(dfUnknown, sName)
{
	lSpecies <- list()
	lSpecies$name <- sName
	dfUnknown <- dfUnknown[order(dfUnknown$start),]
	dfUnknown <- dfUnknown[order(dfUnknown$chrom),]
	dfUnknown <- dfUnknown[order(suppressWarnings(as.numeric(dfUnknown$chrom))),]	
	lGeneIndices <- vector(mode = "list", length = length(dfUnknown$GeneName))
	seq <- paste(dfUnknown$seq, collapse = "", sep = "")
	seqLen <- nchar(seq)
	seqSplit <- strsplit(seq, split = "")
	revSeq <- seqSplit[[1]][seqLen:1]
	
	lSpecies$protSeq <- paste(revSeq, collapse="")
	
	vPepLengths <- nchar(dfUnknown$seq)
	
	vIndex <- 1
	
	vStartIndices <- vector(mode = "double", length = length(vPepLengths))
	
	map <- vector(mode = "integer", length = nchar(lSpecies$protSeq))	
		
	for (i in 1:length(vPepLengths))
	{
		vStartIndices[i] <- vIndex
		vIndex <- vIndex + vPepLengths[i]
	}
	
	for (i in 1:length(vPepLengths))
	{	
		endIndex <- vStartIndices[i] + vPepLengths[i] - 1
		map[vStartIndices[i]:endIndex] <- i
		lGeneIndices[[i]] <- vStartIndices[i]:endIndex 
	}
	rm(vStartIndices)
	lSpecies$map <- map
	lSpecies$genes <- dfUnknown$GeneName
	lSpecies$lGeneIndices <- lGeneIndices
	lSpecies$chromes <- dfUnknown$chrom
	lSpecies$startLoc <- dfUnknown$start
	lSpecies$endLoc <- dfUnknown$end
	lSpecies$lCoincidenceMaps <- list()	
	return(lSpecies)
}

prepareSpeciesMap <- function(dfUnknown, sName)
{
	lSpecies <- list()
	lSpecies$name <- sName
	dfUnknown <- dfUnknown[order(dfUnknown$start),]
	dfUnknown <- dfUnknown[order(dfUnknown$chrom),]
	dfUnknown <- dfUnknown[order(suppressWarnings(as.numeric(dfUnknown$chrom))),]	
	lGeneIndices <- vector(mode = "list", length = length(dfUnknown$GeneName))
	lSpecies$protSeq <- paste(dfUnknown$seq, collapse = '', sep = '')
	vPepLengths <- nchar(dfUnknown$seq)
	
	vIndex <- 1
	vStartIndices <- vector(mode = "double", length = length(vPepLengths))
	map <- vector(mode = "integer", length = nchar(lSpecies$protSeq))
	
	for (i in 1:length(vPepLengths))
	{
		vStartIndices[i] <- vIndex
		vIndex <- vIndex + vPepLengths[i]
	}
	
	for (i in 1:length(vPepLengths))
	{	
		endIndex <- vStartIndices[i] + vPepLengths[i] - 1
		map[vStartIndices[i]:endIndex] <- i
		lGeneIndices[[i]] <- vStartIndices[i]:endIndex 
	}
	rm(vStartIndices)
	lSpecies$map <- map
	lSpecies$genes <- dfUnknown$GeneName
	lSpecies$lGeneIndices <- lGeneIndices
	lSpecies$chromes <- dfUnknown$chrom
	lSpecies$startLoc <- dfUnknown$start
	lSpecies$endLoc <- dfUnknown$end
	lSpecies$lCoincidenceMaps <- list()	
	return(lSpecies)
}

pepMapToSeq <- function(pep, count, seq, map)
{
	pepLen <- nchar(pep)
	vIndices <- gregexpr(pep, seq, fixed = TRUE)[[1]] ## fixed = TRUE for performance

	for (i in 1:pepLen)
	{
		map[vIndices] <- map[vIndices] + count
		vIndices <- vIndices + 1
	}
	return(map)
}

pepMapToSeqMascot <- function(pep, maxMascot, seq, map)
{
	pepLen <- nchar(pep)
	vIndices <- gregexpr(pep, seq, fixed = TRUE)[[1]] ## fixed = TRUE for performance
	count <- rep(maxMascot, length(vIndices))
	
	for (i in 1:pepLen)
	{
		map[vIndices] <- pmax(map[vIndices], count)
		vIndices <- vIndices + 1
	}
	return(map)
}

mapIncidence <- function (pep, seq, genotype, map)
{
	iFound <- vector(mode = "integer", length = length(pep))
	iFound <- 0
	pepLen <- nchar(names(pep))
	
	for (i in 1:length(pep))	
	{
		iFound <- length(which(pep[[i]]$sample == genotype))
		
		if (iFound > 0)
		{	
		
			vIndices <- gregexpr(names(pep)[i], seq, fixed = TRUE)[[1]] ## fixed = TRUE for performance
			
			if (vIndices[1] > 0)	# make sure valid index found by checking first returned index
			{						# make sure it is not -1 (not found) and a valid index
		
				for (j in 1:pepLen[i])
				{
					map[vIndices] <- map[vIndices] + iFound
					vIndices <- vIndices + 1
				}
			}
		}
	}
	
	return(map)
}

mapMascot <- function(pep, seq, genotype, map)
{
	maxMascot <- 0
	pepLen <- nchar(names(pep))
	for (i in 1:length(pep))	
	{
		
		if (length(which(pep[[i]]$sample == genotype)) > 0)
		{	
			vIndices <- gregexpr(names(pep)[i], seq, fixed = TRUE)[[1]] ## fixed = TRUE for performance
			
			if (vIndices[1] > 0)	# make sure valid index found by checking first returned index
			{						# make sure it is not -1 (not found) and a valid index
				maxMascot <- max(pep[[i]]$mScore)		
	
				count <- rep(maxMascot, length(vIndices))
				
				for (j in 1:pepLen[i])
				{
					map[vIndices] <- pmax(map[vIndices], count)
					vIndices <- vIndices + 1
				}
			}
		}
	}
	
	return(map)
}

mapAllPeptides <- function(pep, seq, map)
{
	for (i in 1:length(pep))	
		map <- pepMapToSeq(names(pep)[i], pep[[i]]$count, seq, map)
	
	return(map)
}
##	user  	system 	elapsed 
## 2068.84    1.20 	2070.48 
mapPeptidesPerGenotype <- function(pep, seq, genotype, mapLength)
{
	incidence <- vector(mode = "integer", length = mapLength)
	mascot <- vector(mode = "integer", length = mapLength)
	
	incidence <- mapIncidence(pep, seq, genotype, incidence)
	mascot <- mapMascot(pep, seq, genotype, mascot)

	return(map <- data.frame(incidence, mascot))
}


#   user  system elapsed 
#  2243.72  304.58 2549.42
#mapPeptidesPerGenotype <- function(pep, seq, genotype, mapLength)
#{
#	incidence <- vector(mode = "integer", length = mapLength)
#	mascot <- vector(mode = "integer", length = mapLength)
#	iFound <- 0
#	
#	for (i in 1:length(pep))	
#	{
#		iFound <- length(which(pep[[i]]$sample == genotype))
#		if (iFound > 0)
#		{
#			incidence <- pepMapToSeq(names(pep)[i], iFound, seq, incidence)
#			mascot <- pepMapToSeqMascot(names(pep)[i], max(pep[[i]]$mScore), seq, mascot)
#		}
#	}
#	return(map <- data.frame(incidence, mascot))
#}


mapGenotypes <- function(pep, lSpecies, uGenotypes)
{
	lGenotypes <- vector(mode = "list", length = length(uGenotypes))
	vInitVector <- vector(mode = "integer", length = length(lSpecies$genes))
	vInitVector <- rep.int(0, length(lSpecies$genes))
	iSeqLength <- nchar(lSpecies$protSeq)
	
	for (i in 1:length(uGenotypes))
	{
		lGenotype <- list()
		lGenotype$name <- uGenotypes[i]
		lGenotype$map <- mapPeptidesPerGenotype(pep, lSpecies$protSeq, uGenotypes[i], iSeqLength)
		lGenotype$mascotScore <- vInitVector
		lGenotype$coincidenceScore <- vInitVector
		lGenotypes[[i]]	<- lGenotype
	}
	lSpecies$lCoincidenceMaps <- lGenotypes
	return(lSpecies)
}

getIncidentPeptidesPerSnippet <- function (pep, seq, genotype)
{
	iFound <- vector(mode = "integer", length = length(pep))
	iFound <- 0
	pepLen <- nchar(names(pep))
	dfPeptides <- data.frame(seq = character(), count = integer(), stringsAsFactors = FALSE)
	
	for (i in 1:length(pep))	
	{	
		iFound <- length(which(pep[[i]]$sample == genotype))	## iFound is number of times peptide was found in
																## biological sample
		
		if (iFound > 0)
		{
			vIndices <- gregexpr(names(pep)[i], seq, fixed = TRUE)[[1]] ## fixed = TRUE for performance
																## length of vIndices is number of times that peptide was
																## found along proteomic sequence
			if (vIndices[1] > 0)	# make sure valid index found by checking first returned index
			{						# make sure it is not -1 (not found) and a valid index
				## iFound is number of that peptide to add to list
##				names(pep)[i]) iFound
				
				dfPeptides[nrow(dfPeptides) + 1,] <- c(names(pep[i]), iFound) 
			}			
		}
		
	}
	
	return(dfPeptides)
}

# debug code for getIncidentPeptidesPerSnippet
#iFound <- F
#i <- i + 1
#while (iFound < 1 || (!bStop))
#{
#	iFound <- length(which(pep[[i]]$sample == genotype))
#	if (iFound > 0)
#	{
#		vIndices <- gregexpr(names(pep)[i], seq, fixed = TRUE)[[1]]
#		bStop <- (vIndices[1] > 0)
#	}
#	i <- i + 1
#}
#i <- i - 1
#vIndices <- gregexpr(names(pep)[i], seq, fixed = TRUE)[[1]]
#vIndices
#iFound
#i

## call to generate list of peptides that map to sequence snippet
getSnippetPeptideList <- function(pep, lSpecies, uGenotypes)
{
	lGenotypes <- vector(mode = "list", length = length(uGenotypes))
	
	for (i in 1:length(uGenotypes))
	{
		lGenotype <- list()
		lGenotype$name <- uGenotypes[i]
		lGenotype$dfPeptides <- data.frame()
		lGenotype$dfPeptides <- getIncidentPeptidesPerSnippet(pep, lSpecies$protSeq, uGenotypes[i])
		lGenotypes[[i]]	<- lGenotype
	}
	
	return(lGenotypes)
}



## Timing
##   user  system elapsed 
## 478.11  241.35  719.81 
#calculateScoresPerGene <- function(lSpecies, bIncludeMascot = FALSE)
#{
#	numMaps <- length(lSpecies$lCoincidenceMaps)
#	
#	if (bIncludeMascot)
#	{
#		for (iGene in 1:length(lSpecies$genes))
#		{
#			indices <- which(lSpecies$map == iGene)
#			
#			for (iMap in 1:numMaps)
#			{
#				lSpecies$lCoincidenceMaps[[iMap]][["coincidenceScore"]][iGene] <- mean(lSpecies$lCoincidenceMaps[[iMap]][["map"]]$incidence[indices])
#				lSpecies$lCoincidenceMaps[[iMap]][["mascotScore"]][iGene] <-max(lSpecies$lCoincidenceMaps[[iMap]][["map"]]$mascot[indices])
#			}
#		}
#	}
#	else
#	{
#		for (iGene in 1:length(lSpecies$genes))
#		{
#			indices <- which(lSpecies$map == iGene)
#			
#			for (iMap in 1:numMaps)
#			{
#				lSpecies$lCoincidenceMaps[[iMap]][["coincidenceScore"]][iGene] <- mean(lSpecies$lCoincidenceMaps[[iMap]][["map"]]$incidence[indices])
#			}
#		}		
#	}
#	return(lSpecies)
#}

##   user  system elapsed 
##  28.03    0.05   28.08 
calculateScoresPerGene <- function(lSpecies, bIncludeMascot = FALSE)
{
	numMaps 	<- length(lSpecies$lCoincidenceMaps)
	numGenes 	<- length(lSpecies$genes)

	for (iMap in 1:numMaps)
	{
		for (iGene in 1:numGenes)
		{
			lSpecies$lCoincidenceMaps[[iMap]][["coincidenceScore"]][iGene] <- mean(lSpecies$lCoincidenceMaps[[iMap]][["map"]]$incidence[lSpecies$lGeneIndices[[iGene]]])
			lSpecies$lCoincidenceMaps[[iMap]][["mascotScore"]][iGene] <- max(lSpecies$lCoincidenceMaps[[iMap]][["map"]]$mascot[lSpecies$lGeneIndices[[iGene]]])
		}
	}
	
	return(lSpecies)
}

#calculateScores <- function(lSpecies, bIncludeMascot = FALSE)
#{
#	numMaps <- length(lSpecies$lCoincidenceMaps)
#	
#	for (iGene in 1:length(lSpecies$genes))
#	{
#		indices <- which(lSpecies$map == iGene)
#		offset <- length(indices)
#		cScores <- vector(mode = "double", length = (offset * numMaps))
#		mScores <- vector(mode = "double", length = (offset * numMaps))
#		iStart <- 1
#		
#		for (iMap in 1:numMaps)
#		{
#			cScores[iStart:(iStart+offset-1)] <- lSpecies$lCoincidenceMaps[[iMap]][["map"]]$incidence[indices]
#			mScores[iStart:(iStart+offset-1)] <- lSpecies$lCoincidenceMaps[[iMap]][["map"]]$mascot[indices]
#			iStart <- iStart + offset
#		}
#		lSpecies$lCoincidenceMaps[[1]][["coincidenceScore"]][iGene] <- mean(cScores[1:(iStart-1)])
#		lSpecies$lCoincidenceMaps[[1]][["mascotScore"]][iGene] <- max(mScores[1:(iStart-1)])
#	}
#	return(lSpecies)
#}

plotGraphByGenotypePerChrom <- function(lSpecies, genotype, bMascot = FALSE)
{
	i <- 1
	iMax <- length(lSpecies$lCoincidenceMaps)
	sName <- lSpecies$lCoincidenceMaps[[i]]$name
	
	while ((i <= iMax) && (!is.null(sName)))
	{
		if (lSpecies$lCoincidenceMaps[[i]]$name == genotype)
		{
			if (bMascot)
				return(makePlotPerChrom(lSpecies$lCoincidenceMaps[[i]]$mascotScore, lSpecies$genes, lSpecies$chromes, "Max Mascot Score per Gene", genotype))
			else
				return(makePlotPerChrom(lSpecies$lCoincidenceMaps[[i]]$coincidenceScore, lSpecies$genes, lSpecies$chromes, "Mean Coincidence Score per Gene", genotype))
			
			break
		}
		i <- i + 1	
	}
}
plotPeptideAlignment <- function(lSpecies, dfPeps)
{
	seq <- lSpecies$protSeq
	
	totPeps <- length(dfPeps$count)	
	xLab <- strsplit(seq, "", fixed=TRUE)
	iStartOrder <- NULL
	
	for(i in 1: totPeps)
	{	
		iStartOrder[i] <- gregexpr(dfPeps$seq[i], seq, fixed = TRUE)[[1]][1]
	}

	dfPeps <- dfPeps[order(iStartOrder),]

	seqLengths <- nchar(dfPeps$seq)
	lineW <- 3
	index <- 1
	yLim = c(0,(sum(dfPeps$count)))
	
	plot(1:nchar(seq), lSpecies$lCoincidenceMaps[[1]]$map$incidence, type="l", xlab = "Residue", ylab = "Peptide Frequency", xaxt = "n", ylim = yLim, lwd=lineW)
	

	
	for(i in 1: totPeps)
	{
		xStart <- gregexpr(dfPeps$seq[i], seq, fixed = TRUE)[[1]][1]
		xEnd <- xStart + seqLengths[i] - 1
		for (j in 1:dfPeps$count[i])
		{
			#segments(xStart, index, (xStart + seqLengths[i] - 1), index, lty = "solid", col = "blue", lwd = lineW)
			text(x = xStart:xEnd, y = index, labels = unlist(strsplit(dfPeps$seq[i], "")), col = "black", cex = 0.6)
			index <- index + 1
		}
#		segments(xStart, i, (xStart + seqLengths[i] - 1), i, lty = "solid", col = dfPeps$count[i], lwd = lineW)
	}
#	lines(1:nchar(seq), lSpecies$lCoincidenceMaps[[1]]$map$incidence, type="l", col="black", lwd = lineW)
	axis(side = 1, at = 1:nchar(seq), labels = unlist(xLab), tick = FALSE)
#	axis(side = 2, at =, labels=lSpecies$lCoincidence[[1]]$map$incidence)
#	freq <- as.character(c(1:max(dfPeps$count)))
#	legend("topright", title = "Peptide Frequency", legend = freq, lty=c(1,1), bty = "n", lwd = lineW, col=1:max(dfPeps$count))
	return(dfPeps)
}

plotMultipleGenotypesPerChrom <- function(lSpecies, genotypes, bMascot = FALSE, bLabelGenes = TRUE, yMax = NULL)
{
	i <- 1
	iFound <- 1
	iTargets <- length(genotypes)
	iMax <- length(lSpecies$lCoincidenceMaps)
	sName <- lSpecies$lCoincidenceMaps[[i]]$name
	vIndices <- vector(mode = "integer", length = iTargets) 
	
	while ((i <= iMax) && (!is.null(sName)) && (iFound <= iTargets))
	{
		if (lSpecies$lCoincidenceMaps[[i]]$name %in% genotypes)
		{
			vIndices[iFound] <- i
			iFound <- iFound + 1
		}
		i <- i + 1	
	}
	if (iFound > 1)
	{
		if (bMascot)
			return(makeMultiPlotPerChrom(lSpecies, vIndices, "Max Mascot Score per Gene", genotypes, bMascot, bLabelGenes, yMax))
		else
			return(makeMultiPlotPerChrom(lSpecies, vIndices, "Mean Coincidence Score per Gene", genotypes, bMascot, bLabelGenes, yMax))
		
	}
}

getYLim <- function(lSpecies, vIndices, bMascot)
{
	vMaxY <- vector(mode = "double", length = length(vIndices))
	vMinY <- vector(mode = "double", length = length(vIndices))
	
	if (bMascot)
		for (i in 1:length(vIndices))
		{
			vMaxY[i] <- max(lSpecies$lCoincidenceMaps[[vIndices[i]]]$mascotScore)
			vMinY[i] <- min(lSpecies$lCoincidenceMaps[[vIndices[i]]]$mascotScore)
		}
	else
		for (i in 1:length(vIndices))
		{
			vMaxY[i] <- max(lSpecies$lCoincidenceMaps[[vIndices[i]]]$coincidenceScore)
			vMinY[i] <- min(lSpecies$lCoincidenceMaps[[vIndices[i]]]$coincidenceScore)
		}
	
	## of the Y values of interest, the maximum is held in iMaxY (for scaling Y axis
	iMaxY <- max(vMaxY)
	iMaxY <- iMaxY + iMaxY * 0.05 ## add 5% padding to top of Y axis
	iMinY <- min(vMaxY, 0)
	
	return(c(iMinY, iMaxY))	
}

makeMultiPlotPerChrom <- function(lSpecies, vIndices, yLab, sTitle, bMascot = FALSE, bLabelGenes = TRUE, yMax = NULL)
{
	genes <- lSpecies$genes
	chromes <- lSpecies$chromes
	chromLabels <- unique(chromes)
	lChromIndex  <- vector(mode = "integer", length = length(chromLabels))
	uChromIndex <- vector(mode = "integer", length = length(chromLabels))
	chromLabelPos <- vector(mode = "integer", length = length(chromLabels))
	margins <- c(5,5,4,2) # bottom, left, top, right // default: 5,4,4,2
	sLegendLabels <- NULL
	vYLim <- NULL

	for (i in 1:length(chromLabels))
	{
		lChromIndex[i] <- min(which(chromLabels[i] == chromes))
		uChromIndex[i] <- max(which(chromLabels[i] == chromes))
	}

	chromLabelPos <- ((lChromIndex + uChromIndex) / 2) # position labels at midpoint of chromosome limits

	if (is.null(yMax))
	{
		vYlim <- getYLim(lSpecies, vIndices, bMascot)
	}
	else
	{
		vYlim <-c(0, yMax)
	}
	
	for (i in 1:length(vIndices))
	{
		if (bMascot)
			yValues <- lSpecies$lCoincidenceMaps[[vIndices[i]]]$mascotScore
		else
			yValues <- lSpecies$lCoincidenceMaps[[vIndices[i]]]$coincidenceScore
		
		iMean <- mean(yValues)
		iSd <- sd(yValues)
		geneIndices <- which(yValues > iMean + (4 * iSd))
		geneLabels <- genes[geneIndices]		
		
		if(i == 1)
		{
			par(mar = margins)
			par(bty="n")	
			par(cex.lab=1.5)
			par(cex.axis=1.5)
			
			plot(1:length(yValues), yValues, xlab = "Chromosome", xaxt = "n", ylab = yLab, type="l", ylim = vYlim)
			drawChromeBounds(uChromIndex)
			axis(side = 1, at = chromLabelPos, labels = chromLabels, tick = FALSE)
			
			if (bLabelGenes)
				text(x = geneIndices, y = yValues[geneIndices], geneLabels, pos=3, col=i)
		}
		else
		{
			lines(1:length(yValues), yValues, type="l", col= i)
			
			if (bLabelGenes)
				text(x = geneIndices, y = yValues[geneIndices], geneLabels, pos=3, col=i)
		}
		
		if(is.null(sLegendLabels))
			sLegendLabels <- lSpecies$lCoincidenceMaps[[vIndices[i]]]$name
		
		else
			sLegendLabels <- c(sLegendLabels, lSpecies$lCoincidenceMaps[[vIndices[i]]]$name)
	}
	
	orderedIndices <- order(sLegendLabels)
	legend("topright", sLegendLabels[orderedIndices], lty=c(1,1), col=orderedIndices)	
	return(vYlim)
}

drawChromeBounds <- function(xValues)
{
	yBottomSegment <- par("usr")[3]
	iHeight <- abs((par("usr")[4] - par("usr")[3]) * 0.01)
	yTopSegment <- yBottomSegment + iHeight
	
	yHorizPos <- yTopSegment - (iHeight/2)
	
	xValues <- c(as.integer("0"), xValues)
	segments(xValues, yTopSegment, xValues, yBottomSegment, lty = "solid", col = "black", lwd=2)
	segments(xValues, yHorizPos, xValues[length(xValues)], yHorizPos, lty = "solid", col = "black", lwd=3)	
}

plotMultipleGenotypesPerGene <- function(lSpecies, genes, genotypes, bMascot = FALSE, bLabelGenes = TRUE)
{
	i <- 1
	iFound <- 1
	iTargets <- length(genotypes)
	iMax <- length(lSpecies$lCoincidenceMaps)
	sName <- lSpecies$lCoincidenceMaps[[i]]$name
	vIndices <- vector(mode = "integer", length = iTargets) 
	
	while ((i <= iMax) && (!is.null(sName)) && (iFound <= iTargets))
	{
		if (lSpecies$lCoincidenceMaps[[i]]$name %in% genotypes)
		{
			vIndices[iFound] <- i
			iFound <- iFound + 1
		}
		i <- i + 1	
	}
	if (iFound > 1)
	{
		if (bMascot)
			return(makeMultiPlotPerGene(lSpecies, genes, vIndices, "Mascot Score per Residue", genotypes, bMascot, bLabelGenes))
		else
			return(makeMultiPlotPerGene(lSpecies, genes, vIndices, "Coincidence Score per Residue", genotypes, bMascot, bLabelGenes))
		
	}
}


getGeneStart <- function(lSpecies, gene)
{
	geneIndex <- which(lSpecies$genes == gene)
	return(min(which(lSpecies$map == geneIndex)))
}

getGeneEnd <- function(lSpecies, gene)
{
	geneIndex <- which(lSpecies$genes == gene)
	return(max(which(lSpecies$map == geneIndex)))
}

drawGeneBounds <- function(xValues, iHeight)
{
	yTopSegment <- (par("usr")[4] - par("usr")[3]) * 0.01
	yBottomSegment <- yTopSegment - iHeight
	yHorizPos <- yTopSegment - (iHeight/2)
	
	xValues[1] <- xValues[1] - 1
	xValues[length(xValues)] <- xValues[length(xValues)] + 1
	
	segments(xValues, yTopSegment, xValues, yBottomSegment, lty = "solid", col = "black", lwd=2)
	segments(xValues, yHorizPos, xValues[length(xValues)], yHorizPos, lty = "solid", col = "black", lwd=3)
	return(yHorizPos)
}

savePlotByGenotypePerGene <- function(lSpecies, genes, genotypes, bIncludeMascot = FALSE, iWidth, sName)
{
	sPath <- paste("D:\\Project Source Files\\Malaria_Genotypes_Mascot\\", sName, ".pdf", collapse = '', sep = '')
	pdf(sPath, width=iWidth, height=(0.5*iWidth))
	plotMultipleGenotypesPerGene(lSpecies, genes, genotypes, bIncludeMascot)
	dev.off()
}



makeMultiPlotPerGene <- function(lSpecies, genes, indices, yLab, sTitle, bMascot = FALSE)
{

	vGenesStart <- vector(mode = "integer", length = length(genes))
	vGenesEnd <- vector(mode = "integer", length = length(genes))
	geneLabelPos <- vector(mode = "integer", length = length(genes))
	margins <- c(5,5,4,2) # bottom, left, top, right // default: 5,4,4,2
	vYLim <- vector(mode = "double", length = length(indices))
	vMaxY <- vector(mode = "double", length = length(indices))
	vMinY <- vector(mode = "double", length = length(indices))
	lGenotypeYVals <- vector(mode = "list", length = length(indices))
	vGeneAAVals <- vector(mode = "list", length = length(indices))
	xTicks <- vector(mode = "double", length = length(genes) + 1)
	xTicks[1] <- 0
	temp <- 0
	newY <- NULL
	vAASeq <- NULL
	yOffset <- 0
#	vYLim <- NULL
	
	sLegendLabels <- NULL
	
	## this loop assembles the vectors of start and end locations of genes for plotting on axis
	for (i in 1:length(genes))
	{
		vGenesStart[i] <- getGeneStart(lSpecies, genes[i])
		vGenesEnd[i] <- getGeneEnd(lSpecies, genes[i])
		
		if (i == 1)
			vAASeq <- substr(lSpecies$protSeq, vGenesStart[i], vGenesEnd[i])
		
		else
			vAASeq <- paste(vAASeq, substr(lSpecies$protSeq, vGenesStart[i], vGenesEnd[i]), sep="")
	}
	
	if (bMascot)
		for (i in 1:length(indices))
		{
			vMaxY[indices[i]] <- max(lSpecies$lCoincidenceMaps[[indices[i]]]$map$mascot)
			vMinY[indices[i]] <- min(lSpecies$lCoincidenceMaps[[indices[i]]]$map$mascot)
		}
	else
		for (i in 1:length(indices))
		{
			vMaxY[indices[i]] <- max(lSpecies$lCoincidenceMaps[[indices[i]]]$map$incidence)
			vMinY[indices[i]] <- min(lSpecies$lCoincidenceMaps[[indices[i]]]$map$incidence)
		}
	
	iMaxY <- max(vMaxY) ## of the Y values of interest, the maximum is held in iMaxY (for scaling Y axis
	iMinY <- min(vMaxY)
	
	for (i in 1:length(indices))
	{
		for (j in 1:length(genes))
		{
			if (bMascot)
				newY <- lSpecies$lCoincidenceMaps[[indices[i]]]$map$mascot[vGenesStart[j]:vGenesEnd[j]]
			else
				newY <- lSpecies$lCoincidenceMaps[[indices[i]]]$map$incidence[vGenesStart[j]:vGenesEnd[j]]
				
			if (j == 1)
			{
				yValues <- newY
				geneLabelPos[j] <- length(newY) / 2
				xTicks[j+1] <- length(newY)
			}
			else
			{
				geneLabelPos[j] <- length(yValues) + (length(newY) / 2)
				yValues <- c(yValues, newY)
				xTicks[j+1] <- length(yValues)
			}
		}
		vYLim[i] <- max(yValues)
		lGenotypeYVals[[i]] <- yValues
		
		if(is.null(sLegendLabels))
			sLegendLabels <- lSpecies$lCoincidenceMaps[[indices[i]]]$name
		
		else
			sLegendLabels <- c(sLegendLabels, lSpecies$lCoincidenceMaps[[indices[i]]]$name)			
	}
	
	for (i in 1:length(indices))
	{
		if(i == 1)
		{
			par(mar = margins)
			par(bty="n")	
			par(cex.lab=1.5)
			par(cex.axis=1.5)
			
			plot(1:length(lGenotypeYVals[[i]]), lGenotypeYVals[[i]], ylim=c(as.integer("0"), max(vYLim)), xlab = "Residues", xaxt = "n", ylab = yLab, type="l")
			temp <- drawGeneBounds(xTicks, (0.02*max(vYLim)))
			yOffset <- (par("usr")[4] - par("usr")[3]) * 0.02
			axis(side = 1, at = geneLabelPos, labels = genes, tick = FALSE)
			text(x = 1:nchar(vAASeq), y = -yOffset, unlist(strsplit(vAASeq, "", fixed = TRUE)), col="black", cex=0.5)
		}
		else
		{
			lines(1:length(lGenotypeYVals[[i]]), lGenotypeYVals[[i]], type="l", col= i)
		}
	}

	orderedIndices <- order(sLegendLabels)
	legend("topright", sLegendLabels[orderedIndices], lty=c(1,1), col=orderedIndices)	
	return(vAASeq)
}


makePlotPerChrom <- function(yValues, genes, chromes, yLab, sTitle)
{
	chromLabels <- unique(chromes)
	lChromIndex  <- vector(mode = "integer", length = length(chromLabels))
	uChromIndex <- vector(mode = "integer", length = length(chromLabels))
	chromLabelPos <- vector(mode = "integer", length = length(chromLabels))
	margins <- c(5,5,4,2) # bottom, left, top, right // default: 5,4,4,2
	iMean <- mean(yValues)
	iSd <- sd(yValues)
	geneIndices <- which(yValues > (iMean + 4 * iSd))
	geneLabels <- genes[geneIndices]
	
	for (i in 1:length(chromLabels))
	{
		lChromIndex[i] <- min(which(chromLabels[i] == chromes))
		uChromIndex[i] <- max(which(chromLabels[i] == chromes))
	}
	
	chromLabelPos <- ((lChromIndex + uChromIndex) / 2) # position labels at midpoint of chromosome limits
	par(mar = margins)
	plot(1:length(yValues), yValues, main = sTitle, xlab = "Chromosome", xaxt = "n", ylab = yLab, type="l", cex.lab = 1.5)
	drawChromeBounds(uChromIndex)
	axis(side = 1, at = chromLabelPos, labels = chromLabels, tick = FALSE)
	#axis(side = 1, at = geneIndices, labels = geneLabels, las = 0, tck = 0)
	text(x = geneIndices, y = yValues[geneIndices], geneLabels, pos=3)
	return(geneLabels)
}

getSN <- function(lSpecies, chromes, yValues)
{
	indices <- which(lSpecies$chromes %in% chromes)
	iLen <- length(indices)
	vValues <- vector(mode = "double", length = iLen)
	
	vValues <- yValues[indices]
	
	nom <- max(yValues)
	
	denom <- sqrt(sum(vValues^2)/iLen) 
	
	return(nom/denom)
}

savePlotByGenotypePerChrom <- function(lSpecies, genotypes, bIncludeMascot = FALSE, iWidth)
{
	for (i in 1:length(genotypes))
	{
		sPath <- paste("D:\\Project Source Files\\2016_08_30_Figures\\drugNoDrugStack", ".pdf", collapse = '', sep = '')
		pdf(sPath, width=iWidth, height=(0.5*iWidth))
		plotGraphByGenotypePerChrom(lSpecies, genotypes[i], bIncludeMascot)
		dev.off()
	}
}



mergeGenotypes <- function(lSpecies, genotypes)
{
	i <- 1
	iFound <- 1
	iTargets <- length(genotypes)
	iMax <- length(lSpecies$lCoincidenceMaps)
	sName <- lSpecies$lCoincidenceMaps[[i]]$name
	vIndices <- vector(mode = "integer", length = iTargets)
	gTypes <- list()
	sMergeName <- NULL
	
	lSpecies2 <- list()
	lSpecies2$name 		<- lSpecies$name
	lSpecies2$protSeq 	<- lSpecies$protSeq
	lSpecies2$map 		<- lSpecies$map
	lSpecies2$genes 	<- lSpecies$genes
	lSpecies2$lGeneIndices <- lSpecies$lGeneIndices
	lSpecies2$chromes	<- lSpecies$chromes
	lSpecies2$startLoc 	<- lSpecies$startLoc
	lSpecies2$endLoc 	<- lSpecies$endLoc
	lSpecies2$lCoincidenceMaps <- list()
	
	while ((i <= iMax) && (!is.null(sName)) && (iFound <= iTargets))
	{
		if (lSpecies$lCoincidenceMaps[[i]]$name %in% genotypes)
		{
			vIndices[iFound] <- i
			iFound <- iFound + 1
			
			if(is.null(sMergeName))
				sMergeName <- lSpecies$lCoincidenceMaps[[i]]$name
			else
				sMergeName <- paste(sMergeName, lSpecies$lCoincidenceMaps[[i]]$name, collapse = '', sep = '')
		}
		i <- i + 1	
	}
	
	gTypes$name <- sMergeName
	gTypes$coincidenceScore		<- vector(mode = "double", length = length(lSpecies$lCoincidenceMaps[[1]]$coincidenceScore))
	gTypes$mascotScore			<- vector(mode = "double", length = length(lSpecies$lCoincidenceMaps[[1]]$mascotScore))
	gTypes$map$mascot			<- vector(mode = "double", length = length(lSpecies$lCoincidenceMaps[[1]]$map$mascot))
	gTypes$map$mascot			<- rep(0, length(gTypes$map$mascot))
	
	gTypes$map$incidence		<- vector(mode = "double", length = length(lSpecies$lCoincidenceMaps[[1]]$map$incidence))
	gTypes$map$incidence		<- rep(0, length(gTypes$map$incidence))
	
	for (i in 1:length(vIndices))
	{
		gTypes$map$incidence		<- gTypes$map$incidence + lSpecies$lCoincidenceMaps[[vIndices[i]]]$map$incidence
		gTypes$map$mascot			<- pmax(gTypes$map$mascot, lSpecies$lCoincidenceMaps[[vIndices[i]]]$map$mascot)
	}
	
	lSpecies2$lCoincidenceMaps[[1]] <- gTypes
	return(calculateScoresPerGene(lSpecies2, TRUE))
	
}

doSnippetStuff <- function(lSpecies, pep)
{
	lPepList <- getSnippetPeptideList(pep, lSpecies, c("C2_NoDrug", "C4_NoDrug", "C6_NoDrug"))
	C2 <- lPepList[[1]]$dfPeptides
	C4 <- lPepList[[2]]$dfPeptides
	C6 <- lPepList[[3]]$dfPeptides
	
	dfMerge <- merge(C2, C4, by = "seq", all = TRUE)
	dfMerge <- merge(dfMerge, C6, by = "seq", all = TRUE)
	dfMerge$count.x <- as.numeric(dfMerge$count.x)
	dfMerge$count.y <- as.numeric(dfMerge$count.y)
	dfMerge$count <- as.numeric(dfMerge$count)
	dfMerge <- transform(dfMerge, count = rowSums(dfMerge[,2:4], na.rm = TRUE))
	dfMerge <- dfMerge[c("seq", "count")]
	
	return(plotPeptideAlignment(lSpecies, dfMerge))

}

## sPath <- paste("D:\\Project Source Files\\figures\\", "C2C4C6ChromMasc", ".pdf", collapse = '', sep = '')
## pdf(sPath, width=iWidth, height=(0.5*iWidth))
## plotMultipleGenotypesPerChrom(lHuman, c("C2", "C4", "C6"), TRUE)
## dev.off()

## sPath <- paste("D:\\Project Source Files\\figures\\", "C2C4C6Masc", ".pdf", collapse = '', sep = '')
## pdf(sPath, width=iWidth, height=(0.5*iWidth))
## plotMultipleGenotypesPerGene(lHuman, c("HBB", "HBA1"), c("C2", "C4", "C6"), TRUE)
## dev.off()




## timing...
##  73s to read in MascotMap_Pf_April20_2016_2	~ 0.95 GB file		-> 1.2 million peptides
## 109s to read in MascotMap_Human_April20_2016_2 ~ 1.78 GB file 
## 182s to read in MascotMap_Combined_April20_2016_2 ~ 2.74 GB file   -> 3.5 million peptides
	
#format(object.size(lPep), units = "auto")
#"D:\\Project Source Files\\mascot_map.csv"
#"D:\\Project Source Files\\mascot_plot.pdf"
dfMascot <- loadMascot("D:\\Project Source Files\\mascot_map.csv")
system.time(lPep <- preparePeptideList(dfMascot))
uGenotypes <- unique(dfMascot$sample)

dfHuman <- loadHumanFile("D:\\Project Source Files\\chrHs.csv")
lHuman <- prepareSpeciesMap(dfHuman, "Human")
system.time ( lHuman <- mapGenotypes(lPep, lHuman, uGenotypes) )
system.time ( lHuman <- calculateScoresPerGene(lHuman, T) )

dfPlasmodium <- loadPf("D:\\Project Source Files\\chrPf.csv")
lPf <- prepareSpeciesMap(dfPlasmodium, "Pf")
system.time ( lPf <- mapGenotypes(lPep, lPf, uGenotypes) )
system.time ( lPf <- calculateScoresPerGene(lPf, T) )
#sGenotypes <- c("C2_NoDrug","C4_NoDrug","C6_NoDrug")

dfSnippet <- loadSnippet("D:\\Project Source Files\\chrSnippet.csv")
lSnippet <- prepareSpeciesMap(dfSnippet, "HBBSnip")
lSnippet <- mapGenotypes(lPep, lSnippet, uGenotypes)
lSnippet <- calculateScoresPerGene(lSnippet, T)
 
#plotMultipleGenotypesPerChrom(lHuman, sGenotypes)
#plotMultipleGenotypesPerGene(lHuman, c("HBB", "HBA1"), c("C2", "C4", "C6"), FALSE)
#savePlotByGenotypePerGene(lHuman, c("HBA1", "HBB"), c("C2", "C4", "C6"), FALSE, 20, "zoomNew")

################################
### Merged Data for C2, C4, C6 Across entire Human Genome
#plotMultipleGenotypesPerChrom(lHuman2, c("C2C6C4"), T)

### Data Per Genotoype for C2, C4, C6 Across entire Human Genome
#plotMultipleGenotypesPerChrom(lHuman, c("C2", "C6", "C4"), T)


### Merged Data Data Per Genotoype for C2, C4, C6 For HBB and HBA1 
#plotMultipleGenotypesPerGene(lHuman2, c("HBB", "HBA1"), c("C2C6C4"), FALSE)
#plotMultipleGenotypesPerGene(lHuman, c("HBB", "HBA1"), c("C2", "C6", "C4"), FALSE)
#plotMultipleGenotypesPerGene(lHuman, c("HBB", "HBA1"), uGenotypes, FALSE)

# show old slide

#plotMultipleGenotypesPerGene(lHuman2, c("CA1"), c("C2C6C4"), FALSE)
#plotMultipleGenotypesPerGene(lHuman, c("CA1"), uGenotypes, FALSE)

## after mapGenotypes
#> length(which(lHuman$lCoincidenceMaps[[1]]$map > 0))
#[1] 162560
#> length(which(lHuman$lCoincidenceMaps[[1]]$map$incidence > 0))
#[1] 81280
#> length(which(lHuman$lCoincidenceMaps[[1]]$map$mascot > 0))
#[1] 81280
#		sPath <- paste("D:\\Project Source Files\\2016_08_30_Figures\\C2C4C6GeneCoinc", ".pdf", collapse = '', sep = '')
#		pdf(sPath, width=20, height=10)
#plotMultipleGenotypesPerGene(lHuman, c("HBB", "HBA1"), c("C6_NoDrug", "C2_NoDrug", "C4_NoDrug"), F)
# stacked plot
#plotMultipleGenotypesPerGene(lHuman, c("HBB", "HBA1"), c("C6_NoDrug", "C2_NoDrug", "C4_NoDrug", "C6_Drug", "C2_Drug", "C4_Drug"), F)
#getSN(lHuman2, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), lHuman2$lCoincidenceMaps[[9]]$mascotScore)

##  iHBB <- which(lHumanMerge$genes == "HBB")
##  vHBBIndices <- lHumanMerge$lGeneIndices[iHBB]
##  vSeq <- unlist(strsplit(lHumanMerge$protSeq, "", fixed = TRUE))[unlist(vHBBIndices)]
##  vScore <- lHumanMerge$lCoincidenceMaps[[1]]$map$incidence[unlist(vHBBIndices)]

## snippet code
##	 lSnippetMerge <- mergeGenotypes(lSnippet, c("C2_NoDrug", "C4_NoDrug", "C6_NoDrug"))
##	 temp <- getSnippetPeptideList(lPep, lSnippetMerge, uGenotypes)
##	 temp <- getSnippetPeptideList(lPep, lSnippetMerge, c("C2_NoDrug", "C4_NoDrug", "C6_NoDrug"))
##	C2 <- temp[[1]]$dfPeptides
##	C4 <- temp[[2]]$dfPeptides
##	C6 <- temp[[3]]$dfPeptides
##	tmp <- merge(C2, C4, by = "seq", all = TRUE)
##	tmp <- merge(tmp, C6, by = "seq", all = TRUE)
##	test$count.x <- as.numeric(test$count.x)
##	test$count.y <- as.numeric(test$count.y)
##	test$count <- as.numeric(test$count)
##	tmp1 <- transform(tmp, count = rowSums(test[,1:3], na.rm = TRUE))
##	tmp2 <- tmp1[,1]
##	tmp2 <- cbind(tmp2, tmp1[,4])