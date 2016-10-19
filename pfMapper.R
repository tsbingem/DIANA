## Adapted from IALewis Aug 30, 2016
readEuPath <- function(inPath = file.choose() )
{
	inList <- scan(inPath, sep="\n", what = 'A', allowEscapes = F,  
			strip.white = T, blank.lines.skip = T)
	inList <- paste(inList, collapse = "")
	inList <- strsplit(inList, ">")[[1]][-1]
	inList <- strsplit(inList, " | ", fixed = T)
	return(inList)
}


# TODO: Add comment
# 
# Author: Travis
###############################################################################
generatePfMap <- function()
{
	pf <- readEuPath()
	
	lsDtPf <- NULL
	
	for (i in 1: length(pf))
	{
		vRow <- NULL
		geneName <- pf[[i]][1]
		sSeq <- getSequence(pf[[i]][7])
		loc <- pf[[i]][4]
		
		if (isMitochondrial(pf[[i]][6]))
		{
			sChrom <- "M"
		}
		else
		{
			sChrom <- getChrom(loc)
		}
		
		iStart <-  getStart(loc)
		iEnd <- getEnd(loc)
		vRow <- c(	Index = i, GeneName = geneName, seq = sSeq, 
					chrom = sChrom, start = iStart, end = iEnd)   

		lsDtPf <- rbind(lsDtPf, vRow)
	}
	return(lsDtPf)
}
#write.csv(temp, file="D:\\Project Source Files\\chrPf.csv", row.names = F)
getSequence <- function(sSequence)
{
	return(strsplit(sSequence, "coding")[[1]][2])	
}

isMitochondrial <- function(sString)
{
	return(isTRUE(grep("mitochondrial", sString)>0))
}

getChrom <- function(sString)
{
	return(strsplit(sString, "_")[[1]][2])
}

getStart <- function(sString)
{
	temp <- strsplit(sString, "-")[[1]][1]
	return(strsplit(temp, ":")[[1]][2])
}

getEnd <- function(sString)
{
	temp <- strsplit(sString, "-")[[1]][2]
	return(strsplit(temp, "\\(")[[1]][1])
}