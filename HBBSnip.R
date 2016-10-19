# TODO: Add comment
# 
# Author: Travis
###############################################################################
## snippet of aas of interest
## sAAs <- "LHCDKLHVDPENFRLLGNVLVCV"

generateSnippet <- function(sSnippet)
{
	lsDtPf <- NULL
	sChrom <- "unknown"
	sSeq <- sSnippet
	iStart <- 1
	iEnd <- nchar(sSeq)
	sGeneName <- "sub_gene"
	iIndex <- 1
	vRow <- c(	Index = iIndex, GeneName = sGeneName, seq = sSeq, 
			chrom = sChrom, start = iStart, end = iEnd)
	
	lsDtPf <- rbind(lsDtPf, vRow)
	write.csv(lsDtPf, file="D:\\Project Source Files\\chrSnippet.csv", row.names = F)
}
