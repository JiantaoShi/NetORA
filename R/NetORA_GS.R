#' Network-based enrichment for a list of signatures
#' 
#' This function is wraper of core function NetORA to work on a list of gene sets
#' 
#' @param sgList a vector of gene signature
#' @param queryTable a data frame returned by NetORA_Pre
#' @param PG a vector of background genes
#' @param FDR an FDR cutoff to define gene that are significantly connected a pathway (defaul 0.05)
#' 
#' @author Jiantao Shi
#' @references
#' 
#' Jiantao Shi: NetORA, a package for network-based pathway over-representation analysis.
#' 
#' @return
#' a list
#' @examples
#' library("NetGPA")
#' data(text_2006_12_NetGPA)
#' data(Example_NetGPA)
#' ExE_Hyper  <- Example_NetGPA$ExE_Hyper
#' queryTable <- NetORA_Pre(ExE_Hyper, text_2006_12_NetGPA, progressBar = TRUE)
#' Cancer_GS  <- Example_NetGPA$Cancer_GeneSet
#' PG <- colnames(text_2006_12_NetGPA)
#' mergedT <- NetORA_GS(Cancer_GS, queryTable, PG, FDR = 0.05)
#' @export

NetORA_GS <- function(sgList, queryTable, PG, FDR = 0.05) {

	allGene <- rownames(queryTable)
	PG      <- intersect(PG, allGene)
	netSR   <- allGene[queryTable$queryFDR < FDR]

	Ns <- length(sgList)
	rM <- matrix(NA, Ns, 5)
	
	for(i in 1:Ns){

		SG  <-sgList[[i]]
		rL  <- NetORA_pHyper(SG, netSR, PG)
		rM[i, 1] = rL$nSG
		rM[i, 2] = rL$nSR
		rM[i, 3] = rL$nOverlap
		rM[i, 4] = rL$ES
		rM[i, 5] = rL$pvalue
	}
	sgID   = names(sgList)
	colnames(rM) = c("nSG", "nSR", "nOverlap", "ES", "pvalue")

	mTable = data.frame(sgID, rM) 
	return(mTable)
}


