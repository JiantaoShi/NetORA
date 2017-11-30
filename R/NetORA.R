#' The core function to score a signature
#' 
#' This function is the core function to test whether a signature is statistically enriched in a pathway in network-level.
#' 
#' @param SG a vector of gene signature
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
#' library("NetORA")
#' data(text_2006_12_NetGPA)
#' data(Example_NetGPA)
#' ExE_Hyper  <- Example_NetGPA$ExE_Hyper
#' queryTable <- NetORA_Pre(ExE_Hyper, text_2006_12_NetGPA, progressBar = TRUE)
#' Cancer_GS  <- Example_NetGPA$Cancer_GeneSet
#' PG <- colnames(text_2006_12_NetGPA)
#' mergedT <- NetORA_GS(Cancer_GS, queryTable, PG, FDR = 0.05)
#' @export

NetORA <- function(SG, queryTable, PG, FDR = 0.05) {

	allGene <- rownames(queryTable)
	PG <- intersect(PG, allGene)

	netSR  <- allGene[queryTable$queryFDR < FDR]
	rL     <- NetORA_pHyper(SG, netSR, PG)
	
	return(rL)
}

