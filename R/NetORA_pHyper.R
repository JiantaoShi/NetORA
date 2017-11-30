#' The Hypergeometric Distribution in NetORA
#' 
#' This function implements a test based on Hypergeometric Distribution
#' 
#' @param SG a vector of signature
#' @param SR a vector of pathway
#' @param PG a vector of background genes
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

NetORA_pHyper <- function(SG, SR, PG){

	SG <- intersect(SG, PG)
	SR <- intersect(SR, PG)

	m  <- length(SR)
	n  <- length(PG) - m
	k  <- length(SG)
	q  <- length(intersect(SG, SR))
	ES <- signif((q/k)/(m/(m + n)), 4)

	pvalue <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)

	rL <- list(nSG = length(SG), 
			   nSR = length(SR), 
			   ES  = ES, 
			   pvalue   = pvalue,
			   nOverlap = q,
			   gOverlap = sort(intersect(SG, SR)))
	return(rL)
}

