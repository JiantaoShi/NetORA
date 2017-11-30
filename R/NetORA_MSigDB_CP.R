#' Network-based enrichment for a list of signatures
#' 
#' This function is wraper of core function NetORA to work on a list of gene sets
#' 
#' @param SG a vector of gene signature
#' @param MSigDB_NetORA_GS pre-computed neighbors of Canonical pathways as defined in MSigDB
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
#' data(MSigDB_NetORA_GS)
#' ExE_Hyper  <- Example_NetGPA$ExE_Hyper
#' PG <- colnames(text_2006_12_NetGPA)
#' CancerPathway <- paste0("REACTOME_", Example_NetGPA$CancerPathway)
#' mergedT <- NetORA_MSigDB_CP(ExE_Hyper, MSigDB_NetORA_GS[CancerPathway], PG)
#' @export

NetORA_MSigDB_CP <- function(SG, MSigDB_NetORA_GS, PG){

	Ns = length(MSigDB_NetORA_GS)
	rM = matrix(NA, Ns, 5)

	for(i in 1:Ns){

		SR    = MSigDB_NetORA_GS[[i]]
		rL    = NetORA_pHyper(SG, SR, PG) 
		rM[i, 1] = rL$nSG
		rM[i, 2] = rL$nSR
		rM[i, 3] = rL$nOverlap
		rM[i, 4] = rL$ES
		rM[i, 5] = rL$pvalue
	}
	
	Pathway      = names(MSigDB_NetORA_GS)
	colnames(rM) = c("nSG", "nSR", "nOverlap", "ES", "pvalue")
	FDR = p.adjust(rM[, 5])

	mTable = data.frame(Pathway, rM, FDR)

	return(mTable)
}
