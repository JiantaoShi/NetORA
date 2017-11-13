#' Building sub-network in background
#' 
#' This function implements a preprocessing step that find all neighbor of a given pathway.
#' 
#' @param SR a vector of pathway
#' @param netMatrix a integer matrix which is accepted by NetGPA as a network
#' 
#' @author Jiantao Shi
#' @references
#' 
#' Jiantao Shi: NetORA, a package for network-based pathway over-representation analysis.
#' 
#' @return
#' a data frame
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

NetORA_Pre <- function(SR, netMatrix, progressBar = TRUE){

	check <- require(NetGPA)
	if(!check)
		stop("NetGPA is required.\n")

	allGene <- colnames(netMatrix)

	if(length(SR) < 20)
		stop("Pathway size is too small to build a network.\n")

	queryTable <- NetGPA(as.list(SR), allGene, netMatrix, Pfcutoff = 0.1, progressBar = progressBar)$queryTable

	return(queryTable)
}
