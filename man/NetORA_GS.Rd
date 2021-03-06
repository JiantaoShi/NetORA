% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NetORA_GS.R
\name{NetORA_GS}
\alias{NetORA_GS}
\title{Network-based enrichment for a list of signatures}
\usage{
NetORA_GS(sgList, queryTable, PG, FDR = 0.05)
}
\arguments{
\item{sgList}{a vector of gene signature}

\item{queryTable}{a data frame returned by NetORA_Pre}

\item{PG}{a vector of background genes}

\item{FDR}{an FDR cutoff to define gene that are significantly connected a pathway (defaul 0.05)}
}
\value{
a list
}
\description{
This function is wraper of core function NetORA to work on a list of gene sets
}
\examples{
library("NetORA")
data(text_2006_12_NetGPA)
data(Example_NetGPA)
ExE_Hyper  <- Example_NetGPA$ExE_Hyper
queryTable <- NetORA_Pre(ExE_Hyper, text_2006_12_NetGPA, progressBar = TRUE)
Cancer_GS  <- Example_NetGPA$Cancer_GeneSet
PG <- colnames(text_2006_12_NetGPA)
mergedT <- NetORA_GS(Cancer_GS, queryTable, PG, FDR = 0.05)
}
\references{
Jiantao Shi: NetORA, a package for network-based pathway over-representation analysis.
}
\author{
Jiantao Shi
}
