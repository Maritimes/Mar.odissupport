% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_worrmsTSN.R
\name{do_worrmsTSN}
\alias{do_worrmsTSN}
\title{do_worrmsTSN}
\usage{
do_worrmsTSN(recs = NULL, knownAphias = NULL, logName = logName)
}
\arguments{
\item{recs}{- vector of taxa names for which we already have an APHIAID}

\item{knownAphias}{- dataframe of taxas for which we have the AphiaIDs and a
field called "SCI_COL_CLN" containing and identifying name for the taxa.}

\item{logName}{- this is the name of the logfile in the working directory 
that progress should be appended to.}
}
\description{
This function sends APHIAIDs to worrms to see if it can find a 
corresponding TSN
}
\seealso{
Other speciesCodes: \code{\link{assignDefinitive}},
  \code{\link{cleanPrepareSpecList}},
  \code{\link{do_ritis}}, \code{\link{do_taxize}},
  \code{\link{do_worrmsAphiaID}}, \code{\link{do_worrms}},
  \code{\link{getTaxaIDs}}
}
\author{
Mike McMahon, \email{Mike.McMahon@dfo-mpo.gc.ca}
}
