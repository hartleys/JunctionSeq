\name{fitDispersionFunction}
\docType{methods}
\alias{fitDispersionFunction}
\title{
   Fit Shared Dispersion Function
}
\description{
   Fit dispersion function to share dispersion information between 
   features across the genome.
   
   For most purposes this function should not be needed, but might 
   be useful to advanced users performing non-standard analyses. 
   This and all other standard analysis functions can be run 
   automatically via the \code{\link{runJunctionSeqAnalyses}} 
   function.
}
\usage{
  fitDispersionFunction(jscs, verbose);
}
\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually initially created by 
    \code{\link{readJunctionSeqCounts}}. Size factors must be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}.
  }
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
}
\value{
  A JunctionSeqCountSet, with dispersion results included.
}