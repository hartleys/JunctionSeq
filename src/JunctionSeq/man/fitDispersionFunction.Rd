\name{fitDispersionFunction}
\docType{methods}
\alias{fitDispersionFunction}
\title{
   Fit Shared Dispersion Function
}
\description{
   Fit dispersion function to share dispersion information between 
   features across the genome.
   
   This function is called internally by the \code{\link{runJunctionSeqAnalyses}} 
   function, and thus for most purposes users should not need to call this
   function directly. It may be useful to advanced users performing non-standard
   analyses.
}
\usage{
  fitDispersionFunction(jscs, fitType = c("parametric","local","mean"), 
                        fitDispersionsForExonsAndJunctionsSeparately = TRUE, 
                        advancedMode = TRUE,
                        verbose = TRUE);
}
\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually initially created by 
    \code{\link{readJunctionSeqCounts}}. Size factors must be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}.
  }
  \item{fitType}{
    Which dispersion fit method to use. Not that currently only the parametric (default) method is supported. The other methods are currently nonfunctional and will be added in future releases.
  }
  \item{fitDispersionsForExonsAndJunctionsSeparately}{
    When running a "hybrid" analysis in which both exons and splice junctions are being tested simultaniously, this parameter
    determines whether a single fitted dispersion model should be fitted for both exons and splice junctions, or if separate 
    fitted dispersions should be calculated for each. By default the dispersions are run separately.
  }
  \item{advancedMode}{
    Switches between the current default method and an older variant.
  }
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
}
\value{
  A JunctionSeqCountSet, with dispersion results included.
}


