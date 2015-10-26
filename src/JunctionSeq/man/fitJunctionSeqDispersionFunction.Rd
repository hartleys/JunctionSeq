\name{fitJunctionSeqDispersionFunction}
\docType{methods}
\alias{fitDispersionFunction}
\alias{fitJunctionSeqDispersionFunction}
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
fitJunctionSeqDispersionFunction(jscs, 
    method.GLM = c(c("advanced","DESeq2-style"), c("simpleML","DEXSeq-v1.8.0-style")),
    method.dispFit = c("parametric", "local", "mean"), 
    method.dispFinal = c("shrink","max","fitted","noShare"),
    fitDispersionsForExonsAndJunctionsSeparately = TRUE, 
    verbose = TRUE)

fitDispersionFunction(...) #Equivalent but depreciated.
}
\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually initially created by 
    \code{\link{readJunctionSeqCounts}}. Size factors must be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}.
  }
  \item{method.GLM}{
    Character string. Can be used to apply alternative methodologies or implementations. Intended for advanced users who have strong opinions about the underlying statistical methodologies.
    
    The default is "advanced" or, equivalently, "DESeq2-style". This uses the 
    dispersion estimation methodology used by DESeq2 and DEXSeq v1.12.0 or higher to generate the initial (feature-specific) dispersion estimates. The alternative 
    method is "simpleML" or, equivalently, "DEXSeq-v1.8.0-style". This uses a simpler maximum-likelihood-based method used by the original DESeq and
    by DEXSeq v1.8.0 or less.
  }
  \item{method.dispFit}{
    Character string. Can be used to apply alternative methodologies or implementations. Intended for advanced users who have strong opinions about the underlying statistical methodologies.
    
    Determines the method used to generated "fitted" dispersion estimates. One of "parametric" (the default), "local", or "mean". See the DESeq2 documentation for more information.
  }
  \item{method.dispFinal}{
    Character string. Can be used to apply alternative methodologies or implementations. Intended for advanced users who have strong opinions about the underlying statistical methodologies.

    Determines the method used to arrive at a "final" dispersion estimate. The default, "shrink" uses the maximum a posteriori estimate, combining information from both
    the fitted and feature-specific dispersion estimates. This is the method used by DESeq2 and DEXSeq v1.12.0 and above.
  }
  \item{fitDispersionsForExonsAndJunctionsSeparately}{
    When running a "junctionsAndExons" type analysis in which both exons and splice junctions are being tested simultaniously, this parameter
    determines whether a single fitted dispersion model should be fitted for both exons and splice junctions, or if separate 
    fitted dispersions should be calculated for each. By default the dispersions are run separately.
  }
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
  \item{...}{If using the depreciated fitDispersionFunction command, use the same syntax as above.}
}
\value{
  A JunctionSeqCountSet, with dispersion results included.
}


