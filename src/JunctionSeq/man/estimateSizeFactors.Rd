\name{estimateSizeFactors}
\docType{methods}
\alias{estimateSizeFactors}
\title{
   Estimate Size Factors
}
\description{
   Estimate size factors, which are scaling factors used as "offsets" 
   by the statistical model to make the different samples comparable.
   This is necessary because the different samples may have been
   sequenced to slightly different depths. Additionally, the presence
   of differentially expressed genes may cause the apparent depth of
   many genes to appear different.
   
   This function uses the "geometric" size factor normalization method,
   which is identical to the one used by DESeq, DESeq2, DEXSeq, and 
   the default method used by CuffDiff.
   
   For most purposes this function should not be needed, but might 
   be useful to advanced users performing non-standard analyses. 
   This and all other standard analysis functions can be run 
   automatically via the \code{\link{runJunctionSeqAnalyses}} 
   function.
}
\usage{
  estimateSizeFactors(jscs);
}
\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually initially created by 
    \code{\link{readJunctionSeqCounts}}. Size factors must be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}.
  }
}
\value{
  A JunctionSeqCountSet, with size factors included.
}