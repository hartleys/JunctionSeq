\name{generateAllExpressionEstimates}
\docType{methods}
\alias{generateAllExpressionEstimates}
\title{
   TODO
}
\description{
   ADVANCED USERS ONLY: This function is intended for advanced users only. 
   For general use, using the \code{\link{generateCompleteResults}} function is recommended, 
   which will in turn call this function internally.
}
\usage{
  generateAllExpressionEstimates(jscs,
                                 nCores = 1,
                                 fitExpToVar="condition",
                                 verbose=TRUE)
}
\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually created by \code{\link{runJunctionSeqAnalyses}}.
    
    Alternatively, this can be created manually by \code{\link{readJunctionSeqCounts}}. 
    Dispersions and size factors must then be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}. Hypothesis tests must
    be performed by \code{\link{testJunctionsForDJU}}. 
  }
  \item{nCores}{
    The number of cores to use. Note that multicore functionality may not be available on all platforms.
  }
  \item{fitExpToVar}{
    Character value. The biological condition to fit the model to.
  }
  \item{verbose}{
    Logical value. Determines how much debugging into to message to the console.
  }
}