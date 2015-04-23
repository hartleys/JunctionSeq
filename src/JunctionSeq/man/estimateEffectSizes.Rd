\name{estimateEffectSizes}
\docType{methods}
\alias{estimateEffectSizes}
\title{
   Estimate Effect Factors
}
\description{
   This function runs fits another generalized linear model to 
   the data, this one intended for use in estimating the
   effect sizes and expression estimates.
   
   For most purposes this function should not be needed, but might 
   be useful to advanced users performing non-standard analyses. 
   This and all other standard analysis functions can be run 
   automatically via the \code{\link{runJunctionSeqAnalyses}} 
   function.
}
\usage{
   estimateEffectSizes(jscs, 
        effect.formula = formula(~ condition + countbin + condition : countbin),
        nCores=1, 
        dispColumn="dispersion",
        verbose = TRUE)}
\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually initially created by 
    \code{\link{readJunctionSeqCounts}}. Size factors must be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}.
  }
  
  \item{effect.formula}{
    For advanced users. The base formula for the model used for effect size estimation.
    
    NOTE: the biological condition to be tested must be named "condition".
  }
  \item{nCores}{
    The number of cores to use. Note that multicore functionality may not be available on all platforms. 
    JunctionSeq attempts to use the BiocParallell package if it can be found installed. Otherwise it will
    attempt to fallback to the multicore package. If neither package can be found it will fallback to single
    core execution.
  }
  \item{dispColumn}{
    Character value. The name of the \code{fData(jscs)} column in which the model dispersion is stored.
  }
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
}
\value{
  A JunctionSeqCountSet, with effect size results included.
}