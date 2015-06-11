\name{estimateEffectSizes}
\docType{methods}
\alias{estimateEffectSizes}
\title{
   Estimate Effect Sizes, parameter estimates, etc.
}
\description{
   This function runs fits another generalized linear model to 
   the data, this one intended for use in estimating the
   effect sizes and expression estimates for each analysis.
   
   This function is called internally by the \code{\link{runJunctionSeqAnalyses}} 
   function, and thus for most purposes users should not need to call this
   function directly. It may be useful to advanced users performing non-standard
   analyses.
}
\usage{
   estimateEffectSizes(jscs, 
        effect.formula = formula(~ condition + countbin + condition : countbin),
        geneLevel.formula = formula(~ condition),
        calculate.geneLevel.expression = TRUE,
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
  \item{geneLevel.formula}{
    For advanced users. The base formula for the model used to estimate total gene-level expression.
    
    NOTE: the biological condition to be tested must be named "condition".
  }
  \item{calculate.geneLevel.expression}{
    Logical value. If \code{TRUE}, gene-level expression will be estimated using the same maximum-likelihood method used in other analyses. Default: \code{TRUE}.
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