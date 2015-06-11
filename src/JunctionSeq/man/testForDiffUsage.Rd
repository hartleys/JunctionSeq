\name{testForDiffUsage}
\docType{methods}
\alias{testForDiffUsage}
\title{
   Test Junctions for Differential Junction Usage
}
\description{
   This function runs the hypothesis tests for differential junction
   usage.
   
   This function is called internally by the \code{\link{runJunctionSeqAnalyses}} 
   function, and thus for most purposes users should not need to call this
   function directly. It may be useful to advanced users performing non-standard
   analyses.
}
\usage{
testForDiffUsage( jscs, 
          test.formula0 = formula(~ sample + countbin), 
          test.formula1 = formula(~ sample + countbin + condition : countbin),
          dispColumn="dispersion", nCores=1, 
          use.alternate.method = TRUE, 
          verbose = TRUE)
}

\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually initially created by 
    \code{\link{readJunctionSeqCounts}}. Dispersions and size factors must be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}.
  }
  \item{test.formula0}{
    The formula for the null hypothesis. Note that the condition to be tested must be named "condition".
  }
  \item{test.formula1}{
    The formula for the alternative hypothesis. Note that the condition to be tested must be named "condition".
  }
  \item{dispColumn}{
    Character value. The name of the \code{fData(jscs)} column in which the model dispersion is stored.
  }
  \item{nCores}{
    The number of cores to use. Note that multicore functionality may not be available on all platforms. 
    JunctionSeq attempts to use the BiocParallell package if it can be found installed. Otherwise it will
    attempt to fallback to the multicore package. If neither package can be found it will fallback to single
    core execution.
  }
  \item{use.alternate.method}{
    ADVANCED USERS ONLY: Logical value. Determines whether to use the JunctionSeq standard model framework, or the DEXSeq model framework.
  }
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}

}

\value{
  A JunctionSeqCountSet, with hypothesis test results included.
}
