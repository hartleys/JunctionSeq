\name{estimateJunctionSeqDispersions}
\docType{methods}
\alias{estimateJunctionSeqDispersions}
\title{
   JunctionSeq Dispersion Estimation
}
\description{
   This method estimates the sample dispersion for each 
   counting bin (in other words, each splice junction locus). 
   
   This function is called internally by the \code{\link{runJunctionSeqAnalyses}} 
   function, and thus for most purposes users should not need to call this
   function directly. It may be useful to advanced users performing non-standard
   analyses.
}
\usage{
  estimateJunctionSeqDispersions(jscs, 
          test.formula1 = formula(~ sample + countbin + condition : countbin), 
          meanCountTestableThreshold=7.5, nCores=1,
          test.aggregated.genes = FALSE,
          use.alternate.method  = TRUE,
          verbose = TRUE);
}
\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually initially created by 
    \code{\link{readJunctionSeqCounts}}. Size factors must be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}.
  }
  \item{test.formula1}{
    The model formula. Note that this formula is different from the 
    formula used to calculate parameter estimates and effect size.
    This is because the two noise components (gene-level and 
    countbin-level noise) are folded into the sample term. Since
    we only intend to test the condition-countbin interaction,
    we do not need to model the gene-level differential expression.
    
    NOTE: the biological condition to be tested MUST be named 
    "condition".
  }
  \item{meanCountTestableThreshold}{
    Numeric value. Features with a total mean normalized count of less than this value will be excluded from the analyses.
  }
  \item{nCores}{
    The number of cores to use. Note that multicore functionality may not be available on all platforms.
  }
  \item{test.aggregated.genes}{Logical value. Whether to attempt to test "aggregate genes" which consist of multiple genes that overlap with one another. Note that inclusion of aggregate genes may affect the false discovery rate, since by their very nature aggregate genes will often show differential splice junction usage, as the two genes will often be regulated independently.}
  \item{use.alternate.method}{
    ADVANCED USERS ONLY: Logical value. Determines which model framework to use.
  }
  \item{verbose}{A boolean flag indicating whether or not to print progress information during execution. (Default=FALSE)}

}
\value{
  A JunctionSeqCountSet, with dispersion results included.
}