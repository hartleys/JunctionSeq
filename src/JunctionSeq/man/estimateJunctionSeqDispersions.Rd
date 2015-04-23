\name{estimateJunctionSeqDispersions}
\docType{methods}
\alias{estimateJunctionSeqDispersions}
\title{
   JunctionSeq Dispersion Estimation
}
\description{
   This method estimates the sample dispersion for each 
   counting bin (in other words, each splice junction locus). 
   
   For most purposes this function should not be needed, but might 
   be useful to advanced users performing non-standard analyses. 
   This and all other standard analysis functions can be run 
   automatically via the \code{\link{runJunctionSeqAnalyses}} 
   function.
}
\usage{
  estimateJunctionSeqDispersions(jscs, 
          test.formula1 = formula(~ sample + countbin + condition : countbin), 
          minCount=10, nCores=1,
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
  \item{minCount}{
    Integer value. The minimum mean normalized read count threshold 
    required. Default is 10, which is reasonable for most purposes.
  }
  \item{nCores}{
    The number of cores to use. Note that multicore functionality may not be available on all platforms.
  }
  \item{test.aggregated.genes}{Logical value. Whether to attempt to test "aggregate genes" which consist of multiple genes that overlap with one another. Note that inclusion of aggregate genes may affect the false discovery rate, since by their very nature aggregate genes will often show differential splice junction usage, as the two genes will often be regulated independently.}
  \item{use.alternate.method}{
    ADVANCED USERS ONLY: Logical value. Determines whether to use the JunctionSeq standard model framework, or the DEXSeq model framework.
  }
  \item{verbose}{A boolean flag indicating whether or not to print progress information during execution. (Default=FALSE)}

}
\value{
  A JunctionSeqCountSet, with dispersion results included.
}