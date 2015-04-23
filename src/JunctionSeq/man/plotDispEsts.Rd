\name{plotDispEsts}
\docType{methods}
\alias{plotDispEsts}
\title{
   Plot Fitted and Test-wise Dispersion
}
\description{
    Plots the countbin-specific estimated dispersion and the fitted dispersion curve.
}
\usage{
  plotDispEsts( jscs, ymin, 
                linecol="#ff000080",
                xlab = "mean of normalized counts", 
                ylab = "dispersion",
                log = "xy", cex = 0.45, ... )
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
  \item{ymin}{
    The minimum for the y-axis.
  }
  \item{linecol}{
    The line color to use for the fit line.
  }
  \item{xlab}{
    The label for the x-axis.
  }
  \item{ylab}{
    The label for the y-axis.
  }
  \item{log}{
    Logical value. Determines which axes to log-scale.
  }
  \item{cex}{
    The cex value to be passed to \code{plot}.
  }
  \item{...}{Additional options to pass to plotting functions, particularly graphical parameters.}

}
