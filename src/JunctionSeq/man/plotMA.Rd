\name{plotMA}
\docType{methods}
\alias{plotMA}
\title{
   Generate a MA-Plot
}
\description{
    TODO
}
\usage{
plotMA(jscs, 
       FDR.threshold = 0.05, 
       fc.name = NULL,
       fc.thresh = 1,
       use.pch = 19,
       smooth.nbin = 256, 
       ylim = c( 1 / 1000,1000),
       use.smoothScatter = TRUE,
       label.counts = TRUE,
       label.axes = c(TRUE,TRUE,FALSE,FALSE), 
       show.labels = TRUE,
       par.cex = 1, points.cex = 1, text.cex = 1,
       lines.cex = 8,
       anno.lwd = 2, mar = c(4.1,4.1,3.1,1.1),
       miniTicks = TRUE, ...)
}
\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually created by \code{\link{runJunctionSeqAnalyses}}.
    
    Alternatively, this can be created manually by \code{\link{readJunctionSeqCounts}}. 
    However in this case a number of additional steps will be necessary: 
    Dispersions and size factors must then be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}. Hypothesis tests must
    be performed by \code{\link{testForDiffUsage}}. Effect sizes and parameter
    estimates must be created via \code{\link{estimateEffectSizes}}.
  }
  \item{FDR.threshold}{
    The FDR threshold used to color dots. Tests with an adjusted-p-value more significant than this threshold will be marked in red.
  }
  \item{fc.name}{
    The name of the column to take from fData(jscs).
  }
  \item{fc.thresh}{
    The fold-change threshold required to count a significant locus in the count labels. It will also draw horizontal lines at this threshold.
  }
  
  \item{use.pch}{
    The value of pch to pass to the \code{\link{points}} call.
  }
  \item{use.smoothScatter}{
    Logical. If TRUE, non-significant genes will be ploted with density shading.
  }
  \item{smooth.nbin}{
    The number of bins to smooth, for the density plot, if \code{use.smoothScatter} is TRUE.
  }

  \item{ylim}{
    The y-axis limits.
  } 
  \item{label.counts}{
    Logical. If TRUE, include labels showing the number of loci that pass both the statistical-significance and fold-change threshold in each direction.
  }
  \item{label.axes}{
    Logical vector. Whether to label each axis. Must have length 4; each corresponds to the bottom, left, top, and right axes respectively.
  }  
  \item{show.labels}{
    Logical. If TRUE, include all titles and axes labels.
  }
  \item{par.cex}{
    The cex value to be passed to \code{\link{par}}.
  }
  \item{points.cex}{
    The cex value to be passed to \code{\link{points}}.
  }
  \item{text.cex}{
    The cex value to be passed to \code{\link{text}}.
  }
  \item{lines.cex}{
    The cex value to be passed to \code{\link{lines}}, \code{\link{box}}, and similar.
  }
  \item{anno.lwd}{
    The lwd value to be passed to \code{\link{lines}}, \code{\link{box}}, \code{\link{axis}}, and similar.
  }
  \item{mar}{
    The margin sizes, expressed in lines. see \code{link{par}}.
  }
  \item{miniTicks}{
    Logical. If TRUE, then include "mini tick marks" on the x and y axes.
  }
  
  \item{...}{
    Additional graphical parameters.
  }
}
























