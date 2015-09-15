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
  plotDispEsts( jscs, ylim, xlim, 
               linecol=c("#0000FF","#FF0000"), 
               pointcol = c("#00009980","#99000080"),
               title.main = "Dispersion Estimates", 
               xlab = "mean of normalized counts", 
               ylab = "dispersion",
               miniTicks = TRUE,
               pch = c(1,4),
               par.cex = 1, points.cex = 1, text.cex = 1, lines.cex = 8,
               use.smoothScatter = FALSE, smooth.nbin = 512, nrpoints = 100,
               plot.exon.results = TRUE, 
               plot.junction.results = TRUE, 
               ... )
}

\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually created by \code{\link{runJunctionSeqAnalyses}}.
    
    Alternatively, this can be created manually by \code{\link{readJunctionSeqCounts}}. 
    Dispersions and size factors must then be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}. Hypothesis tests must
    be performed by \code{\link{testForDiffUsage}}. 
  }
  \item{ylim}{
    The plotting range for the y-axis.
  }
  \item{xlim}{
    The plotting range for the x-axis.
  }
  \item{linecol}{
    Character vector of length 2. The line color to use for the fit line. If the fits were performed separately for exons and junctions, the junction line will be drawn with the second color.
  }
  \item{pointcol}{
    Character vector of length 2. The point color to use for the final dispersions. If the fits were performed separately for exons and junctions, the junction points will be drawn with the second color.
  }
  \item{title.main}{
    The main title of the plot.
  }
  \item{xlab}{
    The label for the x-axis.
  }
  \item{ylab}{
    The label for the y-axis.
  }
  \item{miniTicks}{
    Whether or not to plot smaller ticks at the tenth-decades.
  }
  \item{par.cex}{
     The base cex value to be passed to par() immediately before all plots are created. See \code{\link{par}}.
  }
  \item{points.cex}{
     The character expansion value for the plotted points.
  }
  \item{text.cex}{
     The character expansion value for the annotation text (labels, etc).
  }
  \item{lines.cex}{
    The character expansion value for lines. What this means seems to vary depending on the plotting device.
  }
  
  \item{pch}{
    Numeric vector of length 2. Contains the pch parameters for the exon features and the junction features. This determines the character used in plotting. By default the exons are plotted as circles, and junctions are plotted as X's.
  }
  \item{use.smoothScatter}{
    Logical. If TRUE, features will be ploted with density shading rather than having each point plotted.
  }
  \item{smooth.nbin}{
    The number of bins to smooth, for the density plot, if \code{use.smoothScatter} is TRUE.
  }
  \item{nrpoints}{
    The number of extra points to plot, if \code{use.smoothScatter} is TRUE.
  }
  \item{plot.exon.results}{Logical. If \code{TRUE}, plot results for exons. Technically speaking, JunctionSeq can be used to do DEXSeq-style analyses on exon partitions. However this functionality is for advanced users only.}
  \item{plot.junction.results}{Logical. If \code{TRUE}, plot results for splice junctions. For advanced users only.}
  \item{...}{Additional options to pass to plotting functions, particularly graphical parameters.}

}
