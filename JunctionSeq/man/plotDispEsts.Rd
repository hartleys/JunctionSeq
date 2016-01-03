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
               xlab = "Mean Normalized Coverage", 
               ylab = "Dispersion",
               miniTicks = TRUE,
               pch.MLE = 46, pch.MAP = 1, lwd.fitted = 2,
               par.cex = 1, points.cex = 1, text.cex = 1, lines.cex = 8,
               use.smoothScatter = FALSE, smooth.nbin = 512, nrpoints = 100,
               plot.exon.results = TRUE, 
               plot.junction.results = TRUE, 
               anno.lwd = 2,
               mar = c(4.1,4.1,3.1,1.1),
               show.legends = TRUE,
               verbose = TRUE, debug.mode = FALSE,
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
  
  \item{pch.MLE}{
    Numeric. The pch code for the MLE (ie single-feature) dispersion estimate. The default is a small point.
  }
  \item{pch.MAP}{
    Numeric. The pch code for the MAP (ie. final) dispersion estimate. The default is a circle.
  }
  \item{lwd.fitted}{
    Numeric. The width of the dispersion fit line(s).
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
  
  \item{anno.lwd}{
    The lwd value to be passed to \code{\link{lines}}, \code{\link{box}}, \code{\link{axis}}, and similar.
  }
  \item{mar}{
    The margin sizes, expressed in lines. see \code{link{par}}.
  }
  
  \item{show.legends}{
    Logical. If \code{TRUE}, display legends.
  }
  
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
  \item{debug.mode}{if TRUE, send even more debugging and progress messages to the console / stdout.}
  
  \item{...}{Additional options to pass to plotting functions, particularly graphical parameters.}

}
\value{
  This is a side-effecting function, and does not return a value.
}

\examples{
data(exampleDataSet,package="JctSeqData");
plotDispEsts(jscs);

\dontrun{
########################################
#Set up example data:
decoder.file <- system.file(
                  "extdata/annoFiles/decoder.bySample.txt",
                  package="JctSeqData");
decoder <- read.table(decoder.file,
                  header=TRUE,
                  stringsAsFactors=FALSE);
gff.file <- system.file(
            "extdata/cts/withNovel.forJunctionSeq.gff.gz",
            package="JctSeqData");
countFiles <- system.file(paste0("extdata/cts/",
     decoder$sample.ID,
     "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
     package="JctSeqData");
######################
#Run example analysis:
jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
           sample.names = decoder$sample.ID,
           condition=factor(decoder$group.ID),
           flat.gff.file = gff.file,
           analysis.type = "junctionsAndExons"
);
########################################

#Plot dispersions:
plotDispEsts(jscs);

}
}