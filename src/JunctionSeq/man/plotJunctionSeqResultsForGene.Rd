\name{plotJunctionSeqResultsForGene}
\docType{methods}
\alias{plotJunctionSeqResultsForGene}
\title{
   Plot JunctionSeq results
}
\description{
   Creates one results plot for one gene.
}
\usage{
plotJunctionSeqResultsForGene(geneID, jscs, flat.gff.data, 
                              colorRed.FDR.threshold=0.05,
                              plot.type = "expr", 
                              displayTranscripts = FALSE,
                              color = NULL, 
                              use.vst = FALSE, use.log = FALSE,truncateBelowOne = TRUE,
                              exon.rescale.factor = 0.3,
                              label.p.vals = TRUE, 
                              plot.lwd = 3,axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                              par.cex = 1, anno.cex.text = 1,
                              anno.cex.axis=anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                              fit.countbin.names = TRUE,
                              debug.mode = FALSE,
                              plot.exon.results = FALSE, plot.splice.results = TRUE, plot.novel.splice.results = plot.splice.results, 
                              plot.untestable.results = FALSE,
                              show.strand.arrows = 20,
                              sort.features = TRUE,
                              drawCoordinates = TRUE,
                              arrows.length = 0.25,
                              title.main=NULL, title.ylab=NULL,
                              verbose=TRUE,
                              ...)
}
\arguments{
  \item{geneID}{Character string. The gene to the plotted.}
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually created by \code{\link{runJunctionSeqAnalyses}}.
    
    Alternatively, this can be created manually by \code{\link{readJunctionSeqCounts}}. 
    However in this case a number of additional steps will be necessary: 
    Dispersions and size factors must then be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}. Hypothesis tests must
    be performed by \code{\link{testJunctionsForDJU}}. Effect sizes and parameter
    estimates must be created via \code{\link{estimateEffectSizes}}.
  }
  \item{flat.gff.data}{A data frame containing the annotation data. This can be extracted from a flat gff file using the \code{\link{readAnnotationData}} function.}
  
  \item{colorRed.FDR.threshold}{The adjusted-p-value threshold used to determine whether a feature should be marked as "significant" and colored pink. By default this will be the same as the FDR.threshold.}
  
  \item{plot.type}{
    Character string. Determines which plot to produce. Options are: "expr" for "expression", or mean normalized read counts by experimental condition, "rExpr" for 
    "relative" expression relative to gene-level expression, "normCounts" for normalized read counts for each sample, and "rawCounts" for raw read counts for each sample.
  }
  
  \item{displayTranscripts}{Logical. If true, then the full set of annotated transcripts will be displayed below the expression plot (to a maximum of 42 different TX).}
  \item{color}{A vector of R colors, named for each possible value of condition. By default, it will attempt to choose reasonable colors for each condition.}

  \item{use.vst}{Logical. If TRUE, all plots will be scaled via a variance stabilizing transform.}
  \item{use.log}{Logical. If TRUE, all plots will be log-scaled.}
  \item{truncateBelowOne}{Logical. If TRUE, all values between 0 and 1 will be drawn with a linear scale. If use.log is FALSE, this does nothing.}
  \item{exon.rescale.factor}{Numeric. Exons will be proportionately scaled-up so that the exonic regions make up this fraction of the horizontal plotting area. If negative, exons and introns will be plotted to a common scale. }
  \item{label.p.vals}{Logical. If \code{TRUE}, then statistically significant p-values will be labelled.}

  \item{plot.lwd}{the line width for the plotting lines.}
  \item{axes.lwd}{the line width for the axes.}
  \item{anno.lwd}{the line width for the various other annotation lines.}
  \item{par.cex}{
     The base cex value to be passed to par() immediately before all plots are created. See \code{\link{par}}.
  }
  \item{anno.cex.text}{
     The font size multiplier for most annotation text. This will be multiplied by a factor of the par.cex value.
     More specifically: The cex value to be passed to all function calls that take \link{graphical parameters}. See \code{\link{par}}.
  }
  \item{anno.cex.axis}{
     The font size multiplier for the axis text. This will be multiplied by a factor of the par.cex value.
     More specifically: The cex.axis value to be passed to all function calls that take \link{graphical parameters}. See \code{\link{par}}.
  }
  \item{anno.cex.main}{
     The font size multiplier for the main title text. This will be multiplied by a factor of the par.cex value.
     More specifically: The cex.main value to be passed to all function calls that take \link{graphical parameters}. See \code{\link{par}}.
  }
  \item{fit.countbin.names}{Logical. If \code{TRUE}, then splice-junction-locus labels should be rescaled to fit in whatever horizontal space is available.}
  \item{debug.mode}{Logical. If \code{TRUE}, print additional debugging information during execution.}
  
  
  \item{plot.exon.results}{Logical. If TRUE, plot results for exons. Technically speaking, JunctionSeq can be used to do DEXSeq-style analyses on exon partitions. However this functionality is for advanced users only.}
  \item{plot.splice.results}{Logical. If TRUE, plot results for splice junctions. For advanced users only.}
  \item{plot.novel.splice.results}{Logical. If TRUE, plot results for novel splice junctions. If false, novel splice junctions will be ignored.}
  \item{plot.untestable.results}{Logical. If TRUE, plots splice junctions that had coverage that was too low to be tested.}
  
  \item{show.strand.arrows}{The number of strand-direction arrows to display.}
  \item{sort.features}{Logical. If \code{TRUE}, sort features by genomic position.}

  \item{drawCoordinates}{Whether to label the genomic coordinates at the bottom of the plot.}
  
  \item{arrows.length}{The length of the strand-direction arrows, in inches.}
  
  \item{title.main}{
    Character string. Overrides the default main plot title.
  }
  \item{title.ylab}{
    Character string. Overrides the default y-axis label.
  }

  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
  \item{...}{Additional options to pass to plotting functions, particularly graphical parameters.}
}
