\name{plotJunctionSeqResultsForGene}
\docType{methods}
\alias{plotJunctionSeqResultsForGene}
\title{
   Generate a JunctionSeq expression plot.
}
\description{
   Creates one results plot for one gene. Note that this function does not call a plotting device, so it will simply plot to the "current" device. If you want to automatically save images to file, use \code{\link{buildAllPlotsForGene}}, which internally calls this function.
   
   Note that this function has MANY parameters, allowing the user to tweak the appearance of the plots to suit their particular needs and preferences. Don't be daunted: the default parameters are probably fine for most purposes.
}
\usage{
plotJunctionSeqResultsForGene(geneID, jscs, 
                              colorRed.FDR.threshold=0.05,
                              plot.type = "expr", 
                              displayTranscripts = FALSE,
                              color = NULL, 
                              use.vst = FALSE, use.log = TRUE,truncateBelowOne = TRUE,
                              exon.rescale.factor = 0.3,
                              label.p.vals = TRUE, 
                              plot.lwd = 3,axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                              par.cex = 1, anno.cex.text = 1,
                              anno.cex.axis=anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                              fit.countbin.names = TRUE,
                              plot.gene.level.expression = NULL,
                              plot.exon.results, plot.junction.results, plot.novel.junction.results, 
                              plot.untestable.results = FALSE, draw.untestable.annotation = TRUE,
                              show.strand.arrows = 10, arrows.length = 0.125,
                              sort.features = TRUE,
                              drawCoordinates = TRUE,
                              yAxisLabels.inExponentialForm = FALSE,
                              title.main=NULL, title.ylab=NULL,
                              graph.margins = c(2,3,3,2),
                              GENE.annotation.relative.height = 0.2,
                              TX.annotation.relative.height = 0.025,
                              condition.legend.text = NULL, include.TX.names = TRUE, 
                              draw.start.end.sites = TRUE,
                              verbose=TRUE, debug.mode = FALSE,
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
    be performed by \code{\link{testForDiffUsage}}. Effect sizes and parameter
    estimates must be created via \code{\link{estimateEffectSizes}}.
  }
  
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
  
  \item{plot.gene.level.expression}{
    Logical value. If \code{TRUE}, gene-level expression (when applicable) will be plotted beside the sub-element-specific expression in a small seperate plotting box. 
    The gene-level expression will NOT be plotted for the "relative expression" plots (since you can't meaningfully plot gene-level expression relative to itself).
  }
  \item{plot.exon.results}{Logical. If \code{TRUE}, plot results for exons. By default everything that was tested will be plotted.}
  \item{plot.junction.results}{Logical. If \code{TRUE}, plot results for splice junctions. By default everything that was tested will be plotted.}
  \item{plot.novel.junction.results}{Logical. If \code{TRUE}, plot results for novel splice junctions. If false, novel splice junctions will be ignored. By default everything that was tested will be plotted.}
  \item{plot.untestable.results}{Logical. If TRUE, plots the expression of splice junctions that had coverage that was too low to be tested.}
  \item{draw.untestable.annotation}{Logical. If TRUE, draws the annotation for splice junctions that had coverage that was too low to be tested.}

  \item{show.strand.arrows}{The number of strand-direction arrows to display.}
  \item{arrows.length}{The length of the strand-direction arrows, in inches.}

  \item{sort.features}{Logical. If \code{TRUE}, sort features by genomic position.}

  \item{drawCoordinates}{Whether to label the genomic coordinates at the bottom of the plot.}
  \item{yAxisLabels.inExponentialForm}{Logical. If \code{TRUE}, then the y-axis will be labelled in exponential form.}
  
  \item{graph.margins}{
    Numeric vector of length 4. These margins values used (as if for par("mar")) for the main graph. The lower part of the plot uses the same left and right margins. 
  }
  \item{GENE.annotation.relative.height}{
    The height of the "gene track" displayed underneath the main graph, relative to the height of the main graph. By default it is 20 percent.
  }
  \item{TX.annotation.relative.height}{
    For all plots that draw the annotated-transcript set (when the \code{with.TX} parameter is \code{TRUE}), this sets the height of each transcript, as a fraction of the height of the main graph. By default it is 2.5 percent.
  }
  
  \item{title.main}{
    Character string. Overrides the default main plot title.
  }
  \item{title.ylab}{
    Character string. Overrides the default y-axis label.
  }
  \item{condition.legend.text}{
    List or named vector of character strings. This optional parameter can be used to assign labels to each condition variable values. It should be a list or named vector with length equal to \code{factor(condition)}. Each element
    should be named with one of the values from \code{factor(condition)}, and should contain the label. They will be listed in this order in the figure legend.
  }
  
  \item{include.TX.names}{
    Logical value. If \code{TRUE}, then for the plots that include the annotated transcript, the transcript names will be listed. The labels will be drawn at half the size of \code{anno.cex.text}.
  }
  \item{draw.start.end.sites}{
    Logical value. If \code{TRUE}, then transcript start/end sites will be marked on the main gene annotation.
  }
  
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
  \item{debug.mode}{Logical. If \code{TRUE}, print additional debugging information during execution.}
  \item{...}{Additional options to pass to plotting functions, particularly graphical parameters.}
}
