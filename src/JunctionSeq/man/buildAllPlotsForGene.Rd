\name{buildAllPlotsForGene}
\docType{methods}
\alias{buildAllPlotsForGene}
\title{
   Creating result plots
}
\description{
   Saves a large battery of plots displaying the analysis results, for the purposes of data visualization.
}
\usage{

buildAllPlotsForGene(geneID, jscs,
                     outfile.prefix = "./", 
                     flat.gff.data = NULL, flat.gff.file = NULL,
                     use.plotting.device = "png", 
                     use.vst=FALSE, use.log = TRUE, truncateBelowOne = TRUE,
                     exon.rescale.factor = 0.3,
                     with.TX=TRUE,without.TX=TRUE,
                     expr.plot=TRUE,normCounts.plot=TRUE,
                     rExpr.plot=TRUE,rawCounts.plot=TRUE,
                     colorRed.FDR.threshold = 0.01, 
                     color=NULL,
                     plot.exon.results = FALSE, 
                     plot.splice.results = TRUE, 
                     plot.novel.splice.results = plot.splice.results,
                     plot.untestable.results = FALSE,
                     plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                     par.cex = 1, anno.cex.text = 2,
                     anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                     drawCoordinates = TRUE,
                     arrows.length = 0.25,
                     base.plot.height = 22.222, base.plot.width = 22.222, 
                     base.plot.units = "in", 
                     plotting.device.params = list(pointsize = 18, res = 72),
                     openPlottingDeviceFunc = NULL, closePlottingDeviceFunc = NULL,
                     verbose=TRUE,  ...)

}
\arguments{
  \item{geneID}{Character string. Which gene to plot.}

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
  \item{outfile.prefix}{The prefix file path to save the images to.}

  \item{flat.gff.data}{
    A data frame containing the annotation data. This can be extracted from a flat gff file using the \code{\link{readAnnotationData}} function.
    If multiple plots are going to be generated it may be preferable to use this to save runtime.
  }
  \item{flat.gff.file}{A flattened gff file containing the annotation data. Only one of flat.gff.data or flat.gff.file should be set.}
  

  \item{use.plotting.device}{The plotting device to use.}
  \item{use.vst}{Logical. If TRUE, all plots will be scaled via a variance stabilizing transform.}
  \item{use.log}{Logical. If TRUE, all plots will be log-scaled.}
  \item{truncateBelowOne}{Logical. If TRUE, all values below 1 will be linear-scaled. If use.log is FALSE, this does nothing.}
  \item{exon.rescale.factor}{Numeric. Exons will be proportionately scaled-up so that the exonic regions make up this fraction of the horizontal plotting area. If negative, exons and introns will be plotted to a common scale. }

  \item{with.TX}{if TRUE, save expression plots with the full transcripts printed}
  \item{without.TX}{if TRUE, save expression plots with only the compiled exons printed. Note that if this and with.TX.plot are both TRUE, both versions will be saved seperately.}
  \item{expr.plot}{if TRUE, save an expression plot of the expression parameter estimates for each splice site, for each condition.}
  \item{normCounts.plot}{if TRUE, save an expression plot of the normalized mean counts for each splice site, for each sample.}
  \item{rExpr.plot}{if TRUE, save an expression plot of the expression parameter estimates, relative to gene-wide expression, for each splice site, for each condition.}
  \item{rawCounts.plot}{if TRUE, save an expression plot of the raw counts for each splice site, for each sample. Note that these will never be VST-transformed, even when use.vst == TRUE.}
  
  \item{colorRed.FDR.threshold}{The adjusted-p-value threshold used to determine whether a feature should be marked as "significant" and colored pink. By default this will be the same as the FDR.threshold.}
  \item{color}{A vector of R colors, named for each possible value of condition. By default, it will attempt to choose reasonable colors for each condition.}
  
  \item{plot.exon.results}{Logical. If TRUE, plot results for exons. Technically speaking, JunctionSeq can be used to do DEXSeq-style analyses on exon partitions. However this functionality is for advanced users only.}
  \item{plot.splice.results}{Logical. If TRUE, plot results for splice junctions. For advanced users only.}
  \item{plot.novel.splice.results}{Logical. If TRUE, plot results for novel splice junctions. If false, novel splice junctions will be ignored.}
  \item{plot.untestable.results}{Logical. If TRUE, plots splice junctions that had coverage that was too low to be tested.}
  
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
  \item{drawCoordinates}{Whether to label the genomic coordinates at the bottom of the plot.}
  \item{arrows.length}{The length (in inches), of the arrows that indicate gene strandedness.}

  
  \item{base.plot.height}{The base height of the standard-sized plots. The plots that include the transcript set will be 1.5x this height.}
  \item{base.plot.width}{The width of the plots.}
  \item{base.plot.units}{The units of measurement for the plot height and width. Default is px, or pixels.}
  \item{plotting.device.params}{Additional parameters to be passed to the plotting device.}

  \item{openPlottingDeviceFunc}{An R function. This option can be used to use plotting devices other than the ones directly supported by JunctionSeq. This must be a function that must have 3 parameters: filename, heightMult, and widthMult. It should open the desired plotting device. For advanced users only.}
  \item{closePlottingDeviceFunc}{An R function. This must be used in conjunction with openPlottingDeviceFunc. For most devices, you can just use the function "dev.off". For advanced users only.}
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
  \item{...}{Additional options to pass to plotting functions, particularly graphical parameters.}
}
