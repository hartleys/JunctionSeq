\name{buildAllPlotsForGene}
\docType{methods}
\alias{buildAllPlotsForGene}
\title{
   Create and save one or more JunctionSeq expression plots.
}
\description{
   Generates and saves one or more plots, displaying counts or averages for all counting bins across one particular gene. The parameters \code{expr.plot}, \code{normCounts.plot}, \code{rExpr.plot}, and \code{rawCounts.plot} determine
   which plot types are to be generated, and the parameters \code{with.TX} and \code{without.TX} determines whether these plots should include or not include the full transcript information, or if separate plots should
   be generated with and without the full transcript information.
   
   Note that this function has MANY parameters, allowing the user to tweak the behavior and appearance of the plots to suit their particular needs and preferences. Don't be daunted: the default parameters are probably fine for most purposes.
}
\usage{

buildAllPlotsForGene(geneID, jscs,
                     outfile.prefix = "./", 
                     use.plotting.device = "png", 
                     use.vst=FALSE, use.log = TRUE, truncateBelowOne = TRUE,
                     exon.rescale.factor = 0.3,
                     with.TX=TRUE,without.TX=TRUE,
                     expr.plot=TRUE,normCounts.plot=TRUE,
                     rExpr.plot=FALSE,rawCounts.plot=FALSE,
                     colorRed.FDR.threshold = 0.01, 
                     color=NULL,
                     plot.gene.level.expression = NULL,
                     plot.exon.results, plot.junction.results, plot.novel.junction.results,
                     plot.untestable.results = FALSE,
                     plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                     par.cex = 1, anno.cex.text = 1,
                     anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                     drawCoordinates = TRUE,
                     yAxisLabels.inExponentialForm = FALSE,
                     show.strand.arrows = 10, arrows.length = 0.125,
                     graph.margins = c(2,3,3,2),
                     base.plot.height = 22.222, base.plot.width = 22.222, 
                     base.plot.units = "in", 
                     GENE.annotation.relative.height = 0.2, TX.annotation.relative.height = 0.025,
                     autoscale.height.to.fit.TX.annotation = TRUE,
                     autoscale.width.to.fit.bins = 35,
                     plotting.device.params = list(),
                     condition.legend.text, include.TX.names = TRUE,
                     draw.start.end.sites = TRUE,
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
    be performed by \code{\link{testForDiffUsage}}. Effect sizes and parameter
    estimates must be created via \code{\link{estimateEffectSizes}}.
  }
  \item{outfile.prefix}{
    Character string or vector. Sets the prefix file path where image files should be saved. Optionally it can be a vector of strings, assigning a different file prefix to each plot.
  }

  \item{use.plotting.device}{The plotting device to use.}
  \item{use.vst}{Logical. If \code{TRUE}, all plots will be scaled via a variance stabilizing transform.}
  \item{use.log}{Logical. If \code{TRUE}, all plots will be log-scaled.}
  \item{truncateBelowOne}{Logical. If \code{TRUE}, all values below 1 will be linear-scaled. If use.log is FALSE, this does nothing.}
  \item{exon.rescale.factor}{Numeric. Exons will be proportionately scaled-up so that the exonic regions make up this fraction of the horizontal plotting area. If negative, exons and introns will be plotted to a common scale. }
  
  \item{with.TX}{if \code{TRUE}, save expression plots with the full transcripts printed}
  \item{without.TX}{if \code{TRUE}, save expression plots with only the compiled exons printed. Note that if this and with.TX.plot are both TRUE, both versions will be saved seperately.}
  \item{expr.plot}{if \code{TRUE}, save an expression plot of the expression parameter estimates for each splice site, for each condition.}
  \item{normCounts.plot}{if \code{TRUE}, save an expression plot of the normalized mean counts for each splice site, for each sample.}
  \item{rExpr.plot}{if \code{TRUE}, save an expression plot of the expression parameter estimates, relative to gene-wide expression, for each splice site, for each condition.}
  \item{rawCounts.plot}{if \code{TRUE}, save an expression plot of the raw counts for each splice site, for each sample. Note that these will never be VST-transformed, even when use.vst == TRUE.}
  
  \item{colorRed.FDR.threshold}{The adjusted-p-value threshold used to determine whether a feature should be marked as "significant" and colored pink. By default this will be the same as the FDR.threshold.}
  \item{color}{A vector of R colors, named for each possible value of condition. By default, it will attempt to choose reasonable colors for each condition.}
  
  \item{plot.gene.level.expression}{
    Logical value. If \code{TRUE}, gene-level expression (when applicable) will be plotted beside the sub-element-specific expression in a small seperate plotting box. 
    The gene-level expression will NOT be plotted for the "relative expression" plots (since you can't meaningfully plot gene-level expression relative to itself).
  }
  \item{plot.exon.results}{Logical. If \code{TRUE}, plot results for exons. By default everything that was tested will be plotted.}
  \item{plot.junction.results}{Logical. If \code{TRUE}, plot results for splice junctions. By default everything that was tested will be plotted.}
  \item{plot.novel.junction.results}{Logical. If \code{TRUE}, plot results for novel splice junctions. If false, novel splice junctions will be ignored. By default everything that was tested will be plotted.}
  \item{plot.untestable.results}{Logical. If \code{TRUE}, plots splice junctions that had coverage that was too low to be tested.}
  
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
  \item{yAxisLabels.inExponentialForm}{Logical. If \code{TRUE}, then the y-axis will be labelled in exponential form.}

  \item{show.strand.arrows}{The number of strand-direction arrows to display.}
  \item{arrows.length}{The length of the strand-direction arrows, in inches.}
  
  \item{graph.margins}{
    Numeric vector of length 4. These margins values used (as if for par("mar")) for the main graph. The lower part of the plot uses the same left and right margins. 
  }
  
  \item{base.plot.height}{
    The base height of the standard-sized plots. 
    Plots that include the full transcript annotation will be expanded by the height of these additional rows.
    See the withTxPlot.height.multiplier parameter, below.
  }
  \item{base.plot.width}{The base width of the plots (plots with a large number of features may be scaled up, see parameter \code{autoscale.width.to.fit.bins}).}
  \item{base.plot.units}{The units of measurement for the plot height and width. Default is px, or pixels.}
  \item{GENE.annotation.relative.height}{
    The height of the "gene track" displayed underneath the main graph, relative to the height of the main graph. By default it is 20 percent.
  }
  \item{TX.annotation.relative.height}{
    For all plots that draw the annotated-transcript set (when the \code{with.TX} parameter is \code{TRUE}), this sets the height of each transcript, as a fraction of the height of the main graph. By default it is 2.5 percent.
  }
  \item{autoscale.height.to.fit.TX.annotation}{
    Plots that include the full transcript annotation generally need to have a larger height in order to maintain readability. By default, 
    all plots that include transcripts will be expanded vertically by the height of the additional transcripts. This maintains the same appearance and
    aspect ratio of the main graph, but also means that the height of the plot will differ between genes.
    This parameter can be used to override that behavior if a specific figure size is desired. If \code{TRUE}, the \code{base.plot.height} will be used
    as the height of the plot, regardless of how many transcripts are included.
  }
  \item{autoscale.width.to.fit.bins}{
    Integer value. JunctionSeq will automatically go to great lengths to autofit the data in a readable way. By default, any plots that have more than 35 plotting columns
    will be widened linearly to fit the excess columns. This parameter can be used to change that value, or turn it off entirely by setting this parameter to \code{NA}.
  }
  
  \item{plotting.device.params}{Additional parameters to be passed to the plotting device.}


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

  \item{openPlottingDeviceFunc}{An R function. This option can be used to use plotting devices other than the ones directly supported by JunctionSeq. This must be a function that must have 3 parameters: filename, heightMult, and widthMult. It should open the desired plotting device. For advanced users only.}
  \item{closePlottingDeviceFunc}{An R function. This must be used in conjunction with openPlottingDeviceFunc. For most devices, you can just use the function "dev.off". For advanced users only.}
    
  \item{verbose}{if \code{TRUE}, send debugging and progress messages to the console / stdout.}
  \item{...}{Additional options to pass to plotting functions, particularly graphical parameters.}
}
