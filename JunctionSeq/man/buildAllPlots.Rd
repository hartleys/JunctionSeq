\name{buildAllPlots}
\docType{methods}
\alias{buildAllPlots}
\title{
   Create and save a full battery of JunctionSeq expression plots.
}
\description{
   Saves a large battery of plots displaying the analysis results, for the purposes of data visualization. By default it saves a full set of plots for every gene that shows statistical significance and the adjusted-p < 0.01 level.
   Alternatively, it can be supplied with a specific gene list using the \code{gene.list} parameter, and will plot those specific genes.
   
   Note that this function has MANY parameters, allowing the user to tweak the appearance of the plots to suit their particular needs and preferences. 
   Don't be daunted: the default parameters are probably fine for most purposes.
}
\usage{
buildAllPlots(jscs,
    outfile.prefix = "./",
    gene.list = NULL, FDR.threshold = 0.01, max.gene.ct,
    method.selectionCriterion = c("feature-pAdjust", "genewise-pAdjust"),
    use.plotting.device = c("png","CairoPNG","svg",
                            "tiff","cairo_ps","custom"),
    sequencing.type = c("paired-end","single-end"),
    use.vst=FALSE,use.log = TRUE, 
    exon.rescale.factor = 0.3,
    subdirectories.by.type = TRUE,
    ma.plot=TRUE, variance.plot=TRUE,
    with.TX=TRUE,without.TX=TRUE,
    expr.plot=TRUE,normCounts.plot=TRUE,
    rExpr.plot=TRUE,rawCounts.plot=FALSE,
    colorRed.FDR.threshold = FDR.threshold, 
    colorList=list(),
    plot.gene.level.expression = TRUE,
    plot.exon.results, plot.junction.results, plot.novel.junction.results,
    plot.untestable.results = FALSE,
    plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
    gene.lwd = plot.lwd / 2,
    par.cex = 1, anno.cex.text = 1, anno.cex.axis = anno.cex.text, 
    anno.cex.main = anno.cex.text * 1.2,
    drawCoordinates = TRUE,
    yAxisLabels.inExponentialForm = FALSE,
    show.strand.arrows = 1, 
    graph.margins = c(2, 3, 3, 3),
    base.plot.height = 12, base.plot.width = 12, 
    base.plot.units = "in", 
    GENE.annotation.relative.height = 0.1,
    TX.annotation.relative.height = 0.05,
    CONNECTIONS.relative.height = 0.1,
    SPLICE.annotation.relative.height = 0.1,
    TX.margins = c(0,0.5),
    autoscale.height.to.fit.TX.annotation = TRUE,
    autoscale.width.to.fit.bins = 35,
    plotting.device.params = list(), 
    number.plots = FALSE,
    name.files.with.geneID = TRUE,
    condition.legend.text, include.TX.names = TRUE,
    draw.start.end.sites = TRUE,
    openPlottingDeviceFunc, closePlottingDeviceFunc,
    writeHTMLresults = TRUE,
    html.cssFile, html.cssLink, html.imgFileExtension,
    html.plot.height = 90, html.plot.height.units = "vh",
    html.compare.results.list = NULL,
    verbose=TRUE, debug.mode = FALSE,
    INTERNAL.VARS = list(),
    ...)
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
  \item{outfile.prefix}{The prefix file path to save the images to.}

  \item{gene.list}{Character vector. List of genes to plot. Either this variable OR FDR.threshold must be set.}
  \item{FDR.threshold}{If this option is used, genes will be selected for plotting based on the presence of statistically significant junctions. 
        The adjusted-p-value threshold used to determine significance. Only genes containing at least 1 significant feature will be plotted.
  }
  \item{max.gene.ct}{
    Integer or numeric value. This option is usually only used with the FDR.threshold parameter (as opposed to the gene.list parameter).
    This option sets an upper limit to the number of genes to plot. This prevents JunctionSeq from taking too long to complete, or from
    using too much disk space if an enormous number of genes turn out to be significant at the selected significance value. 
    If there are more genes than max.gene.ct, then JunctionSeq will only plot the top max.gene.ct genes.
  }
  
  \item{method.selectionCriterion}{
    Determines the method used to select genes for plotting. If set to "feature-pAdjust", genes will be selected if one or more features show significance.
    If "genewise-pAdjust" is used, gene-wise adjusted p-values will be generated and genes will be filtered on that basis.
  }
  
  \item{use.plotting.device}{The plotting device to use.}
  \item{sequencing.type}{
    The type of sequencing used, either "paired-end" or "single-end". This only affects the labelling of the y-axis, and does not affect the actual plots in any way.
  }
  \item{use.vst}{Logical. If \code{TRUE}, all plots will be scaled via a variance stabilizing transform.}
  \item{use.log}{Logical. If \code{TRUE}, all plots will be log-scaled.}
  \item{exon.rescale.factor}{
    Floating point numeric value. To improve readability the exons drawn in 
    the coordinate annotation are rescaled by default so that they take up 
    30 percent of the x axis. This makes the plots easier to read, as 
    exons are usually much smaller than introns and thus a group of 
    clustered exons can be hard to distinguish when plotted on a simple 
    scale. If this value is set to NA or a value not between 0 and 1 inclusive, 
    then the exons and introns will be drawn on the same scale.
    Note that this function can also take the \code{exonRescaleFunction} parameter, which is
    passed to \code{\link{plotJunctionSeqResultsForGene}}.
  }
  \item{subdirectories.by.type}{
    Logical value. If \code{TRUE}, then subdirectories will be created in the outfile.prefix directory, containing each plot type.
  }

  \item{ma.plot}{if TRUE, generate and save a MA plot. A MA-plot is a plot of fold change versus base mean normalized counts.}
  \item{variance.plot}{if TRUE, generate and save a plot of the dispersion as a function of the base mean.}
  \item{with.TX}{if TRUE, save expression plots with the full transcripts printed}
  \item{without.TX}{if TRUE, save expression plots with only the compiled exons printed. Note that if this and with.TX.plot are both TRUE, both versions will be saved seperately.}
  \item{expr.plot}{if TRUE, save an expression plot of the expression parameter estimates for each splice site, for each condition.}
  \item{normCounts.plot}{if TRUE, save an expression plot of the normalized mean counts for each splice site, for each sample.}
  \item{rExpr.plot}{if TRUE, save an expression plot of the expression parameter estimates, relative to gene-wide expression, for each splice site, for each condition.}
  \item{rawCounts.plot}{if TRUE, save an expression plot of the raw counts for each splice site, for each sample. Note that these will never be VST-transformed, even when use.vst == TRUE.}
  
  \item{colorRed.FDR.threshold}{The adjusted-p-value threshold used to determine whether a feature should be marked as "significant" and colored pink. By default this will be the same as the FDR.threshold.}
  \item{colorList}{A named list of R colors, setting the colors used for various things. See \code{\link{plotJunctionSeqResultsForGene}}}
  
  \item{plot.gene.level.expression}{
    Logical value. If \code{TRUE}, gene-level expression (when applicable) will be plotted beside the sub-element-specific expression in a small seperate plotting box. 
    For the "relative expression" plots the simple mean normalized expression will be plotted (since it doesn't make sense to plot something relative to itself).
  }
  \item{plot.exon.results}{Logical. If \code{TRUE}, plot results for exons. By default everything that was tested will be plotted.}
  \item{plot.junction.results}{Logical. If \code{TRUE}, plot results for splice junctions. By default everything that was tested will be plotted.}
  \item{plot.novel.junction.results}{Logical. If \code{TRUE}, plot results for novel splice junctions. If false, novel splice junctions will be ignored. By default everything that was tested will be plotted.}
  \item{plot.untestable.results}{
    Logical. If \code{TRUE}, plots splice junctions that had coverage that was too low to be tested. 
    Note that, in general, only normCounts and rawCounts plots will have non-\code{NA} values for untestable counting bins.
  }
  
  \item{plot.lwd}{The line width for the plotting lines.}
  \item{axes.lwd}{The line width for the axes.}
  \item{anno.lwd}{The line width for the various other annotation lines.}
  \item{gene.lwd}{The line width used for the gene annotation lines. The default is half the standard line width.}
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

  \item{show.strand.arrows}{The number of strand-direction arrows to display. If equal to \code{1} (the default) then the arrow will extend from the end of the gene drawing, if it is greater than \code{1} then arrows will be drawn
    along the gene length. If it is \code{0} or \code{NA} then arrows will not be drawn.
  }

  \item{graph.margins}{
    Numeric vector of length 4. These margins values used (as if for par("mar")) for the main graph. The lower part of the plot uses the same left and right margins. 
  }
  
  \item{base.plot.height}{
    The base height of the standard-sized plots. 
    Plots that include the full transcript annotation will be expanded by the height of these additional rows.
    See the withTxPlot.height.multiplier parameter, below.
  }
  \item{base.plot.width}{The width of the plots.}
  \item{base.plot.units}{The units of measurement for the plot height and width. Default is px, or pixels.}
  
  
  \item{GENE.annotation.relative.height}{
    The height of the "gene track" displayed underneath the main graph, relative to the height of the main graph. By default it is 20 percent.
  }
  \item{TX.annotation.relative.height}{
    For all plots that draw the annotated-transcript set (when the \code{with.TX} parameter is \code{TRUE}), this sets the height of each transcript, as a fraction of the height of the main graph. By default it is 2.5 percent.
  }
  \item{CONNECTIONS.relative.height}{
    The height of the panel that connects the plotting area to the gene annotation area, relative to the height of the plotting area. This panel has the lines that connects the counting bin columns to their actual loci on the gene. By default it is 10 percent.
  }
  \item{SPLICE.annotation.relative.height}{
    The height of the area that shows the splice junction loci, relative to the size of the plotting area.
  }
  \item{TX.margins}{
    A numeric vector of length 2. The size of the blank space between the gene plot and the transcript list and then beneath the transcript list, relative to the size of each transcript line.
  }
  
  \item{autoscale.height.to.fit.TX.annotation}{
    Logical. Plots that include the full transcript annotation generally need to have a larger height in order to maintain readability. 
    If \code{TRUE}) (the default),
    all plots that include transcripts will be expanded vertically to fit the additional transcripts. This maintains the same appearance and
    aspect ratio of the main graph area, but means that the height of the plot will differ between genes when TX are included.
    This parameter can be used to override that behavior if a specific figure size is desired. 
    If \code{FALSE}, then the height of the entire 
    output image will always be equal to \code{base.plot.height}.
  }
  \item{autoscale.width.to.fit.bins}{
    Integer value. JunctionSeq will automatically go to great lengths to autofit 
    the data in a readable way. By default, any plots that have more than 35 
    plotting columns
    will be widened linearly to fit the excess columns. This parameter can be 
    used to change that value, or turn it off entirely by setting this 
    parameter to \code{NA}.
  }
  \item{condition.legend.text}{
    List or named vector of character strings. This optional parameter can be used to assign labels to each condition variable values. It should be a list or named vector with length equal to \code{factor(condition)}. Each element
    should be named with one of the values from \code{factor(condition)}, and should contain the label. They will be listed in this order in the figure legend.
  }
  
  \item{include.TX.names}{
    Logical value. If \code{TRUE}, then for the plots that include the annotated transcript, the transcript names will be listed. The labels will be drawn at half the size of \code{anno.cex.text}.
  }
  
  \item{plotting.device.params}{Additional parameters to be passed to the plotting device.}
  \item{number.plots}{
     Whether to number each gene in the image names, based on either 
     the order they appear in the input gene.list, or in order of 
     ascending p-values. 
  }
  \item{name.files.with.geneID}{Whether to use the geneID (rather than gene name) for naming the files.}

  \item{draw.start.end.sites}{
    Logical value. If \code{TRUE}, then transcript start/end sites will be marked on the main gene annotation.
  }
  \item{openPlottingDeviceFunc}{An R function. This option can be used to use plotting devices other than the ones directly supported by JunctionSeq. This must be a function that must have 3 parameters: filename, heightMult, and widthMult. It should open the desired plotting device. For advanced users only.}
  \item{closePlottingDeviceFunc}{An R function. This must be used in conjunction with openPlottingDeviceFunc. For most devices, you can just use the function "dev.off". For advanced users only.}

  \item{writeHTMLresults}{If TRUE, write an index html file to present the results in a navigable way.}
  
  \item{html.cssFile}{
    Optional: specify a css file to use. Copies the entire contents of the supplied file into the page directory and links to it with relative links.
  }
  \item{html.cssLink}{
      Optional: specify an external css file to use. This can be an absolute or relative link.
  }
  \item{html.imgFileExtension}{
      The file extension of the image files. This is only needed if you are using a custom device. If you are using one of the default devices, it will autodetect the file extension.
  }
  \item{html.plot.height}{
      Numeric. The base height of the plot, for the plots without TX annotation. The default is 90.
  }
  \item{html.plot.height.units}{
      The units used for the \code{html.plot.height} parameter. The default is "vh", which sets the height relative to the available max height.
  }
  \item{html.compare.results.list}{
      Named list of character strings. (Advanced) Optional parameter that allows you to cross-link multiple analyses for easy navigation 
      between analysis for specific genes of interest.
      In order to create such cross-linking, you will need to run \code{builtAllPlots} separately for each analysis. 
      The \code{outfile.prefix} for each
      run must be a sub-directory of the same parent directory. The \code{html.compare.results.list} must be a named list of these subdirectories.
      \code{names(html.compare.results.list)} must be the title of each analysis as you want it to appear in the navigation links.
      Note: This parameter is incompatible with the \code{number.plots} option.
  }
  
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
  \item{debug.mode}{if TRUE, send even more debugging and progress messages to the console / stdout.}
  
  \item{INTERNAL.VARS}{
    NOT FOR GENERAL USE. Intended only for use by JunctionSeq itself, internally. This is used 
    for passing pre-generated data (when generating many similar plots, for example), and 
    for internally-generated parameters. DO NOT USE.
  }

  \item{...}{Additional options to pass to \code{\link{buildAllPlotsForGene}}, \code{\link{plotJunctionSeqResultsForGene}}, or graphical parameters passed to plotting functions.}
}

\value{
  This is a side-effecting function, and does not return a value.
}

\examples{

data(exampleDataSet,package="JctSeqExData2");
buildAllPlots(jscs);

\dontrun{
########################################
#Set up example data:
decoder.file <- system.file(
                  "extdata/annoFiles/decoder.bySample.txt",
                  package="JctSeqExData2");
decoder <- read.table(decoder.file,
                  header=TRUE,
                  stringsAsFactors=FALSE);
gff.file <- system.file(
            "extdata/cts/withNovel.forJunctionSeq.gff.gz",
            package="JctSeqExData2");
countFiles <- system.file(paste0("extdata/cts/",
     decoder$sample.ID,
     "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
     package="JctSeqExData2");
########################################
#Run example analysis:
jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
           sample.names = decoder$sample.ID,
           condition=factor(decoder$group.ID),
           flat.gff.file = gff.file,
           analysis.type = "junctionsAndExons"
);
########################################

#Generate all plots and the html index
#   Save them as pngs to the current directory:
buildAllPlots(jscs);

}
}