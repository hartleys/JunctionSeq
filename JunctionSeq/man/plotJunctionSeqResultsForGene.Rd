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
    plot.type = c("expr","normCounts","rExpr","rawCounts"), 
    sequencing.type = c("paired-end","single-end"),
    displayTranscripts = FALSE,
    colorList = list(), 
    use.vst = FALSE, use.log = TRUE,
    exon.rescale.factor = 0.3, 
    exonRescaleFunction = c("sqrt","log","linear","34root"),
    label.p.vals = TRUE, 
    plot.lwd = 3, axes.lwd = plot.lwd, 
    anno.lwd = plot.lwd, gene.lwd = plot.lwd / 2,
    par.cex = 1, anno.cex.text = 1,
    anno.cex.axis=anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
    cex.arrows = "auto",
    fit.countbin.names = TRUE, fit.genomic.axis = TRUE, fit.labels = TRUE,
    plot.gene.level.expression = TRUE,
    plot.exon.results, plot.junction.results, plot.novel.junction.results, 
    plot.untestable.results = FALSE, draw.untestable.annotation = TRUE,
    show.strand.arrows = 1, 
    sort.features = TRUE,
    drawCoordinates = TRUE,
    yAxisLabels.inExponentialForm = FALSE,
    title.main, title.ylab, title.ylab.right, 
    graph.margins = c(2, 3, 3, 3),
    GENE.annotation.relative.height = 0.1,
    TX.annotation.relative.height = 0.05,
    CONNECTIONS.relative.height = 0.1,
    SPLICE.annotation.relative.height = 0.1,
    TX.margins = c(0,0.5),
    condition.legend.text = NULL, include.TX.names = TRUE, 
    draw.start.end.sites = TRUE,
    label.chromosome = TRUE,
    splice.junction.drawing.style = c("hyperbola","ellipse",
                                      "triangular","line"),
    draw.nested.SJ = TRUE, merge.exon.parts = TRUE,
    verbose=TRUE, debug.mode = FALSE,
    INTERNAL.VARS = list(),
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
  
  \item{sequencing.type}{
    The type of sequencing used, either "paired-end" or "single-end". This only affects the labelling of the y-axis, and does not affect the actual plots in any way.
  }
  
  \item{displayTranscripts}{Logical. If true, then the full set of annotated transcripts will be displayed below the expression plot (to a maximum of 42 different TX).}
  \item{colorList}{A named list of R colors, setting the colors used for various things. See \code{link{junctionSeqColors}} for more information.}

  \item{use.vst}{Logical. If TRUE, all plots will be scaled via a variance stabilizing transform.}
  \item{use.log}{Logical. If TRUE, all plots will be log-scaled.}
  \item{exon.rescale.factor}{Numeric. Exons will be proportionately scaled-up so that the exonic regions make up this fraction of the horizontal plotting area. If negative, exons and introns will be plotted to a common scale. }
  \item{exonRescaleFunction}{
  Character string. Exonic and intronic regions will be rescaled to be proportional to this transformation of their span. By default the square-root function is used, which shrinks long features and extends short features so that
  they are all still readable and destinguishable against one another. This default option seems to behave well on mammalian genomes. This parameter does nothing if \code{exon.rescale.factor} is negative.
  }
  
  \item{label.p.vals}{Logical. If \code{TRUE}, then statistically significant p-values will be labelled.}

  \item{plot.lwd}{the line width for the plotting lines.}
  \item{axes.lwd}{the line width for the axes.}
  \item{anno.lwd}{the line width for the various other annotation lines.}
  \item{gene.lwd}{the line width used for the gene annotation lines.}
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
  \item{cex.arrows}{
     The font size for the strand-direction arrows in the gene annotation region. 
     The arrows will be sized to equal the dimensions of the letter "M" at this font size.
  }
  
  \item{fit.countbin.names}{Logical. If \code{TRUE}, then splice-junction-locus labels should be rescaled to fit in whatever horizontal space is available.}
  \item{fit.genomic.axis}{Logical. If \code{TRUE}, then the genomic coordinate labels will be auto-scaled down to fit, if needed.}
  \item{fit.labels}{Logical. If \code{TRUE}, then y-axis labels will be auto-scaled down to fit, if needed. Note this only applies to the text labels, not the numeric scales.}

  \item{plot.gene.level.expression}{
    Logical value. If \code{TRUE}, gene-level expression (when applicable) will be plotted beside the sub-element-specific expression in a small seperate plotting box. 
    For the "relative expression" plots the simple mean normalized expression will be plotted (since it doesn't make sense to plot something relative to itself).
  }
  \item{plot.exon.results}{Logical. If \code{TRUE}, plot results for exons. By default everything that was tested will be plotted.}
  \item{plot.junction.results}{Logical. If \code{TRUE}, plot results for splice junctions. By default everything that was tested will be plotted.}
  \item{plot.novel.junction.results}{Logical. If \code{TRUE}, plot results for novel splice junctions. If false, novel splice junctions will be ignored. By default everything that was tested will be plotted.}
  \item{plot.untestable.results}{Logical. If TRUE, plots the expression of splice junctions that had coverage that was too low to be tested.}
  \item{draw.untestable.annotation}{Logical. If TRUE, draws the annotation for splice junctions that had coverage that was too low to be tested.}

  \item{show.strand.arrows}{The number of strand-direction arrows to display. If equal to \code{1} (the default) then the arrow will extend from the end of the gene drawing, if it is greater than \code{1} then arrows will be drawn
    along the gene length. If it is \code{0} or \code{NA} then arrows will not be drawn.
  }

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
  \item{CONNECTIONS.relative.height}{
    The height of the panel that connects the plotting area to the gene annotation area, relative to the height of the plotting area. This panel has the lines that connects the counting bin columns to their actual loci on the gene. By default it is 10 percent.
  }
  \item{SPLICE.annotation.relative.height}{
    The height of the area that shows the splice junction loci, relative to the size of the plotting area.
  }
  \item{TX.margins}{
    A numeric vector of length 2. The size of the blank space between the gene plot and the transcript list and then beneath the transcript list, relative to the size of each transcript line.
  }
  \item{title.main}{
    Character string. Overrides the default main plot title.
  }
  \item{title.ylab}{
    Character string. Overrides the default y-axis label for the left y-axis.
  }
  \item{title.ylab.right}{
    Character string. Overrides the default y-axis label for the right y-axis.
  }
  \item{condition.legend.text}{
    List or named vector of character strings. This optional parameter can be used to assign labels to each condition variable values. It should be a list or named vector with length equal to \code{factor(condition)}. Each element
    should be named with one of the values from \code{factor(condition)}, and should contain the label. They will be listed in this order in the figure legend.
  }
  
  \item{include.TX.names}{
    Logical value. If \code{TRUE}, then for the plots that include the annotated transcript, 
    the transcript names will be listed. The labels will be drawn at half the size 
    of \code{anno.cex.text}.
  }
  \item{draw.start.end.sites}{
    Logical value. If \code{TRUE}, then transcript start/end sites will be marked on the main gene annotation.
  }
  \item{label.chromosome}{
    Logical. If \code{TRUE}, label the chromosome in the left margin. 
    If the text is too long it will be auto-fitted into the available margin.
  }
  \item{splice.junction.drawing.style}{
    The visual style of the splice junctions drawn on the gene annotation. 
    The default uses paired hyperbolas with the ends straightened out. 
    A number of other styles are available.
  }
  \item{draw.nested.SJ}{
    Logical. If \code{TRUE}, overlapping splice junctions will be drawn layered under one another. 
    This can vastly improve readability when there are a large number of overlapping splice junctions. 
    Default is \code{TRUE}.
  }
  \item{merge.exon.parts}{
    Logical. If \code{TRUE}, in the gene annotation plot merge connected exon-fragments and delineate them with dotted lines.
  }
  
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
  \item{debug.mode}{Logical. If \code{TRUE}, print additional debugging information during execution.}
  \item{INTERNAL.VARS}{
    NOT FOR GENERAL USE. Intended only for use by JunctionSeq itself, internally. This is used 
    for passing pre-generated data (when generating many similar plots, for example), and 
    for internally-generated parameters. DO NOT USE.
  }

  \item{...}{Additional options to pass to plotting functions, particularly graphical parameters.}
}
\value{
  This is a side-effecting function, and does not return a value.
}
\examples{
data(exampleDataSet,package="JctSeqData");

plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs);


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

#Make an expression plot for a given gene:
plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs);

#Plot normalized read counts for a given gene:
plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs,
                    plot.type = "normCounts");

#Plot relative expression for a given gene:
plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs,
                    plot.type = "rExpr");

#Plot raw read counts for a given gene:
plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs,
                    plot.type = "rawCounts");

#Same thing, but with isoforms shown:
plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs,
                    plot.type = "rawCounts",
                    displayTranscripts = TRUE);

}
}