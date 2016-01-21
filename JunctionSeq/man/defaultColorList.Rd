\name{defaultColorList}
\docType{data}
\alias{JUNCTIONSEQ.DEFAULT.COLOR.LIST}
\alias{defaultColorList}
\alias{junctionSeqColors}

\title{JunctionSeq Color Parameters}
\description{
  This data set specifies some of the default color and style parameters 
  for JunctionSeq. 
  Any of these parameters can be overridden by using the \code{colorList} 
  parameter in \code{\link{buildAllPlots}}, \code{\link{buildAllPlotsForGene}}, 
  or \code{\link{plotJunctionSeqResultsForGene}}.
}
\usage{JUNCTIONSEQ.DEFAULT.COLOR.LIST}
\format{A list.}
\details{
  The parameter \code{colorList} must be a named list. Any elements with the names listed below will override the default JunctionSeq colors.

    \code{SIG.FEATURE.COLOR} The color of vertical lines in the plotting panel, for significant features.  
    
    \code{NOSIG.FEATURE.COLOR} The color of vertical lines in the plotting panel, for tested but non-significant features.  
    
    \code{UNTESTABLE.FEATURE.COLOR} The color of vertical lines in the plotting panel, for untestable features. 
    
    \code{EXCLUDED.FEATURE.COLOR} The color of vertical lines in the plotting panel, for features that were not part of the analysis (for example: splice junctions in an exon-based analysis). 
    
    \code{SIG.VERTLINE.COLOR, NOSIG.VERTLINE.COLOR, UNTESTABLE.VERTLINE.COLOR}: The colors of the vertical dotted lines in the plotting panel.
   
    \code{SIG.FEATURE.BORDER.COLOR}, \code{NOSIG.FEATURE.BORDER.COLOR}, 
    \code{UNTESTABLE.FEATURE.BORDER.COLOR}, \code{EXCLUDED.FEATURE.BORDER.COLOR}: 
      The color for borders of features in the gene and transcript annotation.
    
    \code{SIG.FEATURE.FILL.COLOR}, \code{NOSIG.FEATURE.FILL.COLOR}, 
    \code{UNTESTABLE.FEATURE.FILL.COLOR}, \code{EXCLUDED.FEATURE.FILL.COLOR}
      The colors for exon shading in the gene and transcript annotation.
    
    \code{KNOWN.SPLICE.LTY}, \code{NOVEL.SPLICE.LTY}: 
    These set the "lty" parameter for the known and novel splice junction diagrams. 
    These must be specified as one of the character specifications. In other words, 
    either a descriptive specification like "solid","dotted",etc., or by an even number
    of nonzero hexadecimal digits. See the section "Line Type Specification" in the
    graphical parameters help (link{par}). You CANNOT use the single-digit index
    specification, as this will fail when used in vector form with the others 
    (technically you could use them, but only if you reset ALL lty parameters).
    
    \code{EXON.CONNECTION.LTY}, \code{NOVEL.SPLICE.CONNECTION.LTY}, 
    \code{KNOWN.SPLICE.CONNECTION.LTY}: The "lty" parameters for the lines
    that connect the features on the gene schematic diagram to the plotting panel.
    These have the same requirements as the previous set of lty parameters.
    
    \code{PLOTTING.LINE.COLORS}
      Character vector. The colors for each condition, in order. Will be shortened to the number of conditions.
}

\value{
  See above.
}

\keyword{datasets}

\examples{
data(exampleDataSet,package="JctSeqData");

#Set a few alternative colors:
buildAllPlotsForGene(geneID = "ENSRNOG00000009281", jscs, 
                     outfile.prefix = "./oddColors.",
                     colorList = list(SIG.FEATURE.COLOR = "red",
                                      SIG.FEATURE.FILL.COLOR = "green",
                                      NOSIG.FEATURE.FILL.COLOR = "blue"
                                      ));
}