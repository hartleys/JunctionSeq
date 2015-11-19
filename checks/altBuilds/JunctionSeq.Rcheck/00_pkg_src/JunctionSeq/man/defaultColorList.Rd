\name{defaultColorList}
\docType{data}
\alias{JUNCTIONSEQ.DEFAULT.COLOR.LIST}
\alias{defaultColorList}
\alias{junctionSeqColors}

\title{JunctionSeq Color Parameters}
\description{
  This data set gives default color parameters for JunctionSeq. 
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
   
    \code{SIG.FEATURE.BORDER.COLOR, NOSIG.FEATURE.BORDER.COLOR, UNTESTABLE.FEATURE.BORDER.COLOR, EXCLUDED.FEATURE.BORDER.COLOR}: 
      The color for borders of features in the gene and transcript annotation.
    
    \code{SIG.FEATURE.FILL.COLOR, NOSIG.FEATURE.FILL.COLOR, UNTESTABLE.FEATURE.FILL.COLOR, EXCLUDED.FEATURE.FILL.COLOR}
      The colors for exon shading in the gene and transcript annotation.
    
    \code{PLOTTING.LINE.COLORS} 
      Character vector. The colors for each condition, in order. Will be shortened to the number of conditions.
}

\keyword{datasets}

\examples{
data(exampleDataSet,package="JctSeqExData2");

#Set a few alternative colors:
buildAllPlotsForGene(geneID = "ENSRNOG00000009281", jscs, 
                     outfile.prefix = "./oddColors.",
                     colorList = list(SIG.FEATURE.COLOR = "red",
                                      SIG.FEATURE.FILL.COLOR = "green",
                                      NOSIG.FEATURE.FILL.COLOR = "blue"
                                      ));
}