\name{writeBedTrack}
\docType{methods}
\alias{writeExprBedTrack}
\alias{writeSigBedTrack}
\title{
   Write splice junction browser tracks
}

\description{
   This function saves the JunctionSeq results in the form of a set of "bed" files designed for use with the UCSC genome browser.
}

\usage{
  writeExprBedTrack(file, jscs, 
                    trackLine,
                    only.with.sig.gene = TRUE,
                    only.sig = FALSE, 
                    only.testable = TRUE, 
                    group.RGB, 
                    use.score = FALSE, 
                    FDR.threshold = 0.05, 
                    count.digits = 1, 
                    includeGeneID = FALSE,
                    includeLocusID = TRUE,
                    includeGroupID = TRUE)

  writeSigBedTrack(file, 
                   jscs, 
                   trackLine,
                   only.sig = TRUE, 
                   only.testable = TRUE, 
                   sig.RGB = "255,0,0", 
                   nonsig.RGB = "0,0,0", 
                   use.score = TRUE, 
                   FDR.threshold = 0.05, 
                   pval.digits = 4, 
                   includeGeneID = FALSE,
                   includeLocusID = TRUE)

}

\arguments{
  \item{file}{Character string. File path for the output bed file.}

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
  \item{trackLine}{
    The "track line" of the bed file. In other words, the first line of the file. 
    By default JunctionSeq will attempt to automatically generate a reasonable track line.
  }
  \item{only.with.sig.gene}{
    Logical. If TRUE, only genes containing statistically significant results will be included.
  }
  \item{only.sig}{
    Logical. If TRUE, only statistically significant splice junctions will be included.
  }
  \item{only.testable}{
    Logical. If TRUE, only junction loci with sufficiently high expression to be tested will be included.
  }
  \item{sig.RGB}{
    Character string. The RGB color for significant genes. Must be in the format "r,g,b", with each value ranging from 0 to 255.
  }
  \item{nonsig.RGB}{
    Character string. The RGB color for non-significant junctions junctions. Must be in the format "r,g,b", with each value ranging from 0 to 255.
  }
  \item{group.RGB}{
    Character string. The RGB color used for each experimental group. Must be in the format "r,g,b", with each value ranging from 0 to 255. Must have a length equal to the number of experimental condition values.
  }
  \item{use.score}{
    Logical. If TRUE, score each junction based on the p-value.
  }
  \item{FDR.threshold}{
    Numeric. The FDR-adjusted p-value threshold to use to assign statistical significance.
  }
  \item{count.digits}{
    Numeric. The number of digits after the decimal point to include for the mean normalized counts.
  }
  \item{pval.digits}{
    Numeric. The number of digits after the decimal point to include for the p-values.
  }  
  \item{includeGeneID}{
    Logical. If TRUE, include the ID of the gene in the "name" field of each line.
  }
  \item{includeLocusID}{
    Logical. If TRUE, include the ID of the junction locus in the "name" field of each line.
  }
  \item{includeGroupID}{
    Logical. If TRUE, include the ID of the group in the "name" field of each line.
  }

}
