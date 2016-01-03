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
     plot.exons = TRUE, plot.junctions = TRUE, plot.novel.junctions = TRUE,
     group.RGB, 
     use.score = FALSE, 
     FDR.threshold = 0.05, 
     count.digits = 1, 
     includeGeneID = FALSE,
     includeLocusID = TRUE,
     includeGroupID = TRUE,
     output.format = c("BED", "GTF", "GFF3"),
     use.gzip = TRUE,
     verbose = TRUE)

writeSigBedTrack(file, 
     jscs, 
     trackLine,
     only.sig = TRUE, 
     only.testable = TRUE, 
     plot.exons = TRUE, plot.junctions = TRUE, plot.novel.junctions = TRUE,
     sig.RGB = "255,0,0", 
     nonsig.RGB = "0,0,0", 
     use.score = TRUE, 
     FDR.threshold = 0.05, 
     pval.digits = 4, 
     includeGeneID = FALSE,
     includeLocusID = TRUE,
     output.format = c("BED", "GTF", "GFF3"),
     use.gzip = TRUE,
     verbose = TRUE)

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
    be performed by \code{\link{testForDiffUsage}}. Effect sizes and parameter
    estimates must be created via \code{\link{estimateEffectSizes}}.
  }
  \item{trackLine}{
    The "track line" of the bed file. In other words, the first line of the file. 
    By default JunctionSeq will attempt to automatically generate a reasonable track line.
  }
  \item{only.with.sig.gene}{
    Logical. If \code{TRUE}, only genes containing statistically significant results will be included.
  }
  \item{only.sig}{
    Logical. If \code{TRUE}, only statistically significant loci will be included.
  }
  \item{only.testable}{
    Logical. If \code{TRUE}, only loci with sufficiently high expression to be tested will be included.
  }
  \item{plot.exons}{
    Logical. If \code{TRUE}, exons will be plotted.
  }
  \item{plot.junctions}{
    Logical. If \code{TRUE}, splice junctions will be plotted.
  }
  \item{plot.novel.junctions}{
    Logical. If \code{TRUE}, novel splice junctions will be plotted (if plot.junctions is also \code{TRUE}).
  }
  \item{sig.RGB}{
    Character string. The RGB color for significant genes. Must be in the format "r,g,b", with each value ranging from 0 to 255.
  }
  \item{nonsig.RGB}{
    Character string. The RGB color for non-significant loci. Must be in the format "r,g,b", with each value ranging from 0 to 255.
  }
  \item{group.RGB}{
    Character string. The RGB color used for each experimental group. Must be in the format "r,g,b", with each value ranging from 0 to 255. Must have a length equal to the number of experimental condition values.
  }
  \item{use.score}{
    Logical. If TRUE, score each locus based on the p-value.
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
    Logical. If TRUE, include the ID of the locus in the "name" field of each line.
  }
  \item{includeGroupID}{
    Logical. If TRUE, include the ID of the group in the "name" field of each line.
  }
  \item{output.format}{
    Character string. The format to use.
  }
  \item{use.gzip}{
    Logical. Whether or not to gzip the bed file.
  }
  \item{verbose}{
    Logical. if TRUE, output debugging/progress information.
  }
}
\value{
  This is a side-effecting function, and does not return a value.
}

\examples{
data(exampleDataSet,package="JctSeqData");
writeExprBedTrack("test.exonCoverage.bed.gz", jscs, 
                  plot.exons = TRUE, plot.junctions = FALSE)

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

#Exon coverage:
writeExprBedTrack("test.exonCoverage.bed.gz", jscs, 
                  plot.exons = TRUE, plot.junctions = FALSE)
#Junction coverage:
writeExprBedTrack("test.jctCoverage.bed.gz", jscs, 
                  plot.exons = FALSE, plot.junctions = TRUE)
#Both Exon and Junction coverage:
writeExprBedTrack("test.featureCoverage.bed.gz", jscs)

#p-values of significant features:
writeSigBedTrack("test.pvals.bed.gz", jscs)

}
}