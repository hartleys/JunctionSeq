\name{readAnnotationData}
\docType{methods}
\alias{readAnnotationData}
\title{
   Read junctionSeq annotation files produced by QoRTs.
}
\description{
   This function reads the "flattened" gff annotation file created by QoRTs. 
   This annotation file contains all the gene, transcript, exon, and junction ID's and their loci.
   In general this function is not used by the end-user, but is called internally by 
   \code{\link{runJunctionSeqAnalyses}} or \code{\link{readJunctionSeqCounts}}.
}
\usage{
   readAnnotationData(flat.gff.file)
}

\arguments{
  \item{flat.gff.file}{
    Character string. The filename of the "flat" gff annotation file. The file may be gzip-compressed. 
    This "flat" gff file must be produced by the QoRTs jar utility using the makeFlatGtf or 
    mergeNovelSplices functions (depending on whether inclusion of novel splice junctions is desired).
  }
}

\value{
  A \code{data.frame} object, containing the annotation information from the flat gff file.
}

\examples{
gff.file <- system.file("extdata/cts/withNovel.forJunctionSeq.gff.gz",
                        package="JctSeqData");

#Parse the GFF file:
annoData <- readAnnotationData(gff.file);
head(annoData);
}


