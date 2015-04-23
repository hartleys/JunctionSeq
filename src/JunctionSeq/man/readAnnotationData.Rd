\name{readAnnotationData}
\docType{methods}
\alias{readAnnotationData}
\title{
   Read junctionSeq annotation files produced by QoRTs.
}
\description{
   This function reads the "flattened" gff annotation file created by QoRTs. This annotation file contains all the gene, transcript, exon, and junction ID's and their loci.
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
