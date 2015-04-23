\docType{class}
\name{JunctionSeqCountSet-class}
\alias{JunctionSeqCountSet}
\alias{JunctionSeqCountSet-class}
\alias{newJunctionSeqCountSet}
\title{JunctionSeqCountSet object and constructors}
\usage{

  newJunctionSeqCountSet( countData, 
                          geneCountData, 
                          design, 
                          geneIDs, 
                          countbinIDs, 
                          featureIntervals=NULL, 
                          transcripts=NULL)

}
\arguments{
  \item{countData}{
    A matrix of junction-level count data of non-negative integer values. 
    The rows correspond to counts for each splice-junction counting bin, 
    the columns correspond to samples. Note that biological 
    replicates should each get their own column, while the 
    counts of technical replicates (i.e., several sequencing 
    runs/lanes from the same sample) should be summed up 
    into a single column.
  }
  \item{geneCountData}{
    A matrix of gene-level count data of non-negative integer values. 
    The rows correspond to counts for each gene, 
    the columns correspond to samples. Note that biological 
    replicates should each get their own column, while the 
    counts of technical replicates (i.e., several sequencing 
    runs/lanes from the same sample) should be summed up 
    into a single column. Must have the same dimensions as countData.
  }
  \item{design}{
    A data frame consisting of all factors to be included in the analysis.
    All columns should be factors. Each column should represent a different variable,
    each row should represent a different sample. The number of rows must equal the
    number of columns in geneCountData and countData.
  }
  \item{geneIDs}{
    A character vector of gene indentifiers for each splice junction. The length must
    equal the number of rows in countData.
  }
  \item{countbinIDs}{
    A character vector of splice-junction-locus indentifiers for each splice junction. The length must
    equal the number of rows in countData.
  }
  
  \item{featureIntervals}{
    Optional. A data.frame with 4 columns: "chr", "start", "end", and "strand". chr and strand should 
    be character vectors or factors, start and end must be integers.  
  }
  \item{transcripts}{
    Optional. Character vector listing the transcripts that each splice junction belongs to. Some junctions
    may belong to more than one transcripts. In this case, transcripts should be separated with the "+" character.
  }
}

\value{
  A JunctionSeqCountSet object. Additional data can be added to the 
  
  The constructor function above SHOULD NOT BE USED in normal operation. Instead you should use the \code{\link{readJunctionSeqCounts}} function.
}

\seealso{
  \code{\link{readJunctionSeqCounts}}
}

\description{
  The \code{JunctionSeqCountSet} is a subclass of the bioconductor class \code{eSet},
  designed to contain all data related to a JunctionSeq analysis.
}