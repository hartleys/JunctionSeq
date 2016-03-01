\name{writeCompleteResults}
\docType{methods}
\alias{writeCompleteResults}
\title{
   Produce output data files, given annotation files and DEXSeq 
   exonCountSet object and DEXSeq results data.
}
\description{
   This function takes the raw DEXSeq results and merges in feature 
   annotations, as well as calculating and merging in a number of 
   different normalized and fitted values for each level of the 
   condition variable.
}
\usage{
writeCompleteResults(jscs, outfile.prefix, 
                        gzip.output = TRUE,
                        FDR.threshold = 0.01,
                        save.allGenes = TRUE, save.sigGenes = TRUE,
                        save.fit = FALSE, save.VST = FALSE,
                        save.bedTracks = TRUE,
                        save.jscs = FALSE,
                        bedtrack.format = c("BED", "GTF", "GFF3"),
                        verbose = TRUE)
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
  \item{outfile.prefix}{A string indicating the filename prefix where output files should be saved.}
  \item{gzip.output}{Logical. If \code{TRUE}, then all ".txt" text files should be gzip-compressed to save space.}
  \item{FDR.threshold}{The adjusted-p-value threshold used to determine statistical significance.}
  \item{save.allGenes}{Logical. Whether to save files containing data for all genes.}
  \item{save.sigGenes}{Logical. Whether to save a separate set of files containing data for only the significant genes. If this and save.allGenes are both true then two sets of files will be generated.}
  \item{save.fit}{Logical. Whether to save model fit data.}
  \item{save.VST}{Logical. Whether to save VST-transformed data.}
  \item{save.bedTracks}{Logical. Whether to save "bed" junction coverage tracks.}
  \item{save.jscs}{Logical. Whether to the entire \code{JunctionSeqCountSet} object. Default is FALSE.}  
  \item{bedtrack.format}{
    Character string. The format to use for the browser tracks.
  }
  \item{verbose}{A boolean flag indicating whether or not to print progress information during execution. (Default=FALSE)}
}
\details{
  Saves a wide variety of data from the analyses.
}

\value{
  This is a side-effecting function, and does not return a value.
}
\examples{
data(exampleDataSet,package="JctSeqData");
#Write results tables and browser track files:
writeCompleteResults(jscs, outfile.prefix = "./results.");

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

#Write results tables and browser track files:
writeCompleteResults(jscs, outfile.prefix = "./results.");

}
}