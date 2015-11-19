\name{estimateJunctionSeqSizeFactors}
\docType{methods}
\alias{estimateJunctionSeqSizeFactors}
\alias{writeSizeFactors}
\title{
   Estimate Size Factors
}
\description{
   Estimate size factors, which are scaling factors used as "offsets" 
   by the statistical model to make the different samples comparable.
   This is necessary because the different samples may have been
   sequenced to slightly different depths. Additionally, the presence
   of differentially expressed genes may cause the apparent depth of
   many genes to appear different.
   
   This function uses the "geometric" size factor normalization method,
   which is identical to the one used by DESeq, DESeq2, DEXSeq, and 
   the default method used by CuffDiff.
   
   This function is called internally by the \code{\link{runJunctionSeqAnalyses}} 
   function, and thus for most purposes users should not need to call this
   function directly. It may be useful to advanced users performing non-standard
   analyses.
}
\usage{
  estimateJunctionSeqSizeFactors(jscs, 
          method.sizeFactors = c("byGenes","byCountbins"), 
          replicateDEXSeqBehavior.useRawBaseMean = FALSE, 
          calcAltSF = TRUE, 
          verbose = FALSE);
  
  writeSizeFactors(jscs, file);
}
\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually initially created by 
    \code{\link{readJunctionSeqCounts}}. Size factors must be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}.
  }
  \item{method.sizeFactors}{
    Character string. Can be used to apply alternative methodologies or implementations. Intended for advanced users who have strong opinions about the underlying statistical methodologies.

    Determines the method used to calculate normalization size factors. By default JunctionSeq uses gene-level expression. As an alternative, feature-level counts can be used as
    they are in DEXSeq. In practice the difference is almost always negligible.
  }
  \item{replicateDEXSeqBehavior.useRawBaseMean}{
    USED ONLY FOR INTERNAL TESTING! NOT INTENDED FOR ACTUAL USE! 
    
    This variable activates an alternative mode in which a (very minor) bug in DEXSeq v1.14.0 and earlier is replicated. If \code{TRUE}, the 
    baseMean and baseVar variables will be computed using raw counts rather than normalized counts. 
    This is used for internal tests in which DEXSeq functionality is replicated precisely and the results are compared against equivalent DEXSeq results.
    Without this option the results would differ slightly (generally by less than 1 hundreth of a percent). 
    
    USED ONLY FOR INTERNAL TESTING! NOT INTENDED FOR ACTUAL USE! 
  }
  \item{calcAltSF}{
    Logical. Determines whether both types of size factor calculations should be generated, and placed in the \code{jscs@altSizeFactors} slot.
  }
  \item{verbose}{if \code{TRUE}, send debugging and progress messages to the console / stdout.}
  \item{file}{
    A file path to write the size factor table.
  }
  \item{...}{If using the (depreciated) estimateSizeFactors command, use the same syntax as above.}
}
\value{
  A JunctionSeqCountSet, with size factors included.
}

\examples{

data(exampleDataSet,package="JctSeqExData2");
jscs <- estimateJunctionSeqSizeFactors(jscs);

\dontrun{
########################################
#Set up example data:
decoder.file <- system.file("extdata/annoFiles/decoder.bySample.txt",package="JctSeqExData2");
decoder <- read.table(decoder.file,
                      header=T,
                      stringsAsFactors=F);
gff.file <- system.file("extdata/cts/withNovel.forJunctionSeq.gff.gz",package="JctSeqExData2");
countFiles <- system.file(paste0("extdata/cts/",
                 decoder$sample.ID,
                 "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
                 package="JctSeqExData2");
########################################
#Advanced Analysis:

#Make a "design" dataframe:
design <- data.frame(condition = factor(decoder$group.ID));
#Read the QoRTs counts.
jscs = readJunctionSeqCounts(countfiles = countFiles,
           samplenames = decoder$sample.ID,
           design = design,
           flat.gff.file = gff.file
);
#Generate the size factors and load them into the JunctionSeqCountSet:
jscs <- estimateJunctionSeqSizeFactors(jscs);
#Estimate feature-specific dispersions:
jscs <- estimateJunctionSeqDispersions(jscs);
#Fit dispersion function and estimate MAP dispersion:
jscs <- fitJunctionSeqDispersionFunction(jscs);
#Test for differential usage:
jscs <- testForDiffUsage(jscs);
#Estimate effect sizes and expression estimates:
jscs <- estimateEffectSizes( jscs);

}
}