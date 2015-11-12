\name{fitJunctionSeqDispersionFunction}
\docType{methods}
\alias{fitJunctionSeqDispersionFunction}
\title{
   Fit Shared Dispersion Function
}
\description{
   Fit dispersion function to share dispersion information between 
   features across the genome.
   
   This function is called internally by the \code{\link{runJunctionSeqAnalyses}} 
   function, and thus for most purposes users should not need to call this
   function directly. It may be useful to advanced users performing non-standard
   analyses.
}
\usage{
fitJunctionSeqDispersionFunction(jscs, 
    method.GLM = c(c("advanced","DESeq2-style"), c("simpleML","DEXSeq-v1.8.0-style")),
    method.dispFit = c("parametric", "local", "mean"), 
    method.dispFinal = c("shrink","max","fitted","noShare"),
    fitDispersionsForExonsAndJunctionsSeparately = TRUE, 
    verbose = TRUE)
}
\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually initially created by 
    \code{\link{readJunctionSeqCounts}}. Size factors must be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}.
  }
  \item{method.GLM}{
    Character string. Can be used to apply alternative methodologies or implementations. Intended for advanced users who have strong opinions about the underlying statistical methodologies.
    
    The default is "advanced" or, equivalently, "DESeq2-style". This uses the 
    dispersion estimation methodology used by DESeq2 and DEXSeq v1.12.0 or higher to generate the initial (feature-specific) dispersion estimates. The alternative 
    method is "simpleML" or, equivalently, "DEXSeq-v1.8.0-style". This uses a simpler maximum-likelihood-based method used by the original DESeq and
    by DEXSeq v1.8.0 or less.
  }
  \item{method.dispFit}{
    Character string. Can be used to apply alternative methodologies or implementations. Intended for advanced users who have strong opinions about the underlying statistical methodologies.
    
    Determines the method used to generated "fitted" dispersion estimates. One of "parametric" (the default), "local", or "mean". See the DESeq2 documentation for more information.
  }
  \item{method.dispFinal}{
    Character string. Can be used to apply alternative methodologies or implementations. Intended for advanced users who have strong opinions about the underlying statistical methodologies.

    Determines the method used to arrive at a "final" dispersion estimate. The default, "shrink" uses the maximum a posteriori estimate, combining information from both
    the fitted and feature-specific dispersion estimates. This is the method used by DESeq2 and DEXSeq v1.12.0 and above.
  }
  \item{fitDispersionsForExonsAndJunctionsSeparately}{
    When running a "junctionsAndExons" type analysis in which both exons and splice junctions are being tested simultaniously, this parameter
    determines whether a single fitted dispersion model should be fitted for both exons and splice junctions, or if separate 
    fitted dispersions should be calculated for each. By default the dispersions are run separately.
  }
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
  \item{...}{If using the depreciated fitDispersionFunction command, use the same syntax as above.}
}
\value{
  A JunctionSeqCountSet, with dispersion results included.
}


\examples{
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