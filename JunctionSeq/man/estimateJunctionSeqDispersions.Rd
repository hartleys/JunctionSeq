\name{estimateJunctionSeqDispersions}
\docType{methods}
\alias{estimateJunctionSeqDispersions}
\title{
   JunctionSeq Dispersion Estimation
}
\description{
   This method estimates the sample dispersion for each 
   counting bin (in other words, each splice junction locus). 
   
   This function is called internally by the \code{\link{runJunctionSeqAnalyses}} 
   function, and thus for most purposes users should not need to call this
   function directly. It may be useful to advanced users performing non-standard
   analyses.
}
\usage{
estimateJunctionSeqDispersions( jscs, 
    method.GLM = c(c("advanced","DESeq2-style"),
                   c("simpleML","DEXSeq-v1.8.0-style")),
    test.formula1 = formula(~ sample + countbin + condition : countbin),
    meanCountTestableThreshold="auto", nCores=1, 
    use.multigene.aggregates = FALSE, 
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
  \item{test.formula1}{
    The model formula. Note that this formula is different from the 
    formula used to calculate parameter estimates and effect size.
    This is because the two noise components (gene-level and 
    countbin-level noise) are folded into the sample term. Since
    we only intend to test the condition-countbin interaction,
    we do not need to model the gene-level differential expression.
    
    NOTE: the biological condition to be tested MUST be named 
    "condition".
  }
  \item{meanCountTestableThreshold}{
    "auto" or Numeric value. Features with a total mean normalized count of less than this value will be excluded from the analyses.
    If left as the default ("auto"), then the cutoff threshold will be determined automatically using the DESeq2 independent filtering method.
  }
  \item{nCores}{
    Either an integer or a BiocParallelParam object. 
    Either way, this determines The number of cores to use. Note that multicore functionality may not be available on all platforms. 
    If parallel execution is not available then JunctionSeq will automatically fallback to single-core execution. See the BiocParallel
    package for more information.
  }
  \item{use.multigene.aggregates}{Logical value. Whether to attempt to test "aggregate genes" which consist of multiple genes that overlap with one another. Note that inclusion of aggregate genes may affect the false discovery rate, since by their very nature aggregate genes will often show differential splice junction usage, as the two genes will often be regulated independently.}

  \item{verbose}{A boolean flag indicating whether or not to print progress information during execution. (Default=FALSE)}

}
\value{
  A JunctionSeqCountSet, with dispersion results included.
}

\examples{
\dontrun{
#Full example (from scratch):
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