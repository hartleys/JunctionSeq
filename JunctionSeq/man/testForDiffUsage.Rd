\name{testForDiffUsage}
\docType{methods}
\alias{testForDiffUsage}
\title{
   Test Junctions for Differential Junction Usage
}
\description{
   This function runs the hypothesis tests for differential junction
   usage.
   
   This function is called internally by the \code{\link{runJunctionSeqAnalyses}} 
   function, and thus for most purposes users should not need to call this
   function directly. It may be useful to advanced users performing non-standard
   analyses.
}
\usage{
testForDiffUsage( jscs, 
    test.formula0 = formula(~ sample + countbin), 
    test.formula1 = formula(~ sample + countbin + condition : countbin),
    method.GLM = c(c("advanced","DESeq2-style"), 
                   c("simpleML","DEXSeq-v1.8.0-style")),
    dispColumn="dispersion", nCores=1, 
    keep.hypothesisTest.fit = FALSE,
    meanCountTestableThreshold = "auto",
    optimizeFilteringForAlpha = 0.01,
    method.cooksFilter = TRUE,
    cooksCutoff,
    pAdjustMethod = "BH",
    verbose = TRUE)
}

\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually initially created by 
    \code{\link{readJunctionSeqCounts}}. Dispersions and size factors must be
    set, usually using functions \code{\link{estimateJunctionSeqSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}.
  }
  \item{test.formula0}{
    The formula for the null hypothesis. Note that the condition to be tested must be named "condition".
  }
  \item{test.formula1}{
    The formula for the alternative hypothesis. Note that the condition to be tested must be named "condition".
  }
  \item{method.GLM}{
    Character string. Can be used to apply alternative methodologies or implementations. Intended for advanced users who have strong opinions about the underlying statistical methodologies.
    
    The default is "advanced" or, equivalently, "DESeq2-style". This uses the 
    model test methodology used by DESeq2 and DEXSeq v1.12.0 or higher. The alternative 
    method is "simpleML" or, equivalently, "DEXSeq-v1.8.0-style". This uses a simpler maximum-likelihood-based method used by the original DESeq and
    by some earlier versions of DEXSeq (v1.8.0 or less).
  }
  \item{dispColumn}{
    Character value. The name of the \code{fData(jscs)} column in which the model dispersion is stored.
  }
  \item{nCores}{
    Either an integer or a BiocParallelParam object. 
    Either way, this determines The number of cores to use. Note that multicore functionality may not be available on all platforms. 
    If parallel execution is not available then JunctionSeq will automatically fallback to single-core execution. See the BiocParallel
    package for more information.
  }
  \item{keep.hypothesisTest.fit}{
    Logical value. If \code{TRUE}, save both complete hypothesis test model fits for every gene. This will require a lot of memory, but may be useful for statistical diagnostics. Default: \code{FALSE}.
  }
  \item{meanCountTestableThreshold}{
    "auto" or Numeric value. Features with a total mean normalized count of less than this value will be excluded from the analyses.
    If left as the default ("auto"), then the cutoff threshold will be determined automatically using the DESeq2 independent filtering method.
  }
  \item{ optimizeFilteringForAlpha }{
    Numeric value between 0 and 1. If \code{meanCountTestableThreshold} is set to "auto" then this sets the adjusted-p-value threshold to optimize against.
  }
  \item{method.cooksFilter}{
    Logical value. if \code{TRUE}, use the cook's filter to detect and remove outliers.
  }
  \item{cooksCutoff}{
    The cook's cutoff threshold to use.
  }
  \item{pAdjustMethod}{
    The p-adjustment method to use with the \code{p.adjust} function.
  }
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}

}

\value{
  A JunctionSeqCountSet, with hypothesis test results included.
}

\examples{
data(exampleDataSet,package="JctSeqData");
jscs <- testForDiffUsage(jscs);

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