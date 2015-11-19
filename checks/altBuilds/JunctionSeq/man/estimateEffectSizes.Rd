\name{estimateEffectSizes}
\docType{methods}
\alias{estimateEffectSizes}
\title{
   Estimate Effect Sizes, parameter estimates, etc.
}
\description{
   This function runs fits another generalized linear model to 
   the data, this one intended for use in estimating the
   effect sizes and expression estimates for each analysis.
   
   This function is called internally by the \code{\link{runJunctionSeqAnalyses}} 
   function, and thus for most purposes users should not need to call this
   function directly. It may be useful to advanced users performing non-standard
   analyses.
}
\usage{
   estimateEffectSizes(jscs, 
        method.expressionEstimation = c("feature-vs-gene","feature-vs-otherFeatures"),
        effect.formula = formula(~ condition + countbin + condition : countbin),
        geneLevel.formula = formula(~ condition),
        calculate.geneLevel.expression = TRUE,
        keep.estimation.fit = FALSE,
        nCores=1, 
        dispColumn="dispersion",
        verbose = TRUE)}
\arguments{
  \item{jscs}{
    A \code{JunctionSeqCountSet}. Usually initially created by 
    \code{\link{readJunctionSeqCounts}}. Size factors must be
    set, usually using functions \code{\link{estimateSizeFactors}} and 
    \code{\link{estimateJunctionSeqDispersions}}.
  }
  
  \item{method.expressionEstimation}{
    Character string. Can be used to apply alternative methodologies or implementations. Intended for advanced users who have strong opinions about the underlying statistical methodologies.

    Determines the methodology used to generate feature expression estimates and relative fold changes. By default each feature is modeled separately. Under the default count-vector method,
    this means that the resultant relative fold changes will be a measure of the relative fold change between the feature and the gene as a whole.
    
    Alternatively, the "feature-vs-otherFeatures" method builds a large, complex model containing all features belonging to the gene. The coefficients for each feature
    are then "balanced" using linear contrasts weighted by the inverse of their variance. In general we have found this method to produce very similar results but 
    less efficiently and less consistently. Additionally, this alternative method "multi-counts" reads that cover more than one feature. This can result in over-weighting
    of exonic regions with a large number of annotated variations in a small genomic area, as each individual read or read-pair may be counted many times in the model.
    
    Under the default option, no read or read-pair is ever counted more than once in a given model.
  }
  \item{effect.formula}{
    For advanced users. The base formula for the model used for effect size estimation.
    
    NOTE: the biological condition to be tested must be named "condition".
  }
  \item{geneLevel.formula}{
    For advanced users. The base formula for the model used to estimate total gene-level expression.
    
    NOTE: the biological condition to be tested must be named "condition".
  }
  \item{calculate.geneLevel.expression}{
    Logical value. If \code{TRUE}, gene-level expression will be estimated using the same maximum-likelihood method used in other analyses. Default: \code{TRUE}.
  }
  \item{keep.estimation.fit}{
    Logical value. If \code{TRUE}, save the complete model fits for every gene. This will require a lot of memory, but may be useful for statistical diagnostics. Default: \code{FALSE}.
  }
  
  \item{nCores}{
    Either an integer or a BiocParallelParam object. 
    Either way, this determines The number of cores to use. Note that multicore functionality may not be available on all platforms. 
    If parallel execution is not available then JunctionSeq will automatically fallback to single-core execution. See the BiocParallel
    package for more information.
  }
  \item{dispColumn}{
    Character value. The name of the \code{fData(jscs)} column in which the model dispersion is stored.
  }
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}
}
\value{
  A JunctionSeqCountSet, with effect size results included.
}

\examples{
data(exampleDataSet,package="JctSeqExData2");
jscs <- estimateEffectSizes(jscs);

\dontrun{
#Full example (from scratch):

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