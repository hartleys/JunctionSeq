\name{runJunctionSeqAnalyses}
\docType{methods}
\alias{runJunctionSeqAnalyses}
\title{
   Run a JunctionSeq analysis.
}
\description{
   This function runs a complete analysis from start to finish. 
   It internally calls functions \code{\link{readAnnotationData}}, 
   \code{\link{readJunctionSeqCounts}},
   \code{\link{estimateSizeFactors}}, \code{\link{estimateJunctionSeqDispersions}}, 
   \code{\link{fitDispersionFunction}}, \code{\link{testJunctionsForDJU}}, and 
   \code{\link{estimateEffectSizes}}.
}
\usage{
  runJunctionSeqAnalyses(sample.files, sample.names, condition, 
                         flat.gff.file, 
                         outfile.prefix, saveState = FALSE,
                         use.splice.sites = TRUE, 
                         use.novel.splice.sites = TRUE, 
                         use.exons = FALSE, 
                         nCores = 1,
                         use.covars = NULL,                          
                         test.formula0 = formula(~ sample + countbin), 
                         test.formula1 = formula(~ sample + countbin + condition : countbin),
                         effect.formula = formula(~ condition + countbin + condition : countbin),
                         gzip.output = TRUE,
                         test.aggregated.genes = FALSE,
                         use.alternate.method = TRUE,
                         verbose = TRUE)

}
\arguments{
  \item{sample.files}{
    A character vector of sample files. Each sample file is a simple tab-delimited file containing two columns. The first column contains the feature name, using the format GENE_ID:SPLICE_SITE_ID, where GENE_ID can be any alphanumeric string, and SPLICE_SITE_ID is a 3-digit number. The 2nd column is the read count for that feature, as a non-negative integer. Do not use normalized read counts, RPKM, FPKM, or anything other than raw read counts, as this will conflict with the DEXSeq normalization.
  }
  \item{sample.names}{
    A character vector of sample names. This must have the same length as sample.files, and should be in the same order.
  }
  \item{condition}{
    A factor vector of condition values. This must have the same length as sample.files and sample.names, and should be listed in the same order.
  }
  \item{flat.gff.file}{
    A flattened gff-formatted annotation file from which the gene counts were generated. 
    Technically optional, but STRONGLY RECOMMENDED, as the annotation data 
    WILL be required by plotting functions.
  }
  \item{outfile.prefix}{
    The prefix of the output files to be written. By default no output files will be created.
  }
  \item{saveState}{
    If TRUE and if outfile.prefix is non-null, then the ecs and res objects will be saved to disk as RData files after all analysis is complete.
  }
  \item{use.splice.sites}{
    FOR ADVANCED USERS ONLY: Logical value. If \code{TRUE}, splice junctions will be included in the analyses.
  }
  \item{use.novel.splice.sites}{
    Logical value. If \code{TRUE}, splice junctions that are marked as novel in the gff annotation will be included in the analysis.
    Note that this is does nothing if there are no novel splice junctions in the count tables. 
  }
  \item{use.exons}{
    FOR ADVANCED USERS ONLY: Logical value. If \code{TRUE}, exon segments will be included in the analyses.
  }
  \item{nCores}{
    The number of cores to use. Note that multicore functionality may not be available on all platforms. 
    JunctionSeq attempts to use the BiocParallell package if it can be found installed. Otherwise it will
    attempt to fallback to the multicore package. If neither package can be found it will fallback to single
    core execution.
  }  
  \item{use.covars}{
    Optional: for advanced users. A data frame containing covariate factors. The names must be included in the model formulas.
  }
  \item{test.formula0}{
    For advanced users. The base formula for the null hypothesis model used in the hypothesis tests.
    
    NOTE: the biological condition to be tested must be named "condition".
  }
  \item{test.formula1}{
    For advanced users. The base formula for the alternate hypothesis model used in the hypothesis tests.
    
    NOTE: the biological condition to be tested must be named "condition".
  }
  \item{effect.formula}{
    For advanced users. The base formula for the model used for effect size estimation.
    
    NOTE: the biological condition to be tested must be named "condition".
  }
  \item{gzip.output}{Determines whether the text output should be gzip compressed or in plaintext.}

  \item{test.aggregated.genes}{
    Logical value. Whether to attempt to test "aggregate genes" which consist of multiple genes that overlap with one another. Note that inclusion of aggregate genes may affect the false discovery rate, since by their very nature aggregate genes will often show differential splice junction usage, as the two genes will often be regulated independently.
  }
  \item{use.alternate.method}{DEPRECIATED: whether to use the recommended JunctionSeq model framework, or to fallback to the DEXSeq-style model framework.}
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}

}