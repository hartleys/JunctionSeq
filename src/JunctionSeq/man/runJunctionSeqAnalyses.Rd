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
   \code{\link{fitDispersionFunction}}, \code{\link{testForDiffUsage}}, and 
   \code{\link{estimateEffectSizes}}.
}
\usage{
  runJunctionSeqAnalyses(sample.files, sample.names, condition, 
                         flat.gff.file, 
                         outfile.prefix, saveState = FALSE,
                         analysis.type = c("junctionsAndExons","junctionsOnly","exonsOnly"),
                         use.exons, use.junctions, 
                         use.novel.junctions = TRUE, 
                         meanCountTestableThreshold = 7.5,
                         nCores = 1,
                         use.covars = NULL,                          
                         test.formula0 = formula(~ sample + countbin), 
                         test.formula1 = formula(~ sample + countbin + condition : countbin),
                         effect.formula = formula(~ condition + countbin + condition : countbin),
                         geneLevel.formula = formula(~ condition),
                         gzip.output = TRUE,
                         test.aggregated.genes = FALSE,
                         fitDispersionsForExonsAndJunctionsSeparately = TRUE,
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
  \item{analysis.type}{
    Character string. One of "junctionsAndExons", "junctionsOnly", or "exonsOnly". This parameter determines 
    what type of analysis is to be performed. By default JunctionSeq tests both splice junction 
    loci and exonic regions for differential usage (a "hybrid" analysis). This parameter
    can be used to limit analyses specifically to either splice junction loci or exonic regions.
  }
  
  \item{use.exons}{
    Logical value. This is an alternate parameterization of the \code{analysis.type} parameter. 
    If \code{TRUE}, then exonic region loci will be included in the analyses and will be tested for 
    differential usage. If this parameter is set, then parameter \code{use.junctions} must also be set.
  }
  \item{use.junctions}{
    Logical value. This is an alternate parameterization of the \code{analysis.type} parameter. 
    If \code{TRUE}, then splice junction loci will be included in the analyses and will be tested for 
    differential usage. If this parameter is set, then parameter \code{use.exons} must also be set.
  }
  \item{use.novel.junctions}{
    Logical value. If \code{TRUE}, then novel splice junctions will not be filtered out prior to analysis.
  }

  \item{meanCountTestableThreshold}{
    Numeric value. Features with a total mean normalized count of less than this value will be excluded from the analyses.
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
  \item{geneLevel.formula}{
    For advanced users. The base formula for the model used to estimate total gene-level expression.
    
    NOTE: the biological condition to be tested must be named "condition".
  }
  \item{gzip.output}{Determines whether the text output should be gzip compressed or in plaintext.}

  \item{test.aggregated.genes}{
    Logical value. Whether to attempt to test "aggregate genes" which consist of multiple genes that overlap with one another. Note that inclusion of aggregate genes may affect the false discovery rate, since by their very nature aggregate genes will often show differential splice junction usage, as the two genes will often be regulated independently.
  }
  \item{fitDispersionsForExonsAndJunctionsSeparately}{
    When running a "hybrid" analysis in which both exons and splice junctions are being tested simultaniously, this parameter
    determines whether a single fitted dispersion model should be fitted for both exons and splice junctions, or if separate 
    fitted dispersions should be calculated for each. By default the dispersions are run separately.
  }
  
  
  \item{use.alternate.method}{DEPRECIATED: whether to use the recommended JunctionSeq model framework, or to fallback to the DEXSeq-style model framework.}
  \item{verbose}{if TRUE, send debugging and progress messages to the console / stdout.}

}