
JS.install <- function(...){
  message("Installing dependencies...");

  message("Installing CRAN package dependencies... ",date());
  message("   This can be called directly using the function JS.install.CRAN.dependencies(...)");
  JS.install.CRAN.dependencies(...);
  message("Done installing CRAN package dependencies. ",date());

  message("Installing Bioconductor package dependencies... ",date());
  message("   This can be called directly using the function JS.install.BIOC.dependencies(...)");
  JS.install.BIOC.dependencies(...);
  message("Done installing Bioconductor package dependencies. ",date());

  message("Done installing package dependencies. ",date());

  message("Beginning JunctionSeq Installation... ",date());
  JS.install.JunctionSeq(...);
  message("Installation complete. ",date());

  message("Testing to see if package can be loaded... ",date());
  library("JunctionSeq");
  message("Done.",date());
}

JS.install.CRAN.dependencies <- function(...){
  if(! suppressWarnings(suppressPackageStartupMessages(require("statmod"))))       install.packages("statmod"        , ...)
  if(! suppressWarnings(suppressPackageStartupMessages(require("plotrix"))))       install.packages("plotrix"        , ...)
  if(! suppressWarnings(suppressPackageStartupMessages(require("stringr"))))       install.packages("stringr"        , ...)
  if(! suppressWarnings(suppressPackageStartupMessages(require("Rcpp"))))          install.packages("Rcpp"           , ...)
  if(! suppressWarnings(suppressPackageStartupMessages(require("RcppArmadillo")))) install.packages("RcppArmadillo"  , ...)
  #if(! suppressWarnings(suppressPackageStartupMessages(require("ggplot2"))))       install.packages("ggplot2"        , ...)
  if(! suppressWarnings(suppressPackageStartupMessages(require("locfit"))))        install.packages("locfit"         , ...)
  #if(! suppressWarnings(suppressPackageStartupMessages(require("Hmisc"))))         install.packages("Hmisc"          , ...)
}

JS.install.BIOC.dependencies <- function(...){
  source("http://bioconductor.org/biocLite.R");
  biocLite();
  if(! suppressWarnings(suppressPackageStartupMessages(require("Biobase"))))              biocLite("Biobase"              , ...);
  if(! suppressWarnings(suppressPackageStartupMessages(require("BiocGenerics"))))         biocLite("BiocGenerics"         , ...);
  if(! suppressWarnings(suppressPackageStartupMessages(require("BiocParallel"))))         biocLite("BiocParallel"         , ...);
  if(! suppressWarnings(suppressPackageStartupMessages(require("GenomicRanges"))))        biocLite("GenomicRanges"        , ...);
  if(! suppressWarnings(suppressPackageStartupMessages(require("IRanges"))))              biocLite("IRanges"              , ...);
  if(! suppressWarnings(suppressPackageStartupMessages(require("S4Vectors"))))            biocLite("S4Vectors"            , ...);
  if(! suppressWarnings(suppressPackageStartupMessages(require("genefilter"))))           biocLite("genefilter"           , ...);
  if(! suppressWarnings(suppressPackageStartupMessages(require("geneplotter"))))          biocLite("geneplotter"          , ...);
  if(! suppressWarnings(suppressPackageStartupMessages(require("SummarizedExperiment")))) biocLite("SummarizedExperiment" , ...);
}

JS.install.JunctionSeq <- function(repos = NULL, type = "source", ...){
  install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz",repos=repos, type=type, ...);
}

JS.install.optional <- function(...){
  if(! suppressWarnings(suppressPackageStartupMessages(require("knitr"))))       install.packages("knitr"        , ...)
  if(! suppressWarnings(suppressPackageStartupMessages(require("Cairo"))))       install.packages("Cairo"        , ...)
  if(! suppressWarnings(suppressPackageStartupMessages(require("pryr"))))        install.packages("pryr"         , ...)
  if(! suppressWarnings(suppressPackageStartupMessages(require("MASS"))))        install.packages("MASS"        , ...)
  
  install.packages("http://hartleys.github.io/JunctionSeq/install/JctSeqExData2_LATEST.tar.gz",repos=NULL, type="source", ...);
}



message("JunctionSeq Installation scripts loaded.");
message("To install CRAN dependencies, use the command:");
message("JS.install.CRAN.dependencies()");
message("   (note: all options will be bassed directly to install.packages())");
message("");
message("To install Bioconductor dependencies, use the command:");
message("JS.install.BIOC.dependencies()");
message("   (note: all options will be bassed directly to biocLite())");
message("");
message("To install JunctionSeq, once all dependencies are installed, use the command:");
message("JS.install.JunctionSeq()");
message("   (note: all options will be bassed directly to install.packages())");
message("");
message("To install everything automatically, use the command:");
message("JS.install()");
message("");
message("Consider installing the optional packages as well, using the command:");
message("JS.install.optional()");




