message("Loading JunctionSeq installer (v0.6.1e)");

attemptLoad <- function(pkg){
  return(suppressWarnings(suppressPackageStartupMessages(require(pkg))));
}

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

JS.install.CRAN.dependencies <- function( ...){
  if(! attemptLoad("statmod"))        install.packages("statmod"        , ...)
  if(! attemptLoad("statmod")) stop(paste0("Installation of package statmod failed"));
  
  if(! attemptLoad("plotrix"))        install.packages("plotrix"        , ...)
  if(! attemptLoad("plotrix")) stop(paste0("Installation of package plotrix failed"));
  
  if(! attemptLoad("stringr"))        install.packages("stringr"        , ...)
  if(! attemptLoad("stringr")) stop(paste0("Installation of package stringr failed"));

  if(! attemptLoad("Rcpp"))           install.packages("Rcpp"           , ...)
  if(! attemptLoad("Rcpp")) stop(paste0("Installation of package Rcpp failed"));

  
  if(! attemptLoad("RcppArmadillo"))  install.packages("RcppArmadillo"  , ...)
  if(! attemptLoad("RcppArmadillo")) stop(paste0("Installation of package RcppArmadillo failed"));

  if(! attemptLoad("locfit"))         install.packages("locfit"         , ...)
  if(! attemptLoad("locfit")) stop(paste0("Installation of package locfit failed"));

  #Deprecated dependencies:
  #if(! attemptLoad("ggplot2"))))       install.packages("ggplot2"        , ...)
  #if(! attemptLoad("Hmisc"))))         install.packages("Hmisc"          , ...)
}

JS.install.BIOC.dependencies <- function(...){
  source("http://bioconductor.org/biocLite.R");
  biocLite();
  if(! attemptLoad("Biobase"))              biocLite("Biobase"              , ...);
  if(! attemptLoad("Biobase")) stop(paste0("Installation of package Biobase failed"));

  
  if(! attemptLoad("BiocGenerics"))         biocLite("BiocGenerics"         , ...);
  if(! attemptLoad("BiocGenerics")) stop(paste0("Installation of package BiocGenerics failed"));

  
  if(! attemptLoad("BiocParallel"))         biocLite("BiocParallel"         , ...);
  if(! attemptLoad("BiocParallel")) stop(paste0("Installation of package BiocParallel failed"));


  if(! attemptLoad("GenomicRanges"))        biocLite("GenomicRanges"        , ...);
  if(! attemptLoad("GenomicRanges")) stop(paste0("Installation of package GenomicRanges failed"));

  if(! attemptLoad("IRanges"))              biocLite("IRanges"              , ...);
  if(! attemptLoad("IRanges")) stop(paste0("Installation of package IRanges failed"));

  if(! attemptLoad("S4Vectors"))            biocLite("S4Vectors"            , ...);
  if(! attemptLoad("S4Vectors")) stop(paste0("Installation of package S4Vectors failed"));

  if(! attemptLoad("genefilter"))           biocLite("genefilter"           , ...);
  if(! attemptLoad("genefilter")) stop(paste0("Installation of package genefilter failed"));

  if(! attemptLoad("geneplotter"))          biocLite("geneplotter"          , ...);
  if(! attemptLoad("geneplotter")) stop(paste0("Installation of package geneplotter failed"));

  if(! attemptLoad("SummarizedExperiment")) biocLite("SummarizedExperiment" , ...);
  if(! attemptLoad("SummarizedExperiment")) stop(paste0("Installation of package SummarizedExperiment failed"));

}

JS.install.JunctionSeq <- function(installFromSource = TRUE, ...){
  if(Sys.info()[["sysname"]] == "Linux"){
    #That's fine. Install from source.
    message("Installing from source (Linux)");
    install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz",repos=NULL, type="source", ...);
  } else if(Sys.info()[["sysname"]] == "Windows"){
    if(installFromSource){
      message("Installing from source (Windows)");
      message("   Note: This requires Rtools. ");
      message("   (https://cran.r-project.org/bin/windows/Rtools/)");
      message("To install without Rtools, you can use the command:");
      message("   JS.install.JunctionSeq(installFromSource=FALSE);");
      message("   (installs a compiled binary)");
      install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz",repos=NULL, type="source", ...);
    } else {
      message("Installing from binary (Windows)");
      install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST-WIN.zip",repos=NULL, type="source", ...);
    }
  } else {
    message("Installing from source (",Sys.info()[["sysname"]],")");
    message("Note: your OS is not one of the supported operating systems. (Linux or Windows)");
    message("You may need to install additional compilers to install from source.");
    message("For OSX you will need Xcode tools and gfortran.");
    install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz",repos=NULL, type="source", ...);
  }
}

JS.install.optional <- function(...){
  if(! attemptLoad("knitr"))       install.packages("knitr"        , ...)
  if(! attemptLoad("knitr")) stop(paste0("Installation of package SummarizedExperiment failed"));

  if(! attemptLoad("Cairo"))       install.packages("Cairo"        , ...)
  if(! attemptLoad("Cairo")) stop(paste0("Installation of package Cairo failed"));

  if(! attemptLoad("pryr"))        install.packages("pryr"         , ...)
  if(! attemptLoad("pryr")) stop(paste0("Installation of package pryr failed"));

  if(! attemptLoad("MASS"))        install.packages("MASS"        , ...)
  if(! attemptLoad("MASS")) stop(paste0("Installation of package MASS failed"));

  source("http://bioconductor.org/biocLite.R");
  biocLite();
  if(! attemptLoad("BiocStyles"))          biocLite("BiocStyles"          , ...);
  if(! attemptLoad("BiocStyles")) stop(paste0("Installation of package BiocStyles failed"));
  
  install.packages("http://hartleys.github.io/JunctionSeq/install/JctSeqExData2_LATEST.tar.gz",repos=NULL, type="source", ...);
}

JS.install.windowsBinary <- function(...){
  install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST-WIN.zip", repos=NULL);
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
message("To install to windows via a binary (does not require Rtools):");
message("JS.install.JunctionSeq(installFromSource = FALSE)");
message("");
message("To install all required packages and JunctionSeq (from source), use the command:");
message("JS.install()");
message("");
message("Consider installing the optional packages as well, using the command:");
message("JS.install.optional()");




