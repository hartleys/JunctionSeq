JUNCTIONSEQVERSION = "0.6.16"

message("Loading JunctionSeq installer (v0.7.0.0)");
message("For JunctionSeq v",JUNCTIONSEQVERSION);


##########
# Internal Functions:
attemptLoad <- function(pkg, verbose = TRUE){
  out <- suppressWarnings(suppressPackageStartupMessages(require(pkg, character.only = TRUE)))
  return(out);
}
attemptLoadAndInstall.CRAN <- function(pkg, verbose = TRUE, ...){
  if(attemptLoad(pkg, verbose = verbose)){
    message("CRAN Package \"",pkg,"\" already installed.");
    return(TRUE);
  } else {
    message("Attempting installation of CRAN package: \"",pkg,"\"");
    install.packages(pkg, ...);

    if(attemptLoad(pkg, verbose = verbose)){
      message("Successfully installed CRAN package: \"",pkg,"\"");
      return(TRUE);
    } else {
      stop(paste0("Installation FAILED for CRAN package: \"",pkg,"\""));
      return(FALSE);
    }
  }
}
attemptLoadAndInstall.BIOC <- function(pkg, verbose = TRUE, ...){
  if(attemptLoad(pkg, verbose = verbose)){
    message("Bioconductor Package \"",pkg,"\" already installed.");
    return(TRUE);
  } else {
    message("Attempting installation of Bioconductor package: \"",pkg,"\"");
    source("https://bioconductor.org/biocLite.R")
    biocLite(pkg);
    
    if(attemptLoad(pkg, verbose = verbose)){
      message("Successfully installed Bioconductor package: \"",pkg,"\"");
      return(TRUE);
    } else {
      stop(paste0("Installation FAILED for Bioconductor package: \"",pkg,"\""));
      return(FALSE);
    }
  }
}

JunctionSeq.CRAN.dep <- c("statmod","plotrix","stringr","Rcpp","RcppArmadillo","locfit");
JunctionSeq.BIOC.dep <- c("Biobase", "BiocGenerics", "BiocParallel", "GenomicRanges", "IRanges","S4Vectors","genefilter","geneplotter","SummarizedExperiment");
JunctionSeq.CRAN.optional <- c("knitr", "Cairo", "pryr","MASS");
JunctionSeq.BIOC.optional <- c("BiocStyles")

##########
#User functions:


JS.install <- function(install.from.source = TRUE, verbose = TRUE){
  message("Installing dependencies...");

  if(verbose) message("Installing CRAN package dependencies... ",date());
  if(verbose) message("   This can be called directly using the function JS.install.CRAN.dependencies(...)");
  JS.install.CRAN.dependencies(verbose = verbose);
  if(verbose) message("Done installing CRAN package dependencies. ",date());

  if(verbose) message("Installing Bioconductor package dependencies... ",date());
  if(verbose) message("   This can be called directly using the function JS.install.BIOC.dependencies(...)");
  JS.install.BIOC.dependencies(verbose = verbose);
  if(verbose) message("Done installing Bioconductor package dependencies. ",date());

  if(verbose) message("Done installing package dependencies. ",date());

  if(verbose) message("Beginning JunctionSeq Installation... ",date());
  JS.install.JunctionSeq(install.from.source = install.from.source, verbose = verbose);
  if(verbose) message("Installation complete. ",date());

  if(verbose) message("Testing to see if package can be loaded... ",date());
  library("JunctionSeq");
  if(verbose) message("Done.",date());
}

JS.install.CRAN.dependencies <- function(verbose = TRUE, ...){
  for(d in JunctionSeq.CRAN.dep){
    attemptLoadAndInstall.CRAN(d, verbose=verbose, ...);
  }
}

JS.install.BIOC.dependencies <- function(verbose = TRUE, ...){
  source("http://bioconductor.org/biocLite.R");
  biocLite();

  for(d in JunctionSeq.BIOC.dep){
    attemptLoadAndInstall.BIOC(d, verbose=verbose, ...);
  }
}

JS.install.JunctionSeq <- function(install.from.source = TRUE, verbose = TRUE, ...){
  if(install.from.source){
    sysname <- Sys.info()[["sysname"]];
    if(sysname == "Windows"){
      if(verbose) message("Installing from source on Windows.");
      if(verbose) message("  NOTE: You MUST have Rtools installed in order to install from source on Windows!");
      if(verbose) message("        If you do not have Rtools, you can either install Rtools or install the ");
      if(verbose) message("        compiled binary using the command:");
      if(verbose) message("  JS.install.JunctionSeq(install.from.source=FALSE)");
      
    } else if(sysname == "Linux"){
      if(verbose) message("Installing from source on Linux.");
    } else if(sysname == "Darwin" || sysname == "Yosemite" || sysname == "Mavericks"){
      if(verbose) message("Installing from source on OSX (",sysname,")");
      if(verbose) message("  NOTE: You MUST have Xcode and gfortran 4.8.2 installed in order to install from source on OSX!");
      if(verbose) message("        You may need other non-standard tools as well, depending on OS version.");
      if(verbose) message("  If you don't want to install these, you can install ");
      if(verbose) message("  using a pre-compiled binary with the command:");
      if(verbose) message("  JS.install.JunctionSeq(install.from.source=FALSE)");
      
    } else {
      if(verbose) message("Installing from source on unknown OS (",sysname,")");
    }
    install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz",repos=NULL, type="source", ...);
  } else {
    sysname <- Sys.info()[["sysname"]];
    if(sysname == "Windows"){
      if(verbose) message("Installing Windows binary");
      install.packages(paste0("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_",JUNCTIONSEQVERSION,".zip"),repos=NULL, ...);
    } else if(sysname == "Linux"){
      stop(paste0("Compiled JunctionSeq binary not available for OS ",sysname));
    } else if(sysname == "Darwin" || sysname == "Yosemite" || sysname == "Mavericks"){
      if(verbose) message("Installing OSX binary. (",sysname,")");
      if(verbose) message("   NOTE: This binary is experimental and may not work on all versions of OSX.");
      if(verbose) message("         If you do encounter an error, PLEASE post it to the JunctionSeq Issues page");
      if(verbose) message("         Found here: http://github.com/hartleys/JunctionSeq/issues");
      install.packages(paste0("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_",JUNCTIONSEQVERSION,".tgz"),repos=NULL, ...);
    } else {
      if(verbose) message("Operating System \"",sysname,"\" not recognized.");
      stop(paste0("Compiled JunctionSeq binary not available for OS ",sysname));
    }
  }
}


JS.install.optional <- function(verbose = TRUE){
  for(d in JunctionSeq.CRAN.optional){
    attemptLoadAndInstall.CRAN(d, verbose=verbose);
  }
  for(d in JunctionSeq.BIOC.optional){
    attemptLoadAndInstall.BIOC(d, verbose=verbose);
  }
  
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




