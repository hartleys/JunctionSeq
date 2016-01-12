
message("Loading JunctionSeq installer (v0.8.4)");

JunctionSeq.CRAN.dep <- c("statmod","plotrix","stringr","locfit","Hmisc");
JunctionSeq.BIOC.dep <- c("Biobase", "BiocGenerics", "BiocParallel", "GenomicRanges", "IRanges","S4Vectors","genefilter","geneplotter","SummarizedExperiment","DESeq2");
JunctionSeq.CRAN.optional <- c("knitr", "Cairo", "pryr","MASS");
JunctionSeq.BIOC.optional <- c("BiocStyle","BiocCheck")


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
    source("http://bioconductor.org/biocLite.R")
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

##########
#User functions:

 
JS.install.CRAN.dependencies <- function(verbose = TRUE, ...){
  for(d in JunctionSeq.CRAN.dep){
    attemptLoadAndInstall.CRAN(d, verbose=verbose, ...);
  }
}

JS.install.BIOC.dependencies <- function(...){
  source("http://bioconductor.org/biocLite.R");
  biocLite();

  for(d in JunctionSeq.BIOC.dep){
    attemptLoadAndInstall.BIOC(d, verbose=verbose, ...);
  }
}

JS.install.JunctionSeq <- function(...){
  install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz",
                   repos=NULL, type="source", ...);
}


JS.install.optional <- function(){
  for(d in JunctionSeq.CRAN.optional){
    attemptLoadAndInstall.CRAN(d, verbose=verbose);
  }
  for(d in JunctionSeq.BIOC.optional){
    attemptLoadAndInstall.BIOC(d, verbose=verbose);
  }
  
  install.packages("http://hartleys.github.io/JunctionSeq/install/JctSeqData_LATEST.tar.gz",repos=NULL, type="source");
}


JS.install <- function(install.from.source = TRUE){
  message("Installing dependencies...");

  message("Installing CRAN package dependencies... ",date());
  message("   This can be called directly using the function JS.install.CRAN.dependencies(...)");
  JS.install.CRAN.dependencies();
  message("Done installing CRAN package dependencies. ",date());

  message("Installing Bioconductor package dependencies... ",date());
  message("   This can be called directly using the function JS.install.BIOC.dependencies(...)");
  JS.install.BIOC.dependencies();
  message("Done installing Bioconductor package dependencies. ",date());

  message("Done installing package dependencies. ",date());

  message("Beginning JunctionSeq Installation... ",date());
  JS.install.JunctionSeq();
  message("Installation complete. ",date());

  message("Testing to see if package can be loaded... ",date());
  library("JunctionSeq");
  message("Done.",date());
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
