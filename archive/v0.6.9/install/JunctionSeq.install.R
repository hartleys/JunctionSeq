

message("this installer is deprecated! use JS.install.R instead!")
message("");
message("Use the R commands:");
message("source(\"http://hartleys.github.io/JunctionSeq/install/JS.install.R\");");
message("JS.install();");


oldInstallerFcn <- function(){
  message("Installing CRAN package dependencies... ",date());
  if(! require("statmod")) install.packages("statmod")
  if(! require("plotrix")) install.packages("plotrix")
  if(! require("stringr")) install.packages("stringr")
  if(! require("Rcpp")) install.packages("Rcpp")
  if(! require("RcppArmadillo")) install.packages("RcppArmadillo")
  if(! require("ggplot")) install.packages("ggplot")
  if(! require("locfit")) install.packages("locfit")
  if(! require("Hmisc")) install.packages("Hmisc")
  message("Done installing CRAN package dependencies. ",date());

  message("Installing Bioconductor package dependencies... ",date());
  source("http://bioconductor.org/biocLite.R");
  biocLite();
  if(! require("Biobase")) biocLite("Biobase");
  if(! require("BiocGenerics")) biocLite("BiocGenerics");
  if(! require("BiocParallel")) biocLite("BiocParallel");
  if(! require("GenomicRanges")) biocLite("GenomicRanges");
  if(! require("IRanges")) biocLite("IRanges");
  if(! require("S4Vectors")) biocLite("S4Vectors");
  if(! require("genefilter")) biocLite("genefilter");
  if(! require("geneplotter")) biocLite("geneplotter");



  message("Done installing Bioconductor package dependencies. ",date());

  message("Done installing package dependencies. ",date());

  message("Beginning JunctionSeq Installation... ",date());

  if(.Platform$OS.type == "unix"){
    message("Auto-detected unix-like system. Attempting to install from source...");
    install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz",repos=NULL, type="source");
  } else {
    message("Auto-detected windows system. Attempting to install from source using Rtools...");
    install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz",repos=NULL, type="source");
  }

  message("Installation complete. ",date());

  message("Testing to see if package can be loaded... ",date());
  library("JunctionSeq");
  message("Done.",date());
}
