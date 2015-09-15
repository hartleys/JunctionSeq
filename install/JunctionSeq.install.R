
message("Installing dependencies...");

message("Installing CRAN package dependencies... ",date());
install.packages("statmod")
install.packages("plotrix")
install.packages("stringr")
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("ggplot")
install.packages("locfit")
install.packages("Hmisc")
message("Done installing CRAN package dependencies. ",date());

message("Installing Bioconductor package dependencies... ",date());
source("http://bioconductor.org/biocLite.R");
biocLite();
biocLite("Biobase");
biocLite("BiocGenerics");
biocLite("BiocParallel");
biocLite("GenomicRanges");
biocLite("IRanges");
biocLite("S4Vectors");
biocLite("genefilter");
biocLite("geneplotter");
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