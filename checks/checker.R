#Script for running additional checks in a Windows environment.

args <- commandArgs(TRUE)
ver <- args[1];
RVER <- args[2];
wd <- args[3];
winVer <- args[4];
message("Starting BiocChecks for v",ver);
#ver <- "0.6.1";

#setwd(paste0("C:\\Users\\hartleys\\work\\nihwork\\home_copy\\projects\\ZZZ-JunctionSeq\\releases\\v",ver,"\\checks\\"));
setwd(wd);
library("BiocCheck");
#js.tar <- paste0("C:\\Users\\hartleys\\work\\nihwork\\home_copy\\projects\\ZZZ-JunctionSeq\\releases\\v",ver,"\\JunctionSeq_",ver,".tar.gz");
js.tar <- paste0("../JunctionSeq_",ver,".tar.gz");

sc <- file(paste0("R-CMD-BiocCheck-",winVer,"-",RVER,".log"), open = 'a');
sink(file = sc , type = "message");
message("Starting BiocCheck... (",date(),")");
message("     R Version: \"",R.Version(),"\"");
test <- BiocCheck(js.tar);
message("Finished BiocCheck. (",date(),")");
sink(file = NULL, type = "message");
close(sc);

#js.dir <- paste0("C:\\Users\\hartleys\\work\\nihwork\\home_copy\\projects\\ZZZ-JunctionSeq\\releases\\v",ver,"\\JunctionSeq\\");
#library("devtools");
##devtools.check():
#sc <- file("WIN64.devtools.check.log", open = 'w');
#sink(file = sc , type = "message");
#message("Starting devtools::check... (",date(),")");
##message("     Sys.info():");
##for(i in 1:length(Sys.info())){
##message("        [",names(Sys.info())[[i]]," = \"",Sys.info()[[i]],"\"]")}
#devtools::check(pkg = js.dir);
#message("Finished devtools::check. (",date(),")");
#sink(file = NULL, type = "message");
#close(sc);




