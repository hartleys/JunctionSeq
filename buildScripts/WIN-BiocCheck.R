#Script for running additional checks in a Windows environment.

args <- commandArgs(TRUE)
js.tar <- args[1];
wd <- args[2]
logFile <- args[3];

setwd(wd);
library("BiocCheck");

sc <- file(logFile, open = 'a');
sink(file = sc , type = "message");
message("Starting BiocCheck... (",date(),")");
message("     R Version: \"",R.Version(),"\"");
test <- BiocCheck(js.tar);
message("Finished BiocCheck. (",date(),")");
sink(file = NULL, type = "message");
close(sc);





