message("#############################################################################################");
message("### Loading packages (",date(),")");
message("#############################################################################################");

library("JunctionSeq");
library("JctSeqExData2");
sessionInfo();

message("#############################################################################################");
message("### Loading data (",date(),")");
message("#############################################################################################");

decoder.file <- system.file("extdata/annoFiles/decoder.bySample.txt",package="JctSeqExData2");
decoder <- read.table(decoder.file,
                      header=TRUE,
                      stringsAsFactors=FALSE);
gff.file <- system.file("extdata/tiny/withNovel.forJunctionSeq.gff.gz",package="JctSeqExData2");
countFiles <- system.file(paste0("extdata/tiny/",
               decoder$sample.ID,
               "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
               package="JctSeqExData2");
des <- data.frame(condition = factor(decoder$group.ID));
sample.files <- countFiles
sample.names <- decoder$sample.ID
condition <-  decoder$group.ID

mv.decoder <- rbind(decoder, c("SAMP7","CASE"), c("SAMP8","CTRL"));
mv.use.covars <- data.frame(smokeStatus = c("Y","Y","N","Y","Y","N","N","N"));
mv.countFiles <- c(countFiles, countFiles[1], countFiles[4]);
mv.des <- data.frame(condition = factor(mv.decoder$group.ID), smokeStatus = mv.use.covars$smokeStatus);

mv.sample.files <<- mv.countFiles
mv.sample.names <<- mv.decoder$sample.ID
mv.condition <<-  mv.decoder$group.ID

g1 <- "ENSRNOG00000056944";
g2 <- "ENSRNOG00000004621";
g3 <- "ENSRNOG00000009281";

message("#############################################################################################");
message("### Loading Functions (",date(),")");
message("#############################################################################################");

runJS <- function(...){
  jscs <<- runJunctionSeqAnalyses(sample.files = countFiles,
           sample.names = decoder$sample.ID,
           condition=factor(decoder$group.ID),
           flat.gff.file = gff.file,
           ...
  );
}
runJS.mv <- function(...){
  jscs <<- runJunctionSeqAnalyses(sample.files = countFiles,
           sample.names = mv.decoder$sample.ID,
           condition = factor(mv.decoder$group.ID),
           use.covars = mv.use.covars,
           flat.gff.file = gff.file,
           ...
  );
}

f.na <- function(x){ ifelse(is.na(x), FALSE,x) }

padInt <- function(x, cols = 4){
  sprintf(paste0("%0",cols,"d"), x);
}

counter <- 0;

outputPrefix <- function(id = "miscTest"){
  counter <<- counter + 1;
  paste0("out/",TEST.ID,".",padInt(counter),".",id,".");
}
plotWrapper <- function(EXPR, id = "miscTest", asp = 1, height = 10, width = 10 * asp, units="in",res=100, pointsize = 12, ...){
  counter <<- counter + 1;
  png(file = paste0("out/",TEST.ID,".",padInt(counter),".",id,".png"), height = height, width = width, units = units, res= res, pointsize =pointsize, ...);
  eval(EXPR);
  dev.off();
}

message("#############################################################################################");
message("### Loading Complete (",date(),")");
message("### For test: ",TEST.ID)
message("#############################################################################################");