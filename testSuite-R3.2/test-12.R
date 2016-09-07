
runJS <- function(...){
  jscs <<- runJunctionSeqAnalyses(sample.files = countFiles,
           sample.names = decoder$sample.ID,
           condition=factor(decoder$group.ID),
           ...
  );
}

runJS();

writeCompleteResults(jscs, outfile.prefix = outputPrefix("NOGFF"), gzip.output = FALSE);
