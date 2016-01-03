
runJS();

writeCompleteResults(jscs, outfile.prefix = outputPrefix("std"), gzip.output = FALSE);

writeCompleteResults(jscs, outfile.prefix = outputPrefix("FDR-1"), gzip.output = FALSE, FDR.threshold = -1);

writeCompleteResults(jscs, outfile.prefix = outputPrefix("noBED"), gzip.output = FALSE, save.bedTracks=FALSE);

writeCompleteResults(jscs, outfile.prefix = outputPrefix("noAll"), gzip.output = FALSE, save.allGenes=FALSE);
writeCompleteResults(jscs, outfile.prefix = outputPrefix("noAll.noBED"), gzip.output = FALSE, save.allGenes=FALSE, save.bedTracks=FALSE);

writeCompleteResults(jscs, outfile.prefix = outputPrefix("noSig"), gzip.output = FALSE, save.sigGenes=FALSE);
writeCompleteResults(jscs, outfile.prefix = outputPrefix("noSig.noBED"), gzip.output = FALSE, save.sigGenes=FALSE, save.bedTracks=FALSE);
