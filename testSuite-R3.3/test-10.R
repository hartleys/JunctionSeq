
getwd();

countFiles <- paste0("legacy-datasets/v0.3.17/tinyData/",
               decoder$sample.ID,
               "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz");
gff.file <- paste0("legacy-datasets/v0.3.17/tinyData/withNovel.forJunctionSeq.gff.gz");

runJS <- function(...){
  jscs <<- runJunctionSeqAnalyses(sample.files = countFiles,
           sample.names = decoder$sample.ID,
           condition=factor(decoder$group.ID),
           flat.gff.file = gff.file,
           ...
  );
}

runJS();


plotWrapper({plotJunctionSeqResultsForGene(geneID=g1, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g2, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g3, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=1)

buildAllPlotsForGene(geneID = g1, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)
buildAllPlotsForGene(geneID = g2, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)
buildAllPlotsForGene(geneID = g3, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)

plotWrapper({plotDispEsts(jscs)});
plotWrapper({plotMA(jscs)});

runJS(meanCountTestableThreshold = 20);

plotWrapper({plotJunctionSeqResultsForGene(geneID=g1, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g2, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g3, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=1)

buildAllPlotsForGene(geneID = g1, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)
buildAllPlotsForGene(geneID = g2, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)
buildAllPlotsForGene(geneID = g3, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)

plotWrapper({plotDispEsts(jscs)});
plotWrapper({plotMA(jscs)});
