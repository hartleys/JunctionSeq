
runJS(method.dispFit="local");

plotWrapper({plotJunctionSeqResultsForGene(geneID=g1, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g2, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g3, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=1)

buildAllPlotsForGene(geneID = g1, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)
buildAllPlotsForGene(geneID = g2, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)
buildAllPlotsForGene(geneID = g3, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)

plotWrapper({plotDispEsts(jscs)});
plotWrapper({plotMA(jscs)});

runJS(meanCountTestableThreshold = 20,method.dispFit="local");

plotWrapper({plotJunctionSeqResultsForGene(geneID=g1, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g2, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g3, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr")}, asp=1)

buildAllPlotsForGene(geneID = g1, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)
buildAllPlotsForGene(geneID = g2, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)
buildAllPlotsForGene(geneID = g3, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)

plotWrapper({plotDispEsts(jscs)});
plotWrapper({plotMA(jscs)});

#Not currently supported:
#runJS(method.dispFit="mean");
#plotWrapper({plotDispEsts(jscs)});
#plotWrapper({plotMA(jscs)});

