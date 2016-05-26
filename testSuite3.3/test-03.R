
runJS(method.dispFit="local");

plotWrapper({plotJunctionSeqResultsForGene(geneID=g1, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", label.p.vals = F)}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g2, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", label.p.vals = F)}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g3, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", label.p.vals = F)}, asp=1)

buildAllPlotsForGene(geneID = g1, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE, label.p.vals = F)
buildAllPlotsForGene(geneID = g2, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE, label.p.vals = F)
buildAllPlotsForGene(geneID = g3, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE, label.p.vals = F)

plotWrapper({plotDispEsts(jscs)});
plotWrapper({plotMA(jscs)});

runJS(meanCountTestableThreshold = 20,method.dispFit="local");

plotWrapper({plotJunctionSeqResultsForGene(geneID=g1, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", splice.junction.drawing.style = "ellipse")}, asp=0.75, id="SJDS-ellipse")
plotWrapper({plotJunctionSeqResultsForGene(geneID=g2, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", splice.junction.drawing.style = "triangular")}, asp=0.75, id="SJDS-triangular")
plotWrapper({plotJunctionSeqResultsForGene(geneID=g3, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", splice.junction.drawing.style = "line")}, asp=1, id="SJDS-line")
plotWrapper({plotJunctionSeqResultsForGene(geneID=g3, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", splice.junction.drawing.style = "line", draw.nested.SJ = F)}, asp=1, id="SJDS-line.nested-F")
plotWrapper({plotJunctionSeqResultsForGene(geneID=g3, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", splice.junction.drawing.style = "ellipse", draw.nested.SJ = F, merge.exon.parts=F)}, asp=1, id="SJDS-ellipse.nested-F.MEP-F")
plotWrapper({plotJunctionSeqResultsForGene(geneID=g3, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", draw.nested.SJ = F, merge.exon.parts=F, label.chromosome=F)}, asp=1, id="nested-F.MEP-F.LC-F")

buildAllPlotsForGene(geneID = g1, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)
buildAllPlotsForGene(geneID = g2, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)
buildAllPlotsForGene(geneID = g3, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE)
buildAllPlotsForGene(geneID = g3, jscs, outfile.prefix = outputPrefix(), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE, label.p.vals = FALSE)

plotWrapper({plotDispEsts(jscs)});
plotWrapper({plotMA(jscs)});

