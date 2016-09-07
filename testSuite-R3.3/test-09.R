
runJS();

buildAllPlots(jscs, outfile.prefix = outputDir("CairoPNG"), use.plotting.device = "CairoPNG", plotting.device.params = list(res = 72), gene.list = c(g1));
buildAllPlots(jscs, outfile.prefix = outputDir("svg"), use.plotting.device = "svg", gene.list = c(g1));



buildAllPlotsForGene(geneID = g1, jscs, outfile.prefix = outputPrefix("tiff"), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE,
use.plotting.device = "tiff", plotting.device.params = list(res = 72)
)


buildAllPlotsForGene(geneID = g1, jscs, outfile.prefix = outputPrefix("cairo_ps"), with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE,
use.plotting.device = "cairo_ps"
)
