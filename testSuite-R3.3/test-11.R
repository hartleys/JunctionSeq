
#gff.file <- system.file("extdata/annoFiles/withNovel.forJunctionSeq.gff.gz",package="JctSeqData");
#countFiles <- system.file(paste0("extdata/cts/",
#               decoder$sample.ID,
#               "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
#               package="JctSeqData");
#runJS <- function(...){
#  jscs <<- runJunctionSeqAnalyses(sample.files = countFiles,
#           sample.names = decoder$sample.ID,
#           condition=factor(decoder$group.ID),
#           flat.gff.file = gff.file,
#           ...
#  );
#}

runJS(gene.names = geneID.to.symbol.file,
      use.multigene.aggregates = TRUE
      );

RES <- 30

buildAllPlots(jscs, outfile.prefix = outputDir("FTT"), 
with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE,
              base.plot.height = 12, base.plot.width = 12, 
              base.plot.units = "in", 
              plotting.device.params = list(res = RES),
              number.plots = FALSE,
              minimalImageFilenames=TRUE,
              name.files.with.geneID = TRUE
)

buildAllPlots(jscs, outfile.prefix = outputDir("TTT"), 
with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE,
              base.plot.height = 12, base.plot.width = 12, 
              base.plot.units = "in", 
              plotting.device.params = list(res = RES),
              number.plots = TRUE,
              minimalImageFilenames=TRUE,
              name.files.with.geneID = TRUE
)

buildAllPlots(jscs, outfile.prefix = outputDir("TFT"), 
with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE,
              base.plot.height = 12, base.plot.width = 12, 
              base.plot.units = "in", 
              plotting.device.params = list(res = RES),
              number.plots = TRUE,
              minimalImageFilenames=FALSE,
              name.files.with.geneID = TRUE
)

buildAllPlots(jscs, outfile.prefix = outputDir("FFT"), 
with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE,
              base.plot.height = 12, base.plot.width = 12, 
              base.plot.units = "in", 
              plotting.device.params = list(res = RES),
              number.plots = FALSE,
              minimalImageFilenames=FALSE,
              name.files.with.geneID = TRUE
)

#######

buildAllPlots(jscs, outfile.prefix = outputDir("FTF"), 
with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE,
              base.plot.height = 12, base.plot.width = 12, 
              base.plot.units = "in", 
              plotting.device.params = list(res = RES),
              number.plots = FALSE,
              minimalImageFilenames=TRUE,
              name.files.with.geneID = FALSE
)

buildAllPlots(jscs, outfile.prefix = outputDir("TTF"), 
with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE,
              base.plot.height = 12, base.plot.width = 12, 
              base.plot.units = "in", 
              plotting.device.params = list(res = RES),
              number.plots = TRUE,
              minimalImageFilenames=TRUE,
              name.files.with.geneID = FALSE
)

buildAllPlots(jscs, outfile.prefix = outputDir("TFF"), 
with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE,
              base.plot.height = 12, base.plot.width = 12, 
              base.plot.units = "in", 
              plotting.device.params = list(res = RES),
              number.plots = TRUE,
              minimalImageFilenames=FALSE,
              name.files.with.geneID = FALSE
)

buildAllPlots(jscs, outfile.prefix = outputDir("FFF"), 
with.TX=TRUE,without.TX=FALSE, expr.plot=FALSE,normCounts.plot=TRUE, rExpr.plot=FALSE,rawCounts.plot=FALSE,
              base.plot.height = 12, base.plot.width = 12, 
              base.plot.units = "in", 
              plotting.device.params = list(res = RES),
              number.plots = FALSE,
              minimalImageFilenames=FALSE,
              name.files.with.geneID = FALSE
)
