
runJS();

plotWrapper({plotJunctionSeqResultsForGene(geneID=g1, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", label.p.vals = F)}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g2, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", label.p.vals = F)}, asp=0.75)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g3, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", label.p.vals = F)}, asp=1)
plotWrapper({plotJunctionSeqResultsForGene(geneID=g3, jscs=jscs, displayTranscripts=TRUE, plot.type = "rExpr", label.p.vals = F)}, asp=1)

plotWrapper({plotMA(jscs)}, id = "std");
plotWrapper({plotMA(jscs,FDR.threshold=0.1)}, id = "FDRthresh0.1");
plotWrapper({plotMA(jscs,FDR.threshold=0.5)}, id = "FDRthresh0.5");
plotWrapper({plotMA(jscs,FDR.threshold=1.0)}, id = "FDRthresh1.0");

plotWrapper({plotMA(jscs,use.pch=1)}, id = "usepch1");
plotWrapper({plotMA(jscs,use.pch=2)}, id = "usepch2");
plotWrapper({plotMA(jscs,use.pch=3)}, id = "usepch3");

plotWrapper({plotMA(jscs,smooth.nbin=512)}, id = "smooth.nbin");
plotWrapper({plotMA(jscs,smooth.nbin=128)}, id = "smooth.nbin");
plotWrapper({plotMA(jscs,smooth.nbin=64)}, id = "smooth.nbin");
plotWrapper({plotMA(jscs,label.counts=FALSE)}, id = "label.counts");
plotWrapper({plotMA(jscs,label.axes = c(FALSE,FALSE,TRUE,TRUE))}, id = "label.axes");
plotWrapper({plotMA(jscs,show.labels = FALSE)}, id = "show.labels");
plotWrapper({plotMA(jscs,par.cex = 2)}, id = "par.cex");
plotWrapper({plotMA(jscs,points.cex = 2)}, id = "points.cex");
plotWrapper({plotMA(jscs,text.cex = 2)}, id = "text.cex");
plotWrapper({plotMA(jscs,lines.cex = 16)}, id = "lines.cex");
plotWrapper({plotMA(jscs,anno.lwd = 4)}, id = "anno.lwd");
plotWrapper({plotMA(jscs,mar=c(2,2,1,1))}, id = "mar");
plotWrapper({plotMA(jscs,miniTicks=FALSE)}, id = "miniTicks");
