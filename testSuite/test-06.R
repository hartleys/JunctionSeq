
runJS();
plotWrapper({plotDispEsts(jscs)}, id = "std");
plotWrapper({plotDispEsts(jscs,linecol=c("purple","orange"))}, id = "linecol");
plotWrapper({plotDispEsts(jscs, pointcol=c("purple","orange"))}, id = "pointcol");

plotWrapper({plotDispEsts(jscs, title.main="Custom Title")}, id = "titlemain");
plotWrapper({plotDispEsts(jscs, xlab="Custom xlab")}, id = "xlab");
plotWrapper({plotDispEsts(jscs, ylab="Custom ylab")}, id = "ylab");
plotWrapper({plotDispEsts(jscs, miniTicks =FALSE)}, id = "miniTicks");
plotWrapper({plotDispEsts(jscs, pch.MLE = 100)}, id = "pch.MLE");
plotWrapper({plotDispEsts(jscs, pch.MAP = 2)}, id = "pch.MAP");
plotWrapper({plotDispEsts(jscs, lwd.fitted = 4)}, id = "lwd.fitted");
plotWrapper({plotDispEsts(jscs, par.cex = 2)}, id = "par.cex");
plotWrapper({plotDispEsts(jscs, points.cex = 2)}, id = "points.cex");
plotWrapper({plotDispEsts(jscs, text.cex = 2)}, id = "text.cex");
plotWrapper({plotDispEsts(jscs, lines.cex = 16)}, id = "lines.cex");

#Smoothscatter dispersion plots are deprecated
#plotWrapper({plotDispEsts(jscs, use.smoothScatter = TRUE)}, id = "use.smoothScatter");
#plotWrapper({plotDispEsts(jscs, use.smoothScatter = TRUE, smooth.nbin=1024)}, id = "use.smoothScatter");
#plotWrapper({plotDispEsts(jscs, use.smoothScatter = TRUE, smooth.nbin=64)}, id = "use.smoothScatter");

plotWrapper({plotDispEsts(jscs,plot.exon.results = FALSE)}, id = "noExon");
plotWrapper({plotDispEsts(jscs,plot.junction.results = FALSE)}, id = "noJunction");
plotWrapper({plotDispEsts(jscs, anno.lwd = 4)}, id = "anno.lwd");
plotWrapper({plotDispEsts(jscs, mar = c(6,3,1,3))}, id = "mar");
plotWrapper({plotDispEsts(jscs, show.legends=FALSE)}, id = "show.legends");
