#These functions were (loosely) based on similar functions created for the DEXSeq package.
#
# Note that DEXSeq is licensed under the GPL v3. Therefore this
#   code packaged together is licensed under the GPL v3, as noted in the LICENSE file.
# Some code snippets are "united states government work" and thus cannot be
#   copyrighted. See the LICENSE file for more information.
#
# The current versions of the original functions upon which these were based can be found
#    here: http://github.com/Bioconductor-mirror/DEXSeq
#
# Updated Authorship and license information can be found here:
#   here: http://github.com/Bioconductor-mirror/DEXSeq/blob/master/DESCRIPTION


##########################################################################
##########################################################################
##########################################################################
######### Plotting Functions:
##########################################################################
##########################################################################
##########################################################################

plotDispEsts <- function( jscs, ylim, xlim, 
  linecol=c("#0000FF","#FF0000"), pointcol = c("#00009980","#99000080"),
  title.main = "Dispersion Estimates", xlab = "Mean Normalized Coverage", ylab = "Dispersion",
  miniTicks = TRUE,
  pch.MLE = 46, pch.MAP = 1,
  lwd.fitted = 2,
  par.cex = 1, points.cex = 1, text.cex = 1,lines.cex = 8,
  use.smoothScatter = FALSE, smooth.nbin = 512, nrpoints = 100,
  plot.exon.results = TRUE, 
  plot.junction.results = TRUE, 
  anno.lwd = 2,
  mar = c(4.1,4.1,3.1,1.1),
  show.legends = TRUE,
  verbose = TRUE, debug.mode = FALSE,
  ... )
{
  
    px = rowMeans( counts( jscs, normalized=TRUE ) );
    #pch <- ifelse(fData(jscs)$featureType == "exonic_part", pch[1], pch[2]);
    pchColor <- ifelse(fData(jscs)$featureType == "exonic_part", pointcol[1], pointcol[2]);

    sel <- rep(TRUE, length(px));
    if(! plot.exon.results){
      sel <- sel & fData(jscs)$featureType != "exonic_part";
    }
    if(! plot.junction.results){
      sel <- sel & (fData(jscs)$featureType != "splice_site" & fData(jscs)$featureType != "novel_splice_site");
    }
    py <- fData(jscs)$dispBeforeSharing
    sel <- sel & f.na(px>0);
    sel <- sel & f.na(py>0);
    sel <- sel & (! is.na(px));
    sel <- sel & (! is.na(py));
    sel <- sel & f.na(py > 1e-6);

    px <- px[sel]
    py <- py[sel]
    
    pyFitted <- fData(jscs)$dispersion[sel];

    ymin <- ((min(py, na.rm=TRUE)));
    ymax <- ((max(py, na.rm=TRUE)));
    xmin <- ((min(px, na.rm=TRUE)));
    xmax <- ((max(px, na.rm=TRUE)));

    if(verbose) message("     abundance ranges from ",xmin, " to ",xmax);
    if(verbose) message("     dispersion ranges from ",ymin, " to ",ymax);
    
    xmin <- max(xmin,1);
    ymin <- quantile(py, 0.025,na.rm=TRUE);
    
    if(missing(ylim)){
       ylim <- c(ymin,ymax);
    }
    if(missing(xlim)){
       xlim <- c(xmin,xmax);
    }
    if(verbose) message("     Plotting dispersions from ",ymin, " to ",ymax);
    ylim <- log10(ylim);
    xlim <- log10(xlim);

    logscaleList <- build.log.scale(min(xlim[1],ylim[1]), max(xlim[2],ylim[2]));
    decade.at <- logscaleList[["decade.at"]];
    decade.labels <- logscaleList[["decade.labels"]];
    ticks.at <- logscaleList[["ticks.at"]];
    
    #decade.min <- min(floor(xlim[1]),floor(ylim[1])) - 1;
    #decade.max <- max(ceiling(ylim[2]),ceiling(xlim[2])) + 1;
    #decade.at <- decade.min:decade.max;
    #decade.labels <- as.expression(sapply(decade.at,function(yda){
    #  substitute(10 ^ x, list(x = yda));
    #}));
    #ticks.at <- unlist(lapply(decade.at, function(D){
    #  log10( (10 ^ D) * 2:9 );
    #}));

    #message("ylim = c(",paste0(ylim,collapse=","),")");
    par(cex = par.cex, mar = mar);
    if(use.smoothScatter){
      if( jscs@dispFunctionType[["finalDispersionMethod"]] == "shrink"){
        stop("Cannot use smoothscatter with method 'shrink'!");
      }
      smoothScatter(log10(px),log10(py),nbin = smooth.nbin, nrpoints = nrpoints,
                     xlim = xlim, ylim=ylim,xlab="",ylab="",axes=FALSE, xaxs="i", yaxs="i",  ...);
    } else {
      if( (! is.null(jscs@dispFunctionType[["finalDispersionMethod"]]))  && jscs@dispFunctionType[["finalDispersionMethod"]] == "shrink"){
        plot(log10(px),log10(py),
                       xlim = xlim, ylim=ylim,xlab="",ylab="",axes=FALSE, xaxs="i", yaxs="i", pch = pch.MLE, cex = points.cex, col = pchColor,
                       ...);
        points(log10(px),log10(pyFitted), pch = pch.MAP, cex = points.cex/2, col = pchColor, ...);
      } else {
        plot(log10(px),log10(py),
                     xlim = xlim, ylim=ylim,xlab="",ylab="",axes=FALSE, xaxs="i", yaxs="i", pch = pch.MAP, cex = points.cex, ...);
      }
    }

    box(lwd = anno.lwd,...);
    axis(1, at = decade.at, labels = decade.labels, tcl = -0.5, las = 1, lwd = anno.lwd, lwd.ticks = anno.lwd, ...);
    if(miniTicks) axis(1, at = ticks.at, labels = FALSE, tcl = -0.25, lwd = anno.lwd, lwd.ticks = anno.lwd / 2,  ...);
    axis(2, at = decade.at, labels = decade.labels, tcl = -0.5, las = 1, lwd = anno.lwd, lwd.ticks = anno.lwd, ...);
    if(miniTicks) axis(2, at = ticks.at, labels = FALSE, tcl = -0.25, lwd = anno.lwd, lwd.ticks = anno.lwd / 2, ...);
    axis(4, at = decade.at, labels = FALSE, tcl = -0.5, las = 1, lwd = anno.lwd, lwd.ticks = anno.lwd, ...);
    if(miniTicks) axis(4, at = ticks.at, labels = FALSE, tcl = -0.25, lwd = anno.lwd, lwd.ticks = anno.lwd / 2, ...);
    
    if(show.legends){
      if((! is.null(jscs@dispFunctionType[["fitDispersionsForExonsAndJunctionsSeparately"]]))  && jscs@dispFunctionType[["fitDispersionsForExonsAndJunctionsSeparately"]]){
           #Exon / Junction Labels (only if both exons and junctions are found and tested:
           legend("bottomright",legend=c("Exons","Junctions"),text.col=linecol, bg="transparent", box.col="transparent")
      }
      legend.pt.cex <- if(pch.MLE == 46){points.cex * 5} else {points.cex}
      if( (! is.null(jscs@dispFunctionType[["finalDispersionMethod"]]))  &&  jscs@dispFunctionType[["finalDispersionMethod"]] == "shrink"){
         legend("topright", bg="transparent", box.col="transparent", seg.len=1, cex = text.cex,
                  legend=c("MLE","Fitted","MAP"),
                  lty =  c(NA,1,NA), 
                  lwd =  c(NA,lwd.fitted,NA),
                  pch =  c(pch.MLE,NA,pch.MAP),
                  pt.cex = c(legend.pt.cex,NA,points.cex/2))
      }
    }
    
    title(main = title.main, xlab = xlab, ylab = ylab , cex.main = text.cex * 1.2, cex.lab = text.cex,...);

    log.xg = seq( min(decade.at), max(decade.at), length.out=200 );
    xg <- 10 ^ log.xg;
    
    if( (! is.null(jscs@dispFunctionType[["fitDispersionsForExonsAndJunctionsSeparately"]])) ){
      if(jscs@dispFunctionType[["fitDispersionsForExonsAndJunctionsSeparately"]]){
        lines(log.xg, log10(jscs@dispFunctionExon(xg)) , col=linecol[1], lwd=lwd.fitted, cex = lines.cex, ...);
        lines(log.xg, log10(jscs@dispFunctionJct(xg)) , col=linecol[2], lwd=lwd.fitted, cex = lines.cex, ...);
      } else {
        lines(log.xg, log10(jscs@dispFunction(xg)) , col=linecol[1], lwd=lwd.fitted, cex = lines.cex, ...);
      }
    }
    
    if(! is.null(attr(jscs, "filterThreshold"))){
      filterThreshold <- attr(jscs, "filterThreshold");
      if(filterThreshold > 0){
        abline(v = log10(filterThreshold),col="gray",lty=3, lwd = lwd.fitted, ...);
      }
    }
}

build.log.scale <- function(logmin, logmax){
    loglim <- c(logmin,logmax);
    decade.min <- floor(loglim[1]) - 1;
    decade.max <- ceiling(loglim[2]) + 1;
    decade.at <- decade.min:decade.max;
    decade.labels <- as.expression(unlist(sapply(decade.at,function(yda){
      if(yda == 0){
        expression(1);
      } else if(yda == 1){
        expression(10);
      } else {
        substitute(10 ^ x, list(x = yda));
      }
    })));
    ticks.at <- unlist(lapply(decade.at, function(D){
      log10( (10 ^ D) * 2:9 );
    }));
    return(list(decade.labels = decade.labels, decade.at = decade.at, ticks.at = ticks.at));
}

plotMA <- function(jscs, 
                           FDR.threshold = 0.05, 
                           fc.name = NULL,
                           fc.thresh = 1,
                           use.pch = 20,
                           smooth.nbin = 256, 
                           ylim = c( 1 / 1000,1000),
                           use.smoothScatter = TRUE,
                           label.counts = TRUE,
                           label.axes = c(TRUE,TRUE,FALSE,FALSE), 
                           show.labels = TRUE,
                           par.cex = 1, points.cex = 1, text.cex = 1,
                           lines.cex = 8,
                           anno.lwd = 2,
                           mar = c(4.1,4.1,3.1,1.1),
                           miniTicks = TRUE, 
                            verbose = TRUE, debug.mode = FALSE,
                            ...){

  sig.thresh <- FDR.threshold;
  color.orange <- "red";
  fdata <- fData(jscs);
  
  if(is.null(fc.name)){
    fc.names <- names(fdata)[which(strStartsWith(names(fdata), "log2FC("))];
    if(length(fc.names) != 1){
      stop("Error: you must specify the name of the column containing the desired fold change variable.");
    }
    fc.name <- fc.names[1];
  }
  if( strStartsWith(fc.name, "log2FC(") ){
    lvlsString <- substr(fc.name, 8, nchar(fc.name) - 1);
    lvls <- strsplit(lvlsString,"/",fixed=TRUE)[[1]];
    MAplot.label <- paste0(" (",lvlsString,")");
  } else {
    lvls <- NULL;
    label.counts <- FALSE;
    MAplot.label <- paste0("");
  }
  
  X <- fData(jscs)$baseMean;
  Y <- fData(jscs)[[fc.name]];
  pvals <- fData(jscs)[["padjust"]];
  
  xlim <- c(log10(1),log10( max(X,na.rm=TRUE) ))
  
  if(is.null(ylim)){
    ymax <- 2 ^ max(abs(Y),na.rm=TRUE);
    ylim <- c(1 / ymax, ymax);
    
    ymax.exp <- 2 ^ ymax;
  }
  ylim.exp <- log2( ylim );
  
  markerlines.lty <- "dashed";
  
  point.color <- sapply(pvals,FUN=function(p){
    if(is.na(p)){
      "grey"
    } else if(p < sig.thresh){
      "red"
    } else {
      "black"
    }
  });
  point.pch <- sapply(Y,function(fc){
    if(is.na(fc)){
      46
    } else if(fc > ylim.exp[2]){
      94
    } else if(fc < ylim.exp[1]){
      118
    } else {
      use.pch
    }
  })
  limited.log2fc <- sapply(Y,function(fc){
    if(is.na(fc)){
      fc
    } else if(fc > ylim.exp[2]){
      ylim.exp[2]
    } else if(fc < ylim.exp[1]){
      ylim.exp[1]
    } else {
      fc
    }
  })

  point.color[point.color == "red"] <- ifelse(f.na((2 ^ abs(limited.log2fc)) < fc.thresh)[point.color == "red"], color.orange,"red");
  is.sig <- f.na(pvals < sig.thresh);
  is.over <- f.na(Y > ylim.exp[2] | Y < ylim.exp[1]);

  par(cex = par.cex, mar = mar);
  if(use.smoothScatter){
    smoothScatter(log10(X),Y,nbin = smooth.nbin, nrpoints = 0,
                 xlim = xlim,ylim=ylim.exp,xlab="",ylab="",axes=FALSE, ...);
    points(log10(X[is.over & (! is.sig)]),limited.log2fc[is.over & (! is.sig)], col= point.color[is.over & (! is.sig)],pch= point.pch[is.over & (! is.sig)],cex= points.cex,  ...);
  } else {
    plot(log10(X[! is.sig]),limited.log2fc[! is.sig], col= point.color[! is.sig], pch= point.pch[! is.sig],cex= points.cex,
                 xlim = xlim,ylim=ylim.exp,xlab="",ylab="",axes=FALSE, ...);
  }
  points(log10(X[is.over & (is.sig)]),limited.log2fc[is.over & (is.sig)], col= point.color[is.over & (is.sig)],pch= point.pch[is.over & (is.sig)],cex= points.cex,  ...);
  points(log10(X[(is.sig & !is.over)]),limited.log2fc[(is.sig & !is.over)], col= point.color[(is.sig & !is.over)],pch= point.pch[(is.sig & !is.over)],cex= points.cex/2, ...);
  
  abline(h=0,col=rgb(255,0,0,100,maxColorValue=255), lwd = anno.lwd, cex = lines.cex)

  box(lwd = anno.lwd);

  if(miniTicks){
    log10.label <- c(-3,-2,-1,0,1,2,3);
    log10.expression.label <- c(expression(10^"-3"), expression(10^"-2"), expression(10^"-1"), 1, expression(10), expression(10^2), expression(10^3));
  } else {
    log10.label <- c(-3,-2,-1,-log10(4),0,log10(4),1,2,3);
    log10.expression.label <- c(expression(10^"-3"), expression(10^"-2"), expression(10^"-1"),expression(1/4), 1, 4, expression(10), expression(10^2), expression(10^3));
  }
  #log10.expression.label <- c(expression(10^-3), expression(10^-2), expression(10^-1),"\u00BC", 1, 4, expression(10), expression(10^2), expression(10^3));
  log10.at <- log10.label * log2(10);

  miniTicks.temp <- c(2,3,4,5,6,7,8,9)
  log10.miniTicks <- log10(c(miniTicks.temp, miniTicks.temp * 10, miniTicks.temp * 100, 1 / miniTicks.temp, 1 / (miniTicks.temp * 10), 1 / (miniTicks.temp * 100) ));
  miniTicks.pts <- log10.miniTicks * log2(10);

  #axis(1);
  if(miniTicks){
     x.at <- log10(c(1,10,100,1000,10000,100000));
     x.label <- c(1, expression(10), expression(10^2), expression(10^3), expression(10^4), expression(10^5));
  } else {
     x.at <- log10(c(1,5,10,100,1000,10000,100000));
     x.label <- c(1,5, expression(10), expression(10^2), expression(10^3), expression(10^4), expression(10^5));
  }

  miniTicks.x <- log10(c(miniTicks.temp, miniTicks.temp * 10, miniTicks.temp * 100, miniTicks.temp * 1000, miniTicks.temp * 10000, miniTicks.temp * 100000, miniTicks.temp * 1000000));

  #axis(2);
  #axis(2,at=log10.at,label=log10.label);

  if(label.axes[1]){
    axis(1,at=x.at,labels=x.label, lwd = anno.lwd, lwd.ticks = anno.lwd,  cex.axis = text.cex);
  } else {
    axis(1,at=x.at,labels=FALSE, lwd = anno.lwd, lwd.ticks = anno.lwd, cex.axis = text.cex);  
  }
  if(miniTicks ) axis(1,at=miniTicks.x,labels=FALSE, lwd = anno.lwd / 2, lwd.ticks = anno.lwd / 2, tcl = -0.25, cex.axis = text.cex);
  if(label.axes[2]){
    axis(2,at=log10.at,labels=log10.expression.label,las=1, lwd = anno.lwd, lwd.ticks = anno.lwd, cex.axis = text.cex);
    if(miniTicks ) axis(2,at = miniTicks.pts, labels = FALSE, lwd = anno.lwd/2, lwd.ticks = anno.lwd / 2, tcl = -0.25)
  } else {
    axis(2,at=log10.at,labels=FALSE,las=1, lwd = anno.lwd, lwd.ticks = anno.lwd, tcl = -0.25, cex.axis = text.cex);
    if(miniTicks ) axis(2,at = miniTicks.pts, labels = FALSE, lwd = anno.lwd/2, lwd.ticks = anno.lwd / 2, tcl = -0.2)
  }
  if(label.axes[4]){
    axis(4,at=log10.at,labels=log10.expression.label,las=1, lwd = anno.lwd, lwd.ticks = anno.lwd, cex.axis = text.cex);
    if(miniTicks ) axis(4,at = miniTicks.pts, labels = FALSE, lwd = anno.lwd/2, lwd.ticks = anno.lwd / 2, tcl = -0.25)
  } else {
    axis(4,at=log10.at,labels=FALSE,las=1, lwd = anno.lwd, lwd.ticks = anno.lwd, tcl = -0.25, cex.axis = text.cex);
    if(miniTicks ) axis(4,at = miniTicks.pts, labels = FALSE, lwd = anno.lwd/2, lwd.ticks = anno.lwd / 2, tcl = -0.2) 
  }

  box(lwd=anno.lwd,cex=lines.cex);
  
  if(label.counts){
    if(is.na(fc.thresh)){
      sig.count.plus <- sum(f.na(Y > 0 & pvals < sig.thresh));
      sig.count.minus <- sum(f.na(Y < 0 & pvals < sig.thresh));
    } else {
      sig.count.plus <- sum(f.na(Y > log2(fc.thresh) & pvals < sig.thresh));
      sig.count.minus <- sum(f.na(Y < -log2(fc.thresh) & pvals < sig.thresh));
      
      if(fc.thresh != 1){
        abline(h=log2(fc.thresh),col="gray",lty=markerlines.lty, lwd = anno.lwd, cex = lines.cex);
        abline(h=-log2(fc.thresh),col="gray",lty=markerlines.lty, lwd = anno.lwd, cex = lines.cex);
      }
    }
    fc.label.1 <- if(fc.thresh == 1){ "" } else { paste0(" (FC > ",fc.thresh,")" )}
    fc.label.2 <- if(fc.thresh == 1){ "" } else { paste0(" (FC < ",fc.thresh,")" )}
    text(par("usr")[2],ylim.exp[2],paste0(sig.count.plus,  " higher in ",lvls[1],fc.label.1), cex = text.cex, adj = c(1.1, 1.5),...);
    text(par("usr")[2],ylim.exp[1],paste0(sig.count.minus, " higher in ",lvls[2],fc.label.2), cex = text.cex, adj = c(1.1,-0.5),...);
  }
  
  if(show.labels){
    title(main = paste0("MA Plot",MAplot.label),
          xlab = "Mean Normalized Coverage",
          ylab = "Fold Change",
          cex.main = 1.2 * text.cex,
          cex.lab = 1 * text.cex
          );
  }
}


plotMA.OLDVER <- function(jscs, FDR.threshold, fc.name = NULL, ylim = c(-5,5), ...){
    fdata <- fData(jscs);
    
    if(is.null(fc.name)){
      fc.names <- names(fdata)[which(strStartsWith(names(fdata), "log2FC("))];
      if(length(fc.names) != 1){
        stop("Error: you must specify the name of the column containing the desired fold change variable.");
      }
      fc.name <- fc.names[1];
    }

    #which.is.fc <- which(substr(names(merged.results.data),1,8) == "log2fold")
    #fc.name <- names(merged.results.data)[fc.col];
    ma.frame <- data.frame(baseMean=fdata$baseMean);
    ma.frame[[fc.name]] <- fdata[[fc.name]];
    ma.frame[["padj"]] <- fdata$padjust;
    
    #, log2FoldChange = merged.results.data[,which.is.fc], padj = merged.results.data$padjust);
    #CairoPNG(paste(outfile.prefix,"ma.plot",".png",sep=""),height=2500,width=1600,pointsize = 18);
    plotMA_HELPER( ma.frame, FDR=FDR.threshold, ylim=ylim, ylab = fc.name, ...)
    #dev.off();
}


plotMA_HELPER <- function(ma.frame, ylim,FDR=0.05,
  col = ifelse(ma.frame$padj>=FDR, "gray32", "red3"),
  linecol = "#ff000080",
  xlab = "mean of normalized counts", ylab = expression(log[2]~fold~change),
  log = "x", cex=0.45, ...)
{
  #if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x))))
  #  stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")

  ma.frame = subset(ma.frame, ma.frame$baseMean != 0)
  py = ma.frame[,2];
  if(missing(ylim)){
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  }
  plot(ma.frame$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  abline(h=0, lwd=4, col=linecol);
}