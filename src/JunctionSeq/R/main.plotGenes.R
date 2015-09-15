
##########################################################################
######### Generating Plots From Results:
##########################################################################

buildAllPlots <- function(jscs,
                          outfile.prefix = "./",
                          #flat.gff.data = NULL, flat.gff.file = NULL, 
                          gene.list = NULL, FDR.threshold = 0.01, 
                          use.plotting.device = c("png","CairoPNG","svg","tiff","cairo_ps","custom"),
                          sequencing.type = c("paired-end","single-end"),
                          use.vst=FALSE,use.log = TRUE, truncateBelowOne = TRUE, 
                          exon.rescale.factor = 0.3,
                          subdirectories.by.type = TRUE,
                          ma.plot=TRUE, variance.plot=TRUE,
                          with.TX=TRUE,without.TX=TRUE,
                          expr.plot=TRUE,normCounts.plot=TRUE,
                          rExpr.plot=TRUE,rawCounts.plot=FALSE,
                          colorRed.FDR.threshold = FDR.threshold, 
                          color=NULL,
                          plot.gene.level.expression = TRUE,
                          plot.exon.results = NULL, plot.junction.results = NULL, plot.novel.junction.results = NULL,
                          plot.untestable.results = FALSE,
                          plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                          par.cex = 1, anno.cex.text = 1, anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                          drawCoordinates = TRUE, yAxisLabels.inExponentialForm = FALSE,
                          show.strand.arrows = 10, arrows.length = 0.125,
                          graph.margins = c(2,3,3,3),
                          base.plot.height = 12, base.plot.width = 12, 
                          base.plot.units = "in", 
                          GENE.annotation.relative.height = 0.2, TX.annotation.relative.height = 0.05, 
                          autoscale.height.to.fit.TX.annotation = TRUE,
                          autoscale.width.to.fit.bins = 35,
                          plotting.device.params = list(), 
                          number.plots = FALSE, name.files.with.geneID = TRUE,
                          condition.legend.text = NULL, include.TX.names = TRUE, draw.start.end.sites = TRUE,
                          openPlottingDeviceFunc = NULL, closePlottingDeviceFunc = NULL,
                          writeHTMLresults = TRUE,
                          html.cssFile = NULL, html.cssLink = NULL, html.imgFileExtension = NULL,
                          html.plot.height = 90, html.plot.height.units = "vh",
                          verbose=TRUE, debug.mode = FALSE,
                          ...){
  
  
  condition <- jscs@phenoData$condition;
  flat.gff.data <- jscs@flatGffData;
  
  use.plotting.device <- match.arg(use.plotting.device);
  sequencing.type <- match.arg(sequencing.type);
  
  if(is.null(openPlottingDeviceFunc) & is.null(closePlottingDeviceFunc)){
    devFunctions <- getPlottingDeviceFunc(use.plotting.device = use.plotting.device,
                                            base.plot.height = base.plot.height,
                                            base.plot.width = base.plot.width,
                                            base.plot.units = base.plot.units,
                                            plotting.device.params = plotting.device.params);
    openPlottingDeviceFunc <- devFunctions[[1]];
    closePlottingDeviceFunc <- devFunctions[[2]];
    
    if(is.null(html.imgFileExtension)){
      html.imgFileExtension <- getPlottingDeviceFileExtension(use.plotting.device);
    }
  } else if(use.plotting.device == "custom"){
      stop("Fatal error: custom plotting device is selected, but openPlottingDeviceFunc or closePlottingDeviceFunc is not set.");
  }
  
  
  gtf.format <- TRUE;
  
  if(is.null(gene.list)){
    sig.features <- which(fData(jscs)$padjust < FDR.threshold);
    gene.list <- unique(as.character(fData(jscs)$geneID[sig.features]));
    if(verbose) message("> buildAllPlots: Found ", length(gene.list), " significant genes to plot, at adjusted-p-value threshold ", FDR.threshold);
  } else {
    if(verbose) message("> buildAllPlots: Found ", length(gene.list), " genes to plot.");
  }
  
  if(length(gene.list) > 0){
    #if(is.null(flat.gff.data)){
    #  if(is.null(flat.gff.file)){
    #    stop("FATAL ERROR: buildAllPlots: either flat.gff.file or flat.gff.data must be specified!");
    #  }
    #  flat.gff.data <- readAnnotationData(flat.gff.file);
    #}

    #numcond <- length(levels(condition));

    if(verbose) message("> buildAllPlots: Starting plotting...");

    FDR <- colorRed.FDR.threshold

    total.gene.ct <- length(gene.list);
    number.width <- floor(log10(total.gene.ct)) + 1;
    if( number.plots ){
      geneNum.strings <- paste0(formatC(1:total.gene.ct, width = number.width, format = 'd', flag = '0'), "-");
    } else {
      geneNum.strings <- rep("",total.gene.ct);
    }

    #if(verbose){
    #   message(paste0("total gene ct: ",total.gene.ct,", column width: ",number.width));
    #}
    
    if(subdirectories.by.type){
      if(! file.exists(outfile.prefix)){
        dir.create(outfile.prefix);
      }
      
      if(expr.plot & with.TX          & (! file.exists(paste0(outfile.prefix,"/exprTX"))))       dir.create(paste0(outfile.prefix,"/exprTX"));
      if(expr.plot & without.TX       & (! file.exists(paste0(outfile.prefix,"/expr"))))         dir.create(paste0(outfile.prefix,"/expr"));
      if(normCounts.plot & with.TX    & (! file.exists(paste0(outfile.prefix,"/normCountsTX")))) dir.create(paste0(outfile.prefix,"/normCountsTX"));
      if(normCounts.plot & without.TX & (! file.exists(paste0(outfile.prefix,"/normCounts"))))   dir.create(paste0(outfile.prefix,"/normCounts"));
      if(rExpr.plot & with.TX         & (! file.exists(paste0(outfile.prefix,"/rExprTX"))))      dir.create(paste0(outfile.prefix,"/rExprTX"));
      if(rExpr.plot & without.TX      & (! file.exists(paste0(outfile.prefix,"/rExpr"))))        dir.create(paste0(outfile.prefix,"/rExpr"));
      if(rawCounts.plot & with.TX     & (! file.exists(paste0(outfile.prefix,"/rawCountsTX"))))  dir.create(paste0(outfile.prefix,"/rawCountsTX"));
      if(rawCounts.plot & without.TX  & (! file.exists(paste0(outfile.prefix,"/rawCounts"))))    dir.create(paste0(outfile.prefix,"/rawCounts"));
    }
    
    if(variance.plot){
      openPlottingDeviceFunc(paste(outfile.prefix,"dispersion-plot","",sep=""),heightMult=0.6,widthMult=0.6);
      plotDispEsts( jscs, par.cex = par.cex, points.cex = anno.cex.text, text.cex = anno.cex.text, ... );
      dev.off();
    }

    if(ma.plot){
      #which.is.fc <- which(substr(names(merged.results.data),1,8) == "log2fold")
      ##print("which is fc:");
      ##print(substr(names(merged.results.data),1,8));
      #ma.frame <- data.frame(baseMean=merged.results.data$baseMean, log2FoldChange = merged.results.data[,which.is.fc], padj = merged.results.data$padjust);

      fc.cols <- which(strStartsWith(names(fData(jscs)), "log2FC("));

      for(fc.name in names(fData(jscs))[fc.cols]){
        if(verbose) message("> buildAllPlots: Generating MA-Plot (",fc.name,")");

        fc.title <- sub("/","vs",fc.name, fixed = TRUE);
        openPlottingDeviceFunc(paste(outfile.prefix,"ma-plot-",fc.title,sep=""),heightMult=0.6,widthMult=0.6);
        plotMA( jscs, FDR.threshold=colorRed.FDR.threshold, fc.name = fc.name, par.cex = par.cex, text.cex = anno.cex.text, points.cex = anno.cex.text, ... );
        dev.off();
      }
    }
    

    if(writeHTMLresults){
      if(verbose) message("> buildAllPlots: Writing HTML results index.");
      if((! file.exists(paste0(outfile.prefix,"/htmlFiles")))) dir.create(paste0(outfile.prefix,"/htmlFiles"));
      
      JunctionSeqHTML(jscs = jscs,
                            outfile.dir=outfile.prefix, mainFile="testForDU.html", 
                            gene.list = gene.list,
                            FDR.threshold = FDR.threshold,
                            colorRed.FDR.threshold = colorRed.FDR.threshold,
                            plotting.device.ext = html.imgFileExtension, 
                            use.vst= use.vst,use.log = use.log,
                            subdirectories.by.type = subdirectories.by.type,
                            ma.plot=ma.plot, variance.plot=variance.plot,
                            with.TX=with.TX,without.TX=without.TX,
                            expr.plot=expr.plot,normCounts.plot=normCounts.plot,
                            rExpr.plot=rExpr.plot,rawCounts.plot=rawCounts.plot,
                            plot.exon.results = plot.exon.results, plot.junction.results = plot.junction.results, plot.novel.junction.results = plot.novel.junction.results,
                            plot.untestable.results = plot.untestable.results,
                            base.html.height = html.plot.height, base.html.height.units = html.plot.height.units,
                            GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, 
                            autoscale.height.to.fit.TX.annotation = autoscale.height.to.fit.TX.annotation,
                            autoscale.width.to.fit.bins = autoscale.width.to.fit.bins,
                            number.plots = geneNum.strings,
                            css.file = html.cssFile, css.link = html.cssLink);
     if(verbose) message("> buildAllPlots: Finished writing HTML results index.");
    }
    
    geneNum <- 1;
    for(geneID in gene.list){
      message(paste("> buildAllPlots: starting geneID:",geneID));
      geneNum.string <- geneNum.strings[geneNum];
      
      if(subdirectories.by.type){
        outfile.prefixes <- paste0(outfile.prefix, c("/exprTX/","/expr/",
                                                    "/normCountsTX/","/normCounts/",
                                                    "/rExprTX/","/rExpr/",
                                                    "/rawCountsTX/","/rawCounts/"),geneNum.string);
      } else {
        outfile.prefixes <- paste0(outfile.prefix, geneNum.string);
      }
      
      buildAllPlotsForGene(geneID = geneID, jscs = jscs,
                            outfile.prefix = outfile.prefixes,
                            #flat.gff.data = flat.gff.data, 
                            use.plotting.device = use.plotting.device,
                            use.vst=use.vst, use.log = use.log, truncateBelowOne = truncateBelowOne,
                            exon.rescale.factor = exon.rescale.factor, plot.gene.level.expression = plot.gene.level.expression,
                            with.TX=with.TX,without.TX=without.TX,
                            expr.plot=expr.plot,normCounts.plot=normCounts.plot,
                            rExpr.plot=rExpr.plot,rawCounts.plot=rawCounts.plot,
                            colorRed.FDR.threshold = colorRed.FDR.threshold, 
                            color= color,
                            plot.exon.results = plot.exon.results, 
                            plot.junction.results = plot.junction.results,
                            plot.novel.junction.results = plot.novel.junction.results,
                            plot.untestable.results = plot.untestable.results,
                            plot.lwd = plot.lwd, axes.lwd = axes.lwd, anno.lwd = anno.lwd, 
                            drawCoordinates = drawCoordinates,
                            arrows.length = arrows.length, show.strand.arrows = show.strand.arrows,
                            graph.margins = graph.margins,
                            par.cex = par.cex,
                            anno.cex.text = anno.cex.text,
                            anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,
                            base.plot.height = base.plot.height, base.plot.width = base.plot.width,
                            base.plot.units = base.plot.units,
                            GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height,
                            autoscale.height.to.fit.TX.annotation = autoscale.height.to.fit.TX.annotation,
                            autoscale.width.to.fit.bins = autoscale.width.to.fit.bins,
                            name.files.with.geneID = name.files.with.geneID,
                            plotting.device.params = plotting.device.params,
                            yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm,
                            condition.legend.text = condition.legend.text, include.TX.names = include.TX.names, draw.start.end.sites=draw.start.end.sites,
                            verbose=verbose, debug.mode = debug.mode,
                            sequencing.type=sequencing.type,
                            ...);

      geneNum <- geneNum + 1;
    }
    if(verbose) message("> buildAllPlots: Plotting complete.");
  }
  if(verbose) message("> buildAllPlots: Plotting and data writing complete.");
  
  
}
#html.cssFile = NULL, html.cssLink = NULL, html.imgFileExtension = NULL,

##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
######### Build one gene's plot:
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################

#pointsize = 18, res = 150
buildAllPlotsForGene <- function(geneID,jscs,
                          outfile.prefix = "./", 
                          #flat.gff.data = NULL, flat.gff.file = NULL,
                          use.plotting.device = c("png","CairoPNG","svg","tiff","cairo_ps","custom"),
                          sequencing.type = c("paired-end","single-end"),
                          use.vst=FALSE, use.log = TRUE, truncateBelowOne = TRUE,
                          exon.rescale.factor = 0.3,   
                          with.TX=TRUE,without.TX=TRUE,
                          expr.plot=TRUE, normCounts.plot=TRUE,
                          rExpr.plot=TRUE, rawCounts.plot=FALSE,
                          colorRed.FDR.threshold = 0.01, 
                          color=NULL,
                          plot.gene.level.expression = TRUE,
                          plot.exon.results = NULL, plot.junction.results = NULL, plot.novel.junction.results = NULL,
                          plot.untestable.results = FALSE,
                          plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                          par.cex = 1, 
                          name.files.with.geneID = TRUE,
                          anno.cex.text = 1,
                          anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                          drawCoordinates = TRUE, yAxisLabels.inExponentialForm = FALSE,
                          show.strand.arrows = 10, arrows.length = 0.125,
                          graph.margins = c(2,3,3,3),
                          base.plot.height = 12, base.plot.width = 12, 
                          base.plot.units = "in", 
                          GENE.annotation.relative.height = 0.2, TX.annotation.relative.height = 0.05,
                          autoscale.height.to.fit.TX.annotation = TRUE,
                          autoscale.width.to.fit.bins = 35,
                          plotting.device.params = list(),
                          condition.legend.text = NULL, include.TX.names = TRUE, draw.start.end.sites = TRUE,
                          openPlottingDeviceFunc = NULL, closePlottingDeviceFunc = NULL,
                          verbose=TRUE, debug.mode = FALSE,  ...){
    message(paste("starting buildAllPlotsForGene() for geneID:",geneID));
    condition <- jscs@phenoData$condition
    flat.gff.data <- jscs@flatGffData;
    
    sequencing.type <- match.arg(sequencing.type);
    use.plotting.device <- match.arg(use.plotting.device);
    
    if(length(outfile.prefix) == 1){
      outfile.prefix <- rep(outfile.prefix, 8);
    }
    
    #if("HIDDENPARAM_FIX_SVG_RHEL5" %in% names(plotting.device.params)){
    #   closeDeviceFunc <- function(filename){
    #     dev.off();
    #     #Hack (doesn't work), attempting to workaround a flaw in one particular linux distribution where svg files are rendered incorrectly.
    #     #svgfixfile <- system.file("extdata/FIX_SVG.pl", package="JunctionSeq", mustWork=TRUE);
    #     #svgfixcommand <- paste("perl ",svgfixfile," ", filename, " > ", filename);
    #     #system(svgfixcommand);
    #   };
    #   plotting.device.params <- plotting.device.params[- which(names(plotting.device.params) == "HIDDENPARAM_FIX_SVG_RHEL5")];
    #} else {
    #   closeDeviceFunc <- function(filename){
    #     dev.off();
    #   }
    #}
    #if(is.null(openPlottingDeviceFunc)){
    #    openPlottingDeviceFunc <- getPlottingDeviceFunc(use.plotting.device = use.plotting.device,
    #                                          base.plot.height = base.plot.height,
    #                                          base.plot.width = base.plot.width,
    #                                          base.plot.units = base.plot.units,
    #                                          plotting.device.params = plotting.device.params);
    #}
    
    if(is.null(openPlottingDeviceFunc) & is.null(closePlottingDeviceFunc)){
      devFunctions <- getPlottingDeviceFunc(use.plotting.device = use.plotting.device,
                                              base.plot.height = base.plot.height,
                                              base.plot.width = base.plot.width,
                                              base.plot.units = base.plot.units,
                                              plotting.device.params = plotting.device.params);
      openPlottingDeviceFunc <- devFunctions[[1]];
      closePlottingDeviceFunc <- devFunctions[[2]];
    } else if(use.plotting.device == "custom"){
      stop("Fatal error: custom plotting device is selected, but openPlottingDeviceFunc or closePlottingDeviceFunc is not set.");
    }
    
    
    
    FDR <- colorRed.FDR.threshold;
    #if(is.null(flat.gff.data)){
    #   if(is.null(flat.gff.file)){
    #      stop("FATAL ERROR: buildAllPlotsForGene: either flat.gff.file or flat.gff.data must be specified!");
    #   }
    #   flat.gff.data <- readAnnotationData(flat.gff.file);
    #}
    
    transcripts <- sapply(sapply(flat.gff.data$transcripts[which(flat.gff.data$gene_id==geneID & flat.gff.data$featureType == "exonic_part")],toString), function(x){strsplit(x, "+",fixed=TRUE)})
    trans <- Reduce(union, transcripts)
    tx.ct <- length(trans);
    
    if(autoscale.height.to.fit.TX.annotation){
      GENE.annotation.height <- GENE.annotation.relative.height * 10;
      TX.annotation.height <- TX.annotation.relative.height * 10;
      withTxPlot.height.multiplier <- (10+1+GENE.annotation.height+TX.annotation.height*tx.ct) / (10+1+GENE.annotation.height);
    } else {
      withTxPlot.height.multiplier <- 1;
    }
    
    if(is.na(autoscale.width.to.fit.bins) | autoscale.width.to.fit.bins == 0){
      width.multiplier <- 1;
    } else {
      num.cols <- getColCt(geneID = geneID, merged.data = fData(jscs), 
                                plot.exon.results = plot.exon.results, plot.junction.results = plot.junction.results, plot.novel.junction.results=plot.novel.junction.results, 
                                plot.untestable.results=plot.untestable.results);
      if(num.cols > autoscale.width.to.fit.bins){
        width.multiplier <- num.cols / autoscale.width.to.fit.bins;
      } else {
        width.multiplier <- 1;
      }
    }
    
    if(verbose & debug.mode){
      message("Found ",tx.ct," TX");
      message("Final Dimensions,   No TX:",base.plot.width * width.multiplier,"x",base.plot.height);
      message("Final Dimensions, With TX:",base.plot.width * width.multiplier,"x",base.plot.height * withTxPlot.height.multiplier);
    }
    
    geneName <- if(name.files.with.geneID){
      geneID;
    } else {
      jscs@flatGffGeneData$gene_name[jscs@flatGffGeneData$geneID == geneID];
    }
    
    
    if(expr.plot){
      plot.type <- "expr"
      if(with.TX){
        outfile <- paste(outfile.prefix[1],geneName,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs, colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd , par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results , plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[2],geneName,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,...)
        closePlottingDeviceFunc();  
      }
    }
    if(normCounts.plot){
      plot.type <- "normCounts"
      if(with.TX){
        outfile <- paste(outfile.prefix[3],geneName,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[4],geneName,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,...)
        closePlottingDeviceFunc();  
      }
    }
    if(rExpr.plot){
      plot.type <- "rExpr"
      if(with.TX){
        outfile <- paste(outfile.prefix[5],geneName,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[6],geneName,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,...)
        closePlottingDeviceFunc();  
      }
    }
    if(rawCounts.plot){
      plot.type <- "rawCounts"
      if(with.TX){
        outfile <- paste(outfile.prefix[7],geneName,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=FALSE,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[8],geneName,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=FALSE,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,...)
        closePlottingDeviceFunc();  
      }
    }
}

##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
##############################################################################################################################################################################################################################
######### Gene Plot Engine:

USE.MARGIN.MEX <- FALSE;

plotJunctionSeqResultsForGene <- function(geneID, jscs, 
                              colorRed.FDR.threshold=0.05,
                              plot.type = "expr", 
                              sequencing.type = c("paired-end","single-end"),
                              displayTranscripts = FALSE,
                              color = NULL, 
                              use.vst = FALSE, use.log = TRUE,truncateBelowOne = TRUE,
                              exon.rescale.factor = 0.3,
                              label.p.vals = TRUE, 
                              plot.lwd = 3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                              par.cex = 1, 
                              #cex.master = 1,
                              #cex.axis = cex.master,
                              #cex.ylabel = cex.master,
                              #cex.main = cex.master * 1.2,
                              #cex.featureLabels = "auto",
                              #cex.pvalues = "auto",
                              #cex.TxLabels = "auto",
                              anno.cex.text = 1,
                              anno.cex.axis=anno.cex.text, anno.cex.main = anno.cex.text * 1.2, cex.arrows = 1,
                              fit.countbin.names = TRUE,
                              plot.gene.level.expression = TRUE,
                              plot.exon.results = NULL, plot.junction.results = NULL, plot.novel.junction.results = NULL, 
                              plot.untestable.results = FALSE, draw.untestable.annotation = TRUE,
                              show.strand.arrows = 10, arrows.length = 0.125,
                              sort.features = TRUE,
                              drawCoordinates = TRUE,
                              yAxisLabels.inExponentialForm = FALSE,
                              title.main=NULL, title.ylab=NULL, title.ylab.right=NULL, 
                              graph.margins = c(2,3.5,3,3),
                              GENE.annotation.relative.height = 0.2, TX.annotation.relative.height = 0.05,
                              condition.legend.text = NULL, include.TX.names = TRUE, draw.start.end.sites = TRUE,
                              label.chromosome=TRUE, 
                              verbose=TRUE, debug.mode = FALSE, 
                              ...)
{
    tryCatch({
      GENE.annotation.height <- GENE.annotation.relative.height * 10;
      TX.annotation.height <- TX.annotation.relative.height * 10;
      flat.gff.data <- jscs@flatGffData;
      
      if(is.null(plot.exon.results)){
        plot.exon.results <- any( fData(jscs)$featureType == "exonic_part" );
      }
      if(is.null(plot.junction.results)){
        plot.junction.results <- any( fData(jscs)$featureType == "splice_site" );
      }
      if(is.null(plot.novel.junction.results)){
        plot.novel.junction.results <- any( fData(jscs)$featureType == "novel_splice_site" );
      }
      
      geneName <- jscs@flatGffGeneData$gene_name[jscs@flatGffGeneData$geneID == geneID];

     if(is.null(condition.legend.text)){
       condition.legend.text <- levels(jscs@phenoData$condition);
       names(condition.legend.text) <- condition.legend.text;
     } else {
       if(is.null(names(condition.legend.text))){
         warning("names(condition.legend.text is NULL! condition.legend.text mis-formatted. Must be a list or character vector, with element names equal to the levels of pData(jscs)$condition. Falling back.");
         condition.legend.text <- levels(jscs@phenoData$condition);
         names(condition.legend.text) <- condition.legend.text;
       } else if(! all(levels(jscs@phenoData$condition) %in% names(condition.legend.text))){
         warning("Not all levels contained in names(condition.legend.text)! condition.legend.text mis-formatted. Must be a list or character vector, with element names equal to the levels of pData(jscs)$condition. Falling back.");
         condition.legend.text <- levels(jscs@phenoData$condition);
         names(condition.legend.text) <- condition.legend.text;
       }
     }
     
     
     geneStrand <-  as.character(jscs@flatGffGeneData[["aggregateGeneStrand"]][ jscs@flatGffGeneData[["geneID"]] == geneID ])
     txSetString <- as.character(jscs@flatGffGeneData[["tx_set"]][ jscs@flatGffGeneData[["geneID"]] == geneID ])
     txStrandString <-  as.character(jscs@flatGffGeneData[["tx_strands"]][ jscs@flatGffGeneData[["geneID"]] == geneID ])
     txSet <-  strsplit(as.character(txSetString),"+",fixed=TRUE)[[1]];
     txStrand <-  strsplit(as.character(txStrandString),",",fixed=TRUE)[[1]];
     txStrandMap <- as.list(txStrand);
     names(txStrandMap) <- txSet;
     
     gene.level.buffer <- 0.5;

     if(is.null(plot.gene.level.expression)){
       #rExpr plot now shows gene-level expr:
       #if(plot.type == "rExpr"){
       #  plot.gene.level.expression <- FALSE;
       #} else 
       if(use.vst){
         plot.gene.level.expression <- FALSE;
       } else {
         plot.gene.level.expression <- TRUE;
       }
     }
     #if(plot.gene.level.expression & plot.type == "rExpr"){
     #  warning("WARNING: plotting of gene-level expression is not supported for relative expression plots (it doesn't make sense to do so). Errors are likely to follow.");
     #}
     if(plot.gene.level.expression & use.vst){
       warning("WARNING: plotting of gene-level expression is not supported for vst-transformed plots. Errors are likely to follow.");
     }
     
     FDR <- colorRed.FDR.threshold;
     merged.data <- fData(jscs);
     condition <- jscs@phenoData$condition;

     if(verbose){
        displayTXstring <- if(displayTranscripts) " (with TX)" else "";
        message("> pJSRfG(): ", geneID, ", plot.type: ", plot.type,displayTXstring);
     }
     
     chrom.label <- as.character(merged.data$chr[1]);
     
     rt <- merged.data$geneID == geneID;
     if(! plot.exon.results){
       if(debug.mode) message(">     Removing ",sum(rt & merged.data$featureType == "exonic_part")," exonic_part features. ",sum(rt & merged.data$featureType != "exonic_part")," features remaining");
       rt <- rt & merged.data$featureType != "exonic_part" ;
     }
     if(! plot.junction.results){
       if(debug.mode) message(">     Removing ",sum(rt & merged.data$featureType == "splice_site")," splice_site features. ",sum(rt & merged.data$featureType != "splice_site")," features remaining");
       rt <- rt & merged.data$featureType != "splice_site" ;
     }
     if(! plot.novel.junction.results){
       if(debug.mode) message(">     Removing ",sum(rt & merged.data$featureType == "novel_splice_site")," novel_splice_site features. ",sum(rt & merged.data$featureType != "novel_splice_site")," features remaining");
       rt <- rt & merged.data$featureType != "novel_splice_site" ;
     }
     untestable.rt <- which(rt & (! merged.data$testable));
     if(! plot.untestable.results){
       if(debug.mode) message(">     Removing ",sum(rt & ! merged.data$testable)," untestable features. ",sum(rt & merged.data$testable)," features remaining");
       rt <- rt & merged.data$testable;
     }
     if(debug.mode)   message(">     Plotting ",sum(rt), " features.");

     rt <- which(rt);
     if(length(rt) == 0){
       message("NO FEATURES TO PLOT!");
     } else {

       if(sort.features){
         #sort by start, then end (depreciated):
         #rt <- rt[order(merged.data$start[rt],merged.data$end[rt])];
         #New method: sort by midpoint. Looks better?
         rt <- rt[order(((merged.data$end[rt] - merged.data$start[rt])/2) + merged.data$start[rt] )];
       }

       rango <- 1:length(rt)
       draw.legend <- TRUE;
       condition.names <- levels(condition);
       sample.names <- sampleNames(jscs@phenoData);
       colorcode.title <- TRUE;
       numcond <- length(condition.names);

       if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached step 2.");

       rt.allExon <- which(flat.gff.data$gene_id==geneID & flat.gff.data$featureType == "exonic_part");
       rango.allExon <- 1:length(rt.allExon);

       rt.allJunction <- which(flat.gff.data$gene_id==geneID & (flat.gff.data$featureType == "splice_site" | flat.gff.data$featureType == "novel_splice_site"));
       rango.allJunction <- 1:length(rt.allJunction);

       ####### DETERMINE COLORS, IF THE USER DOES NOT PROVIDE ONE PER SAMPLE THE COUNT WILL OBTAIN THEM CORRESPONDING TO THEIR DESIGN ####
       ##### determine colors if not provided by user ######
       if(is.null(color)){
          #color<-rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond)), maxColorValue=255, alpha=175)
          default.color.list <- c("red","blue","orange","green3","purple","cyan1", "magenta","yellow3","tan4")

          if(numcond > length(default.color.list)){
             message("Too many condition values, the default color selection may not look good! Set your own colors by setting the \"color\" parameter.");
             color<-rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond)), maxColorValue=255, alpha=175)
          } else {
             color <- color2transparentVector(default.color.list[1:numcond], t = 175);
          }
       }
    if(debug.mode & verbose) message(">    pJSRforGene(): ","color = ", paste0(color, collapse=","));
    if(debug.mode & verbose) message(">    pJSRforGene(): ","length(color) = ", length(color));
    if(debug.mode & verbose) message(">    pJSRforGene(): ","length(condition.names) = ", length(condition.names));
    if(debug.mode & verbose) message(">    pJSRforGene(): ","condition.names = ", paste0(condition.names,collapse=","));

       names(color) <- condition.names;

       y.axis.title <- "";
       main.title <- "";

    if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached step 3.");

       #For reference:
       #jscs@plottingEstimates    <- list(exprEstimate = exprEstimate,
       #                                  geneExprEstimate = otherExprEstimate,
       #                                  relExprEstimate = relExprEstimate, 
       #                                  normCounts = normCounts);
       #jscs@plottingEstimatesVST <- list(exprEstimateVST = exprEstimateVST, 
       #                                   geneExprEstimateVST = otherExprEstimateVST, 
       #                                  relExprEstimateVST = relExprEstimateVST, 
       #                                  normCountsVST = normCountsVST);

       if(truncateBelowOne){
         convertY <- function(y){
           ifelse(y < 0, ifelse(is.infinite(y), INTERNAL.NINF.VALUE, (1-exp(y)) * INTERNAL.NINF.VALUE), y);
         }
       } else {
         convertY <- function(y){
           y
         }
       }

       if(plot.type == "rExpr"){
         count <- if(use.vst) { 
                    vst( jscs@plottingEstimates[["relExprEstimate"]][rt,, drop=FALSE], jscs);
                  } else  if(use.log) {
                    apply(log10(jscs@plottingEstimates[["relExprEstimate"]][rt,, drop=FALSE]),c(1,2),FUN=convertY);
                  } else {
                    jscs@plottingEstimates[["relExprEstimate"]][rt,, drop=FALSE];
                  }
         color.count <- rep(color[condition.names],each=nrow(count))
         y.axis.title <- "Relative Coverage";
         y.axis.title.right <- "Gene-Level Mean Normalized Counts";
         main.title <- paste0("Relative Coverage (",geneName,")");
         
         if(plot.gene.level.expression) {
           geneCount <- if(use.vst){ 
                    plot.gene.level.expression <- FALSE;
                    NULL;
                  } else if(use.log){
                    sapply(log10(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]][rownames(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]]) == geneID ,]),FUN=convertY);
                  } else {
                    jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]][rownames(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]]) == geneID ,];
           }
           color.geneCount <- color[condition.names];
         }
         
       } else if(plot.type == "normCounts"){
         count <- if(use.vst){ 
                    vst(jscs@plottingEstimates[["normCounts"]][rt,1:length(sample.names), drop=FALSE], jscs) 
                  } else if(use.log) {
                    apply(log10(jscs@plottingEstimates[["normCounts"]][rt,1:length(sample.names), drop=FALSE]),c(1,2),FUN=convertY)
                  } else {
                    jscs@plottingEstimates[["normCounts"]][rt,1:length(sample.names), drop=FALSE]
                  }
         if(plot.gene.level.expression) { 
           geneCount <- if(use.vst){ 
                    plot.gene.level.expression <- FALSE;
                    NULL;
                  } else if(use.log){
                    sapply(log10(jscs@geneCountData[rownames(jscs@geneCountData) == geneID,] / sizeFactors(jscs)),FUN=convertY)
                  } else {
                    jscs@geneCountData[rownames(jscs@geneCountData) == geneID,] / sizeFactors(jscs);
           }
           color.geneCount <- color[as.character(condition)];
         }
         color.count <- rep(color[as.character(condition)], each=nrow(count));

         y.axis.title <- "Normalized Counts";
         y.axis.title.right <- "Gene-Level Normalized Counts";
         main.title <- paste0("Normalized Counts (",geneName,")");
       } else if(plot.type == "rawCounts"){
         count <- if(use.vst){ 
                    #jscs@plottingEstimatesVST[["normCountsVST"]][rt,1:length(sample.names), drop=FALSE];
                    vst(jscs@countVectors[rt,1:length(sample.names), drop=FALSE], jscs);
                  } else if(use.log) {
                    apply(log10(jscs@countVectors[rt,1:length(sample.names), drop=FALSE]),c(1,2),FUN=function(y){max(INTERNAL.NINF.VALUE,y)})
                  } else {
                    jscs@countVectors[rt,1:length(sample.names), drop=FALSE]
                  }
         if(plot.gene.level.expression) {
           geneCount <- if(use.vst){ 
                    plot.gene.level.expression <- FALSE;
                    NULL;
                  } else if(use.log){
                    sapply(log10(jscs@geneCountData[rownames(jscs@geneCountData) == geneID,]),FUN=convertY)
                  } else {
                    jscs@geneCountData[rownames(jscs@geneCountData) == geneID,]
           }
           color.geneCount <- color[as.character(condition)];
         }
         color.count <- rep(color[as.character(condition)], each=nrow(count))
         y.axis.title <- "Raw Counts";
         y.axis.title.right <- "Raw Gene-Level Counts";
         main.title <- paste0("Raw Counts (",geneName,")");
       } else if(plot.type == "expr"){
         count <- if(use.vst){ 
                    vst(jscs@plottingEstimates[["exprEstimate"]][rt,, drop=FALSE], jscs);
                  } else if(use.log) {
                    apply(log10(jscs@plottingEstimates[["exprEstimate"]][rt,, drop=FALSE]),c(1,2),FUN=convertY);
                  } else {
                    jscs@plottingEstimates[["exprEstimate"]][rt,, drop=FALSE];
                  }
         if(plot.gene.level.expression) {
           geneCount <- if(use.vst){ 
                    plot.gene.level.expression <- FALSE;
                    NULL;
                  } else if(use.log){
                    sapply(log10(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]][rownames(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]]) == geneID ,]),FUN=convertY);
                  } else {
                    jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]][rownames(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]]) == geneID ,];
           }
           color.geneCount <- color[condition.names];
         }
         color.count <- rep(color[condition.names],each=nrow(count))
         y.axis.title <- "Mean Normalized Counts";
         y.axis.title.right <- "Gene-Level Mean Normalized Counts";
         main.title <- paste0("Mean Normalized Coverage (",geneName,")");
       } else {
         stop(paste0("FATAL ERROR: Unknown plot type! plot.type = \"",plot.type,"\""));
       }
       
       ###########################################################
       #Default titles and axes labels:
       fragment.label <- if(sequencing.type == "paired-end") "Read-Pair" else "Read";
       
       if(plot.type == "rExpr"){
         main.title <- paste0("Relative Coverage (",geneName,")");
         y.axis.title <- paste0("Relative Coverage");
         y.axis.title.right <- paste0(fragment.label,"s per Sample, Gene-Level");
         
       } else if(plot.type == "normCounts"){
         main.title <- paste0("Normalized Counts (",geneName,")");
         y.axis.title <- paste0("Normalized Counts");
         y.axis.title.right <- paste0("Normalized Counts, Gene-Level");
         
       } else if(plot.type == "rawCounts"){
         main.title <- paste0("Raw Counts (",geneName,")");
         y.axis.title <- paste0("Raw ",fragment.label," Counts");
         y.axis.title.right <- paste0("Raw ",fragment.label," Counts, Gene-Level");
         
       } else if(plot.type == "expr"){
         main.title <- paste0("Mean Normalized Coverage (",geneName,")");
         y.axis.title <- paste0(fragment.label,"s per Sample");
         y.axis.title.right <- paste0(fragment.label,"s per Sample, Gene-Level");
         
       } else {
         stop(paste0("FATAL ERROR: Unknown plot type! plot.type = \"",plot.type,"\""));
       }
       # Option: Allow users to explicitly override labels and titles:
       if(! is.null(title.main)) main.title <- title.main;
       if(! is.null(title.ylab)) y.axis.title <- title.ylab;
       if(! is.null(title.ylab.right)) y.axis.title.right <- title.ylab.right;
       
       if(! plot.gene.level.expression) {
         geneCount <- NULL;
         color.geneCount <- NULL;
       }
       count <- as.matrix(count);

    if(debug.mode & verbose) message(">    pJSRforGene(): ", "Reached step 4.");


       #SET COLORS FOR ANNOTATION:
       SIG.VERTLINE.COLOR <- "#F219ED60";
       NOSIG.VERTLINE.COLOR <- "#AAAAAA60";
       UNTESTABLE.VERTLINE.COLOR <- "#22222260";

       SIG.FEATURE.COLOR <- "#F219ED";
       NOSIG.FEATURE.COLOR <- "#111111";
       UNTESTABLE.FEATURE.COLOR <- "#DDDDDD";
       EXCLUDED.FEATURE.COLOR <- "#DDDDDD";

       SIG.FEATURE.BORDER.COLOR <- "#000000";
       NOSIG.FEATURE.BORDER.COLOR <- "#000000";
       UNTESTABLE.FEATURE.BORDER.COLOR <- "#AAAAAA";
       EXCLUDED.FEATURE.BORDER.COLOR <- "#000000";

       SIG.FEATURE.FILL.COLOR <- "#F219ED";
       NOSIG.FEATURE.FILL.COLOR <- "#CCCCCC";
       UNTESTABLE.FEATURE.FILL.COLOR <- "#EEEEEE";
       EXCLUDED.FEATURE.FILL.COLOR <- "#CCCCCC";

       intervals<-(0:nrow(count))/nrow(count)

       numexons<-nrow(count)
       each <- merged.data$padjust[rt]
       vertline.col <- ifelse(merged.data$testable[rt], ifelse(f.na(merged.data$padjust[rt] <= FDR), SIG.VERTLINE.COLOR, NOSIG.VERTLINE.COLOR), UNTESTABLE.VERTLINE.COLOR);
       annolink.col <- ifelse(merged.data$testable[rt], ifelse(f.na(merged.data$padjust[rt] <= FDR), SIG.FEATURE.COLOR,  NOSIG.FEATURE.COLOR),  UNTESTABLE.FEATURE.COLOR);
       exonlty <- rep(1,length(vertline.col));
       exonlty[as.character(merged.data$featureType[rt]) == "novel_splice_site"] <- 2 

       is.sig.feature <- f.na(each <= FDR);
       sig.feature <- which(is.sig.feature);

    if(debug.mode & verbose) message(">    pJSRforGene(): ", "Reached step 5.");

       ################## DETERMINE THE LAYOUT OF THE PLOT DEPENDING OF THE OPTIONS THE USER PROVIDES ###########
          sub <- data.frame(start=merged.data$start[rt], 
                            end=merged.data$end[rt], 
                            chr=merged.data$chr[rt], 
                            strand=merged.data$strand[rt], 
                            is.exon = (merged.data$featureType[rt] == "exonic_part"), 
                            is.testable = merged.data$testable[rt],
                            is.sig = is.sig.feature,
                            col = vertline.col,
                            featureType = merged.data$featureType[rt], 
                            featureID = rownames(merged.data)[rt],
                            stringsAsFactors=F);
          sub.allExon <- data.frame(start=flat.gff.data$start[rt.allExon], 
                                    end=flat.gff.data$end[rt.allExon], 
                                    chr=flat.gff.data$chrom[rt.allExon], 
                                    strand=flat.gff.data$strand[rt.allExon], 
                                    is.exon = (flat.gff.data$featureType[rt.allExon] == "exonic_part"),
                                    featureID = as.character(flat.gff.data$featureName[rt.allExon]), 
                                    stringsAsFactors=F);
          sub.allJunction <- data.frame(start=flat.gff.data$start[rt.allJunction], 
                                    end=flat.gff.data$end[rt.allJunction], 
                                    chr=flat.gff.data$chrom[rt.allJunction], 
                                    strand=flat.gff.data$strand[rt.allJunction], 
                                    is.exon = (flat.gff.data$featureType[rt.allJunction] == "exonic_part"), 
                                    feature.type = flat.gff.data$featureType[rt.allJunction] ,
                                    is.novel = (flat.gff.data$featureType[rt.allJunction] == "novel_splice_site"), 
                                    featureID = as.character(flat.gff.data$featureName[rt.allJunction]), 
                                    stringsAsFactors=F);

       testable.featureIDs <- sub$featureID[sub$is.testable];
       sig.featureIDs <- sub$featureID[sub$is.sig];
       untestable.featureIDs <- rownames(merged.data)[untestable.rt];

       if(debug.mode & verbose){ message("testable.featureIDs: ",paste0(testable.featureIDs,collapse=",")) }
       if(debug.mode & verbose){ message("sig.featureIDs: ",paste0(sig.featureIDs,collapse=",")) }
       if(debug.mode & verbose){ message("untestable.featureIDs: ",paste0(untestable.featureIDs,collapse=",")) }


       sub.allExon$is.testable   <- sub.allExon$featureID %in% testable.featureIDs;
       sub.allExon$is.sig        <- sub.allExon$featureID %in% sig.featureIDs;
       sub.allExon$is.untestable <- sub.allExon$featureID %in% untestable.featureIDs;
       sub.allJunction$is.testable <- sub.allJunction$featureID %in% testable.featureIDs;
       sub.allJunction$is.sig      <- sub.allJunction$featureID %in% sig.featureIDs;
       sub.allJunction$is.untestable <- sub.allJunction$featureID %in% untestable.featureIDs;

       sub.allExon$lineColor     <- ifelse(sub.allExon$is.testable,
                                           ifelse(sub.allExon$is.sig,
                                                  SIG.FEATURE.COLOR,
                                                  NOSIG.FEATURE.COLOR),
                                           ifelse(sub.allExon$is.untestable, 
                                                  UNTESTABLE.FEATURE.COLOR, 
                                                  EXCLUDED.FEATURE.COLOR));
       sub.allExon$fillColor     <- ifelse(sub.allExon$is.testable,ifelse(sub.allExon$is.sig,SIG.FEATURE.FILL.COLOR,NOSIG.FEATURE.FILL.COLOR),ifelse(sub.allExon$is.untestable, UNTESTABLE.FEATURE.FILL.COLOR, EXCLUDED.FEATURE.FILL.COLOR));
       sub.allExon$borderColor   <- ifelse(sub.allExon$is.testable,ifelse(sub.allExon$is.sig,SIG.FEATURE.BORDER.COLOR,NOSIG.FEATURE.COLOR),ifelse(sub.allExon$is.untestable, UNTESTABLE.FEATURE.BORDER.COLOR, EXCLUDED.FEATURE.BORDER.COLOR));
       sub.allJunction$lineColor <- ifelse(sub.allJunction$is.testable,ifelse(sub.allJunction$is.sig,SIG.FEATURE.COLOR,NOSIG.FEATURE.COLOR),ifelse(sub.allJunction$is.untestable, UNTESTABLE.FEATURE.COLOR, EXCLUDED.FEATURE.COLOR));
       sub.allJunction$lty <- ifelse(sub.allJunction$is.novel,2,1); 


       sig.feature.names <- sub$featureID[is.sig.feature & sub$is.exon];
       allExon.isSig <- sub.allExon$featureID %in% sig.featureIDs;
       allExon.exonCol <- ifelse(allExon.isSig, "#F219ED", "#CCCCCC");

       #TEMPORARY DEBUGGING MESSAGES:
    if(debug.mode == 2 & verbose){
         message(">    Debugging Info:");
         message("        Exons:");
         message("           ", paste0(names(sub.allExon),collapse='\t'));
         for(i in 1:length(sub.allExon$start)){
           message("           ", paste0(sub.allExon[i,],collapse='\t'));
         }
         message(">       dim(count) = ", paste0(dim(count), collapse=","));
         message(">       rt: ", paste0(rt, collapse=","));
         message(">       rango: ", paste0(rango, collapse=","));
         message(">       ncol(count): ", ncol(count));
    }

          sub.sig <- sub[sig.feature,, drop=FALSE];
           #print(sub.sig);
           #print("sub.allExon");
           #print(sub.allExon);
          rel.calc.min <- min(sub.allJunction$start, sub.allExon$start)
          rel.calc.max <- max(sub.allJunction$end,   sub.allExon$end)

          if(is.na(exon.rescale.factor) | exon.rescale.factor <= 0 | exon.rescale.factor >= 1){
            rel <- (data.frame(sub.allExon$start, sub.allExon$end))-rel.calc.min;
            rel <- rel/(rel.calc.max - rel.calc.min);
            rescale.iv <- NULL;
          } else {
            rescale.iv <- generate.interval.scale(
              data.frame(
                start = c(sub.allExon$start, sub.allJunction$start),
                end = c(sub.allExon$end, sub.allJunction$end),
                is.exon = c(sub.allExon$is.exon, sub.allJunction$is.exon)
              ),exon.rescale.factor);
            rel <- data.frame(start = rescale.coords(sub$start,rescale.iv), 
                              end   = rescale.coords(sub$end,  rescale.iv));
          }

          transcripts <- sapply(sapply(flat.gff.data$transcripts[rt.allExon],toString), function(x){strsplit(x, "+",fixed=TRUE)})
          trans <- Reduce(union, transcripts)
          if(displayTranscripts==TRUE){
             mat <- 1:4
             hei<-c(10,1, GENE.annotation.height, TX.annotation.height * length(trans));
          }else{
             mat<-1:3
             hei<-c(10,1, GENE.annotation.height)
          }
          #Add another transcript height to the bottom, to serve as a lower margin.
          hei <- c(hei, TX.annotation.height)
          mat <- c(mat, length(mat)+1)
          layout(matrix(mat), heights=hei)
          
          if(debug.mode & verbose){
            message(">   length(trans) = ",length(trans));
            message(">   FINAL LAYOUT:");
            message(">      heights = [", paste0(hei, collapse=","),"]");
          }
          
    if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached step 6.");

          ylimn <- c(min(min(count,na.rm=TRUE),0), max(count, na.rm=TRUE));
          if((! use.vst) & use.log ) ylimn[1] <- INTERNAL.NINF.VALUE;

          p.values.labels <- ifelse(each<=FDR, format(each,digits=3), "");

          if(any(sub$is.exon)){
            italicize.label <- ! sub$is.exon;
          } else {
            italicize.label <- NULL;
          }

          intervals <- drawPlot(matr=count, ylimn,jscs, 
                   intervals, rango, textAxis=y.axis.title, geneLevelAxisTitle = y.axis.title.right,
                   rt=rt, color.count=color.count, 
                   colorlines=vertline.col, 
                   countbinIDs = merged.data$countbinID[rt],
                   use.vst=use.vst, use.log = use.log,plot.type=plot.type,
                   main.title=main.title,draw.legend=draw.legend,
                   color.key=color,condition.names=condition.names,
                   p.values=p.values.labels,draw.p.values=TRUE, 
                   plot.lwd=plot.lwd, axes.lwd = axes.lwd, 
                   anno.lwd = anno.lwd, par.cex = par.cex, 
                   anno.cex.text = anno.cex.text,
                   anno.cex.axis = anno.cex.axis, 
                   anno.cex.main = anno.cex.main,
                   fit.countbin.names = fit.countbin.names,
                   debug.mode = debug.mode, plot.gene.level.expression = plot.gene.level.expression, geneCount = geneCount, color.geneCount = color.geneCount,
                   yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, italicize.label = italicize.label, condition.legend.text = condition.legend.text,
                   rel = rel, annolink.col = annolink.col, exonlty = exonlty, graph.margins = graph.margins,
                   ...);
    if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached end of step 6.");


    if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached step 7.");

          #########PLOT THE GENE MODEL:
          if(USE.MARGIN.MEX){
            par(mar=c(0, graph.margins[2], 0, graph.margins[4]), cex = par.cex, mex = anno.cex.text); 
          } else {
            par(mar=c(0, graph.margins[2], 0, graph.margins[4]), cex = par.cex); 
          }

          plot.new();
          plot.window(xlim=c(0, 1.04), ylim=c(0,1), xaxs = "i");
          
          # lines linking exons / splices to their column:
          segments(
                        apply((rbind(rel[rango,2], rel[rango, 1])), 2, median), 
                        0, #par("usr")[3], (old version: lines are connected.)
                        apply(rbind(intervals[rango], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2)), 2, median), 
                        1, col=annolink.col, lty = exonlty, lwd = plot.lwd, cex = anno.cex.text,cex.axis=anno.cex.main, cex.main=anno.cex.main, xpd=NA, ...) #col=colorlinesB,...)
          

          #axes.lwd = axes.lwd, anno.lwd = anno.lwd,
          par(mar=c(2, graph.margins[2], 0, graph.margins[4]), cex = par.cex);

          #Get start/end sites:
          startSites <- c();
          endSites <- c();
          for(i in 1:length(trans)){
             logicexons <- sapply(transcripts, function(x){any(x == trans[i])})
             startSites <- c(startSites, sub.allExon$start[logicexons][1]);
             endSites <- c(endSites, sub.allExon$end[logicexons][sum(logicexons)]);
          }
          startSites <- unique(startSites);
          endSites <- unique(endSites);

          drawGene(rel.calc.min, 
                   rel.calc.max, 
                   tr=sub, 
                   tr.allExon=sub.allExon, 
                   tr.allJunction = sub.allJunction,
                   #tr.untestable=sub.untestable,
                   rango, 
                   rescale.iv = rescale.iv, 
                   exoncol=annolink.col,
                   allExon.exonCol=allExon.exonCol, 
                   names, 
                   trName="Gene model", 
                   anno.cex.text=anno.cex.text, par.cex = par.cex, 
                   exonlty = exonlty, 
                   plot.lwd = plot.lwd, anno.lwd=anno.lwd, 
                   show.strand.arrows = show.strand.arrows, 
                   geneStrand = geneStrand,
                   cex.axis=anno.cex.axis, cex.main=anno.cex.main, 
                   arrows.length = arrows.length, 
                   draw.untestable.annotation = draw.untestable.annotation,
                   draw.start.end.sites = draw.start.end.sites, startSites = startSites, endSites = endSites,
                   cex.arrows = cex.arrows, chrom.label = chrom.label, label.chromosome = label.chromosome,
                   ...)
     # graph.margins = graph.margins,
           #Maybe make these options external at some point?
           include.endpoints.on.coordinates <- FALSE;
           num.coord.miniticks.per.tick <- 10;

       if(drawCoordinates){
             if(! is.null(rescale.iv)){
               pretty.x <- pretty(c(rel.calc.min,rel.calc.max),n=5);
               pretty.interval <- pretty.x[2] - pretty.x[1];
               pretty.x <- pretty.x[pretty.x > rel.calc.min & pretty.x < rel.calc.max];
               rescaled.pretty.x <- rescale.coords(pretty.x,rescale.iv);

               if(num.coord.miniticks.per.tick > 0){
                 rel.coord.miniticks <- pretty.interval * (1:(num.coord.miniticks.per.tick-1)) / num.coord.miniticks.per.tick;
                 unscaled.coord.miniticks <- unlist(lapply(pretty.x, function(a){ a + rel.coord.miniticks }));
                 unscaled.coord.miniticks <- c(pretty.x[1] - rel.coord.miniticks, unscaled.coord.miniticks);
                 rescaled.coord.miniticks <- rescale.coords(unscaled.coord.miniticks, rescale.iv) * (rel.calc.max-rel.calc.min) + rel.calc.min;
               } else {
                 rescaled.coord.miniticks <- FALSE;
               }

               if(include.endpoints.on.coordinates){
                 if(min(rescaled.pretty.x) > 0.05){
                   pretty.x <- c(rel.calc.min, pretty.x);
                   rescaled.pretty.x <- c(0,rescaled.pretty.x);
                 }
                 if(max(rescaled.pretty.x) < 0.95){
                   pretty.x <- c(pretty.x, rel.calc.max);
                   rescaled.pretty.x <- c(rescaled.pretty.x,1);
                 }
               }
               rescaled.pretty.x <- rescaled.pretty.x * (rel.calc.max-rel.calc.min) + rel.calc.min;
             } else {
               pretty.x <- pretty(c(rel.calc.min,rel.calc.max),n=5);
               coord.miniticks <- FALSE;
             }
             
             usr <- par("usr");
             cxy <- par("cxy");
             #anno.cex.coordAxis <- anno.cex.axis;
             pretty.x <- sprintf("%0.f",pretty.x);
             smallest.width.coordAxis <- min(abs(rescaled.pretty.x[-1] - rescaled.pretty.x[-length(rescaled.pretty.x)])) * 0.8;
             anno.cex.coordAxis <- shrink.character.vector(pretty.x, curr.cex = anno.cex.axis, max.width = smallest.width.coordAxis);
             
             segments(x0 = rescaled.pretty.x,            y0 = usr[3], x1 = rescaled.pretty.x,            y1 = usr[3] - (cxy[2] / 2), xpd=NA, lwd = anno.lwd, ...);
             lines(c(rel.calc.min,rel.calc.max), c(par("usr")[3], par("usr")[3]) , lwd = axes.lwd, xpd = NA, ...);
             segments(x0 = rescaled.coord.miniticks,     y0 = usr[3], x1 = rescaled.coord.miniticks,     y1 = usr[3] - (cxy[2] / 4), xpd=NA, lwd = anno.lwd, ...);
             text(rescaled.pretty.x, usr[3] - (cxy[2] * (3/4)), pretty.x, cex.axis = anno.cex.coordAxis, xpd = NA, adj = c(0.5,1), ...);
             
             #JS.axis(1, at = rescaled.pretty.x, labels = pretty.x, cex.axis=anno.cex.axis, tcl = -0.5, lwd = anno.lwd, xpd = NA, ...);
             #lines(c(rel.calc.min,rel.calc.max), c(par("usr")[3], par("usr")[3]) , lwd = axes.lwd, xpd = NA, ...);
             #JS.axis(1, at = rescaled.coord.miniticks, labels = FALSE, cex.axis=anno.cex.axis, tcl = -0.25, lwd = anno.lwd, xpd = NA, ...);
             #axis(1, at = rescaled.pretty.x,labels=pretty.x,cex.axis=anno.cex.axis, tcl = -0.5, lwd = axes.lwd, ...);
             #axis(1, at = c(rel.calc.min,rel.calc.max),labels=FALSE, cex.axis = anno.cex.axis, tcl = 0, lwd = axes.lwd,...);
             #axis(1, at = rescaled.coord.miniticks,labels=FALSE,cex.axis=anno.cex.axis, tcl = -0.25, lwd = axes.lwd, ...);
       }
       if(label.chromosome){
                chrom.label.width.max <- abs(par("usr")[1] - device.limits()[1]) * 0.9;
                chrom.label.cex <- shrink.character.vector(chrom.label, curr.cex = anno.cex.text, max.width = chrom.label.width.max);
                text(par("usr")[1],par("usr")[3],chrom.label, cex = chrom.label.cex, adj = c(1.1,0.5),xpd=NA, font = 2, ...);
       }
       if(displayTranscripts){
         if(USE.MARGIN.MEX){
           par(cex = par.cex, mar=c(0, graph.margins[2], 0, graph.margins[4]), mex = anno.cex.text);
         } else {
           par(cex = par.cex, mar=c(0, graph.margins[2], 0, graph.margins[4]));
         }
         plot.new();
         plot.window(xlim=c(rel.calc.min, rel.calc.max + (rel.calc.max - rel.calc.min) * 0.04), ylim=c(0,length(trans)), xaxs = "i");
         
         for(i in 1:length(trans)){
            ymin <- length(trans) - i;
            if(include.TX.names){
              trName = trans[i];
            } else {
              trName = NULL;
            } 

            logicexons <- sapply(transcripts, function(x){any(x == trans[i])})
            tr <- sub.allExon[logicexons,]  #   data.frame(start = sub.allExon$start[logicexons==1], end = sub.allExon$end[logicexons==1], featureType = sub.allExon$featureType[logicexons==1], stringsAsFactors = F);
            curr.exoncol <- ifelse(allExon.isSig[logicexons],"#F219ED", "#CCCCCC");
            drawTranscript(rel.calc.min, 
                           rel.calc.max, 
                           ymin = ymin,
                           tr=tr, 
                           tr.allJunction = sub.allJunction,
                           rango=1:nrow(tr), 
                           rescale.iv=rescale.iv, 
                           #exoncol=curr.exoncol, 
                           names=c(), 
                           trName=trName, 
                           trStrand = txStrandMap[[trName]],
                           draw.strand = geneStrand == ".",
                           par.cex = par.cex, 
                           anno.cex.text = anno.cex.text,
                           sub.sig = sub.sig,
                           anno.lwd=anno.lwd, 
                           cex.axis=anno.cex.axis, 
                           cex.main=anno.cex.main, 
                           cex.arrows = cex.arrows,
                           ...)
         }
      }
      par(mar = c(0,0,0,0));
      plot.new();
      par(mar = c(0,0,0,0));
      if(debug.mode & verbose) message("> pJSRfG(): "," Done.");
    }
  }, error = function(e){
    message("Error caught while attempting plotJunctionSeqResultsForGene");
    message("---------------------");
    message("     Error text:");
    message("     ",e);
    message("");
    message("---------------------");
    message("     Input parameters:")
    message("     geneID = ",geneID);
    message("     colorRed.FDR.threshold = ", colorRed.FDR.threshold);
    message("     plot.type = ",plot.type);
    message("     displayTranscripts = ",displayTranscripts);
    message("     color = ",color);
    message("     use.vst = ",use.vst);
    message("     use.log = ",use.log);
    message("     truncateBelowOne = ",truncateBelowOne);
    message("     exon.rescale.factor = ",exon.rescale.factor);
    message("     label.p.vals = ",label.p.vals);
    message("     plot.lwd = ",plot.lwd);
    message("     axes.lwd = ",axes.lwd);
    message("     anno.lwd = ",anno.lwd);
    message("     fit.countbin.names = ",fit.countbin.names);
    message("     plot.gene.level.expression = ",plot.gene.level.expression);
    message("     plot.exon.results = ",plot.exon.results);
    message("     plot.junction.results = ",plot.junction.results);
    message("     plot.novel.junction.results = ",plot.novel.junction.results);
    message("     plot.untestable.results = ",plot.untestable.results);
    message("     draw.untestable.annotation = ",draw.untestable.annotation);
    message("     show.strand.arrows = ",show.strand.arrows);
    message("     sort.features = ",sort.features);
    message("     drawCoordinates = ",drawCoordinates);
    message("     yAxisLabels.inExponentialForm = ",yAxisLabels.inExponentialForm);
    message("---------------------");
    
    warning("Error caught while attempting plotJunctionSeqResultsForGene");
    warning("Error text:");
    warning(e);
  })
}

