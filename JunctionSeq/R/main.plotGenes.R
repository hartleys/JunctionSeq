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
######### Generating Plots From Results:
##########################################################################

buildAllPlots <- function(jscs,
                          outfile.prefix = "./",
                          #flat.gff.data = NULL, flat.gff.file = NULL, 
                          gene.list = NULL, FDR.threshold = 0.01, max.gene.ct,
                          method.selectionCriterion = c("feature-pAdjust", "genewise-pAdjust"),
                          use.plotting.device = c("png","CairoPNG","svg","tiff","cairo_ps","custom"),
                          sequencing.type = c("paired-end","single-end"),
                          use.vst=FALSE,use.log = TRUE, #truncateBelowOne = TRUE, 
                          exon.rescale.factor = 0.3,
                          subdirectories.by.type = TRUE,
                          ma.plot=TRUE, variance.plot=TRUE,
                          with.TX=TRUE,without.TX=TRUE,
                          expr.plot=TRUE,normCounts.plot=TRUE,
                          rExpr.plot=TRUE,rawCounts.plot=FALSE,
                          colorRed.FDR.threshold = FDR.threshold, 
                          colorList=list(),
                          plot.gene.level.expression = TRUE,
                          plot.exon.results = NULL, plot.junction.results = NULL, plot.novel.junction.results = NULL,
                          plot.untestable.results = FALSE,
                          plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, gene.lwd = plot.lwd / 2,
                          par.cex = 1, anno.cex.text = 1, anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                          drawCoordinates = TRUE, yAxisLabels.inExponentialForm = FALSE,
                          show.strand.arrows = 1, #arrows.length = 0.125,
                          graph.margins = c(2,3,3,3),
                          base.plot.height = 12, base.plot.width = 12, 
                          base.plot.units = "in", 
                          GENE.annotation.relative.height = 0.15, TX.annotation.relative.height = 0.05, CONNECTIONS.relative.height = 0.1,
                          SPLICE.annotation.relative.height = 0.1,
                          TX.margins = c(0,0.5),
                          autoscale.height.to.fit.TX.annotation = TRUE,
                          autoscale.width.to.fit.bins = 35,
                          plotting.device.params = list(), 
                          number.plots = FALSE, name.files.with.geneID = TRUE,
                          condition.legend.text = NULL, include.TX.names = TRUE, draw.start.end.sites = TRUE,
                          openPlottingDeviceFunc = NULL, closePlottingDeviceFunc = NULL,
                          writeHTMLresults = TRUE,
                          html.cssFile = NULL, html.cssLink = NULL, html.imgFileExtension = NULL,
                          html.plot.height = 90, html.plot.height.units = "vh",
                          html.compare.results.list = NULL,
                          minimalImageFilenames = writeHTMLresults,
                          verbose=TRUE, debug.mode = FALSE,
                          INTERNAL.VARS = list(),
                          ...){
  condition <- jscs@phenoData$condition
  flat.gff.data <- jscs@flatGffData
  
  method.selectionCriterion <- match.arg(method.selectionCriterion)
  use.plotting.device <- match.arg(use.plotting.device)
  sequencing.type <- match.arg(sequencing.type)
  
  if(is.null(openPlottingDeviceFunc) || is.null(closePlottingDeviceFunc)){
    if(use.plotting.device == "custom"){
        stop("Fatal error: custom plotting device is selected, but openPlottingDeviceFunc or closePlottingDeviceFunc is not set.")
    }
    devFunctions <- getPlottingDeviceFunc(use.plotting.device = use.plotting.device,
                                            base.plot.height = base.plot.height,
                                            base.plot.width = base.plot.width,
                                            base.plot.units = base.plot.units,
                                            plotting.device.params = plotting.device.params)
    openPlottingDeviceFunc <- devFunctions[[1]]
    closePlottingDeviceFunc <- devFunctions[[2]]
  }
  if(is.null(html.imgFileExtension)){
    html.imgFileExtension <- getPlottingDeviceFileExtension(use.plotting.device)
  }
  
  gtf.format <- TRUE
  
  if(is.null(gene.list)){
    if(method.selectionCriterion == "genewise-pAdjust"){
      genewise.sig <- JS.perGeneQValue(pvals = fData(jscs)$pvalue, wTest = fData(jscs)$testable, fData(jscs)$geneID)
      gene.list <- as.character(names(genewise.sig)[genewise.sig < FDR.threshold])
      
      if(is.null(fData(jscs)$geneWisePadj)){
        fData(jscs)$geneWisePadj <- sapply(as.character(fData(jscs)$geneID), function(g){ if(any(g == names(genewise.sig))){ genewise.sig[[g]]; } else { NA;} })
      }
      
      gene.list <- gene.list[order(genewise.sig[gene.list])]
      

      
      if(verbose) message("> buildAllPlots: Found ", length(gene.list), " significant genes, with gene-wise adjusted-p-value threshold ", FDR.threshold)
    } else if(method.selectionCriterion == "feature-pAdjust"){
      sig.features <- which(fData(jscs)$padjust < FDR.threshold)
      gene.list <- unique(as.character(fData(jscs)$geneID[sig.features]))

      gene.list.pval <- sapply(gene.list, function(g){
        min(fData(jscs)$padjust[fData(jscs)$geneID == g], na.rm=TRUE)
      })
      gene.list <- gene.list[order(gene.list.pval)]
      
      if(verbose) message("> buildAllPlots: Found ", length(gene.list), " genes with at least one significant exon, at adjusted-p-value threshold ", FDR.threshold)
    }

  } else {
    if(verbose) message("> buildAllPlots: Found ", length(gene.list), " genes to plot.")
    
    if(! all(gene.list %in% unique(fData(jscs)$geneID) )){
      notFound <- ! (gene.list %in% unique(fData(jscs)$geneID))
      stop(paste0("FATAL ERROR: Not all GeneID's found in dataset! Example: geneID \"", gene.list[notFound][1],"\" is missing!" ))
    }
  }
  
  #If there are too many genes, cut it down!
  if(! missing(max.gene.ct)){
    if(! is.numeric(max.gene.ct)){
      stop("Parameter max.gene.ct must be NUMERIC.")
    }
    if(length(gene.list) > max.gene.ct){
      gene.list <- gene.list[1:max.gene.ct]
      if(verbose) message("> buildAllPlots: Too many genes found. Only plotting the first ",max.gene.ct, " genes.")
    }
  }

    if(verbose) message("> buildAllPlots: Starting plotting...")

    FDR <- colorRed.FDR.threshold

    total.gene.ct <- length(gene.list)
    number.width <- max(1, floor(log10(total.gene.ct)) + 1)
    plot.nums <- paste0(formatC(1:total.gene.ct, width = number.width, format = 'd', flag = '0'));
    if( number.plots || minimalImageFilenames){
      geneNum.strings <- paste0(plot.nums,"-");
    } else {
      geneNum.strings <- rep("",total.gene.ct)
    }
    
    if(subdirectories.by.type){
      if(! file.exists(outfile.prefix)){
        dir.create(outfile.prefix)
      }
      
      if(expr.plot && with.TX          && (! file.exists(paste0(outfile.prefix,"/exprTX"))))       dir.create(paste0(outfile.prefix,"/exprTX"))
      if(expr.plot && without.TX       && (! file.exists(paste0(outfile.prefix,"/expr"))))         dir.create(paste0(outfile.prefix,"/expr"))
      if(normCounts.plot && with.TX    && (! file.exists(paste0(outfile.prefix,"/normCountsTX")))) dir.create(paste0(outfile.prefix,"/normCountsTX"))
      if(normCounts.plot && without.TX && (! file.exists(paste0(outfile.prefix,"/normCounts"))))   dir.create(paste0(outfile.prefix,"/normCounts"))
      if(rExpr.plot && with.TX         && (! file.exists(paste0(outfile.prefix,"/rExprTX"))))      dir.create(paste0(outfile.prefix,"/rExprTX"))
      if(rExpr.plot && without.TX      && (! file.exists(paste0(outfile.prefix,"/rExpr"))))        dir.create(paste0(outfile.prefix,"/rExpr"))
      if(rawCounts.plot && with.TX     && (! file.exists(paste0(outfile.prefix,"/rawCountsTX"))))  dir.create(paste0(outfile.prefix,"/rawCountsTX"))
      if(rawCounts.plot && without.TX  && (! file.exists(paste0(outfile.prefix,"/rawCounts"))))    dir.create(paste0(outfile.prefix,"/rawCounts"))
    }
    
    if(variance.plot){
      if(verbose) message("> buildAllPlots: Generating Dispersion Plot")
      
      openPlottingDeviceFunc(paste(outfile.prefix,"dispersion-plot","",sep=""),heightMult=0.6,widthMult=0.6)
      plotDispEsts( jscs, par.cex = par.cex, points.cex = anno.cex.text, text.cex = anno.cex.text, verbose = verbose, debug.mode = debug.mode, ... )
      dev.off()
    }

    if(ma.plot){
      fc.cols <- which(strStartsWith(names(fData(jscs)), "log2FC("))

      for(fc.name in names(fData(jscs))[fc.cols]){
        if(verbose) message("> buildAllPlots: Generating MA-Plot (",fc.name,")")

        fc.title <- sub("/","vs",fc.name, fixed = TRUE)
        openPlottingDeviceFunc(paste(outfile.prefix,"ma-plot-",fc.title,sep=""),heightMult=0.6,widthMult=0.6)
        plotMA( jscs, FDR.threshold=colorRed.FDR.threshold, fc.name = fc.name, par.cex = par.cex, text.cex = anno.cex.text, points.cex = anno.cex.text, verbose = verbose, debug.mode = debug.mode, ... )
        dev.off()
      }
    }
    

    if(writeHTMLresults){
      if(verbose) message("> buildAllPlots: Writing HTML results index.")
      if((! file.exists(paste0(outfile.prefix,"/htmlFiles")))) dir.create(paste0(outfile.prefix,"/htmlFiles"))
      
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
                            GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, CONNECTIONS.relative.height = CONNECTIONS.relative.height,TX.margins=TX.margins,
                            autoscale.height.to.fit.TX.annotation = autoscale.height.to.fit.TX.annotation,
                            autoscale.width.to.fit.bins = autoscale.width.to.fit.bins,
                            number.plots = geneNum.strings,
                            css.file = html.cssFile, css.link = html.cssLink,
                            minimalImageFilenames = minimalImageFilenames,
                            compare.analysis.list = html.compare.results.list, verbose = verbose, debug.mode = debug.mode,
                            INTERNAL.VARS = INTERNAL.VARS)
     if(verbose) message("> buildAllPlots: Finished writing HTML results index.")
    }
    
      if(is.null(plot.exon.results)){
        plot.exon.results <- any( fData(jscs)$featureType == "exonic_part" )
      }
      if(is.null(plot.junction.results)){
        plot.junction.results <- any( fData(jscs)$featureType == "splice_site" | fData(jscs)$featureType == "novel_splice_site" )
      }
      if(is.null(plot.novel.junction.results)){
        if(plot.junction.results){
          plot.novel.junction.results <- any( fData(jscs)$featureType == "novel_splice_site" )
        } else {
          plot.novel.junction.results <- FALSE
        }
      }
    
    geneNum <- 1
    for(geneID in gene.list){
      if(verbose) message(paste0("> buildAllPlots: starting geneID: ",geneID," (",geneNum," of ",length(gene.list),")"))
      geneNum.string <- geneNum.strings[geneNum]
      
      if(subdirectories.by.type){
        outfile.prefixes <- paste0(outfile.prefix, c("/exprTX/","/expr/",
                                                    "/normCountsTX/","/normCounts/",
                                                    "/rExprTX/","/rExpr/",
                                                    "/rawCountsTX/","/rawCounts/"),geneNum.string)
      } else {
        outfile.prefixes <- paste0(outfile.prefix, geneNum.string)
      }
      
      buildAllPlotsForGene(geneID = geneID, jscs = jscs,
                            outfile.prefix = outfile.prefixes,
                            #flat.gff.data = flat.gff.data, 
                            use.plotting.device = use.plotting.device,
                            use.vst=use.vst, use.log = use.log, #truncateBelowOne = truncateBelowOne,
                            exon.rescale.factor = exon.rescale.factor, plot.gene.level.expression = plot.gene.level.expression,
                            with.TX=with.TX,without.TX=without.TX,
                            expr.plot=expr.plot,normCounts.plot=normCounts.plot,
                            rExpr.plot=rExpr.plot,rawCounts.plot=rawCounts.plot,
                            colorRed.FDR.threshold = colorRed.FDR.threshold, 
                            colorList= colorList,
                            plot.exon.results = plot.exon.results, 
                            plot.junction.results = plot.junction.results,
                            plot.novel.junction.results = plot.novel.junction.results,
                            plot.untestable.results = plot.untestable.results,
                            plot.lwd = plot.lwd, axes.lwd = axes.lwd, anno.lwd = anno.lwd, gene.lwd= gene.lwd,
                            drawCoordinates = drawCoordinates,
                            show.strand.arrows = show.strand.arrows,
                            graph.margins = graph.margins,
                            par.cex = par.cex,
                            anno.cex.text = anno.cex.text,
                            anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,
                            base.plot.height = base.plot.height, base.plot.width = base.plot.width,
                            base.plot.units = base.plot.units,
                            GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, CONNECTIONS.relative.height = CONNECTIONS.relative.height,TX.margins=TX.margins,
                            SPLICE.annotation.relative.height=SPLICE.annotation.relative.height,
                            autoscale.height.to.fit.TX.annotation = autoscale.height.to.fit.TX.annotation,
                            autoscale.width.to.fit.bins = autoscale.width.to.fit.bins,
                            name.files.with.geneID = name.files.with.geneID,
                            plotting.device.params = plotting.device.params,
                            yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm,
                            condition.legend.text = condition.legend.text, include.TX.names = include.TX.names, draw.start.end.sites=draw.start.end.sites,
                            verbose=verbose, debug.mode = debug.mode,
                            sequencing.type=sequencing.type,
                            minimalImageFilenames = minimalImageFilenames,
                            INTERNAL.VARS = INTERNAL.VARS,
                            ...)

      geneNum <- geneNum + 1
    }
    if(verbose) message("> buildAllPlots: Plotting complete.")
  if(verbose) message("> buildAllPlots: Plotting and data writing complete.")
  
  
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
                          use.vst=FALSE, use.log = TRUE, #truncateBelowOne = TRUE,
                          exon.rescale.factor = 0.3,   
                          with.TX=TRUE,without.TX=TRUE,
                          expr.plot=TRUE, normCounts.plot=TRUE,
                          rExpr.plot=TRUE, rawCounts.plot=FALSE,
                          colorRed.FDR.threshold = 0.01, 
                          colorList=list(),
                          plot.gene.level.expression = TRUE,
                          plot.exon.results = NULL, plot.junction.results = NULL, plot.novel.junction.results = NULL,
                          plot.untestable.results = FALSE,
                          plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, gene.lwd = plot.lwd / 2,
                          par.cex = 1, 
                          name.files.with.geneID = TRUE,
                          anno.cex.text = 1,
                          anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                          drawCoordinates = TRUE, yAxisLabels.inExponentialForm = FALSE,
                          show.strand.arrows = 1, #arrows.length = 0.125,
                          graph.margins = c(2,3,3,3),
                          base.plot.height = 12, base.plot.width = 12, 
                          base.plot.units = "in", 
                          GENE.annotation.relative.height = 0.15, TX.annotation.relative.height = 0.05, CONNECTIONS.relative.height = 0.1,
                          SPLICE.annotation.relative.height = 0.1,
                          TX.margins = c(0,0.5),
                          autoscale.height.to.fit.TX.annotation = TRUE,
                          autoscale.width.to.fit.bins = 35,
                          plotting.device.params = list(),
                          condition.legend.text = NULL, include.TX.names = TRUE, draw.start.end.sites = TRUE,
                          draw.nested.SJ = TRUE,
                          openPlottingDeviceFunc = NULL, closePlottingDeviceFunc = NULL,
                          minimalImageFilenames = FALSE,
                          verbose=TRUE, debug.mode = FALSE, 
                          INTERNAL.VARS = list(),
                          ...){
    message(paste("starting buildAllPlotsForGene() for geneID:",geneID))
    condition <- jscs@phenoData$condition
    flat.gff.data <- jscs@flatGffData
    
    sequencing.type <- match.arg(sequencing.type)
    use.plotting.device <- match.arg(use.plotting.device)
    
    if(length(outfile.prefix) == 1){
      outfile.prefix <- rep(outfile.prefix, 8)
    }
    
    if(! any( geneID == fData(jscs)$geneID )){
      stop(paste0("FATAL ERROR: GeneID \"",geneID,"\" not found in dataset!" ))
    }
    
  if(is.null(openPlottingDeviceFunc) || is.null(closePlottingDeviceFunc)){
    if(use.plotting.device == "custom"){
        stop("Fatal error: custom plotting device is selected, but openPlottingDeviceFunc or closePlottingDeviceFunc is not set.")
    }
    devFunctions <- getPlottingDeviceFunc(use.plotting.device = use.plotting.device,
                                            base.plot.height = base.plot.height,
                                            base.plot.width = base.plot.width,
                                            base.plot.units = base.plot.units,
                                            verbose = verbose, debug.mode = debug.mode,
                                            plotting.device.params = plotting.device.params)
    openPlottingDeviceFunc <- devFunctions[[1]]
    closePlottingDeviceFunc <- devFunctions[[2]]
  }
  
    #INTERNAL.VARS <- list()
    
    if(draw.nested.SJ){
      rt.allJunction <- which(jscs@flatGffData$gene_id == geneID & (jscs@flatGffData$featureType != "exonic_part"))
      
      if(length(rt.allJunction) > 0){
        if(debug.mode) message("rt.allJunction [",paste0(rt.allJunction,collapse=","),"]")
        tr.splice <- data.frame(featureID = as.character(jscs@flatGffData$featureName[rt.allJunction]), start = jscs@flatGffData$start[rt.allJunction], end = jscs@flatGffData$end[rt.allJunction])
        tr.splice$span <- tr.splice$end - tr.splice$start
        if(debug.mode) message("tr.splice:")
        if(debug.mode) print(tr.splice)
        INTERNAL.VARS <- c(INTERNAL.VARS, get.nested.heights(tr.splice, 0, 1, verbose = verbose, debug.mode = debug.mode))
      } else {
        INTERNAL.VARS <- c(INTERNAL.VARS, maxDepth = 0)
      }
    } else {
      INTERNAL.VARS <- c(INTERNAL.VARS, maxDepth = 0)
    }
    
    FDR <- colorRed.FDR.threshold
    
    transcripts <- sapply(sapply(flat.gff.data$transcripts[which(flat.gff.data$gene_id==geneID & flat.gff.data$featureType == "exonic_part")],toString), function(x){strsplit(x, "+",fixed=TRUE)})
    trans <- Reduce(union, transcripts)
    tx.ct <- length(trans)
    
    withTxPlot.height.multiplier <- getAutofitTxRelativeHeight(tx.ct, autoscale.height.to.fit.TX.annotation = autoscale.height.to.fit.TX.annotation, GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, CONNECTIONS.relative.height = CONNECTIONS.relative.height, TX.margins = TX.margins,SPLICE.annotation.relative.height=SPLICE.annotation.relative.height)

    
    if(is.na(autoscale.width.to.fit.bins) || autoscale.width.to.fit.bins == 0){
      width.multiplier <- 1
    } else {
      num.cols <- getColCt(geneID = geneID, merged.data = fData(jscs), 
                                plot.exon.results = plot.exon.results, plot.junction.results = plot.junction.results, plot.novel.junction.results=plot.novel.junction.results, 
                                plot.untestable.results=plot.untestable.results)
      if(num.cols > autoscale.width.to.fit.bins){
        width.multiplier <- num.cols / autoscale.width.to.fit.bins
      } else {
        width.multiplier <- 1
      }
    }
    
    if(verbose && debug.mode){
      message("Found ",tx.ct," TX")
      message("Final Dimensions,   No TX:",base.plot.width * width.multiplier,"x",base.plot.height)
      message("Final Dimensions, With TX:",base.plot.width * width.multiplier,"x",base.plot.height * withTxPlot.height.multiplier)
    }
    
    geneName <- if(minimalImageFilenames){
      ""
    } else if(name.files.with.geneID){
      paste0(geneID,"-");
    } else {
      paste0(jscs@flatGffGeneData$gene_name[jscs@flatGffGeneData$geneID == geneID],"-");
    }
    
      if(is.null(plot.exon.results)){
        plot.exon.results <- any( fData(jscs)$featureType == "exonic_part" )
      }
      if(is.null(plot.junction.results)){
        plot.junction.results <- any( fData(jscs)$featureType == "splice_site" | fData(jscs)$featureType == "novel_splice_site" )
      }
      if(is.null(plot.novel.junction.results)){
        if(plot.junction.results){
          plot.novel.junction.results <- any( fData(jscs)$featureType == "novel_splice_site" )
        } else {
          plot.novel.junction.results <- FALSE
        }
      }
    
    
    if(expr.plot){
      plot.type <- "expr"
      plot.type.title <- IMAGE.NAME.MAP[[plot.type]];
      if(with.TX){
        outfile <- paste(outfile.prefix[1],geneName,"",plot.type.title,IMAGE.NAME.TX,"",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier)
        plotJunctionSeqResultsForGene(geneID, jscs, colorRed.FDR.threshold = FDR,colorList = colorList, plot.type = plot.type, use.vst=use.vst,use.log=use.log,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd , par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results , plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,gene.lwd= gene.lwd,CONNECTIONS.relative.height = CONNECTIONS.relative.height,TX.margins=TX.margins,SPLICE.annotation.relative.height=SPLICE.annotation.relative.height,draw.nested.SJ=draw.nested.SJ,INTERNAL.VARS = INTERNAL.VARS,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[2],geneName,"",plot.type.title,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier)
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,colorList = colorList, plot.type = plot.type, use.vst=use.vst,use.log=use.log ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,gene.lwd= gene.lwd,CONNECTIONS.relative.height = CONNECTIONS.relative.height,TX.margins=TX.margins,SPLICE.annotation.relative.height=SPLICE.annotation.relative.height,draw.nested.SJ=draw.nested.SJ,INTERNAL.VARS = INTERNAL.VARS,...)
        closePlottingDeviceFunc();  
      }
    }
    if(normCounts.plot){
      plot.type <- "normCounts"
      plot.type.title <- IMAGE.NAME.MAP[[plot.type]];
      if(with.TX){
        outfile <- paste(outfile.prefix[3],geneName,"",plot.type.title,IMAGE.NAME.TX,"",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier)
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,colorList = colorList, plot.type = plot.type, use.vst=use.vst,use.log=use.log ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,gene.lwd= gene.lwd,CONNECTIONS.relative.height = CONNECTIONS.relative.height,TX.margins=TX.margins,SPLICE.annotation.relative.height=SPLICE.annotation.relative.height,draw.nested.SJ=draw.nested.SJ,INTERNAL.VARS = INTERNAL.VARS,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[4],geneName,"",plot.type.title,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier)
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,colorList = colorList, plot.type = plot.type, use.vst=use.vst,use.log=use.log ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,gene.lwd= gene.lwd,CONNECTIONS.relative.height = CONNECTIONS.relative.height,TX.margins=TX.margins,SPLICE.annotation.relative.height=SPLICE.annotation.relative.height,draw.nested.SJ=draw.nested.SJ,INTERNAL.VARS = INTERNAL.VARS,...)
        closePlottingDeviceFunc();  
      }
    }
    if(rExpr.plot){
      plot.type <- "rExpr"
      plot.type.title <- IMAGE.NAME.MAP[[plot.type]];
      if(with.TX){
        outfile <- paste(outfile.prefix[5],geneName,"",plot.type.title,IMAGE.NAME.TX,"",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier)
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,colorList = colorList, plot.type = plot.type, use.vst=use.vst,use.log=use.log ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,gene.lwd= gene.lwd,CONNECTIONS.relative.height = CONNECTIONS.relative.height,TX.margins=TX.margins,SPLICE.annotation.relative.height=SPLICE.annotation.relative.height,draw.nested.SJ=draw.nested.SJ,INTERNAL.VARS = INTERNAL.VARS,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[6],geneName,"",plot.type.title,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier)
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,colorList = colorList, plot.type = plot.type, use.vst=use.vst,use.log=use.log ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,gene.lwd= gene.lwd,CONNECTIONS.relative.height = CONNECTIONS.relative.height,TX.margins=TX.margins,SPLICE.annotation.relative.height=SPLICE.annotation.relative.height,draw.nested.SJ=draw.nested.SJ,INTERNAL.VARS = INTERNAL.VARS,...)
        closePlottingDeviceFunc();  
      }
    }
    if(rawCounts.plot){
      plot.type <- "rawCounts"
      plot.type.title <- IMAGE.NAME.MAP[[plot.type]];
      if(with.TX){
        outfile <- paste(outfile.prefix[7],geneName,"",plot.type.title,IMAGE.NAME.TX,"",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier)
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,colorList = colorList, plot.type = plot.type, use.vst=FALSE,use.log=use.log ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,gene.lwd= gene.lwd,CONNECTIONS.relative.height = CONNECTIONS.relative.height,TX.margins=TX.margins,SPLICE.annotation.relative.height=SPLICE.annotation.relative.height,draw.nested.SJ=draw.nested.SJ,INTERNAL.VARS = INTERNAL.VARS,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[8],geneName,"",plot.type.title,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier)
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,colorList = colorList, plot.type = plot.type, use.vst=FALSE,use.log=use.log,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows, debug.mode = debug.mode,sequencing.type=sequencing.type,gene.lwd= gene.lwd,CONNECTIONS.relative.height = CONNECTIONS.relative.height,TX.margins=TX.margins,SPLICE.annotation.relative.height=SPLICE.annotation.relative.height,draw.nested.SJ=draw.nested.SJ,INTERNAL.VARS = INTERNAL.VARS,...)
        closePlottingDeviceFunc();  
      }
    }
}

getAutofitTxRelativeHeight <- function(TX.ct, autoscale.height.to.fit.TX.annotation = TRUE,
                              GENE.annotation.relative.height = 0.15, 
                              TX.annotation.relative.height = 0.05, 
                              CONNECTIONS.relative.height = 0.1,
                              SPLICE.annotation.relative.height=0.1,
                              TX.margins = c(0,0.5)){
      
    if(autoscale.height.to.fit.TX.annotation){
      GENE.annotation.height <- GENE.annotation.relative.height * 10
      TX.annotation.height <- TX.annotation.relative.height * 10
      CONNECTIONS.height <- CONNECTIONS.relative.height * 10
      withTxPlot.height.multiplier <- (10+CONNECTIONS.height + GENE.annotation.height + SPLICE.annotation.relative.height + (TX.annotation.height * (TX.ct + sum(TX.margins)))  ) / (10+CONNECTIONS.height + GENE.annotation.height + SPLICE.annotation.relative.height + (TX.annotation.height * TX.margins[2]))
    } else {
      withTxPlot.height.multiplier <- 1
    }
    return(withTxPlot.height.multiplier)
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

JUNCTIONSEQ.DEFAULT.COLOR.LIST <- list(
    SIG.VERTLINE.COLOR = "#F219ED60",
    NOSIG.VERTLINE.COLOR = "#99999960",
    UNTESTABLE.VERTLINE.COLOR = "#CCCCCC60",

    SIG.FEATURE.COLOR = "#F219ED",
    NOSIG.FEATURE.COLOR = "#111111",
    UNTESTABLE.FEATURE.COLOR = "#CCCCCC",
    EXCLUDED.FEATURE.COLOR = "#111111",

    SIG.FEATURE.BORDER.COLOR = "#000000",
    NOSIG.FEATURE.BORDER.COLOR = "#000000",
    UNTESTABLE.FEATURE.BORDER.COLOR = "#AAAAAA",
    EXCLUDED.FEATURE.BORDER.COLOR = "#000000",

    SIG.FEATURE.FILL.COLOR = "#F219ED",
    NOSIG.FEATURE.FILL.COLOR = "#CCCCCC",
    UNTESTABLE.FEATURE.FILL.COLOR = "#F5F5F5",
    EXCLUDED.FEATURE.FILL.COLOR = "#CCCCCC",
    
    KNOWN.SPLICE.LTY = "solid",
    NOVEL.SPLICE.LTY = "32",
    EXON.CONNECTION.LTY = "solid",
    NOVEL.SPLICE.CONNECTION.LTY = "32",
    KNOWN.SPLICE.CONNECTION.LTY = "solid",
    
    PLOTTING.LINE.COLORS = color2transparentVector(c("red","blue","orange","green3","purple","cyan1", "magenta","yellow3","tan4"), t = 175)
)

USE.MARGIN.MEX <- FALSE

plotJunctionSeqResultsForGene <- function(geneID, jscs, 
                              colorRed.FDR.threshold=0.01,
                              plot.type = c("expr","normCounts","rExpr","rawCounts"), 
                              sequencing.type = c("paired-end","single-end"),
                              displayTranscripts = FALSE,
                              colorList = list(), 
                              use.vst = FALSE, use.log = TRUE,
                              exon.rescale.factor = 0.3, exonRescaleFunction = c("sqrt","log","linear","34root"), #intron.break.threshold = 0.2, 
                              label.p.vals = TRUE, 
                              plot.lwd = 3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, gene.lwd = plot.lwd / 2,
                              par.cex = 1, 
                              anno.cex.text = 1,
                              anno.cex.axis=anno.cex.text, anno.cex.main = anno.cex.text * 1.2, cex.arrows = "auto",
                              fit.countbin.names = TRUE, fit.genomic.axis = TRUE, fit.labels = TRUE,
                              plot.gene.level.expression = TRUE,
                              plot.exon.results = NULL, plot.junction.results = NULL, plot.novel.junction.results = NULL, 
                              plot.untestable.results = FALSE, draw.untestable.annotation = TRUE,
                              show.strand.arrows = 1, #arrows.length = 0.125,
                              sort.features = TRUE,
                              drawCoordinates = TRUE,
                              yAxisLabels.inExponentialForm = FALSE,
                              title.main=NULL, title.ylab=NULL, title.ylab.right=NULL, 
                              graph.margins = c(2,3,3,3),
                              GENE.annotation.relative.height = 0.15, TX.annotation.relative.height = 0.05, CONNECTIONS.relative.height = 0.1,
                              SPLICE.annotation.relative.height = 0.1,
                              TX.margins = c(0,0.5),
                              condition.legend.text = NULL, include.TX.names = TRUE, draw.start.end.sites = TRUE,
                              label.chromosome=TRUE, 
                              splice.junction.drawing.style = c("hyperbola","ellipse","triangular","line"),
                              draw.nested.SJ = TRUE, merge.exon.parts = TRUE,
                              verbose=TRUE, debug.mode = FALSE, 
                              INTERNAL.VARS = list(),
                              ...)
{
    tryCatch({
      GENE.annotation.height <- GENE.annotation.relative.height * 10
      TX.annotation.height <- TX.annotation.relative.height * 10
      CONNECTIONS.height <- CONNECTIONS.relative.height * 10
      SPLICE.annotation.height <- SPLICE.annotation.relative.height * 10
      
      if(! any( geneID == fData(jscs)$geneID )){
        stop(paste0("FATAL ERROR: GeneID \"",geneID,"\" not found in dataset!" ))
      }
      
      #flat.gff.data <- jscs@flatGffData[jscs@flatGffData$gene_id==geneID, ];
      #merged.data <- fData(jscs)[fData(jscs)$geneID == geneID,];
      flat.gff.data <- jscs@flatGffData;
      merged.data <- fData(jscs);
      
      condition <- jscs@phenoData$condition
      
      exonRescaleFunction <- match.arg(exonRescaleFunction)
      sequencing.type <- match.arg(sequencing.type)
      splice.junction.drawing.style <- match.arg(splice.junction.drawing.style)
      plot.type <- match.arg(plot.type)
      
      truncateBelowOne <- TRUE
      
      if(is.null(plot.exon.results)){
        plot.exon.results <- any( merged.data$featureType == "exonic_part" )
      }
      if(is.null(plot.junction.results)){
        plot.junction.results <- any( merged.data$featureType == "splice_site" | merged.data$featureType == "novel_splice_site" )
      }
      if(is.null(plot.novel.junction.results)){
        if(plot.junction.results){
          plot.novel.junction.results <- any( merged.data$featureType == "novel_splice_site" )
        } else {
          plot.novel.junction.results <- FALSE
        }
      }
      
      flip.splicing <- if( plot.junction.results ){
        FALSE
      } else {
        TRUE
      }
      
      geneName <- jscs@flatGffGeneData$gene_name[jscs@flatGffGeneData$geneID == geneID]

     if(is.null(condition.legend.text)){
       condition.legend.text <- levels(jscs@phenoData$condition)
       names(condition.legend.text) <- condition.legend.text
     } else {
       if(is.null(names(condition.legend.text))){
         warning("names(condition.legend.text is NULL! condition.legend.text mis-formatted. Must be a list or character vector, with element names equal to the levels of pData(jscs)$condition. Falling back.")
         condition.legend.text <- levels(jscs@phenoData$condition)
         names(condition.legend.text) <- condition.legend.text
       } else if(! all(levels(jscs@phenoData$condition) %in% names(condition.legend.text))){
         warning("Not all levels contained in names(condition.legend.text)! condition.legend.text mis-formatted. Must be a list or character vector, with element names equal to the levels of pData(jscs)$condition. Falling back.")
         condition.legend.text <- levels(jscs@phenoData$condition)
         names(condition.legend.text) <- condition.legend.text
       }
     }
     
     if(! is.null( jscs@flatGffGeneData[["aggregateGeneStrand"]] )){
       geneStrand <-  as.character(jscs@flatGffGeneData[["aggregateGeneStrand"]][ jscs@flatGffGeneData[["geneID"]] == geneID ])
       txSetString <- as.character(jscs@flatGffGeneData[["tx_set"]][ jscs@flatGffGeneData[["geneID"]] == geneID ])
       txStrandString <-  as.character(jscs@flatGffGeneData[["tx_strands"]][ jscs@flatGffGeneData[["geneID"]] == geneID ])
       txSet <-  strsplit(as.character(txSetString),"+",fixed=TRUE)[[1]]
       txStrand <-  strsplit(as.character(txStrandString),",",fixed=TRUE)[[1]]
       txStrandMap <- as.list(txStrand)
       names(txStrandMap) <- txSet
     } else {
       txStrandMap <- list()
       geneStrand <- "."
     }

     #SET COLORS FOR ANNOTATION:
     final.color.list <- overmerge.list(JUNCTIONSEQ.DEFAULT.COLOR.LIST, colorList)
      SIG.VERTLINE.COLOR = final.color.list[["SIG.VERTLINE.COLOR"]]
      NOSIG.VERTLINE.COLOR = final.color.list[["NOSIG.VERTLINE.COLOR"]]
      UNTESTABLE.VERTLINE.COLOR = final.color.list[["UNTESTABLE.VERTLINE.COLOR"]]

      SIG.FEATURE.COLOR = final.color.list[["SIG.FEATURE.COLOR"]]
      NOSIG.FEATURE.COLOR = final.color.list[["NOSIG.FEATURE.COLOR"]]
      UNTESTABLE.FEATURE.COLOR = final.color.list[["UNTESTABLE.FEATURE.COLOR"]]
      EXCLUDED.FEATURE.COLOR = final.color.list[["EXCLUDED.FEATURE.COLOR"]]

      SIG.FEATURE.BORDER.COLOR = final.color.list[["SIG.FEATURE.BORDER.COLOR"]]
      NOSIG.FEATURE.BORDER.COLOR = final.color.list[["NOSIG.FEATURE.BORDER.COLOR"]]
      UNTESTABLE.FEATURE.BORDER.COLOR = final.color.list[["UNTESTABLE.FEATURE.BORDER.COLOR"]]
      EXCLUDED.FEATURE.BORDER.COLOR = final.color.list[["EXCLUDED.FEATURE.BORDER.COLOR"]]

      SIG.FEATURE.FILL.COLOR = final.color.list[["SIG.FEATURE.FILL.COLOR"]]
      NOSIG.FEATURE.FILL.COLOR = final.color.list[["NOSIG.FEATURE.FILL.COLOR"]]
      UNTESTABLE.FEATURE.FILL.COLOR = final.color.list[["UNTESTABLE.FEATURE.FILL.COLOR"]]
      EXCLUDED.FEATURE.FILL.COLOR = final.color.list[["EXCLUDED.FEATURE.FILL.COLOR"]]

      PLOTTING.LINE.COLORS = final.color.list[["PLOTTING.LINE.COLORS"]]
     
     gene.level.buffer <- 0.5

     if(is.null(plot.gene.level.expression)){
       if(use.vst){
         plot.gene.level.expression <- FALSE
       } else {
         plot.gene.level.expression <- TRUE
       }
     }
     if(plot.gene.level.expression && use.vst){
       warning("WARNING: plotting of gene-level expression is not supported for vst-transformed plots. Errors are likely to follow.")
     }
     
     FDR <- colorRed.FDR.threshold
     
     

     if(verbose){
        displayTXstring <- if(displayTranscripts) " (with TX)" else ""
        message("> pJSRfG(): ", geneID, ", plot.type: ", plot.type,displayTXstring)
     }
     
     chrom.label <- as.character(merged.data$chr[merged.data$geneID == geneID][1])
     
     rt <- merged.data$geneID == geneID
     if(! plot.exon.results){
       if(debug.mode) message(">     Removing ",sum(rt & merged.data$featureType == "exonic_part")," exonic_part features. ",sum(rt & merged.data$featureType != "exonic_part")," features remaining")
       rt <- rt & merged.data$featureType != "exonic_part" 
     }
     if(! plot.junction.results){
       if(debug.mode) message(">     Removing ",sum(rt & merged.data$featureType == "splice_site")," splice_site features. ",sum(rt & merged.data$featureType != "splice_site")," features remaining")
       rt <- rt & merged.data$featureType != "splice_site" 
     }
     if(! plot.novel.junction.results){
       if(debug.mode) message(">     Removing ",sum(rt & merged.data$featureType == "novel_splice_site")," novel_splice_site features. ",sum(rt & merged.data$featureType != "novel_splice_site")," features remaining")
       rt <- rt & merged.data$featureType != "novel_splice_site" 
     }
     untestable.rt <- which(rt & (! merged.data$testable))
     if(! plot.untestable.results){
       if(debug.mode) message(">     Removing ",sum(rt & ! merged.data$testable)," untestable features. ",sum(rt & merged.data$testable)," features remaining")
       rt <- rt & merged.data$testable
     }
     if(debug.mode)   message(">     Plotting ",sum(rt), " features.")

     rt <- which(rt)
     if(length(rt) == 0){
       message("NO FEATURES TO PLOT!")
     } else {

       rt.allExon <- which(flat.gff.data$gene_id==geneID & flat.gff.data$featureType == "exonic_part")
       rango.allExon <- 1:length(rt.allExon)

       rt.allJunction <- which(flat.gff.data$gene_id==geneID & (flat.gff.data$featureType == "splice_site" | flat.gff.data$featureType == "novel_splice_site"))
       rango.allJunction <- 1:length(rt.allJunction)

          rescale.iv <- generate.interval.scale(
              data.frame(
                start = c(flat.gff.data$start[rt.allExon], flat.gff.data$start[rt.allJunction]),
                end = c(flat.gff.data$end[rt.allExon], flat.gff.data$end[rt.allJunction]),
                is.exon = c(rep(TRUE,length(rt.allExon)), rep(FALSE,length(rt.allJunction)))
              ),exon.rescale.factor, exonRescaleFunction, debug.mode = debug.mode)
          rel <- data.frame(start = rescale.coords(merged.data$start[rt],rescale.iv), 
                            end   = rescale.coords(merged.data$end[rt],  rescale.iv))
          
       if(sort.features){
         #Reorder based on rescaled values:
         rt <- rt[order(((rel$end - rel$start)/2) + rel$start )]
       }

       rango <- 1:length(rt)
       draw.legend <- TRUE
       condition.names <- levels(condition)
       sample.names <- sampleNames(jscs@phenoData)
       colorcode.title <- TRUE
       numcond <- length(condition.names)

       if(debug.mode && verbose) message(">    pJSRforGene(): ","Reached step 2.")

       ####### DETERMINE COLORS, IF THE USER DOES NOT PROVIDE ONE PER SAMPLE THE COUNT WILL OBTAIN THEM CORRESPONDING TO THEIR DESIGN ####
       ##### determine colors if not provided by user ######
       
       
          default.color.list <- PLOTTING.LINE.COLORS

          if(numcond > length(default.color.list)){
             message("Too many condition values, the default color selection may not look good! Set your own colors by setting the \"color\" parameter.")
             color <- rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond)), maxColorValue=255, alpha=175)
          } else {
             color <- color2transparentVector(default.color.list[1:numcond], t = 175)
          }
          
    if(debug.mode && verbose) message(">    pJSRforGene(): ","color = ", paste0(color, collapse=","))
    if(debug.mode && verbose) message(">    pJSRforGene(): ","length(color) = ", length(color))
    if(debug.mode && verbose) message(">    pJSRforGene(): ","length(condition.names) = ", length(condition.names))
    if(debug.mode && verbose) message(">    pJSRforGene(): ","condition.names = ", paste0(condition.names,collapse=","))

       names(color) <- condition.names

       y.axis.title <- ""
       main.title <- ""

    if(debug.mode && verbose) message(">    pJSRforGene(): ","Reached step 3.")

       if(truncateBelowOne){
         convertY <- function(y){
           ifelse(y < 0, ifelse(is.infinite(y), INTERNAL.NINF.VALUE, (1-exp(y)) * INTERNAL.NINF.VALUE), y)
         }
       } else {
         convertY <- function(y){
           y
         }
       }

       if(plot.type == "rExpr"){
         count <- if(use.vst) { 
                    vst( jscs@plottingEstimates[["relExprEstimate"]][rt,, drop=FALSE], jscs)
                  } else  if(use.log) {
                    apply(log10(jscs@plottingEstimates[["relExprEstimate"]][rt,, drop=FALSE]),c(1,2),FUN=convertY)
                  } else {
                    jscs@plottingEstimates[["relExprEstimate"]][rt,, drop=FALSE]
                  }
         color.count <- rep(color[condition.names],each=nrow(count))
         y.axis.title <- "Relative Coverage"
         y.axis.title.right <- "Gene-Level Mean Normalized Counts"
         main.title <- paste0("Relative Coverage (",geneName,")")
         
         if(plot.gene.level.expression) {
           geneCount <- if(use.vst){ 
                    plot.gene.level.expression <- FALSE
                    NULL
                  } else if(use.log){
                    sapply(log10(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]][rownames(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]]) == geneID ,]),FUN=convertY)
                  } else {
                    jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]][rownames(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]]) == geneID ,]
           }
           color.geneCount <- color[condition.names]
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
                    plot.gene.level.expression <- FALSE
                    NULL
                  } else if(use.log){
                    sapply(log10(jscs@geneCountData[rownames(jscs@geneCountData) == geneID,] / sizeFactors(jscs)),FUN=convertY)
                  } else {
                    jscs@geneCountData[rownames(jscs@geneCountData) == geneID,] / sizeFactors(jscs)
           }
           color.geneCount <- color[as.character(condition)]
         }
         color.count <- rep(color[as.character(condition)], each=nrow(count))

         y.axis.title <- "Normalized Counts"
         y.axis.title.right <- "Gene-Level Normalized Counts"
         main.title <- paste0("Normalized Counts (",geneName,")")
       } else if(plot.type == "rawCounts"){
         count <- if(use.vst){ 
                    vst(jscs@countVectors[rt,1:length(sample.names), drop=FALSE], jscs)
                  } else if(use.log) {
                    apply(log10(jscs@countVectors[rt,1:length(sample.names), drop=FALSE]),c(1,2),FUN=function(y){max(INTERNAL.NINF.VALUE,y)})
                  } else {
                    jscs@countVectors[rt,1:length(sample.names), drop=FALSE]
                  }
         if(plot.gene.level.expression) {
           geneCount <- if(use.vst){ 
                    plot.gene.level.expression <- FALSE
                    NULL
                  } else if(use.log){
                    sapply(log10(jscs@geneCountData[rownames(jscs@geneCountData) == geneID,]),FUN=convertY)
                  } else {
                    jscs@geneCountData[rownames(jscs@geneCountData) == geneID,]
           }
           color.geneCount <- color[as.character(condition)]
         }
         color.count <- rep(color[as.character(condition)], each=nrow(count))
         y.axis.title <- "Raw Counts"
         y.axis.title.right <- "Raw Gene-Level Counts"
         main.title <- paste0("Raw Counts (",geneName,")")
       } else if(plot.type == "expr"){
         count <- if(use.vst){ 
                    vst(jscs@plottingEstimates[["exprEstimate"]][rt,, drop=FALSE], jscs)
                  } else if(use.log) {
                    apply(log10(jscs@plottingEstimates[["exprEstimate"]][rt,, drop=FALSE]),c(1,2),FUN=convertY)
                  } else {
                    jscs@plottingEstimates[["exprEstimate"]][rt,, drop=FALSE]
                  }
         if(plot.gene.level.expression) {
           geneCount <- if(use.vst){ 
                    plot.gene.level.expression <- FALSE
                    NULL
                  } else if(use.log){
                    sapply(log10(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]][rownames(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]]) == geneID ,]),FUN=convertY)
                  } else {
                    jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]][rownames(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]]) == geneID ,]
           }
           color.geneCount <- color[condition.names]
         }
         color.count <- rep(color[condition.names],each=nrow(count))
         y.axis.title <- "Mean Normalized Counts"
         y.axis.title.right <- "Gene-Level Mean Normalized Counts"
         main.title <- paste0("Mean Normalized Coverage (",geneName,")")
       } else {
         stop(paste0("FATAL ERROR: Unknown plot type! plot.type = \"",plot.type,"\""))
       }
       
       ###########################################################
       #Default titles and axes labels:
       fragment.label <- if(sequencing.type == "paired-end") "Read-Pair" else "Read"
       
       if(plot.type == "rExpr"){
         main.title <- paste0("Relative Coverage (",geneName,")")
         y.axis.title <- paste0("Relative Coverage")
         y.axis.title.right <- paste0(fragment.label,"s per Sample, Gene-Level")
         
       } else if(plot.type == "normCounts"){
         main.title <- paste0("Normalized Counts (",geneName,")")
         y.axis.title <- paste0("Normalized Counts")
         y.axis.title.right <- paste0("Normalized Counts, Gene-Level")
         
       } else if(plot.type == "rawCounts"){
         main.title <- paste0("Raw Counts (",geneName,")")
         y.axis.title <- paste0("Raw ",fragment.label," Counts")
         y.axis.title.right <- paste0("Raw ",fragment.label," Counts, Gene-Level")
         
       } else if(plot.type == "expr"){
         main.title <- paste0("Mean Normalized Coverage (",geneName,")")
         y.axis.title <- paste0(fragment.label,"s per Sample")
         y.axis.title.right <- paste0(fragment.label,"s per Sample, Gene-Level")
         
       } else {
         stop(paste0("FATAL ERROR: Unknown plot type! plot.type = \"",plot.type,"\""))
       }
       # Option: Allow users to explicitly override labels and titles:
       if(! is.null(title.main)) main.title <- title.main
       if(! is.null(title.ylab)) y.axis.title <- title.ylab
       if(! is.null(title.ylab.right)) y.axis.title.right <- title.ylab.right
       
       if(! plot.gene.level.expression) {
         geneCount <- NULL
         color.geneCount <- NULL
       }
       count <- as.matrix(count)

    if(debug.mode && verbose) message(">    pJSRforGene(): ", "Reached step 4.")

       intervals<-(0:nrow(count))/nrow(count)

       numexons<-nrow(count)
       each <- merged.data$padjust[rt]
       vertline.col <- ifelse(merged.data$testable[rt], ifelse(f.na(merged.data$padjust[rt] <= FDR), SIG.VERTLINE.COLOR, NOSIG.VERTLINE.COLOR), UNTESTABLE.VERTLINE.COLOR)
       annolink.col <- ifelse(merged.data$testable[rt], ifelse(f.na(merged.data$padjust[rt] <= FDR), SIG.FEATURE.COLOR,  NOSIG.FEATURE.COLOR),  UNTESTABLE.FEATURE.COLOR)

       exonlty <- rep(final.color.list[["EXON.CONNECTION.LTY"]],length(vertline.col))
       exonlty[as.character(merged.data$featureType[rt]) == "novel_splice_site"] <- final.color.list[["NOVEL.SPLICE.CONNECTION.LTY"]]
       exonlty[as.character(merged.data$featureType[rt]) == "splice_site"] <- final.color.list[["KNOWN.SPLICE.CONNECTION.LTY"]]

       is.sig.feature <- f.na(each <= FDR)
       sig.feature <- which(is.sig.feature)

    if(debug.mode && verbose) message(">    pJSRforGene(): ", "Reached step 5.")

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
                            stringsAsFactors=FALSE)
          sub.allExon <- data.frame(start=flat.gff.data$start[rt.allExon], 
                                    end=flat.gff.data$end[rt.allExon], 
                                    chr=flat.gff.data$chrom[rt.allExon], 
                                    strand=flat.gff.data$strand[rt.allExon], 
                                    is.exon = (flat.gff.data$featureType[rt.allExon] == "exonic_part"),
                                    featureID = as.character(flat.gff.data$featureName[rt.allExon]), 
                                    stringsAsFactors=FALSE)
          sub.allJunction <- data.frame(start=flat.gff.data$start[rt.allJunction], 
                                    end=flat.gff.data$end[rt.allJunction], 
                                    chr=flat.gff.data$chrom[rt.allJunction], 
                                    strand=flat.gff.data$strand[rt.allJunction], 
                                    is.exon = (flat.gff.data$featureType[rt.allJunction] == "exonic_part"), 
                                    feature.type = flat.gff.data$featureType[rt.allJunction] ,
                                    is.novel = (flat.gff.data$featureType[rt.allJunction] == "novel_splice_site"), 
                                    featureID = as.character(flat.gff.data$featureName[rt.allJunction]), 
                                    stringsAsFactors=FALSE)
          sub.allJunction$is.plotted <- sub.allJunction$featureID %in% sub$featureID
          
          #print(sub.allJunction)
          
       testable.featureIDs <- sub$featureID[sub$is.testable]
       sig.featureIDs <- sub$featureID[sub$is.sig]
       untestable.featureIDs <- rownames(merged.data)[untestable.rt]

       if(debug.mode && verbose){ message("testable.featureIDs: ",paste0(testable.featureIDs,collapse=",")) }
       if(debug.mode && verbose){ message("sig.featureIDs: ",paste0(sig.featureIDs,collapse=",")) }
       if(debug.mode && verbose){ message("untestable.featureIDs: ",paste0(untestable.featureIDs,collapse=",")) }


       sub.allExon$is.testable   <- sub.allExon$featureID %in% testable.featureIDs
       sub.allExon$is.sig        <- sub.allExon$featureID %in% sig.featureIDs
       sub.allExon$is.untestable <- sub.allExon$featureID %in% untestable.featureIDs
       sub.allJunction$is.testable <- sub.allJunction$featureID %in% testable.featureIDs
       sub.allJunction$is.sig      <- sub.allJunction$featureID %in% sig.featureIDs
       sub.allJunction$is.untestable <- sub.allJunction$featureID %in% untestable.featureIDs

       sub.allExon$lineColor     <- ifelse(sub.allExon$is.testable,
                                           ifelse(sub.allExon$is.sig,
                                                  SIG.FEATURE.COLOR,
                                                  NOSIG.FEATURE.COLOR),
                                           ifelse(sub.allExon$is.untestable, 
                                                  UNTESTABLE.FEATURE.COLOR, 
                                                  EXCLUDED.FEATURE.COLOR))
       sub.allExon$fillColor     <- ifelse(sub.allExon$is.testable,ifelse(sub.allExon$is.sig,SIG.FEATURE.FILL.COLOR,NOSIG.FEATURE.FILL.COLOR),ifelse(sub.allExon$is.untestable, UNTESTABLE.FEATURE.FILL.COLOR, EXCLUDED.FEATURE.FILL.COLOR))
       sub.allExon$borderColor   <- ifelse(sub.allExon$is.testable,ifelse(sub.allExon$is.sig,SIG.FEATURE.BORDER.COLOR,NOSIG.FEATURE.COLOR),ifelse(sub.allExon$is.untestable, UNTESTABLE.FEATURE.BORDER.COLOR, EXCLUDED.FEATURE.BORDER.COLOR))
       sub.allJunction$lineColor <- ifelse(sub.allJunction$is.testable,ifelse(sub.allJunction$is.sig,SIG.FEATURE.COLOR,NOSIG.FEATURE.COLOR),ifelse(sub.allJunction$is.untestable, UNTESTABLE.FEATURE.COLOR, EXCLUDED.FEATURE.COLOR))
       sub.allJunction$lty <- ifelse(sub.allJunction$is.novel,final.color.list[["NOVEL.SPLICE.LTY"]],final.color.list[["KNOWN.SPLICE.LTY"]])

       sig.feature.names <- sub$featureID[is.sig.feature & sub$is.exon]
       allExon.isSig <- sub.allExon$featureID %in% sig.featureIDs
       allExon.exonCol <- ifelse(allExon.isSig, "#F219ED", "#CCCCCC")

       #TEMPORARY DEBUGGING MESSAGES:
    if(debug.mode == 2 && verbose){
         message(">    Debugging Info:")
         message("        Exons:")
         message("           ", paste0(names(sub.allExon),collapse='\t'))
         for(i in 1:length(sub.allExon$start)){
           message("           ", paste0(sub.allExon[i,],collapse='\t'))
         }
         message(">       dim(count) = ", paste0(dim(count), collapse=","))
         message(">       rt: ", paste0(rt, collapse=","))
         message(">       rango: ", paste0(rango, collapse=","))
         message(">       ncol(count): ", ncol(count))
    }

          sub.sig <- sub[sig.feature,, drop=FALSE]
          rel.calc.min <- min(sub.allJunction$start, sub.allExon$start)
          rel.calc.max <- max(sub.allJunction$end,   sub.allExon$end)



          transcripts <- sapply(sapply(flat.gff.data$transcripts[rt.allExon],toString), function(x){strsplit(x, "+",fixed=TRUE)})
          trans <- Reduce(union, transcripts)
          if(displayTranscripts==TRUE){
             mat <- 1:4
             hei<-c(10,CONNECTIONS.height, GENE.annotation.height + SPLICE.annotation.height, TX.annotation.height * ( length(trans) + TX.margins[1] + TX.margins[2] ))
          }else{
             mat<-1:4
             hei<-c(10,CONNECTIONS.height, GENE.annotation.height + SPLICE.annotation.height, TX.annotation.height * TX.margins[2] )
          }
          
          layout(matrix(mat), heights=hei)
          
          if(debug.mode && verbose){
            message(">   length(trans) = ",length(trans))
            message(">   FINAL LAYOUT:")
            message(">      heights = [", paste0(hei, collapse=","),"]")
          }
          
    if(debug.mode && verbose) message(">    pJSRforGene(): ","Reached step 6.")

          ylimn <- c(min(min(count,na.rm=TRUE),0), max(count, na.rm=TRUE))
          if((! use.vst) && use.log ) ylimn[1] <- INTERNAL.NINF.VALUE

          p.values.labels <- ifelse(f.na(each<=FDR), format(each,digits=3), "")

          if(any(sub$is.exon)){
            italicize.label <- ! sub$is.exon
          } else {
            italicize.label <- NULL
          }
          
          plotWindowXmax <- if(plot.gene.level.expression){  (length(intervals) + 1) / length(intervals)  } else { 1.00 }
          
          intervals <- drawPlot(matr=count, ylimn,jscs, 
                   intervals, rango, textAxis=y.axis.title, geneLevelAxisTitle = y.axis.title.right,
                   rt=rt, color.count=color.count, 
                   colorlines=vertline.col, 
                   countbinIDs = merged.data$countbinID[rt],
                   use.vst=use.vst, use.log = use.log,plot.type=plot.type,
                   main.title=main.title,draw.legend=draw.legend,
                   color.key=color,condition.names=condition.names,
                   p.values=p.values.labels,draw.p.values=label.p.vals, 
                   plot.lwd=plot.lwd, axes.lwd = axes.lwd, 
                   anno.lwd = anno.lwd, par.cex = par.cex, 
                   anno.cex.text = anno.cex.text,
                   anno.cex.axis = anno.cex.axis, 
                   anno.cex.main = anno.cex.main,
                   fit.countbin.names = fit.countbin.names,
                   debug.mode = debug.mode, plot.gene.level.expression = plot.gene.level.expression, geneCount = geneCount, color.geneCount = color.geneCount,
                   yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, italicize.label = italicize.label, condition.legend.text = condition.legend.text,
                   annolink.col = annolink.col, exonlty = exonlty, graph.margins = graph.margins, plotWindowXmax = plotWindowXmax, fit.labels =fit.labels,
                   ...)
    if(debug.mode && verbose) message(">    pJSRforGene(): ","Reached end of step 6.")


    if(debug.mode && verbose) message(">    pJSRforGene(): ","Reached step 7.")

          #########PLOT THE GENE MODEL:
          if(USE.MARGIN.MEX){
            par(mar=c(0, graph.margins[2], 0, graph.margins[4]), cex = par.cex, mex = anno.cex.text); 
          } else {
            par(mar=c(0, graph.margins[2], 0, graph.margins[4]), cex = par.cex); 
          }
          
          plot.new()
          plot.window(xlim=c(0, plotWindowXmax), ylim=c(0,1), xaxs = "i")
          rel <- data.frame(start = rescale.coords(sub$start,rescale.iv), 
                            end   = rescale.coords(sub$end,  rescale.iv))
            connection.lines.bottom <- apply((rbind(rel[rango,2], rel[rango, 1])), 2, median) * plotWindowXmax
            connection.lines.top    <- apply(rbind(intervals[rango], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2)), 2, median)
          
          segments(     
                        connection.lines.bottom, 
                        0, #par("usr")[3], (old version: lines are connected.)
                        connection.lines.top, 
                        1, col=annolink.col, lty = exonlty, lwd = gene.lwd, cex = anno.cex.text,cex.axis=anno.cex.main, cex.main=anno.cex.main, xpd=NA, ...) #col=colorlinesB,...)
          
          par(mar=c(1.5, graph.margins[2], 0, graph.margins[4]), cex = par.cex)

          #Get start/end sites:
          startSites <- c()
          endSites <- c()
          for(i in 1:length(trans)){
             logicexons <- sapply(transcripts, function(x){any(x == trans[i])})
             startSites <- c(startSites, sub.allExon$start[logicexons][1])
             endSites <- c(endSites, sub.allExon$end[logicexons][sum(logicexons)])
          }
          startSites <- unique(startSites)
          endSites <- unique(endSites)

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
                   plot.lwd = gene.lwd, anno.lwd=anno.lwd, 
                   show.strand.arrows = show.strand.arrows, 
                   geneStrand = geneStrand,
                   cex.axis=anno.cex.axis, cex.main=anno.cex.main, 
                   draw.untestable.annotation = draw.untestable.annotation,
                   draw.start.end.sites = draw.start.end.sites, startSites = startSites, endSites = endSites,
                   cex.arrows = cex.arrows, chrom.label = chrom.label, label.chromosome = label.chromosome,
                   splice.junction.drawing.style = splice.junction.drawing.style,
                   draw.nested.SJ = draw.nested.SJ, merge.exon.parts= merge.exon.parts,
                   plot.untestable.results = plot.untestable.results,
                   exon.height = GENE.annotation.height / (SPLICE.annotation.height + GENE.annotation.height),
                   INTERNAL.VARS = INTERNAL.VARS,
                   flip.splicing = flip.splicing,
                   ...)
           #Maybe make these options external at some point?
           include.endpoints.on.coordinates <- FALSE
           num.coord.miniticks.per.tick <- 10

       if(drawCoordinates){
             if(! is.null(rescale.iv)){
               pretty.x <- pretty(c(rel.calc.min,rel.calc.max),n=5)
               pretty.interval <- pretty.x[2] - pretty.x[1]
               pretty.x <- pretty.x[pretty.x > rel.calc.min & pretty.x < rel.calc.max]
               
               rescaled.pretty.x <- rescale.coords(pretty.x,rescale.iv)

               if(num.coord.miniticks.per.tick > 0){
                 rel.coord.miniticks <- pretty.interval * (1:(num.coord.miniticks.per.tick-1)) / num.coord.miniticks.per.tick
                 unscaled.coord.miniticks <- unlist(lapply(pretty.x, function(a){ a + rel.coord.miniticks }))
                 unscaled.coord.miniticks <- c(pretty.x[1] - rel.coord.miniticks, unscaled.coord.miniticks)
                 rescaled.coord.miniticks <- rescale.coords(unscaled.coord.miniticks, rescale.iv) * (rel.calc.max-rel.calc.min) + rel.calc.min
               } else {
                 rescaled.coord.miniticks <- FALSE
               }
               
               if(include.endpoints.on.coordinates){
                 if(min(rescaled.pretty.x) > 0.05){
                   pretty.x <- c(rel.calc.min, pretty.x)
                   rescaled.pretty.x <- c(0,rescaled.pretty.x)
                 }
                 if(max(rescaled.pretty.x) < 0.95){
                   pretty.x <- c(pretty.x, rel.calc.max)
                   rescaled.pretty.x <- c(rescaled.pretty.x,1)
                 }
               }
               rescaled.pretty.x <- rescaled.pretty.x * (rel.calc.max-rel.calc.min) + rel.calc.min
             } else {
               pretty.x <- pretty(c(rel.calc.min,rel.calc.max),n=5)
               coord.miniticks <- FALSE
             }
             
             usr <- par("usr")
             cxy <- par("cxy")
             pretty.x <- sprintf("%0.f",pretty.x)
             smallest.width.coordAxis <- min(abs(rescaled.pretty.x[-1] - rescaled.pretty.x[-length(rescaled.pretty.x)]))
             #Fit the genomic labels
             if(fit.genomic.axis){ 
               anno.cex.coordAxis <- shrink.character.vector(paste0(pretty.x,"0"), curr.cex = anno.cex.axis, max.width = smallest.width.coordAxis)
               #If it won't fit, remove labels until it does fit:
               if(anno.cex.coordAxis * 2 < anno.cex.axis){
                 anno.cex.coordAxis <- anno.cex.axis / 2
                 coordAxis.widths <- abs(rescaled.pretty.x[-1] - rescaled.pretty.x[-length(rescaled.pretty.x)])
                 for(i in 2:length(pretty.x)){
                   if((strwidth(pretty.x[i-1], cex = anno.cex.coordAxis) / 2) + (strwidth(pretty.x[i], cex = anno.cex.coordAxis)/2) > abs(rescaled.pretty.x[i-1] - rescaled.pretty.x[i])){
                     pretty.x[i] <- ""
                   }
                 }
               }
             } else {
               anno.cex.coordAxis <- anno.cex.axis
             }
             devlim <- device.limits()
             coord.ticks.top <- usr[3]
             coord.mainTicks.bottom <- usr[3] - (cxy[2] / 2)
             coord.text.top <- usr[3] - (cxy[2] * (3/4))
             
             coord.miniTicks.bottom <- coord.ticks.top - abs(coord.mainTicks.bottom - coord.ticks.top)/2
             
             segments(x0 = rescaled.pretty.x,            y0 = usr[3], x1 = rescaled.pretty.x,            y1 = usr[3] - (cxy[2] / 2), xpd=NA, lwd = anno.lwd, ...)
             lines(c(rel.calc.min,rel.calc.max), c(par("usr")[3], par("usr")[3]) , lwd = axes.lwd, xpd = NA, ...)
             segments(x0 = rescaled.coord.miniticks,     y0 = usr[3], x1 = rescaled.coord.miniticks,     y1 = usr[3] - (cxy[2] / 4), xpd=NA, lwd = anno.lwd, ...)
             text(rescaled.pretty.x, usr[3] - (cxy[2] * (3/4)), pretty.x, cex = anno.cex.coordAxis, xpd = NA, adj = c(0.5,1), ...)
             
       }
       if(label.chromosome){
                chrom.label.width.max <- abs(par("usr")[1] - device.limits()[1]) * 0.9
                chrom.label.cex <- shrink.character.vector(chrom.label, curr.cex = anno.cex.text, max.width = chrom.label.width.max)
                text(par("usr")[1],par("usr")[3],chrom.label, cex = chrom.label.cex, adj = c(1.1,0.5),xpd=NA, font = 2, ...)
       }
       if(displayTranscripts){
         if(USE.MARGIN.MEX){
           par(cex = par.cex, mar=c(0, graph.margins[2], 0, graph.margins[4]), mex = anno.cex.text)
         } else {
           par(cex = par.cex, mar=c(0, graph.margins[2], 0, graph.margins[4]))
         }
         plot.new()
         plot.window(xlim=c(rel.calc.min, rel.calc.max), ylim=c(-TX.margins[2],length(trans) + TX.margins[1]), xaxs = "i", yaxs = "i")
         
         for(i in 1:length(trans)){
            ymin <- length(trans) - i
            if(include.TX.names){
              trName = trans[i]
            } else {
              trName = NULL
            }
            logicexons <- sapply(transcripts, function(x){any(x == trans[i])})
            tr <- sub.allExon[logicexons,]  #   data.frame(start = sub.allExon$start[logicexons==1], end = sub.allExon$end[logicexons==1], featureType = sub.allExon$featureType[logicexons==1], stringsAsFactors = F)
            curr.exoncol <- ifelse(allExon.isSig[logicexons],"#F219ED", "#CCCCCC")
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
                           anno.lwd=gene.lwd, 
                           cex.axis=anno.cex.axis, 
                           cex.main=anno.cex.main, 
                           cex.arrows = cex.arrows,
                           ...)
         }
      } else {
        #TEMP CHANGE: CHANGE ME BACK!
        par(mar = c(0,0,0,0))
        plot.new()
      }

      par(mar = c(0,0,0,0))
      if(debug.mode && verbose) message("> pJSRfG(): "," Done.")
    }
  }, error = function(e){
    message("Error caught while attempting plotJunctionSeqResultsForGene")
    message("---------------------")
    message("     Error text:")
    message("     ",e)
    message("")
    message("---------------------")
    message("     Input parameters:")
    message("     geneID = ",geneID)
    message("     plot.type = ",plot.type)
    message("     displayTranscripts = ",displayTranscripts)
    message("---------------------")
    
    warning("Error caught while attempting plotJunctionSeqResultsForGene")
    warning("Error text:")
    warning(e)
  })
}

