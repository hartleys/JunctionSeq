

########################################################################################
### INPUT:



########################################################################################
### OUTPUT:
#buildAllPlots <- function(jscs,
#                          outfile.prefix = "./",
#                          #flat.gff.data = NULL, flat.gff.file = NULL, 
#                          gene.list = NULL, FDR.threshold = 0.01, 
#                          use.plotting.device = "png", 
#                          use.vst=FALSE,use.log = TRUE, truncateBelowOne = TRUE, 
#                          exon.rescale.factor = 0.3,  
#                          subdirectories.by.type = TRUE,
#                          ma.plot=FALSE, variance.plot=FALSE,
#                          with.TX=TRUE,without.TX=TRUE,
#                          expr.plot=TRUE,normCounts.plot=TRUE,
#                          rExpr.plot=FALSE,rawCounts.plot=FALSE,
#                          colorRed.FDR.threshold = FDR.threshold, 
#                          color=NULL,
#                          plot.gene.level.expression = NULL,
#                          plot.exon.results = NULL, plot.junction.results = NULL, plot.novel.junction.results = NULL,
#                          plot.untestable.results = FALSE,
#                          plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
#                          par.cex = 1, anno.cex.text = 1, anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
#                          drawCoordinates = TRUE, yAxisLabels.inExponentialForm = FALSE,
#                          show.strand.arrows = 10, arrows.length = 0.125,
#                          graph.margins = c(2,3,3,2),
#                          base.plot.height = 12, base.plot.width = 12, 
#                          base.plot.units = "in", 
#                          GENE.annotation.relative.height = 0.2, TX.annotation.relative.height = 0.05, 
#                          autoscale.height.to.fit.TX.annotation = TRUE,
#                          autoscale.width.to.fit.bins = 35,
#                          plotting.device.params = list(), 
#                          number.plots = FALSE,
#                          condition.legend.text = NULL, include.TX.names = TRUE, draw.start.end.sites = TRUE,
#                          openPlottingDeviceFunc = NULL, closePlottingDeviceFunc = NULL,
#                          verbose=TRUE, 
#                          ...){
#expr.plot=TRUE
#normCounts.plot=TRUE
#rExpr.plot=FALSE
#rawCounts.plot=FALSE
#
#HTML output:

JunctionSeqHTML <- function(jscs, 
                            outfile.dir="./", mainFile="testForDU.html", 
                            gene.list = NULL, 
                            FDR.threshold=0.01,
                            colorRed.FDR.threshold = FDR.threshold,
                            plotting.device.ext = ".png", 
                            use.vst=FALSE,use.log = TRUE,
                            subdirectories.by.type = TRUE,
                            ma.plot=FALSE, variance.plot=FALSE,
                            with.TX=TRUE,without.TX=TRUE,
                            expr.plot=TRUE,normCounts.plot=TRUE,
                            rExpr.plot=FALSE,rawCounts.plot=FALSE,
                            plot.exon.results = NULL, plot.junction.results = NULL, plot.novel.junction.results = NULL,
                            plot.untestable.results = FALSE,
                            base.html.height = 90, base.html.height.units = "vh",
                            GENE.annotation.relative.height = 0.2, TX.annotation.relative.height = 0.05, 
                            autoscale.height.to.fit.TX.annotation = TRUE,
                            autoscale.width.to.fit.bins = 35,
                            number.plots = NULL,
                            css.file = NULL, css.link = NULL,
                            verbose = TRUE, debug.mode = FALSE){
                            
   if(is.null(css.link)){
     if(is.null(css.file)){
       if(verbose) message("Copying default css stylesheet.");
       default.css.filepath <- system.file("extdata/styles.css", 
                               package="JunctionSeq",
                               mustWork=TRUE);
       cssLines <- readLines(default.css.filepath);
       cssOut <- file(paste0(outfile.dir,"/styles.css"), "w");
       for(line in cssLines){
         writeLines(line, cssOut);
       }
       close(cssOut);
     } else {
       if(verbose) message("Copying css stylesheet: \"",css.file,"\"");
       cssLines <- readLines(css.file);
       cssOut <- file(paste0(outfile.dir,"/styles.css"), "w");
       for(line in cssLines){
         writeLines(line, cssOut);
       }
       close(cssOut);
     }
     css.path.MAIN <- "styles.css";
     if(subdirectories.by.type){
       css.path.SUB <- "../styles.css";
     } else {
       css.path.SUB <- "styles.css";
     }
   } else {
     if(verbose) message("Linking to external css stylesheet: \"",css.link,"\"");
     css.path.MAIN <- css.link;
     css.path.SUB <- css.link;
   }
   
   fitExpToVar="condition";
   options(stringsAsFactors = F);
   #if(require(hwriter)){
     pf <- file(paste0(outfile.dir,"/", mainFile), "w");
     #p <- hwriter::openPage(file.path(path, file));
     
     writeLines(paste0("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">"), pf);
     writeLines(paste0("<html><head><title> JunctionSeq Results </title>"), pf);
     writeLines(paste0("<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">"), pf);
     writeLines(paste0("<link rel=\"stylesheet\" type=\"text/css\" href=\"",css.path.MAIN,"\">"), pf);
     writeLines(paste0("</head> <body>"), pf);
     
     writeLines(paste0("<h1> JunctionSeq Results </h1>"),pf);
     writeLines(paste0("<h4> Generated using JunctionSeq v",packageVersion("JunctionSeq")," at ",date()," </h4>"),pf);
     
     if(ma.plot | variance.plot){
       writeLines(paste0("<table><tr><th> Summary Plots </th></tr>"),pf);
       writeLines("<tr>",pf);
       if(variance.plot){
         disp.file <- paste0("dispersion-plot",plotting.device.ext);
         writeLines(paste0("<td><a href=\"",disp.file,"\"> ",
                           "<img src=\"",disp.file,"\" "," alt=\"Dispersion Plot\" style=\"width: 45vw;\">",
                           "</a></td>"),pf);
       }
       if(ma.plot){
         fc.cols <- which(strStartsWith(names(fData(jscs)), "log2FC("));
         fc.names <- names(fData(jscs))[fc.cols];
         if(length(fc.names) == 1){
           fc.file <- paste0("ma-plot-",sub("/","vs",fc.names[1], fixed = TRUE), plotting.device.ext);
           writeLines(paste0("<td><a href=\"",fc.file,"\"> ",
                             "<img src=\"",fc.file,"\" "," alt=\"MA Plot\" style=\"width: 45vw;\">",
                             "</a></td>"),pf);
         } else {
           for(j in 1:length(fc.names)){
             fc.file <- paste0("ma-plot-",sub("/","vs",fc.names[j], fixed = TRUE), plotting.device.ext);
             writeLines(paste0("<td><a href=\"",fc.file,"\"> ",
                               "<img src=\"",fc.file,"\" "," alt=\"MA Plot (",fc.names[j],")\" style=\"width: 45vw;\">",
                               "</a></td>"),pf);
             if(j %% 2 == 1 & j < length(fc.names)){
               writeLines("</tr><tr>",pf);
             }
           }
         }
       }

       writeLines(paste0("</tr></table>"),pf);
     }
     
     if(! is.null(gene.list)){
       #hwriter("Selected Genes:", p, heading=2);
       writeLines("<h2>Selected Genes:</h2><br>", pf)
     } else {
       if(all(! f.na(fData(jscs)$padjust < FDR.threshold))){
         #hwriter(paste0("No genes found with significant features (FDR < ",FDR.threshold,")"), p, heading=2);
         writeLines(paste0("<h2>No genes found with significant features (FDR < ",FDR.threshold,")</h2><br>"), pf)
         close(pf, splash=TRUE)
         return(FALSE);
       } else {
         gene.list <- unique(as.character(fData(jscs)$geneID[ f.na(fData(jscs)$padjust < FDR.threshold) ]));
         #hwriter(paste0("Genes with significant features (FDR < ",FDR.threshold,")"), p, heading=2);
         writeLines(paste0("<h2>Genes with significant features (FDR < ",FDR.threshold,")</h2><br>"), pf)
       }
     }
     mainTable <- data.frame(geneID = as.character(gene.list), stringsAsFactors=F);
     geneAnno <- as.data.frame(t(sapply(gene.list, function(g){
       geneRows <- which(fData(jscs)$geneID == g);
       c(as.character(fData(jscs)$chr[geneRows[1]]), 
         as.numeric(min(fData(jscs)$start[geneRows])), 
         as.numeric(max(fData(jscs)$end[geneRows])),
         as.character(fData(jscs)$strand[geneRows[1]])
       );
     })));
     if(is.null(number.plots)){
          number.plots <- rep("",length(gene.list));
     }
     #message("Gene List: [",paste0(gene.list,collapse=","),"]");
     
     #message("1");
     #print(geneAnno);
     colnames(geneAnno) <- c("chr","start","end","strand");
     mainTable <- cbind.data.frame(mainTable, geneAnno);
     mainTable$chr <- as.character(mainTable$chr);
     #message("2");
     geneBaseMeans <- rowMeans(jscs@geneCountData[match(gene.list,rownames(jscs@geneCountData)),] / sizeFactors(jscs));
     mainTable$baseMean <- sprintf("%.1f",geneBaseMeans);

     mainTable$mostSigID <- sapply(gene.list, function(g){
       geneRows <- which(fData(jscs)$geneID == g);
       fData(jscs)$countbinID[ geneRows[which.min( fData(jscs)$padjust[geneRows])] ];
     })
     
     mainTable$mostSigPadjust <- sapply(gene.list, function(g){
       geneRows <- which(fData(jscs)$geneID == g);
       fData(jscs)$padjust[ geneRows[which.min( fData(jscs)$padjust[geneRows])] ];
     })
     mainTable$mostSigPadjust <- sprintf("%.3g",mainTable$mostSigPadjust);
     
     titleLine <- paste0("<tr><td rowspan=2>", paste0( colnames(mainTable), collapse="</td> <td rowspan=2>") )
     
     plotColSpan <- (expr.plot + normCounts.plot + rExpr.plot + rawCounts.plot) * (with.TX + without.TX);
     
     if(plotColSpan > 0){
       titleLine <- paste0(titleLine,"<td colspan=",plotColSpan,">Plots and Tables</td>");
       titleLine <- paste0(titleLine,"</tr><tr>");
     }
     
     if(expr.plot){
       if(without.TX){
         mainTable$expr.plot <- paste0("<a href=\"htmlFiles/",gene.list,"-expr.plot.html\">","PLOT","</a><br>",
                                       "<a href=\"htmlFiles/",gene.list,"-expr.table.html\">","TABLE","</a>");
         titleLine <- paste0(titleLine,"<td>expr</td>")
       }
       if(with.TX){
         mainTable$expr.plot.TX <- paste0("<a href=\"htmlFiles/",gene.list,"-expr-withTX.plot.html\">","PLOT","</a><br>",
                                       "<a href=\"htmlFiles/",gene.list,"-expr.table.html\">","TABLE","</a>");
         titleLine <- paste0(titleLine,"<td>exprTX</td>")
       }
     }
     if(normCounts.plot){
       if(without.TX){
       mainTable$normCt.plot <- paste0("<a href=\"htmlFiles/",gene.list,"-normCt.plot.html\">","PLOT","</a><br>",
                                       "<a href=\"htmlFiles/",gene.list,"-normCt.table.html\">","TABLE","</a>");
       titleLine <- paste0(titleLine,"<td>normCt</td>")
       }
       if(with.TX){
       mainTable$normCt.plot.TX <- paste0("<a href=\"htmlFiles/",gene.list,"-normCt-withTX.plot.html\">","PLOT","</a><br>",
                                       "<a href=\"htmlFiles/",gene.list,"-normCt.table.html\">","TABLE","</a>");
       titleLine <- paste0(titleLine,"<td>normCtTX</td>")
       }
     }
     if(rExpr.plot){
       if(without.TX){
       mainTable$rExpr.plot <- paste0("<a href=\"htmlFiles/",gene.list,"-rExpr.plot.html\">","PLOT","</a><br>",
                                     "<a href=\"htmlFiles/",gene.list,"-rExpr.table.html\">","TABLE","</a>");
       titleLine <- paste0(titleLine,"<td>relExpr</td>")
       }
       if(with.TX){
       mainTable$rExpr.plot.TX <- paste0("<a href=\"htmlFiles/",gene.list,"-rExpr-withTX.plot.html\">","PLOT","</a><br>",
                                     "<a href=\"htmlFiles/",gene.list,"-rExpr.table.html\">","TABLE","</a>");
       titleLine <- paste0(titleLine,"<td>relExprTX</td>")
       }
     }
     if(rawCounts.plot){
       if(without.TX){
       mainTable$rawCounts.plot <- paste0("<a href=\"htmlFiles/",gene.list,"-rawCt.plot.html\">","PLOT","</a><br>",
                                     "<a href=\"htmlFiles/",gene.list,"-rawCt.table.html\">","TABLE","</a>");
       titleLine <- paste0(titleLine,"<td>rawCt</td>")
       }
       if(with.TX){
       mainTable$rawCounts.plot.TX <- paste0("<a href=\"htmlFiles/",gene.list,"-rawCt-withTX.plot.html\">","PLOT","</a><br>",
                                     "<a href=\"htmlFiles/",gene.list,"-rawCt.table.html\">","TABLE","</a>");
       titleLine <- paste0(titleLine,"<td>rawCtTX</td>")
       }
     }
     titleLines <- paste0(titleLine,"</tr>");
     
     
     writeLines("<table BORDER=1>",pf);
     #writeLines(paste0("<tr><td>", paste0( colnames(mainTable), collapse=" </td><td> " ) ,"</td></tr>" ), pf);
     writeLines(titleLine,pf);
     for(i in 1:nrow(mainTable)){
       writeLines("   <tr>", pf);
       #for(j in 1:ncol(mainTable)){
         line <- paste0("<td>", paste0( as.character(mainTable[i,]), collapse=" </td><td> " ) ,"</td> </tr>" );
       #}
       writeLines(line,pf);
       #writeLines("   </tr>", pf);
     }
     writeLines("</table>",pf);
     
     
     #hwriter::hwrite();
     #close(p, splash=TRUE)
   #} else {
   #  message("package hwriter not found. Skipping HTML results pages.");
   #}
   #
   writeLines(paste0("</body> </html>"), pf);
   close(pf);
   
   flat.gff.data <- jscs@flatGffData;

   #TEMP VERSION, needs to be upgraded?

   
   for(i in 1:length(gene.list)){
     if(i == 1){
       prev.g <- NULL;
     } else {
       prev.g <- paste0(gene.list[i - 1]);
     }
     if(i == length(gene.list)){
       next.g <- NULL;  
     } else {
       next.g <- paste0(gene.list[i + 1]);
     }
     g <- gene.list[i];
     
     navTable <- makeNavTable(g, prev.g = prev.g, next.g = next.g, expr.plot=expr.plot, normCounts.plot=normCounts.plot, rExpr.plot=rExpr.plot, rawCounts.plot=rawCounts.plot, with.TX = with.TX, without.TX = without.TX, mainFile = mainFile);
     
     #taken from func buildAllPlotsForGene
     
      transcripts <- sapply(sapply(flat.gff.data$transcripts[which(flat.gff.data$gene_id==g & flat.gff.data$featureType == "exonic_part")],
                                   toString), 
                            function(x){strsplit(x, "+",fixed=TRUE)}
                            );
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
        num.cols <- getColCt(geneID = g, merged.data = fData(jscs), 
                                  plot.exon.results = plot.exon.results, plot.junction.results = plot.junction.results, plot.novel.junction.results=plot.novel.junction.results, 
                                  plot.untestable.results=plot.untestable.results);
        if(num.cols > autoscale.width.to.fit.bins){
          width.multiplier <- num.cols / autoscale.width.to.fit.bins;
        } else {
          width.multiplier <- 1;
        }
      }
      
     #makePlotPage(plotfile = paste0("../expr/",g,"-expr.png"),
     #             htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-expr.plot.html"),
     #             navTable = navTable,
     #             width = floor(base.html.height.px * width.multiplier),
     #             height = floor(base.html.height.px)  );
     
     makeAllPlotPages(g = g, prev.g = prev.g, next.g = next.g,
                      gene.number = number.plots[i],
                      navTable = navTable,
                      outfile.dir = outfile.dir, jscs = jscs,
                      plotting.device.ext = plotting.device.ext,
                      width.multiplier = width.multiplier,
                      withTxPlot.height.multiplier = withTxPlot.height.multiplier,
                      subdirectories.by.type = subdirectories.by.type,
                      mainFile = mainFile,
                      expr.plot = expr.plot, normCounts.plot = normCounts.plot,  rExpr.plot =  rExpr.plot, rawCounts.plot = rawCounts.plot,
                      with.TX=with.TX, without.TX=without.TX,
                      base.html.height = base.html.height, base.html.height.units = base.html.height.units,
                      number.plots = number.plots,
                      FDR.threshold = colorRed.FDR.threshold,
                      css.path = css.path.SUB);
   }
   
   return(TRUE);
}


formatArray <- function( fmt, a){
   apply(a,MARGIN=c(1,2), function(x){
      sprintf(fmt,x);
   });
}


makeAllPlotPages <- function(g, prev.g, next.g, gene.number,
                            navTable,
                            outfile.dir,
                            jscs, 
                            plotting.device.ext,
                            width.multiplier,
                            withTxPlot.height.multiplier,
                            subdirectories.by.type = TRUE,
                            mainFile="testForDU.html", 
                            expr.plot=TRUE,normCounts.plot=TRUE,
                            rExpr.plot=FALSE,rawCounts.plot=FALSE,
                            with.TX=TRUE,without.TX=TRUE,
                            base.html.height = 90, base.html.height.units = "vh",
                            number.plots = FALSE,
                            FDR.threshold = 0.01,
                            css.path = NULL
                              ){
  geneRows <- which(fData(jscs)$geneID == g);    
  td.base <- fData(jscs)[geneRows, c("countbinID","chr","start","end","testable","status","baseMean","dispersion","pvalue","padjust")]
  pval.numeric <- td.base$padjust;
  isTestable <- td.base$testable;
  names(td.base)[1] <- "binID";
  td.base$status <- as.character(td.base$status);
  td.base$chr <- as.character(td.base$chr);
  td.base$baseMean <- sprintf("%.1f",td.base$baseMean);
  td.base$dispersion <- sprintf("%.3g",td.base$dispersion);
  td.base$pvalue <- sprintf("%.4g",td.base$pvalue);
  td.base$padjust <- sprintf("%.4g",td.base$padjust);
  
  LFC.idx <- which(grepl("^log2FC",colnames(fData(jscs))));
    for(i in LFC.idx){
      td.base[[colnames(fData(jscs))[i]]] <- sprintf("%.3f",fData(jscs)[geneRows, i]);
    }
  #if( length(LFC.idx) == 1){
  #  td.base$Log2FC <- fData(jscs)[geneRows, LFC.idx];
  #  td.base$Log2FC <- sprintf("%.1f",td.base$Log2FC);
  #} else {
#
  #}
  
  titleline.base <- colnames(td.base);
  for(i in which(grepl("^log2FC",titleline.base))){
    titleline.base[i] <- sub("\\(","<br>\\(",titleline.base[i]);
  }
  titleline.base <- paste0("<tr><td rowspan=2>",  paste0(titleline.base,collapse="</td><td rowspan=2>"),"</td>");

  
  dirPath <- if(subdirectories.by.type){
    c("../expr/",
      "../exprTX/",
      "../normCounts/",
      "../normCountsTX/",
      "../rExpr/",
      "../rExprTX/",
      "../rawCounts/",
      "../rawCountsTX/");
  } else {
    rep("../",8);
  }
  
  if(expr.plot){
     if(without.TX){
       html.suffix <- "-expr.plot.html";
       makePlotPage(plotfile = paste0(dirPath[1],gene.number,g,"-expr",plotting.device.ext),
                  pageTitle = paste0(g, " expression"),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-expr.plot.html"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height),
                  base.html.height.units = base.html.height.units,
                  css.path = css.path);
     }
     if(with.TX){
       html.suffix <-  "-expr-withTX.plot.html";
       makePlotPage(plotfile = paste0(dirPath[2],gene.number,g,"-expr-withTX",plotting.device.ext),
                  pageTitle = paste0(g, " expression (with TX)"),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-expr-withTX.plot.html"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height * withTxPlot.height.multiplier),
                  base.html.height.units = base.html.height.units,
                  css.path = css.path);
     }


     td <- formatArray("%.1f",jscs@plottingEstimates[["exprEstimate"]][geneRows,])
     colnames(td) <- sub("expr_","", colnames(td),fixed=TRUE)
     titleline <- paste0(titleline.base,
                         "<td colspan=",ncol(td),"> Mean Normalized Coverage Counts </td> </tr> ",
                         "<tr> <td>",paste0(colnames(td),collapse="</td><td>"),"</td></tr>");
     td <- cbind(td.base, td);

     html.suffix <- "-expr.table.html";
     makeTablePage(pvals = pval.numeric, isTestable = isTestable, td = td,
                  titleline = titleline,
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-expr.table.html"),
                  pageTitle = paste0("Feature Coverage/Expression Table (",g,")"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  css.path = css.path,
                  FDR.threshold = FDR.threshold);
  }
  if(normCounts.plot){
     if(without.TX){
     html.suffix <-  "-normCt.plot.html";
     makePlotPage(plotfile = paste0(dirPath[3],gene.number,g,"-normCounts",plotting.device.ext),
                  pageTitle = paste0(g, " norm counts"),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-normCt.plot.html"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height)  ,
                  base.html.height.units = base.html.height.units,
                  css.path = css.path);
     }
     if(with.TX){
     html.suffix <-  "-normCt-withTX.plot.html";
     makePlotPage(plotfile = paste0(dirPath[4],gene.number,g,"-normCounts-withTX",plotting.device.ext),
                  pageTitle = paste0(g, " norm counts"),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-normCt-withTX.plot.html"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height  * withTxPlot.height.multiplier)  ,
                  base.html.height.units = base.html.height.units,
                  css.path = css.path);
     }
     td <- formatArray("%.1f",jscs@plottingEstimates[["normCounts"]][geneRows,])
     td <- td[,1:(ncol(td)/2)];
     colnames(td) <- sub("normCount_","", colnames(td),fixed=TRUE)
     titleline <- paste0(titleline.base,
                         "<td colspan=",ncol(td),"> Normalized Coverage Counts </td> </tr> ",
                         "<tr> <td>",paste0(colnames(td),collapse="</td><td>"),"</td></tr>");
     td <- cbind(td.base, td);
     html.suffix <- "-normCt.table.html";
     makeTablePage(pvals = pval.numeric, isTestable = isTestable, td = td,
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-normCt.table.html"),
                  pageTitle = paste0("Normalized Counts Table (",g,")"),
                  titleline = titleline,
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  css.path = css.path,
                  FDR.threshold = FDR.threshold);
  }
  if(rExpr.plot ){
     if(without.TX){
     html.suffix <-  "-rExpr.plot.html";
     makePlotPage(plotfile = paste0(dirPath[5],gene.number,g,"-rExpr",plotting.device.ext),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rExpr.plot.html"),
                  pageTitle = paste0(g, " relative expression"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height) ,
                  base.html.height.units = base.html.height.units,
                  css.path = css.path );
     }
     if(with.TX){
     html.suffix <-  "-rExpr-withTX.plot.html";
     makePlotPage(plotfile = paste0(dirPath[6],gene.number,g,"-rExpr-withTX",plotting.device.ext),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rExpr-withTX.plot.html"),
                  pageTitle = paste0(g, " relative expression (with TX)"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height * withTxPlot.height.multiplier) ,
                  base.html.height.units = base.html.height.units,
                  css.path = css.path );
     }
     td <- formatArray("%.1f",jscs@plottingEstimates[["relExprEstimate"]][geneRows,])
     colnames(td) <- sub("relExpr_","", colnames(td),fixed=TRUE)
     titleline <- paste0(titleline.base,
                         "<td colspan=",ncol(td),"> Relative Coverage </td> </tr> ",
                         "<tr> <td>",paste0(colnames(td),collapse="</td><td>"),"</td></tr>");
     td <- cbind(td.base, td);
     html.suffix <- "-rExpr.table.html";
     makeTablePage(pvals = pval.numeric, isTestable = isTestable, td = td,
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rExpr.table.html"),
                  pageTitle = paste0("Relative Coverage Table, relative to gene-wide expression (",g,")"),
                  titleline = titleline,
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  css.path = css.path,
                  FDR.threshold = FDR.threshold);
  }
  if(rawCounts.plot ){
     
     if(without.TX){
     html.suffix <- "-rawCt.plot.html";
     makePlotPage(plotfile = paste0(dirPath[7],gene.number,g,"-rawCounts",plotting.device.ext),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rawCt.plot.html"),
                  pageTitle = paste0(g, " Raw Counts"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height) ,
                  base.html.height.units = base.html.height.units,
                  css.path = css.path );
     }
     if(with.TX){
     html.suffix <- "-rawCt-withTX.plot.html";
     makePlotPage(plotfile = paste0(dirPath[8],gene.number,g,"-rawCounts-withTX",plotting.device.ext),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rawCt-withTX.plot.html"),
                  pageTitle = paste0(g, " Raw Counts (with TX)"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height * withTxPlot.height.multiplier) ,
                  base.html.height.units = base.html.height.units,
                  css.path = css.path );
     }
     td <- formatArray("%.1f",counts(jscs)[geneRows,])
     titleline <- paste0(titleline.base,
                         "<td colspan=",ncol(td),"> Raw Counts </td> </tr> ",
                         "<tr> <td>",paste0(colnames(td),collapse="</td><td>"),"</td></tr>");
     td <- cbind(td.base, td);
     html.suffix <- "-rawCt.table.html";
     makeTablePage(pvals = pval.numeric, isTestable = isTestable, td = td,
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rawCt.table.html"),
                  pageTitle = paste0("Raw Counts Table (",g,")"),
                  titleline = titleline,
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g,mainFile,prefix="./",suffix=html.suffix),
                  css.path = css.path,
                  FDR.threshold = FDR.threshold);     
  }

}

makeGeneNavTable <- function(g, prev.g, next.g, mainFile, prefix, suffix){
  prev.link <- if(is.null(prev.g)){
    "&lt&lt Prev Gene"
  } else {
    paste0("<a href=\"",prefix,prev.g,suffix,"\"> &lt&lt Prev Gene </a>")
  }
  next.link <- if(is.null(next.g)){
    "Next Gene &gt&gt"
  } else {
    paste0("<a href=\"",prefix,next.g,suffix,"\"> Next Gene &gt&gt </a>")
  }
  out <- paste0(
          "<table style=\"","margin-left:auto; margin-right:auto;","\">",
            "<tr>",
              "<td  style=\"font-size:x-large;\">",
                 prev.link,
              "</td>",
              "<td style=\"font-size:x-large;\" COLSPAN=2>",
                "<a href=\"../",mainFile,"\"> TOP </a>",
              "</td>",
              "<td  style=\"font-size:x-large;\">",
                 next.link,
              "</td>",
            "</tr>",
          "</table>"
  );
  return(out);
  #out <- paste0(
  #        "<table style=\"","margin-left:auto; margin-right:auto;","\">",
  #          "<tr>",
  #            "<td style=\"font-size:xx-large;\" COLSPAN=2>",
  #              "<a href=\"../",mainFile,"\"> TOP </a>",
  #            "</td>",
  #          "</tr>",
  #          "<tr>",
  #            "<td  style=\"font-size:xx-large;\">",
  #               prev.link,
  #            "</td>",
  #            "<td  style=\"font-size:xx-large;\">",
  #               next.link,
  #            "</td>",
  #          "</tr>",
  #        "</table> <br>"
  #);
}

makeNavTable <- function(g, prev.g, next.g, expr.plot, normCounts.plot, rExpr.plot, rawCounts.plot, with.TX, without.TX, mainFile){
     
     navTable.1 <- ""; #paste0("<b><td rowspan = 3> <a href=\"../",mainFile,"\">BACK</a> </td>");
     navTable.2 <- "";
     navTable.3 <- "";
     
     nav.colSpan <- 0;
     if(expr.plot){
       navTable.1 <- paste0(navTable.1, "<td>expr</td>");
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-expr.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-expr.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1;
     }
     if(expr.plot){
       navTable.1 <- paste0(navTable.1, "<td>exprTX</td>");
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-expr-withTX.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-expr.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1;
     }
     
     if(normCounts.plot){
       navTable.1 <- paste0(navTable.1, "<td>normCt</td>");
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-normCt.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-normCt.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1;
     }
     if(normCounts.plot){
       navTable.1 <- paste0(navTable.1, "<td>normCtTX</td>");
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-normCt-withTX.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-normCt.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1;
     }
     
     if(rExpr.plot){
       navTable.1 <- paste0(navTable.1, "<td>rExpr</td>");
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-rExpr.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-rExpr.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1;
     }
     if(rExpr.plot){
       navTable.1 <- paste0(navTable.1, "<td>rExprTX</td>");
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-rExpr-withTX.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-rExpr.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1;
     }
     
     if(rawCounts.plot){
       navTable.1 <- paste0(navTable.1, "<td>rawCt</td>");
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-rawCt.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-rawCt.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1;
     }
     if(rawCounts.plot){
       navTable.1 <- paste0(navTable.1, "<td>rawCtTX</td>");
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-rawCt-withTX.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-rawCt.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1;
     }
     
     #navTable.1 <- paste0(navTable.1,"<td rowspan=3><a href=\"",g,"-results.html","\"> RESULTS </a> </td></b>");
     navTable <- c("<table style=\"","margin-left:auto;margin-right:auto;","\">",
                       "<th COLSPAN=",nav.colSpan,">",
                          "Gene ",g,
                       "</th>",
                       paste0("<tr>",c(navTable.1,navTable.2,navTable.3),"</tr>"),
                   "</table>");
     
     #table2   <- paste0("<table>",
     #                      "<tr> <td COLSPAN=2>",
     #                         "<a href=\"../",mainFile,"\">> TOP </a>",
     #                      "</td> </tr>",
     #                      "<tr>",
     #                         "<td>",
     #                           "<a href=\">",,"</a>",
     #                         "</td>",
     #                      "</tr>",
     #                   "</table>");
     
     return(navTable);
}

makePlotPage <- function(plotfile, pageTitle = "", htmlfile, navTable, geneNavTable, width, height, css.path = NULL, base.html.height.units = "vh"){
  f <- file(htmlfile, "w");
  
     writeLines(paste0("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">"), f);
     writeLines(paste0("<html><head><title> JunctionSeq Results (",pageTitle,") </title>"), f);
     writeLines(paste0("<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">"), f);
     if(! is.null(css.path)){
       writeLines(paste0("<link rel=\"stylesheet\" type=\"text/css\" href=\"",css.path,"\">"), f);
     }
     writeLines(paste0("</head> <body>"), f);
  
  writeLines(geneNavTable,f);
  #writeLines(paste0("<h1>",pageTitle,"</h1>"),f);
  #writeLines("<br>",f);
  writeLines(
            paste0("<a href=\"",plotfile,"\"><img src=\"",plotfile,"\" "," alt=\"JunctionSeq Plot\" style=\"display:block; margin-left:auto; margin-right:auto; height:",height,base.html.height.units,";\"> </a>"),f);
  
  writeLines("<br>",f);
  writeLines(navTable,f);
  writeLines(paste0("</body></html>"), f);
  close(f);
}

makeStyleKeyTable <- function(pval.thresholds){
  return(
    c("<table class=\"dataTable\">",
      paste0(" <tr><th COLSPAN=",length(pval.thresholds) + 2,"> Color Codes </th></tr>"),
      " <tr>",
      "   <td class=\"untestable\">",
      "      NO TEST",
      "   </td>",
      "   <td class=\"sigNo\">",
      paste0("      p &gt ",pval.thresholds[1],""),
      "   </td>",
      "   <td class=\"sig1\">",
      paste0("      p &lt ",pval.thresholds[1],""),
      "   </td>",
      "   <td class=\"sig2\">",
      paste0("      p &lt ",pval.thresholds[2],""),
      "   </td>",
      "   <td class=\"sig3\">",
      paste0("      p &lt ",pval.thresholds[3],""),
      "   </td>",
      "   <td class=\"sig4\">",
      paste0("      p &lt ",pval.thresholds[4],""),
      "   </td>",
      " </tr>",
      "</table>")
  );
}

makeTablePage <- function(pvals, isTestable, td, pageTitle, htmlfile, navTable, geneNavTable, titleline, css.path = NULL,
                          FDR.threshold = 0.01, pval.thresholds = NULL){
  
  if(is.null(pval.thresholds)){
    pval.thresholds <- c(FDR.threshold, 
                         FDR.threshold / 10, 
                         FDR.threshold / 100, 
                         FDR.threshold / 1000);
  }
  
  f <- file(htmlfile, "w");
  
     writeLines(paste0("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">"), f);
     writeLines(paste0("<html><head><title> JunctionSeq Results (",pageTitle,") </title>"), f);
     writeLines(paste0("<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">"), f);
     if(! is.null(css.path)){
       writeLines(paste0("<link rel=\"stylesheet\" type=\"text/css\" href=\"",css.path,"\">"), f);
     }
     writeLines(paste0("</head> <body>"), f);
  
  
  writeLines(geneNavTable,f);
  writeLines(paste0("<h1>",pageTitle,"</h1>"),f);
  writeLines("<br>",f);
  
  writeLines("<table border=1 class=\"dataTable\">",f);
  #writeLines(paste0("<tr><td>", paste0(colnames(td),collapse="</td><td>") ,"</td></tr>"),f);
  #writeLines("<div class=\"tableNumbers\">",f);
  writeLines(titleline,f);
  for(i in 1:nrow(td)){
    
    if(! isTestable[i]){
      classText <- paste0("untestable");
    } else if(pvals[i] < pval.thresholds[4]){
      classText <- paste0("sig4");
    } else if(pvals[i] < pval.thresholds[3]){
      classText <- paste0("sig3");
    } else if(pvals[i] < pval.thresholds[2]){
      classText <- paste0("sig2");
    } else if(pvals[i] < pval.thresholds[1]){
      classText <- paste0("sig1");
    } else {
      classText <- paste0("sigNo");
    }
    
    #writeLines(classText,f);
    #rowLine <- as.character(td[i,]);
    
    rowLine <- if(td[i, colnames(td) == "status"] != "OK"){
      RL <- sapply(1:length(td[i,]),function(j){
        if(colnames(td)[j] != "status"){
          paste0("<td>",as.character(td[i,j]),"</td>");
        } else {
          paste0("<td style=\"font-size:xx-small\">",as.character(td[i,j]),"</td>");
        }
      });
      paste0(RL,collapse=" ");
    } else {
      paste0("<td>", paste0(as.character(td[i,]),collapse="</td><td>"),"</td>");
    }
    
    writeLines(paste0("<tr class=\"",classText,"\">",rowLine,"</tr></div>"),f);
    #writeLines("</div>",f);
  }
  #writeLines("</div>",f);
  writeLines("</table>",f);
  
  writeLines("<br>",f);
  writeLines(makeStyleKeyTable(pval.thresholds),f);
  
  writeLines("<br>",f);
  writeLines(navTable,f);
  
  writeLines(paste0("</body></html>"), f);
  close(f);
}

#######################################################################
#DEXSEQ.for.reference.makePlotPage <- function(plotFileName, ptowrite, gene, whichtag, links, width, height){
#   allopts <- c("expression", "splicing", "counts", "transcripts")
#   opts <- allopts %in% whichtag
#   onlytag <- allopts[max(which(opts))]
#   pagename <- sapply(strsplit(as.character(gene), "\\+"), "[[", 1)
#   genpage <- openPage(paste(ptowrite, pagename, onlytag, ".html", sep=""))
#   hwrite(links, table=TRUE, border=0, genpage)
#   
#   hwrite(hmakeTag("iframe", src=plotFileName, width=width, height=height, border=0), page=genpage )
#   close(genpage, splash=TRUE)
#}
#######################################################################
#backupDEXSeqHTML <- function(){   
#   stopifnot( is( object, "DEXSeqResults" ) )
#   ######## GET THE RESULT TABLE READY ##########
#   genomicData <- as.data.frame( object$genomicData )
#   results <- data.frame( object[, c( "groupID", "featureID", "exonBaseMean", "dispersion", "pvalue", "padj" )], stringsAsFactors=TRUE)
#   results <- cbind( results, genomicData )
#   
#   results[,c("dispersion", "pvalue", "padj")] <- round(results[,c("dispersion", "pvalue", "padj")], 3)
#
#   dexseqR <- elementMetadata( object )$type == "DEXSeq results"
#
#   if(sum(dexseqR, na.rm=TRUE) > 0){
#      results <-
#          cbind(
#              results,
#              round(
#                  as.data.frame( object[,which(dexseqR)] ), 3
#                  )
#              )
#   }
#
#   rownames(results) <- NULL
#   sampleData <- object@sampleData
#   
#   numcond <- length(unique(sampleData[[fitExpToVar]]))
#   if(is.null(color)){
#      color<-rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond)), maxColorValue=255, alpha=175)
#   }
#   names(color) <- sort(levels(sampleData[[fitExpToVar]]))
#   dir.create(file.path(path, "files"), recursive=TRUE)
#
#   if(is.null(genes)){
#      gns <- as.character(unique(results$groupID[which(results$padj < FDR)]))
#   }else{
#      gns <- genes
#   }
#
#   
#   if(!all(gns %in% object$groupID)){
#      stop("The geneIDs provided are not in the ecs object")}
#   if(length(gns)==0){ 
#      stop("There are no significant results in the test... nothing to report")}
#
#   p<-openPage(file.path(path, file))
#   hwrite('DEXSeq differential exon usage test', p, heading=1)
#   hwrite('Experimental design', p, heading=2)
#   cond<-as.matrix( as.data.frame( sampleData[,!colnames(sampleData) %in% "sizeFactor"] ) )
#   rownames(cond) <- NULL
#   condcolor <- matrix(rep("white", nrow(cond)*ncol(cond)), nrow(cond))
#   condcolor[,which(colnames(cond) %in% fitExpToVar)] <- color[as.character(sampleData[[fitExpToVar]])]
#   if(!is.null(color.samples)){
#      condcolor[,1] <- color.samples}
#   hwrite(cond, bgcolor=condcolor, p)
#
#   formulas <- elementMetadata(object)[colnames(object) == "pvalue","description"]
#   formulas <- sapply( strsplit(formulas, "vs|p-value:" ), "[", c(2, 3))
#   formulas <- as.vector( gsub("'", "", formulas) )
#   
#   hwrite(paste("\n\nformulaDispersion = ", formulas[1], sep=""), p, heading=3)
#   hwrite(paste("\nformula0 = ", formulas[2], sep=""), p, heading=3)
#   hwrite(paste("\nformula1 = ", formulas[1], sep=""), p, heading=3)
#   hwrite('testForDEU result table', p, heading=2)
#   ptowrite <- file.path(path, "files/")
#   ######### prepare colors for table results
#  
#   m2col <- colorRampPalette(c("#FF706B", "#FEE08B", "white"), space="rgb")(5)
#   matcol <- matrix(rep(m2col[5], ncol(results)*nrow(results)), nrow(results))
#   ######### COLOR LEGEND = I BET THERE IS A SMARTER WAY OF DOING THIS:
#   j <- 4
#   for(i in c(0.25, 0.1, 0.05, 0.01)){
#      matcol[which(results$padjust <= i),] <- m2col[j]
#      j <- j-1
#   }
#   legend <- hwrite(c("<= 0.01", "<= 0.05", "<= 0.1", "<= 0.25", "> 0.25"), bgcolor=m2col)
#
#  makePagesForGene <- function(gene){
##      print( gene )
#      back <- hwrite("back", link=file.path("..", file))
#      nameforlinks <- sapply(strsplit(gene, "\\+"), "[[", 1)
#      otherlinks <- hwrite(c("counts", "expression", "splicing", "transcripts", "results"), link=c(paste(nameforlinks, "counts.html", sep=""), paste(nameforlinks, "expression.html", sep=""), paste(nameforlinks, "splicing.html", sep=""), paste(nameforlinks, "transcripts.html", sep=""), paste(nameforlinks, "results.html", sep="")), table=FALSE)
#      loc <- as.character(results$groupID) %in% as.character(gene)
#      ### this makes the page where to explore the pvalues ###
#      subres <- results[loc,]
#      submatcol <- matcol[loc,]
#      rownames(subres) <- NULL
#      genpage <- openPage(paste(ptowrite, nameforlinks, "results.html", sep=""))
#      hwrite(c(back, otherlinks), table=TRUE, border=0, genpage)
#      hwrite(legend, page=genpage)
#      hwrite(as.matrix(subres), bgcolor=submatcol, table.class="sortable", style='margin:16px; border:0px solid black; border-width:0px; width:200px', table=TRUE, page=genpage)
#      close(genpage, splash=TRUE)
#      ### MAKE THE PLOT HTML PAGES FOR expression, counts and splicing
#      makePlotPage( object, ptowrite=ptowrite, gene=gene, whichtag="expression", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=700, h=7)
#      makePlotPage( object, ptowrite=ptowrite, gene=gene, whichtag="counts", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=700, h=7)
#      makePlotPage( object, ptowrite=ptowrite, gene=gene, whichtag="splicing", links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=700, h=7)
#
#      transcripts <- object$transcripts[loc]
#      if(length( unlist( transcripts ) ) > 0 ) {
#         trans <- Reduce(union, transcripts)
#         h <- ifelse(length(trans) > 10, 7+(length(trans)*0.3), 7)
#      if(sum(loc) > 30){   ############# if there are more than 30 exons, increase the size of the plotting region
#         h <- h + (sum(loc)*0.1)
#      }
#      makePlotPage( object, ptowrite=ptowrite, gene=gene, whichtag=c("expression", "transcripts"), links=c(back, otherlinks), color=color, color.samples=color.samples, FDR=FDR, fitExpToVar=fitExpToVar, width=1200, height=h*100, h=h)
#      }
#      return()
#   }
#
#
#   bplapply( gns, makePagesForGene, BPPARAM=BPPARAM )
#   
#   results <- results[as.character(results$groupID) %in% gns,]
#
#   splitCols <- split( seq_len(nrow( results ) ), results$groupID )
#   
#   genetable <- lapply( splitCols, function(x){
#       data.frame(
#           chr=unique( results$seqnames[x] ),
#           start=min( results$start[x] ),
#           end=max( results$end[x] ),
#           total_exons = length(x),
#           exon_changes = sum( results$padj[x] < FDR, na.rm=TRUE) )
#   })
#   genetable <- do.call(rbind, genetable)
#   genetable <- cbind( geneID=rownames(genetable), genetable )
#
#   if(class(mart) == "Mart"){
#      if(attributes(mart)$dataset != ""){
#      forvalues <- strsplit(as.character(genetable$geneID), "\\+")
#      names(forvalues) <- genetable$geneID
#      if(length(filter) > 1){
#         warning("length(filter) > 2, only first element will be taken")
#         filter <- filter[1]
#      }
#      extra <- getBM(attributes=c(filter, attributes), filters=filter, values=forvalues, mart=mart)
#      fromart <- lapply(genetable$geneID, function(x){
#         sep <- do.call(c, strsplit(as.character(x), "\\+"))
#         extra[which(extra[,filter] %in% sep),]
#      })
#
#      extra <- sapply(attributes, 
#         function(r){
#           unlist(
#              lapply(fromart, 
#                 function(x){
#                    paste(x[,r], collapse=" ")
#                  }
#              )
#           )
#         }
#      )
#      genetable <- cbind(geneID=genetable$geneID, extra, genetable[,2:length(genetable)])
#      }else{
#         warning("No dataset in biomart specified")
#      }
#   }else if( mart != ""){
#      warning("Please provide a Mart class object for parameter mart")
#   }  
#
#   if( !is.null( extraCols )){
#      genetable <- cbind(extraCols[match(genetable$geneID, rownames(extraCols)),], genetable)
#   }
#   
#  
#   genetable$geneID <- sapply(as.character(genetable$geneID), function(m){w <- strsplit(m, "\\+");ns <- sapply(w, "[[", 1);hwrite(paste(unlist(w), collapse=" "), link=paste("files/", ns, "expression.html", sep=""))})
#   rownames(genetable) <- NULL
#   hwrite(genetable, page=p, table=TRUE, table.class="table-layout:fixed", style='margin:16px; border:0px solid black; border-width:1px; width:20%')
#   close(p, splash=TRUE)
#}



########################################################################################

#This trivial function fixes a minor bug in the "svg" device found in certain older versions of redhat linux.
fixRHEL5_svg_bug <- function(infile,outfile){
  inlines <- readLines(infile);
  outlines <- sub("symbol id = \"glyph", "symbol overflow=\"visible\" id=\"glyph",inlines, fixed=TRUE);
  writeLines(outlines,outfile);
}
#s/symbol id=\"glyph/symbol overflow=\"visible\" id=\"glyph/;

#testing:
# expr.data <- read.table("Ctrl_DvN.junctSeq.resultssigGenes.results.txt.gz",header=T,stringsAsFactors=F);
# write.sig.bed.file(file = "test.3.bed.gz", jscs = jscs);
# write.expr.bed.file(file = "test.4.bed.gz", jscs=jscs);

#file = "test.2.bed.gz";

writeExprBedTrack <- function(file, jscs, 
                                      trackLine = "track name='JctExpr' description='Junction Coverage Estimates, by group' itemRgb='On' visibility=3",
                                      only.with.sig.gene = TRUE,
                                      only.sig = FALSE, 
                                      only.testable = TRUE, 
                                      plot.exons = TRUE, plot.junctions = TRUE,
                                      plot.novel.junctions = TRUE,
                                      group.RGB = NULL, 
                                      #nongroup.RGB = NULL, 
                                      use.score = FALSE, 
                                      FDR.threshold = 0.05, 
                                      count.digits = 1, 
                                      includeGeneID = FALSE,
                                      includeLocusID = TRUE,
                                      includeGroupID = TRUE,
                                      output.format = c("BED","GTF","GFF3"),
                                      verbose = TRUE){
  
  estimates <- jscs@plottingEstimates[["exprEstimate"]];
  estimate.names <- sapply(colnames(estimates), function(s){substring(s,6,nchar(s))});
  numcond <- length(estimate.names);
  
  if(is.null(group.RGB)){
    default.color.list <- c("red","blue","orange","green3","purple","cyan1", "magenta","yellow3","tan4")
    if(numcond > length(default.color.list)){
      warning(paste0("WARNING: more than ",default.color.list," possible values for condition variable! Cannot select contrasting colors! Recommending setting colors manually using the group.RGB variable!"));
      color<- colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond))
    } else {
      color <- t(col2rgb(default.color.list[1:numcond]));
    }
    group.RGB <- sapply(1:nrow(color), function(i){
      paste0(floor(color[i,]),collapse=",");
    })
  }
  
  keep.features <- rep(TRUE,length(nrow(estimates)));
  res.data <- fData(jscs);
  
  #filtering:
  if(only.testable){
    keep.features <- keep.features & res.data$testable;
  }
  sig.features <- which(fData(jscs)$padjust < FDR.threshold);
  gene.list <- unique(as.character(fData(jscs)$geneID[sig.features]));
  if(only.with.sig.gene){
    keep.features <- keep.features & fData(jscs)$geneID %in% gene.list;
  }
  if(only.sig){
    keep.features <- keep.features & res.data$padjust < FDR.threshold;
  }
  if(! plot.exons){
    keep.features <- keep.features & res.data$featureType != "exonic_part";
  }
  if(! plot.junctions){
    keep.features <- keep.features & (res.data$featureType != "splice_site" & res.data$featureType != "novel_splice_site");
  }
  if(! plot.novel.junctions){
    keep.features <- keep.features & res.data$featureType != "novel_splice_site";
  }
  
  res.data <- res.data[keep.features,];
  estimates <- estimates[keep.features,];
  
  chrom <- res.data$chr;
  chromStart <- res.data$start;
  chromEnd <- res.data$end;
  strand <- res.data$strand;
  featureType <- res.data$featureType;
  
  if(use.score){
    score <- sapply(res.data$padjust, function(p){
      if(p > 0.05){
        0;
      } else {
        floor((0.05 - p) * 18000) + 100;
      }
    });
  } else {
    score <- rep(1000,length(chrom));
  }
  
  featureID <- rep("",nrow(estimates));
  if(includeGeneID & includeLocusID){
    featureID <- paste0("", res.data$featureID,"_");
  } else if(includeGeneID){
    featureID <- paste0("", res.data$geneID,"_");
  } else if(includeLocusID){
    featureID <- paste0("", res.data$countbinID,"");
  }
  
  idx <- order( rep(1:nrow(estimates), ncol(estimates)) );
  #counts <- unlist(lapply(1:ncol(estimates), function(i){ estimates[,i] }))[idx];
  
  if(includeGroupID){
    out.featureNames <- unlist(lapply(1:ncol(estimates), function(i){
      paste0(estimate.names[i],":", featureID, "(", 
         sprintf(paste0("%.",count.digits,"f"), estimates[,i])
         ,")");
    }))[idx];
  } else {
    out.featureNames <- unlist(lapply(1:ncol(estimates), function(i){
      paste0(featureID, "(", 
         sprintf(paste0("%.",count.digits,"f"), estimates[,i]) ,
         ")");
    }))[idx];
  }
  
  #message("(1) Starting write.junction.bed.file: ",file);
  write.junction.bed.file(
    file = file,
    trackLine = trackLine,
    chrom = rep(chrom, each=numcond),
    chromStart = rep(chromStart,each=numcond),
    chromEnd = rep(chromEnd, each = numcond),
    strand = rep(strand, each = numcond),
    featureName = out.featureNames,
    featureScore = rep(score,each = numcond),
    featureRGB = rep(group.RGB, nrow(estimates)),
    featureType = rep(featureType, each = numcond),
    output.format = output.format,
    verbose = verbose
  );
}

writeSigBedTrack <- function(file, 
                               jscs, 
                               trackLine = "track name='sigJct' description='Significant Splice Junction Loci' useScore=1 visibility=3",
                               only.sig = TRUE, only.testable = TRUE, 
                               plot.exons = TRUE, plot.junctions = TRUE,
                               plot.novel.junctions = TRUE,
                               sig.RGB = "255,0,0", nonsig.RGB = "0,0,0", 
                               use.score = TRUE, FDR.threshold = 0.05, 
                               pval.digits = 4, 
                               includeGeneID = FALSE,
                               includeLocusID = TRUE,
                               output.format = c("BED","GTF","GFF3"),
                               verbose = TRUE){ #fullFeatureID includeFeatureID
  res.data <- fData(jscs);
  
  keep.features <- rep(TRUE, nrow(res.data));
  if(only.testable){
    keep.features <- keep.features & res.data$testable;
  }
  if(only.sig){
    keep.features <- keep.features & res.data$padjust < FDR.threshold;
  }
  if(! plot.exons){
    keep.features <- keep.features & res.data$featureType != "exonic_part";
  }
  if(! plot.junctions){
    keep.features <- keep.features & (res.data$featureType != "splice_site" & res.data$featureType != "novel_splice_site");
  }
  if(! plot.novel.junctions){
    keep.features <- keep.features & res.data$featureType != "novel_splice_site";
  }
  res.data <- res.data[keep.features,,drop=FALSE];
  
  chrom <- res.data$chr;
  chromStart <- res.data$start;
  chromEnd <- res.data$end;
  strand <- res.data$strand;
  featureType <- res.data$featureType;
  
  if(use.score){
    score <- sapply(res.data$padjust, function(p){
      if(p > 0.05){
        0;
      } else {
        floor((0.05 - p) * 18000) + 100;
      }
    });
  } else {
    score <- rep(1000,length(chrom));
  }
  
  if(only.sig){
    featureRGB <- rep(sig.RGB, length(chrom));
  } else {
    featureRGB <- ifelse(res.data$padjust < FDR.threshold, sig.RGB, nonsig.RGB);
  }
  
  if(includeLocusID & includeGeneID){
    featureID <- res.data$featureID;
  } else if(includeGeneID){
    featureID <- res.data$geneID;
  } else if(includeLocusID) {
    featureID <- res.data$countbinID;
  } else {
    featureID <- rep("",nrow(res.data))
  }
  
  featureName <- paste0(
    featureID,
    "(",sprintf(paste0("%.",pval.digits,"f"), res.data$padjust),")"
  );
  
  
  #message("(2) Starting write.junction.bed.file: ",file);
  write.junction.bed.file(
    file = file,
    trackLine = trackLine,
    chrom = chrom,
    chromStart = chromStart,
    chromEnd = chromEnd,
    strand = strand,
    featureName = featureName,
    featureScore = score,
    featureRGB = featureRGB,
    featureType = featureType,
    output.format = output.format,
    verbose = verbose);
}

write.junction.bed.file <- function(file, trackLine = NULL, chrom, chromStart, chromEnd, strand, featureName, featureScore, featureRGB, featureType,
                                    output.format = c("BED","GTF","GFF3"),
                                    verbose = TRUE){
  
  #message("write.junction.bed.file: ",file);
  #message("is.null(featureType): ",is.null(featureType));
  output.format <- match.arg(output.format);
  
  if(output.format == "BED"){
    out.bed <- data.frame(
      chrom = chrom,
      chromStart = ifelse(featureType == "exonic_part" , chromStart,chromStart - 1),
      chromEnd = ifelse(featureType == "exonic_part" , chromEnd,chromEnd + 1),
      featureName = featureName,
      featureScore = featureScore,
      strand = strand,
      thickStart = ifelse(featureType == "exonic_part" , chromStart,chromStart - 1),
      thickEnd = ifelse(featureType == "exonic_part" , chromEnd,chromEnd + 1),
      itemRGB = featureRGB,
      blockCount  = ifelse(featureType == "exonic_part", 1, 2),
      blockSizes  = ifelse(featureType == "exonic_part", as.character(chromEnd - chromStart), "1,1"),
      blockStarts = ifelse(featureType == "exonic_part", "0",  paste0("0,",chromEnd - (chromStart - 1)))
      #blockCount = rep(2,length(chrom)),
      #blockSizes = rep("1,1",length(chrom)),
      #blockStarts = paste0("0,",chromEnd - (chromStart - 1))
    );

    gzf <- gzfile(paste0(file,""),"w");
    if(! is.null(trackLine)){
      write.table(trackLine,gzf, quote=F,col.names=F,row.names=F);
    }

    write.table(out.bed, gzf, quote=F, col.names=F,row.names=F,sep='\t');
    close(gzf);
  } else {
    write.junction.bed.file.GTFCONVERT(file = file, trackLine = trackLine, output.format = output.format, 
                                       chrom = chrom, chromStart = chromStart, chromEnd = chromEnd, strand = strand, featureName = featureName, featureRGB = featureRGB, featureType = featureType,
                                       verbose = verbose);
  }
}


write.junction.bed.file.GTFCONVERT <- function(file, trackLine = NULL, 
                                               output.format = c("GTF","GFF3"),
                                               chrom, chromStart, chromEnd, strand, featureName, featureScore, featureRGB, featureType,
                                               verbose = TRUE){
  
  output.format <- match.arg(output.format);
  #message("write.junction.bed.file: ",file);
  #message("is.null(featureType): ",is.null(featureType));
  in.data <- data.frame(chrom = chrom, chromStart = chromStart, chromEnd = chromEnd, strand = strand, 
                        featureName=featureName, featureScore = featureScore, featureRGB = featureRGB, featureType = featureType);
  
  
  out.GTF <- do.call(rbind.data.frame, lapply( 1:nrow(in.data), function(i){
    group <- if(output.format == "GTF"){
      paste0("gene_id ",featureName[i],"; transcript_id ",featureName[i]);
    } else if(output.format == "GFF3"){
      paste0("ID=",featureName[i],";Name=",featureName[i]);
    }
    if(featureType[i] == "exonic_part"){
      data.frame(chrom = chrom[i], source = "JunctionSeq", feature = "exon", 
                 start = chromStart[i]+1, end = chromEnd[i], score = featureScore[i], strand = strand[i],
                 frame = ".", group = group);
    } else {
      data.frame(chrom = rep(chrom[i],2), source = rep("JunctionSeq",2), feature = rep("exon",2), 
                 start = c(chromStart[i], chromEnd[i]+1), end = c(chromStart[i],chromEnd[i]+1), score = rep(featureScore[i],2), strand = rep(strand[i],2),
                 frame = rep(".",2), group = rep(group,2));
    }
  }))
  
  gzf <- gzfile(paste0(file,""),"w");
  if(! is.null(trackLine)){
    write.table(trackLine,gzf, quote=F,col.names=F,row.names=F);
  }
  write.table(out.GTF, gzf, quote=F, col.names=F,row.names=F,sep='\t');
  close(gzf);
}

write.simple.table.gz <- function(d, file, use.gzip = TRUE, ...){
   if(use.gzip){
      gzf <- gzfile(paste0(file,".gz"),"w");
      write.table(d, gzf, ...);
      close(gzf);
   } else {
      write.table(d, file, ...);
   }
}



write.table.gz <- function(write.data, file, use.gzip = TRUE, sep = "	", quote=FALSE, row.names = FALSE, ...){
   if(use.gzip){
      gzf <- gzfile(paste0(file,".gz"),"w");
      write.fmt.table(write.data, gzf, row.names = row.names, quote = quote, sep = sep, ...);
      close(gzf);
   } else {
      write.fmt.table(write.data, file, row.names = row.names, quote = quote, sep = sep, ...);
   }
}

write.fmt.table <- function(d, file, row.names = FALSE, quote = FALSE, sep = sep, ...){
  if(row.names == FALSE){
    out.data <- d;
  } else {
    if(row.names == TRUE){
      row.names.name <- "ROWNAME";
    } else {
      row.names.name <- row.names;
    }
    
    out.data <- cbind.data.frame( row.names(d), d );
    
    names(out.data) <- c(row.names.name, names(d));
  }
  
  write.table(out.data, file = file, row.names = FALSE, quote = quote, sep = sep, ...);
}
