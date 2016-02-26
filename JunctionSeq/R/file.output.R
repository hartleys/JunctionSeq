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
                            html.fixed.dim = c("height","autofit","width"),
                            base.html.height = 90, base.html.height.units = "vh",
                            html.width = 90, html.width.units = "vw",
                            GENE.annotation.relative.height = 0.15, TX.annotation.relative.height = 0.05, CONNECTIONS.relative.height = 0.1, SPLICE.annotation.relative.height = 0.1, TX.margins = c(0,0.5),
                            autoscale.height.to.fit.TX.annotation = TRUE,
                            autoscale.width.to.fit.bins = 35,
                            number.plots = NULL,
                            css.file = NULL, css.link = NULL,
                            compare.analysis.list = NULL,
                            minimalImageFilenames = TRUE,
                            INTERNAL.VARS = INTERNAL.VARS,
                            verbose = TRUE, debug.mode = FALSE){
   
   html.fixed.dim <- match.arg(html.fixed.dim)
   
   if(is.null(css.link)){
     if(is.null(css.file)){
       if(verbose) message("   Copying default css stylesheet.")
       default.css.filepath <- system.file("extdata/styles.css", 
                               package="JunctionSeq",
                               mustWork=TRUE)
       cssLines <- readLines(default.css.filepath)
       cssOut <- file(paste0(outfile.dir,"/styles.css"), "w")
       for(line in cssLines){
         writeLines(line, cssOut)
       }
       close(cssOut)
     } else {
       if(verbose) message("   Copying css stylesheet: \"",css.file,"\"")
       cssLines <- readLines(css.file)
       cssOut <- file(paste0(outfile.dir,"/styles.css"), "w")
       for(line in cssLines){
         writeLines(line, cssOut)
       }
       close(cssOut)
     }
     css.path.MAIN <- "styles.css"
     if(subdirectories.by.type){
       css.path.SUB <- "../styles.css"
     } else {
       css.path.SUB <- "styles.css"
     }
   } else {
     if(verbose) message("   Linking to external css stylesheet: \"",css.link,"\"")
     css.path.MAIN <- css.link
     css.path.SUB <- css.link
   }
   if(verbose) message("   Writing html index. ",date())
   fitExpToVar="condition"
   options(stringsAsFactors = FALSE)
     pf <- file(paste0(outfile.dir,"/", mainFile), "w")
     
     writeLines(paste0("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">"), pf)
     writeLines(paste0("<html><head><title> JunctionSeq Results </title>"), pf)
     writeLines(paste0("<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">"), pf)
     writeLines(paste0("<link rel=\"stylesheet\" type=\"text/css\" href=\"",css.path.MAIN,"\">"), pf)
     writeLines(paste0("</head> <body> <div style=\"margin:25px\">"), pf)
     
     writeLines(paste0("<h1> JunctionSeq Results </h1>"),pf)
     writeLines(paste0("<h4> Generated using JunctionSeq v",packageVersion("JunctionSeq")," at ",date()," </h4>"),pf)
     
     analysisNavTable <- makeAnalysisNavTable(compare.analysis.list = compare.analysis.list, htmlPrefix = paste0("../"), htmlSuffix = paste0("/",mainFile), style = "", verbose = verbose)
     writeLines(analysisNavTable,pf)
     
     if(verbose) message("   Writing formula data. ",date())
     writeLines(paste0("<h3>","Formulas:","</h3>"),pf)
     for(i in 1:length(jscs@formulas)){
       writeLines(paste0("<code>",names(jscs@formulas)[i],":",as.character(jscs@formulas[i]),"</code><br>"),pf)
     }
     
     if(verbose) message("   Writing methods data. ",date())
     if(! is.null(attr(jscs,"AltMethods"))){
       writeLines(paste0("<h3>","Methods:","</h3>"),pf)
       writeLines("<table>",pf)
       altMethods <- attr(jscs,"AltMethods")
       writeLines(paste0("<tr><td><code>jscs@analysis.type</code></td><td> <code>&quot;",as.character(jscs@analysisType),"&quot;</code></td></tr>"),pf)
       for(i in 1:length(altMethods)){
         writeLines(paste0("<tr><td><code>",names(altMethods)[i],"</code></td><td> <code>&quot;",altMethods[i],"&quot; </code></td></tr>"),pf)
       }
       writeLines("</table>",pf)
     }
     
     if(verbose) message("   Writing sample data. ",date())
     writeLines(paste0("<h3>","Sample Data:","</h3>"),pf)
     sampData <- as.matrix(pData(jscs)[colnames(pData(jscs)) != "countfiles"])
     writeLines(paste0("<table>"),pf)
     writeLines(paste0("<tr><th>",paste0("sampleID",colnames(sampData),collapse="</th><th>"),"</th></tr>"),pf)
     for(i in 1:(nrow(sampData))){
       writeLines(paste0("<tr><td>",paste0(rownames(sampData)[i],sampData[i,],collapse="</td><td>"),"</td></tr>"),pf)
     }
     writeLines(paste0("</table>"),pf)
     
     if(verbose) message("   Writing dispersion data. ",date())
     writeLines(paste0("<h3>","Dispersion Method:","</h3>"),pf)
     writeLines(paste0("<table><tr>"),pf)
     for(dft in names(jscs@dispFunctionType)){
        writeLines(paste0("<th>",dft,"</th>"),pf)
     }
     writeLines(paste0("</tr><tr>"),pf)
     for(dft in names(jscs@dispFunctionType)){
        writeLines(paste0("<td>",jscs@dispFunctionType[[dft]],"</td>"),pf)
     }
     writeLines(paste0("</tr></table><br>"),pf)
     
     if(jscs@dispFunctionType[["fitType"]] == "parametric"){
       
       if(is.null(attr(jscs@dispFunction,"coefficients")) || is.null(attr(jscs@dispFunctionExon,"coefficients")) || is.null(attr(jscs@dispFunctionJct,"coefficients"))){
         #do nothing?
       } else {
       writeLines(paste0("<h3>","Dispersion Fits:","</h3>"),pf)
       writeLines(paste0("<table>"),pf)
       writeLines(paste0("<tr><th>Fit</th><th>asymptDisp</th> <th>extraPois</th></tr>"),pf)
       writeLines(paste0("<tr>","<th>Overall</th>","<td>",paste0(attr(jscs@dispFunction,"coefficients"),collapse="</td><td>"),"</td>","</tr>"),pf)
       writeLines(paste0("<tr>","<th>Exonic regions</th>","<td>",paste0(attr(jscs@dispFunctionExon,"coefficients"),collapse="</td><td>"),"</td>","</tr>"),pf)
       writeLines(paste0("<tr>","<th>Splice Junctions</th>","<td>",paste0(attr(jscs@dispFunctionJct,"coefficients"),collapse="</td><td>"),"</td>","</tr>"),pf)
       writeLines(paste0("</table><br>"),pf)
       }
     }
     
     if(verbose) message("   Writing summary plots. ",date())
     writeLines(paste0("<h3>","Summary Plots:","</h3>"),pf)
     if(ma.plot | variance.plot){
       writeLines(paste0("<table>"),pf)
       writeLines("<tr>",pf)
       if(variance.plot){
         disp.file <- paste0("dispersion-plot",plotting.device.ext)
         writeLines(paste0("<td><a href=\"",disp.file,"\"> ",
                           "<img src=\"",disp.file,"\" "," alt=\"Dispersion Plot\" style=\"width: 45vw;\">",
                           "</a></td>"),pf)
       }
       if(ma.plot){
         fc.cols <- which(strStartsWith(names(fData(jscs)), "log2FC("))
         fc.names <- names(fData(jscs))[fc.cols]
         if(length(fc.names) == 1){
           fc.file <- paste0("ma-plot-",sub("/","vs",fc.names[1], fixed = TRUE), plotting.device.ext)
           writeLines(paste0("<td><a href=\"",fc.file,"\"> ",
                             "<img src=\"",fc.file,"\" "," alt=\"MA Plot\" style=\"width: 45vw;\">",
                             "</a></td>"),pf)
         } else {
           for(j in 1:length(fc.names)){
             fc.file <- paste0("ma-plot-",sub("/","vs",fc.names[j], fixed = TRUE), plotting.device.ext)
             writeLines(paste0("<td><a href=\"",fc.file,"\"> ",
                               "<img src=\"",fc.file,"\" "," alt=\"MA Plot (",fc.names[j],")\" style=\"width: 45vw;\">",
                               "</a></td>"),pf)
             if(j %% 2 == 1 & j < length(fc.names)){
               writeLines("</tr><tr>",pf)
             }
           }
         }
       }

       writeLines(paste0("</tr></table>"),pf)
     }
     
     if(length(gene.list) == 0){
       writeLines("<h2>No genes plotted!</h2><br>", pf)
       writeLines(paste0("</div> </body> </html>"), pf)
       close(pf)
       if(verbose) message("   Html index complete. ",date())
       if(verbose) message("   Finished all html files. ",date())
       return(TRUE)
     } else {
       writeLines("<h2>Gene Table:</h2><br>", pf)

       mainTable <- data.frame(geneID = as.character(gene.list), stringsAsFactors=FALSE)
       geneData <- jscs@flatGffGeneData
       geneAnno <- as.data.frame(t(sapply(gene.list, function(g){
         geneRows <- which(fData(jscs)$geneID == g)
         nameRow <- which(geneData$geneID == g)
         c(as.character(geneData[nameRow,"gene_name"]),
           as.character(fData(jscs)$chr[geneRows[1]]), 
           as.numeric(min(fData(jscs)$start[geneRows])), 
           as.numeric(max(fData(jscs)$end[geneRows])),
           as.character(fData(jscs)$strand[geneRows[1]])
         )
       })))
       if(is.null(number.plots)){
            number.plots <- rep("",length(gene.list))
       }

       if(verbose) message("   Compiling data table. ",date())
       colnames(geneAnno) <- c("name","chr","start","end","strand")
       mainTable <- cbind.data.frame(mainTable, geneAnno)
       mainTable$chr <- as.character(mainTable$chr)
       geneBaseMeans <- rowMeans(jscs@geneCountData[match(gene.list,rownames(jscs@geneCountData)),, drop=FALSE] / sizeFactors(jscs))
       mainTable$baseMean <- sprintf("%.1f",geneBaseMeans)

       if(! is.null(fData(jscs)$geneWisePadj)){
         mainTable$geneWisePadj <- sapply(gene.list, function(g){
           geneRows <- which(fData(jscs)$geneID == g)
           min( fData(jscs)$geneWisePadj[geneRows] , na.rm = TRUE)
         })
       }

       mainTable$mostSigID <- sapply(gene.list, function(g){
         geneRows <- which(fData(jscs)$geneID == g)
         fData(jscs)$countbinID[ geneRows[which.min( fData(jscs)$padjust[geneRows])] ]
       })


       mainTable$mostSigPadjust <- sapply(gene.list, function(g){
         geneRows <- which(fData(jscs)$geneID == g)
         fData(jscs)$padjust[ geneRows[which.min( fData(jscs)$padjust[geneRows])] ]
       })
       mainTable$mostSigPadjust <- sprintf("%.3g",mainTable$mostSigPadjust)

       gene.row.list <- lapply(gene.list, function(g){  which(fData(jscs)$geneID == g) })

       numExons <- sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$featureType[geneRows] == "exonic_part", na.rm = TRUE)
       })
       numKnown <- sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$featureType[geneRows] == "splice_site", na.rm = TRUE)
       })
       numNovel <- sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$featureType[geneRows] == "novel_splice_site", na.rm = TRUE)
       })

       exonsSig <- sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$padjust[geneRows] < FDR.threshold & fData(jscs)$featureType[geneRows] == "exonic_part", na.rm = TRUE)
       })
       knownSig <- sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$padjust[geneRows] < FDR.threshold & fData(jscs)$featureType[geneRows] == "splice_site", na.rm = TRUE)
       })
       novelSig <- sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$padjust[geneRows] < FDR.threshold & fData(jscs)$featureType[geneRows] == "novel_splice_site", na.rm = TRUE)
       })

       mainTable$numFeatures = paste0(numExons,"/",numKnown,"/",numNovel)
       mainTable$numSig = paste0(exonsSig,"/",knownSig,"/",novelSig)

       tableTitles <- colnames(mainTable)
       tableTitles[tableTitles == "numFeatures"] <- "&#35; Features<br><div style=\"font-size: smaller\">(Exon/Known/Novel)</div>"
       tableTitles[tableTitles == "numSig"] <- "&#35; Sig"

       titleLine <- paste0("<tr><td rowspan=2>", paste0(tableTitles, collapse="</td> <td rowspan=2>") )

       plotColSpan <- (expr.plot + normCounts.plot + rExpr.plot + rawCounts.plot) * (with.TX + without.TX)

       if(plotColSpan > 0){
         titleLine <- paste0(titleLine,"<td colspan=",plotColSpan,">Plots and Tables</td>")
         titleLine <- paste0(titleLine,"</tr><tr>")
       }

       ##############




       ##############

       if(expr.plot){
         if(without.TX){
           mainTable$expr.plot <- paste0("<a href=\"htmlFiles/",gene.list,"-expr.plot.html\">","PLOT","</a><br>",
                                         "<a href=\"htmlFiles/",gene.list,"-expr.table.html\">","TABLE","</a>")
           titleLine <- paste0(titleLine,"<td>expr</td>")
         }
         if(with.TX){
           mainTable$expr.plot.TX <- paste0("<a href=\"htmlFiles/",gene.list,"-expr-withTX.plot.html\">","PLOT","</a><br>",
                                         "<a href=\"htmlFiles/",gene.list,"-expr.table.html\">","TABLE","</a>")
           titleLine <- paste0(titleLine,"<td>exprTX</td>")
         }
       }
       if(normCounts.plot){
         if(without.TX){
         mainTable$normCt.plot <- paste0("<a href=\"htmlFiles/",gene.list,"-normCounts.plot.html\">","PLOT","</a><br>",
                                         "<a href=\"htmlFiles/",gene.list,"-normCounts.table.html\">","TABLE","</a>")
         titleLine <- paste0(titleLine,"<td>normCt</td>")
         }
         if(with.TX){
         mainTable$normCt.plot.TX <- paste0("<a href=\"htmlFiles/",gene.list,"-normCounts-withTX.plot.html\">","PLOT","</a><br>",
                                         "<a href=\"htmlFiles/",gene.list,"-normCounts.table.html\">","TABLE","</a>")
         titleLine <- paste0(titleLine,"<td>normCtTX</td>")
         }
       }
       if(rExpr.plot){
         if(without.TX){
         mainTable$rExpr.plot <- paste0("<a href=\"htmlFiles/",gene.list,"-rExpr.plot.html\">","PLOT","</a><br>",
                                       "<a href=\"htmlFiles/",gene.list,"-rExpr.table.html\">","TABLE","</a>")
         titleLine <- paste0(titleLine,"<td>relExpr</td>")
         }
         if(with.TX){
         mainTable$rExpr.plot.TX <- paste0("<a href=\"htmlFiles/",gene.list,"-rExpr-withTX.plot.html\">","PLOT","</a><br>",
                                       "<a href=\"htmlFiles/",gene.list,"-rExpr.table.html\">","TABLE","</a>")
         titleLine <- paste0(titleLine,"<td>relExprTX</td>")
         }
       }
       if(rawCounts.plot){
         if(without.TX){
         mainTable$rawCounts.plot <- paste0("<a href=\"htmlFiles/",gene.list,"-rawCounts.plot.html\">","PLOT","</a><br>",
                                       "<a href=\"htmlFiles/",gene.list,"-rawCounts.table.html\">","TABLE","</a>")
         titleLine <- paste0(titleLine,"<td>rawCt</td>")
         }
         if(with.TX){
         mainTable$rawCounts.plot.TX <- paste0("<a href=\"htmlFiles/",gene.list,"-rawCounts-withTX.plot.html\">","PLOT","</a><br>",
                                       "<a href=\"htmlFiles/",gene.list,"-rawCounts.table.html\">","TABLE","</a>")
         titleLine <- paste0(titleLine,"<td>rawCtTX</td>")
         }
       }
       titleLines <- paste0(titleLine,"</tr>")



       if(verbose) message("   Writing data table. ",date())
       writeLines("<table BORDER=1>",pf)
       writeLines(titleLine,pf)
       for(i in 1:nrow(mainTable)){
         writeLines("   <tr>", pf)
           line <- paste0("<td>", paste0( as.character(mainTable[i,]), collapse=" </td><td> " ) ,"</td> </tr>" )
         writeLines(line,pf)
       }
       writeLines("</table>",pf)

   writeLines(paste0("</div> </body> </html>"), pf)
   close(pf)
   if(verbose) message("   Html index complete. ",date())
   flat.gff.data <- jscs@flatGffData

   #TEMP VERSION, needs to be upgraded?

   if(verbose) message("   Writing pages. ",date())
       for(i in 1:length(gene.list)){

         if(i == 1){
           prev.g <- NULL
         } else {
           prev.g <- paste0(gene.list[i - 1])
         }
         if(i == length(gene.list)){
           next.g <- NULL;  
         } else {
           next.g <- paste0(gene.list[i + 1])
         }
         g <- gene.list[i]

         navTable <- makeNavTable(g, expr.plot=expr.plot, normCounts.plot=normCounts.plot, rExpr.plot=rExpr.plot, rawCounts.plot=rawCounts.plot, with.TX = with.TX, without.TX = without.TX, mainFile = mainFile)
         #taken from func buildAllPlotsForGene

          transcripts <- sapply(sapply(flat.gff.data$transcripts[which(flat.gff.data$gene_id==g & flat.gff.data$featureType == "exonic_part")],
                                       toString), 
                                function(x){strsplit(x, "+",fixed=TRUE)}
                                )
          trans <- Reduce(union, transcripts)
          tx.ct <- length(trans)

          withTxPlot.height.multiplier <- getAutofitTxRelativeHeight(tx.ct, autoscale.height.to.fit.TX.annotation = autoscale.height.to.fit.TX.annotation, GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, CONNECTIONS.relative.height = CONNECTIONS.relative.height, TX.margins = TX.margins,SPLICE.annotation.relative.height=SPLICE.annotation.relative.height)

          if(is.na(autoscale.width.to.fit.bins) || autoscale.width.to.fit.bins == 0){
            width.multiplier <- 1
          } else {
            num.cols <- getColCt(geneID = g, merged.data = fData(jscs), 
                                      plot.exon.results = plot.exon.results, plot.junction.results = plot.junction.results, plot.novel.junction.results=plot.novel.junction.results, 
                                      plot.untestable.results=plot.untestable.results)
            if(num.cols > autoscale.width.to.fit.bins){
              width.multiplier <- num.cols / autoscale.width.to.fit.bins
            } else {
              width.multiplier <- 1
            }
          }

         makeAllPlotPages(g = g, prev.g = prev.g, next.g = next.g, first.g = gene.list[1], last.g = gene.list[length(gene.list)],
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
                          html.fixed.dim = html.fixed.dim,
                          html.width = html.width, html.width.units = html.width.units,
                          number.plots = number.plots,
                          FDR.threshold = colorRed.FDR.threshold,
                          css.path = css.path.SUB,
                          compare.analysis.list = compare.analysis.list,
                          minimalImageFilenames = minimalImageFilenames,
                          verbose = verbose, debug.mode = debug.mode)
       }

     if(verbose) message("   Finished all html files. ",date())
     return(TRUE)
   }
}


formatArray <- function( fmt, a){
   as.matrix(apply(a,MARGIN=c(1,2), function(x){
      sprintf(fmt,x)
   }))
}


makeGeneWiseTable <- function(jscs, gene.list, FDR.threshold = 0.05, verbose = TRUE, debug.mode = FALSE){
       if(verbose) message("   Compiling data table. ",date())
       
       mainTable <- data.frame(geneID = as.character(gene.list), stringsAsFactors=FALSE)
       row.names(mainTable) <- gene.list
       mainTable <- AnnotatedDataFrame(mainTable)
       varMetadata(mainTable)["geneID", "labelDescription"] <- "Gene Unique Identifier"
       
       noGenes <- (length(gene.list) == 0)
       
       if(! noGenes){
         geneData <- jscs@flatGffGeneData
         geneAnno <- as.data.frame(t(sapply(gene.list, function(g){
           geneRows <- which(fData(jscs)$geneID == g)
	   nameRow <- which(geneData$geneID == g)
           c(as.character(geneData[nameRow,"gene_name"]),
             as.character(fData(jscs)$chr[geneRows[1]]), 
             as.numeric(min(fData(jscs)$start[geneRows])), 
             as.numeric(max(fData(jscs)$end[geneRows])),
             as.character(fData(jscs)$strand[geneRows[1]])
           )
         })))
         colnames(geneAnno) <- c("name","chr","start","end","strand")
       } else {
         geneAnno <- data.frame(name = character(), chr = character(), start = numeric(), end = numeric(), strand = character())
       }
       
       mainTable$name <- as.character(geneAnno$name)
       mainTable$chr <- as.character(geneAnno$chr)
       mainTable$start <- geneAnno$start
       mainTable$end <- geneAnno$end
       mainTable$strand <- geneAnno$strand
       varMetadata(mainTable)[c("name","chr","start","end","strand"), "labelDescription"] <- 
                              c("Gene name(s)",
                                "Gene chromosome",
                                "Gene start",
                                "Gene end",
                                "Gene strand")
       
       #message("2")
       geneBaseMeans <- if(noGenes){numeric()} else { rowMeans(jscs@geneCountData[match(gene.list,rownames(jscs@geneCountData)),, drop=FALSE] / sizeFactors(jscs))}
       mainTable$baseMean <- if(noGenes){character()} else {sprintf("%.1f",geneBaseMeans)}
       varMetadata(mainTable)["baseMean", "labelDescription"] <- "Gene BaseMean (simple normalized mean read or read-pair count per sample)"
       
       if(! is.null(fData(jscs)$geneWisePadj)){
         mainTable$geneWisePadj <- if(noGenes){ numeric()} else {sapply(gene.list, function(g){
           geneRows <- which(fData(jscs)$geneID == g)
           min( fData(jscs)$geneWisePadj[geneRows] , na.rm = TRUE)
         })}
         varMetadata(mainTable)["geneWisePadj", "labelDescription"] <- "Gene-level adjusted p-value. P-value for the hypothesis that one or more features are DU."
       }

       mainTable$mostSigID <- if(noGenes){character()} else {sapply(gene.list, function(g){
         geneRows <- which(fData(jscs)$geneID == g)
         fData(jscs)$countbinID[ geneRows[which.min( fData(jscs)$padjust[geneRows])] ]
       })}
       varMetadata(mainTable)["mostSigID", "labelDescription"] <- "Feature ID of the most singificant feature."
       
       mainTable$mostSigPadjust <- if(noGenes){numeric()} else {sapply(gene.list, function(g){
         geneRows <- which(fData(jscs)$geneID == g)
         fData(jscs)$padjust[ geneRows[which.min( fData(jscs)$padjust[geneRows])] ]
       })}
       mainTable$mostSigPadjust <- sprintf("%.3g",mainTable$mostSigPadjust)
       varMetadata(mainTable)["mostSigPadjust", "labelDescription"] <- "Adjusted p-value of the most singificant feature."
       
       gene.row.list <- if(noGenes){ integer() } else {lapply(gene.list, function(g){  which(fData(jscs)$geneID == g) })}

       mainTable$numExons <- if(noGenes){ integer() } else { sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$featureType[geneRows] == "exonic_part", na.rm = TRUE)
       })}
       varMetadata(mainTable)["numExons", "labelDescription"] <- "Number of distinct exonic regions belonging to the gene."

       mainTable$numKnown <- if(noGenes){ integer() } else {sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$featureType[geneRows] == "splice_site", na.rm = TRUE)
       })}
       varMetadata(mainTable)["numKnown", "labelDescription"] <- "Number of distinct known splice sites belonging to the gene."
       mainTable$numNovel <- if(noGenes){ integer() } else {sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$featureType[geneRows] == "novel_splice_site", na.rm = TRUE)
       })}
       varMetadata(mainTable)["numNovel", "labelDescription"] <- "Number of distinct novel splice sites belonging to the gene."

       mainTable$exonsSig <- if(noGenes){ integer() } else {sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$padjust[geneRows] < FDR.threshold & fData(jscs)$featureType[geneRows] == "exonic_part", na.rm = TRUE)
       })}
       varMetadata(mainTable)["exonsSig", "labelDescription"] <- paste0("Number of signficant exonic regions at p-adjust < ", FDR.threshold)
       mainTable$knownSig <- if(noGenes){ integer() } else {sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$padjust[geneRows] < FDR.threshold & fData(jscs)$featureType[geneRows] == "splice_site", na.rm = TRUE)
       })}
       varMetadata(mainTable)["knownSig", "labelDescription"] <- paste0("Number of signficant known splice junctions at p-adjust < ", FDR.threshold)
       mainTable$novelSig <- if(noGenes){ integer() } else {sapply(gene.row.list, function(geneRows){
         sum(fData(jscs)$padjust[geneRows] < FDR.threshold & fData(jscs)$featureType[geneRows] == "novel_splice_site", na.rm = TRUE)
       })}
       varMetadata(mainTable)["novelSig", "labelDescription"] <- paste0("Number of signficant novel splice junctions at p-adjust < ", FDR.threshold)

       mainTable$numFeatures = if(noGenes){ character()} else {paste0(mainTable$numExons,"/",mainTable$numKnown,"/",mainTable$numNovel)}
       varMetadata(mainTable)["numFeatures", "labelDescription"] <- "Number exonic regions / num known SJ / num novel SJ"
       mainTable$numSig = if(noGenes){ character()} else {paste0(mainTable$exonsSig,"/",mainTable$knownSig,"/",mainTable$novelSig)}
       varMetadata(mainTable)["numSig", "labelDescription"] <- "Number sig exonic regions / num sig known SJ / num sig novel SJ"
       
       return(mainTable)
}

makeAllPlotPages <- function(g, prev.g, next.g, first.g, last.g, gene.number,
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
                            html.fixed.dim = "height",
                            base.html.height = 90, base.html.height.units = "vh",
                            html.width = 90, html.width.units = "vw",
                            number.plots = FALSE,
                            FDR.threshold = 0.01,
                            css.path = NULL,
                            compare.analysis.list = NULL,
                            minimalImageFilenames = TRUE,
                            verbose = TRUE, debug.mode = FALSE
                              ){
  geneRows <- which(fData(jscs)$geneID == g);    
  td.base <- fData(jscs)[geneRows, c("countbinID","chr","start","end","testable","status","baseMean","dispersion","pvalue","padjust"), drop = FALSE]
  pval.numeric <- td.base$padjust
  isTestable <- td.base$testable
  names(td.base)[1] <- "binID"
  td.base$status <- as.character(td.base$status)
  td.base$chr <- as.character(td.base$chr)
  td.base$baseMean <- sprintf("%.1f",td.base$baseMean)
  td.base$dispersion <- sprintf("%.3g",td.base$dispersion)
  td.base$pvalue <- sprintf("%.4g",td.base$pvalue)
  td.base$padjust <- sprintf("%.4g",td.base$padjust)
  

  
  LFC.idx <- which(grepl("^log2FC",colnames(fData(jscs))))
    for(i in LFC.idx){
      td.base[[colnames(fData(jscs))[i]]] <- sprintf("%.3f",fData(jscs)[geneRows, i])
    }
  
  titleline.base <- colnames(td.base)
  for(i in which(grepl("^log2FC",titleline.base))){
    titleline.base[i] <- sub("\\(","<br>\\(",titleline.base[i])
  }
  titleline.base <- paste0("<tr><td rowspan=2>",  paste0(titleline.base,collapse="</td><td rowspan=2>"),"</td>")

  dirPath <- if(subdirectories.by.type){
    c("../expr/",
      "../exprTX/",
      "../normCounts/",
      "../normCountsTX/",
      "../rExpr/",
      "../rExprTX/",
      "../rawCounts/",
      "../rawCountsTX/")
  } else {
    rep("../",8)
  }
  htmlSuffixes <- c(
      "-expr.plot.html",
      "-expr-withTX.plot.html",
      "-normCounts.plot.html",
      "-normCounts-withTX.plot.html",
      "-rExpr.plot.html",
      "-rExpr-withTX.plot.html",
      "-rawCounts.plot.html",
      "-rawCounts-withTX.plot.html"
  )
  tableSuffixes <- c(
      "-expr.table.html",
      "-expr.table.html",
      "-normCounts.table.html",
      "-normCounts.table.html",
      "-rExpr.table.html",
      "-rExpr.table.html",
      "-rawCounts.table.html",
      "-rawCounts.table.html"
  )
  suffixes <- paste0("-",IMAGE.NAMES); #see 00.minor.utils.R
  
  pageTitles <- c(
      "expression",
      "expression (with TX)",
      "norm counts",
      "norm counts (with TX)",
      "relative expression",
      "relative exprssion (with TX)",
      "raw counts",
      "raw counts (with TX)"
  )
  pageTables <- c(
     
  )

  pages.idx <- 1:length(suffixes)
  names(pages.idx) <- suffixes
  
  if(! without.TX){
    pages.idx <- pages.idx[ grepl(IMAGE.NAME.TX,names(pages.idx), fixed=TRUE) ]
  }
  if(! with.TX){
    pages.idx <- pages.idx[ ! grepl(IMAGE.NAME.TX,names(pages.idx), fixed=TRUE) ]
  }
  
  if(! expr.plot){
    pages.idx <- pages.idx[ ! grepl(IMAGE.NAME.MAP[["expr"]],names(pages.idx), fixed=TRUE) ]
  }
  if(! normCounts.plot){
    pages.idx <- pages.idx[ ! grepl(IMAGE.NAME.MAP[["normCounts"]],names(pages.idx), fixed=TRUE) ]
  }
  if(! rExpr.plot){
    pages.idx <- pages.idx[ ! grepl(IMAGE.NAME.MAP[["rExpr"]],names(pages.idx), fixed=TRUE) ]
  }
  if(! rawCounts.plot){
    pages.idx <- pages.idx[ ! grepl(IMAGE.NAME.MAP[["rawCounts"]],names(pages.idx), fixed=TRUE) ]
  }
  
  if(minimalImageFilenames){
    imgGeneName <- ""
  } else {
    imgGeneName <- paste0(g,"-");
  }
  
  if(expr.plot){
     if(without.TX){
       html.suffix <- "-expr.plot.html"
       makePlotPage(plotfile = paste0(dirPath[1],gene.number,imgGeneName,"expr",plotting.device.ext),
                  pageTitle = paste0(g, " expression"),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-expr.plot.html"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height),
                  base.html.height.units = base.html.height.units,
                  html.fixed.dim = html.fixed.dim,
                  html.width = html.width, html.width.units = html.width.units,
                  css.path = css.path)
     }
     if(with.TX){
       html.suffix <-  "-expr-withTX.plot.html"
       makePlotPage(plotfile = paste0(dirPath[2],gene.number,imgGeneName,"expr-TX",plotting.device.ext),
                  pageTitle = paste0(g, " expression (with TX)"),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-expr-withTX.plot.html"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height * withTxPlot.height.multiplier),
                  base.html.height.units = base.html.height.units,
                  html.fixed.dim = html.fixed.dim,
                  html.width = html.width, html.width.units = html.width.units,
                  css.path = css.path)
     }
     
     
     td <- formatArray("%.1f",jscs@plottingEstimates[["exprEstimate"]][geneRows,, drop=FALSE])
     colnames(td) <- sub("expr_","", colnames(td),fixed=TRUE)
     titleline <- paste0(titleline.base,
                         "<td colspan=",ncol(td),"> Mean Normalized Coverage Counts </td> </tr> ",
                         "<tr> <td>",paste0(colnames(td),collapse="</td><td>"),"</td></tr>")


     td <- cbind(td.base, td)
     

     
     html.suffix <- "-expr.table.html"
     makeTablePage(pvals = pval.numeric, isTestable = isTestable, td = td,
                  titleline = titleline,
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-expr.table.html"),
                  pageTitle = paste0("Feature Coverage/Expression Table (",g,")"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  css.path = css.path,
                  FDR.threshold = FDR.threshold)
  }
  if(normCounts.plot){
     if(without.TX){
     html.suffix <-  "-normCounts.plot.html"
     makePlotPage(plotfile = paste0(dirPath[3],gene.number,imgGeneName,"normCts",plotting.device.ext),
                  pageTitle = paste0(g, " norm counts"),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-normCounts.plot.html"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height)  ,
                  base.html.height.units = base.html.height.units,
                  html.fixed.dim = html.fixed.dim,
                  html.width = html.width, html.width.units = html.width.units,
                  css.path = css.path)
     }
     if(with.TX){
     html.suffix <-  "-normCounts-withTX.plot.html"
     makePlotPage(plotfile = paste0(dirPath[4],gene.number,imgGeneName,"normCts-TX",plotting.device.ext),
                  pageTitle = paste0(g, " norm counts"),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-normCounts-withTX.plot.html"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height  * withTxPlot.height.multiplier)  ,
                  base.html.height.units = base.html.height.units,
                  html.fixed.dim = html.fixed.dim,
                  html.width = html.width, html.width.units = html.width.units,
                  css.path = css.path)
     }
     td <- formatArray("%.1f",jscs@plottingEstimates[["normCounts"]][geneRows,, drop=FALSE])
     td <- td[,1:(ncol(td)/2), drop = FALSE]
     colnames(td) <- sub("normCount_","", colnames(td),fixed=TRUE)
     titleline <- paste0(titleline.base,
                         "<td colspan=",ncol(td),"> Normalized Coverage Counts </td> </tr> ",
                         "<tr> <td>",paste0(colnames(td),collapse="</td><td>"),"</td></tr>")
     td <- cbind(td.base, td)
     html.suffix <- "-normCounts.table.html"
     makeTablePage(pvals = pval.numeric, isTestable = isTestable, td = td,
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-normCounts.table.html"),
                  pageTitle = paste0("Normalized Counts Table (",g,")"),
                  titleline = titleline,
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  css.path = css.path,
                  FDR.threshold = FDR.threshold)
  }
  if(rExpr.plot ){
     if(without.TX){
     html.suffix <-  "-rExpr.plot.html"
     makePlotPage(plotfile = paste0(dirPath[5],gene.number,imgGeneName,"rExpr",plotting.device.ext),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rExpr.plot.html"),
                  pageTitle = paste0(g, " relative expression"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height) ,
                  base.html.height.units = base.html.height.units,
                  html.fixed.dim = html.fixed.dim,
                  html.width = html.width, html.width.units = html.width.units,
                  css.path = css.path )
     }
     if(with.TX){
     html.suffix <-  "-rExpr-withTX.plot.html"
     makePlotPage(plotfile = paste0(dirPath[6],gene.number,imgGeneName,"rExpr-TX",plotting.device.ext),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rExpr-withTX.plot.html"),
                  pageTitle = paste0(g, " relative expression (with TX)"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height * withTxPlot.height.multiplier) ,
                  base.html.height.units = base.html.height.units,
                  html.fixed.dim = html.fixed.dim,
                  html.width = html.width, html.width.units = html.width.units,
                  css.path = css.path )
     }
     td <- formatArray("%.1f",jscs@plottingEstimates[["relExprEstimate"]][geneRows,, drop=FALSE])
     colnames(td) <- sub("relExpr_","", colnames(td),fixed=TRUE)
     titleline <- paste0(titleline.base,
                         "<td colspan=",ncol(td),"> Relative Coverage </td> </tr> ",
                         "<tr> <td>",paste0(colnames(td),collapse="</td><td>"),"</td></tr>")
     td <- cbind(td.base, td)
     html.suffix <- "-rExpr.table.html"
     makeTablePage(pvals = pval.numeric, isTestable = isTestable, td = td,
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rExpr.table.html"),
                  pageTitle = paste0("Relative Coverage Table, relative to gene-wide expression (",g,")"),
                  titleline = titleline,
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  css.path = css.path,
                  FDR.threshold = FDR.threshold)
  }
  if(rawCounts.plot ){
     
     if(without.TX){
     html.suffix <- "-rawCounts.plot.html"
     makePlotPage(plotfile = paste0(dirPath[7],gene.number,imgGeneName,"rawCts",plotting.device.ext),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rawCounts.plot.html"),
                  pageTitle = paste0(g, " Raw Counts"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height) ,
                  base.html.height.units = base.html.height.units,
                  html.fixed.dim = html.fixed.dim,
                  html.width = html.width, html.width.units = html.width.units,
                  css.path = css.path )
     }
     if(with.TX){
     html.suffix <- "-rawCounts-withTX.plot.html"
     makePlotPage(plotfile = paste0(dirPath[8],gene.number,imgGeneName,"rawCts-TX",plotting.device.ext),
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rawCounts-withTX.plot.html"),
                  pageTitle = paste0(g, " Raw Counts (with TX)"),
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  width = floor(base.html.height * width.multiplier),
                  height = floor(base.html.height * withTxPlot.height.multiplier) ,
                  base.html.height.units = base.html.height.units,
                  html.fixed.dim = html.fixed.dim,
                  html.width = html.width, html.width.units = html.width.units,
                  css.path = css.path )
     }
     td <- formatArray("%.1f",counts(jscs)[geneRows,, drop=FALSE])
     titleline <- paste0(titleline.base,
                         "<td colspan=",ncol(td),"> Raw Counts </td> </tr> ",
                         "<tr> <td>",paste0(colnames(td),collapse="</td><td>"),"</td></tr>")
     td <- cbind(td.base, td)
     html.suffix <- "-rawCounts.table.html"
     makeTablePage(pvals = pval.numeric, isTestable = isTestable, td = td,
                  htmlfile = paste0(outfile.dir,"/htmlFiles/",g,"-rawCounts.table.html"),
                  pageTitle = paste0("Raw Counts Table (",g,")"),
                  titleline = titleline,
                  navTable = navTable,
                  geneNavTable = makeGeneNavTable(g,prev.g,next.g, first.g, last.g,mainFile,prefix="./",suffix=html.suffix, compare.analysis.list = compare.analysis.list),
                  css.path = css.path,
                  FDR.threshold = FDR.threshold);     
  }

}

makeAnalysisNavTable <- function(compare.analysis.list = NULL, htmlPrefix, htmlSuffix, style = "style=\"margin-left:auto; margin-right:auto;\"", verbose = FALSE){
  out <- ""
  if(! is.null(compare.analysis.list)){
     if(verbose) message("Making cross-analysis links.")
     
     out <- paste0(out,"<table ",style,">\n")
     
     if(! is.list(compare.analysis.list)){
        stop("ERROR: html.compare.results.list must be NULL or of type \"list\"!")
     }

     if(is.list(compare.analysis.list[[1]])){
       #Nested list: multi-row table!
       if(verbose) message("Writing multi-level analysis list.")
       
       for(i in 1:length(compare.analysis.list)){
          curr <- compare.analysis.list[[i]]
          if(! is.list(curr)){
            stop("ERROR: Attempting multi-level list, encountered non-list elemend!")
          }
          if(  is.null(names(curr))  ){
             names(curr) <- curr
          }
          
          out <- paste0(out,"   <tr>\n")
          for(j in 1:length(curr)){
            analysisName <- names(curr)[j]
            analysisDir <- curr[j]
            out <- paste0(out,"      <td><a href=\"",htmlPrefix,analysisDir,htmlSuffix,"\"> ",analysisName," </a> </td>\n")
          }
          out <- paste0(out,"   </tr>\n")
       }
     } else {
       if(is.null(names(compare.analysis.list))){
          names(compare.analysis.list) <- compare.analysis.list
       }
       out <- paste0(out,"   <tr>\n")
       for(i in 1:length(compare.analysis.list)){
         analysisName <- names(compare.analysis.list)[i]
         analysisDir <- compare.analysis.list[i]
         out <- paste0(out,"      <td><a href=\"",htmlPrefix,analysisDir,htmlSuffix,"\"> ",analysisName," </a> </td>\n")
       }
       out <- paste0(out,"   </tr>\n")
     }
     out <- paste0(out,"</table>\n",
                       "<br>\n")
  } else {
     # Do nothing.
  }
  return(out)
}

makeGeneNavTable <- function(g, prev.g, next.g, first.g, last.g, compare.analysis.list = NULL, mainFile, prefix, suffix){
  out <- makeAnalysisNavTable(compare.analysis.list, 
                              htmlPrefix = paste0("../../"), 
                              htmlSuffix = paste0("/htmlFiles/",g,suffix))
  
  
  first.link <- paste0("<a href=\"",prefix,first.g,suffix,"\"> &lt&lt </a>")
  last.link <-  paste0("<a href=\"",prefix,last.g,suffix,"\"> &gt&gt </a>")
  
  prev.link <- if(is.null(prev.g)){
    "&lt Prev Gene"
  } else {
    paste0("<a href=\"",prefix,prev.g,suffix,"\"> &lt Prev Gene </a>")
  }
  next.link <- if(is.null(next.g)){
    "Next Gene &gt"
  } else {
    paste0("<a href=\"",prefix,next.g,suffix,"\"> Next Gene &gt </a>")
  }
  out <- paste0(
          out,
          "<table style=\"","margin-left:auto; margin-right:auto;","\"> \n",
            "   <tr>",
              "<td  style=\"font-size:x-large;\">",
                first.link,
              "</td>",
              "<td  style=\"font-size:x-large;\">",
                 prev.link,
              "</td>",
              "<td style=\"font-size:x-large;\" COLSPAN=2>",
                "<a href=\"../",mainFile,"\"> TOP </a>",
              "</td>",
              "<td  style=\"font-size:x-large;\">",
                 next.link,
              "</td>",
              "<td  style=\"font-size:x-large;\">",
                 last.link,
              "</td>",
            "</tr> \n",
          "</table> \n"
  )
  return(out)
}

makeNavTable <- function(g, expr.plot, normCounts.plot, rExpr.plot, rawCounts.plot, with.TX, without.TX, mainFile){
     
     navTable.1 <- ""; #paste0("<b><td rowspan = 3> <a href=\"../",mainFile,"\">BACK</a> </td>")
     navTable.2 <- ""
     navTable.3 <- ""
     
     nav.colSpan <- 0
     if(expr.plot & without.TX){
       navTable.1 <- paste0(navTable.1, "<td>expr</td>")
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-expr.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-expr.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1
     }
     if(expr.plot & with.TX){
       navTable.1 <- paste0(navTable.1, "<td>exprTX</td>")
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-expr-withTX.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-expr.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1
     }
     
     if(normCounts.plot & without.TX){
       navTable.1 <- paste0(navTable.1, "<td>normCt</td>")
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-normCounts.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-normCounts.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1
     }
     if(normCounts.plot & with.TX){
       navTable.1 <- paste0(navTable.1, "<td>normCtTX</td>")
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-normCounts-withTX.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-normCounts.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1
     }
     
     if(rExpr.plot & without.TX){
       navTable.1 <- paste0(navTable.1, "<td>rExpr</td>")
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-rExpr.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-rExpr.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1
     }
     if(rExpr.plot & with.TX){
       navTable.1 <- paste0(navTable.1, "<td>rExprTX</td>")
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-rExpr-withTX.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-rExpr.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1
     }
     
     if(rawCounts.plot & without.TX){
       navTable.1 <- paste0(navTable.1, "<td>rawCt</td>")
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-rawCounts.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-rawCounts.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1
     }
     if(rawCounts.plot & with.TX){
       navTable.1 <- paste0(navTable.1, "<td>rawCtTX</td>")
       navTable.2 <- paste0(navTable.2, paste0(
                                  "<td><a href=\"",g,"-rawCounts-withTX.plot.html\">","PLOT","</a></td>"
                              ))
       navTable.3 <- paste0(navTable.3, paste0(
                                  "<td><a href=\"",g,"-rawCounts.table.html\">","TABLE","</a></td>"
                              ))
       nav.colSpan <- nav.colSpan + 1
     }
     
     navTable <- c("<table style=\"","margin-left:auto;margin-right:auto;","\">",
                       "<th COLSPAN=",nav.colSpan,">",
                          "Gene ",g,
                       "</th>",
                       paste0("<tr>",c(navTable.1,navTable.2,navTable.3),"</tr>"),
                   "</table><br><br>")
     
     
     return(navTable)
}

makePlotPage <- function(plotfile, pageTitle = "", htmlfile, navTable, geneNavTable, width, height, css.path = NULL, 
                         base.html.height.units = "vh", 
                         html.fixed.dim = "height", html.width = 90, html.width.units = "vw"){
  f <- file(htmlfile, "w")
  
     writeLines(paste0("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">"), f)
     writeLines(paste0("<html><head><title> JunctionSeq Results (",pageTitle,") </title>"), f)
     writeLines(paste0("<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">"), f)
     if(! is.null(css.path)){
       writeLines(paste0("<link rel=\"stylesheet\" type=\"text/css\" href=\"",css.path,"\">"), f)
     }
     writeLines(paste0("</head> <body>"), f)
  
  writeLines(geneNavTable,f)
  imgStyle <- if(html.fixed.dim == "autofit"){ 
     paste0("display:block; margin-left:auto; margin-right:auto; max-height: 98vh; max-width:95vw;")
  } else if(html.fixed.dim == "width"){
     paste0("display:block; margin-left:auto; margin-right:auto; width:",html.width,html.width.units,"")
  } else {
     paste0("display:block; margin-left:auto; margin-right:auto; height:",height,base.html.height.units,"")
  }
  
  writeLines(
            paste0("<table style=\"margin-left:auto; margin-right:auto;\"> <tr><td>",
                   "<a href=\"",plotfile,"\"><img src=\"",plotfile,"\" "," alt=\"JunctionSeq Plot\" style=\"",imgStyle,"\"> </a>",
                   "</td></tr></table>"),f)
  
  writeLines("<br>",f)
  writeLines(navTable,f)
  writeLines(paste0("</body></html>"), f)
  close(f)
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
  )
}

makeTablePage <- function(pvals, isTestable, td, pageTitle, htmlfile, navTable, geneNavTable, titleline, css.path = NULL,
                          FDR.threshold = 0.01, pval.thresholds = NULL){
  
  if(is.null(pval.thresholds)){
    pval.thresholds <- c(FDR.threshold, 
                         FDR.threshold / 10, 
                         FDR.threshold / 100, 
                         FDR.threshold / 1000)
  }
  
  f <- file(htmlfile, "w")
  
     writeLines(paste0("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">"), f)
     writeLines(paste0("<html><head><title> JunctionSeq Results (",pageTitle,") </title>"), f)
     writeLines(paste0("<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">"), f)
     if(! is.null(css.path)){
       writeLines(paste0("<link rel=\"stylesheet\" type=\"text/css\" href=\"",css.path,"\">"), f)
     }
     writeLines(paste0("</head> <body> <div style=\"margin:25px\">"), f)
  
  
  writeLines(geneNavTable,f)
  writeLines(paste0("<h1>",pageTitle,"</h1>"),f)
  writeLines("<br>",f)
  
  writeLines("<table border=1 class=\"dataTable\">",f)
  writeLines(titleline,f)
  for(i in 1:nrow(td)){
    
    if(! isTestable[i]){
      classText <- paste0("untestable")
    } else if(pvals[i] < pval.thresholds[4]){
      classText <- paste0("sig4")
    } else if(pvals[i] < pval.thresholds[3]){
      classText <- paste0("sig3")
    } else if(pvals[i] < pval.thresholds[2]){
      classText <- paste0("sig2")
    } else if(pvals[i] < pval.thresholds[1]){
      classText <- paste0("sig1")
    } else {
      classText <- paste0("sigNo")
    }
    
    rowLine <- if(td[i, colnames(td) == "status"] != "OK"){
      RL <- sapply(1:length(td[i,]),function(j){
        if(colnames(td)[j] != "status"){
          paste0("<td>",as.character(td[i,j]),"</td>")
        } else {
          paste0("<td style=\"font-size:xx-small\">",as.character(td[i,j]),"</td>")
        }
      })
      paste0(RL,collapse=" ")
    } else {
      paste0("<td>", paste0(as.character(td[i,]),collapse="</td><td>"),"</td>")
    }
    
    writeLines(paste0("<tr class=\"",classText,"\">",rowLine,"</tr></div>"),f)
  }
  writeLines("</table>",f)
  
  writeLines("<br>",f)
  writeLines(makeStyleKeyTable(pval.thresholds),f)
  
  writeLines("<br>",f)
  writeLines(navTable,f)
  
  writeLines(paste0("</div></body></html>"), f)
  close(f)
}

########################################################################################

#This trivial function fixes a minor bug in the "svg" device found in certain older versions of redhat linux.
fixRHEL5_svg_bug <- function(infile,outfile){
  inlines <- readLines(infile)
  outlines <- sub("symbol id = \"glyph", "symbol overflow=\"visible\" id=\"glyph",inlines, fixed=TRUE)
  writeLines(outlines,outfile)
}
#s/symbol id=\"glyph/symbol overflow=\"visible\" id=\"glyph/

writeExprBedTrack <- function(file, jscs, 
                                      trackLine = "track name='JctExpr' description='Junction Coverage Estimates, by group' itemRgb='On' visibility=3",
                                      only.with.sig.gene = FALSE,
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
                                      use.gzip = TRUE,
                                      verbose = TRUE){
  
  estimates <- jscs@plottingEstimates[["exprEstimate"]]
  estimate.names <- sapply(colnames(estimates), function(s){substring(s,6,nchar(s))})
  numcond <- length(estimate.names)
  
  if(is.null(group.RGB)){
    default.color.list <- c("red","blue","orange","green3","purple","cyan1", "magenta","yellow3","tan4")
    if(numcond > length(default.color.list)){
      warning(paste0("WARNING: more than ",default.color.list," possible values for condition variable! Cannot select contrasting colors! Recommending setting colors manually using the group.RGB variable!"))
      color<- colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond))
    } else {
      color <- t(col2rgb(default.color.list[1:numcond]))
    }
    group.RGB <- sapply(1:nrow(color), function(i){
      paste0(floor(color[i,]),collapse=",")
    })
  }
  
  keep.features <- rep(TRUE,length(nrow(estimates)))
  res.data <- fData(jscs)
  
  #filtering:
  if(only.testable){
    keep.features <- keep.features & res.data$testable
  }
  sig.features <- which(fData(jscs)$padjust < FDR.threshold)
  gene.list <- unique(as.character(fData(jscs)$geneID[sig.features]))
  if(only.with.sig.gene){
    keep.features <- keep.features & fData(jscs)$geneID %in% gene.list
  }
  if(only.sig){
    keep.features <- keep.features & res.data$padjust < FDR.threshold
  }
  if(! plot.exons){
    keep.features <- keep.features & res.data$featureType != "exonic_part"
  }
  if(! plot.junctions){
    keep.features <- keep.features & (res.data$featureType != "splice_site" & res.data$featureType != "novel_splice_site")
  }
  if(! plot.novel.junctions){
    keep.features <- keep.features & res.data$featureType != "novel_splice_site"
  }
  
  res.data <- res.data[keep.features,,drop=FALSE]
  estimates <- estimates[keep.features,,drop=FALSE]
  
  chrom <- res.data$chr
  chromStart <- res.data$start
  chromEnd <- res.data$end
  strand <- res.data$strand
  featureType <- res.data$featureType
  
  if(use.score){
    score <- sapply(res.data$padjust, function(p){
      if(p > 0.05){
        0
      } else {
        floor((0.05 - p) * 18000) + 100
      }
    })
  } else {
    score <- rep(1000,length(chrom))
  }
  
  featureID <- rep("",nrow(estimates))
  if(includeGeneID & includeLocusID){
    featureID <- paste0("", res.data$featureID,"_")
  } else if(includeGeneID){
    featureID <- paste0("", res.data$geneID,"_")
  } else if(includeLocusID){
    featureID <- paste0("", res.data$countbinID,"")
  }
  
  idx <- order( rep(1:nrow(estimates), ncol(estimates)) )
  
  if(includeGroupID){
    out.featureNames <- unlist(lapply(1:ncol(estimates), function(i){
      paste0(estimate.names[i],":", featureID, "(", 
         sprintf(paste0("%.",count.digits,"f"), estimates[,i])
         ,")")
    }))[idx]
  } else {
    out.featureNames <- unlist(lapply(1:ncol(estimates), function(i){
      paste0(featureID, "(", 
         sprintf(paste0("%.",count.digits,"f"), estimates[,i]) ,
         ")")
    }))[idx]
  }
  
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
    verbose = verbose,
    use.gzip = use.gzip
  )
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
                               use.gzip = TRUE,
                               verbose = TRUE){ #fullFeatureID includeFeatureID
  res.data <- fData(jscs)
  
  keep.features <- rep(TRUE, nrow(res.data))
  if(only.testable){
    keep.features <- keep.features & res.data$testable
  }
  if(only.sig){
    keep.features <- keep.features & res.data$padjust < FDR.threshold
  }
  if(! plot.exons){
    keep.features <- keep.features & res.data$featureType != "exonic_part"
  }
  if(! plot.junctions){
    keep.features <- keep.features & (res.data$featureType != "splice_site" & res.data$featureType != "novel_splice_site")
  }
  if(! plot.novel.junctions){
    keep.features <- keep.features & res.data$featureType != "novel_splice_site"
  }
  res.data <- res.data[keep.features,,drop=FALSE]
  
  chrom <- res.data$chr
  chromStart <- res.data$start
  chromEnd <- res.data$end
  strand <- res.data$strand
  featureType <- res.data$featureType
  
  if(use.score){
    score <- sapply(res.data$padjust, function(p){
      if(p > 0.05){
        0
      } else {
        floor((0.05 - p) * 18000) + 100
      }
    })
  } else {
    score <- rep(1000,length(chrom))
  }
  
  if(only.sig){
    featureRGB <- rep(sig.RGB, length(chrom))
  } else {
    featureRGB <- ifelse(res.data$padjust < FDR.threshold, sig.RGB, nonsig.RGB)
  }
  
  if(includeLocusID & includeGeneID){
    featureID <- res.data$featureID
  } else if(includeGeneID){
    featureID <- res.data$geneID
  } else if(includeLocusID) {
    featureID <- res.data$countbinID
  } else {
    featureID <- rep("",nrow(res.data))
  }
  
  featureName <- paste0(
    featureID,
    "(",sprintf(paste0("%.",pval.digits,"f"), res.data$padjust),")"
  )
  
  
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
    verbose = verbose,
    use.gzip = use.gzip)
}

write.junction.bed.file <- function(file, trackLine = NULL, chrom, chromStart, chromEnd, strand, featureName, featureScore, featureRGB, featureType,
                                    output.format = c("BED","GTF","GFF3"), use.gzip = TRUE,
                                    verbose = TRUE){
  
  output.format <- match.arg(output.format)
  
  if(output.format == "BED"){
    if(use.gzip){
      gzf <- gzfile(paste0(file,""),"w")
    } else {
      gzf <- file(paste0(file,""),"w")
    }
    if(! is.null(trackLine)){
      write.table(trackLine,gzf, quote=FALSE,col.names=FALSE,row.names=FALSE)
    }
    if(length(chrom) > 0){
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
      )
      write.table(out.bed, gzf, quote=FALSE, col.names=FALSE,row.names=FALSE,sep='\t')
    }
    close(gzf)
  } else {
    write.junction.bed.file.GTFCONVERT(file = file, trackLine = trackLine, output.format = output.format, 
                                       chrom = chrom, chromStart = chromStart, chromEnd = chromEnd, strand = strand, featureName = featureName, featureRGB = featureRGB, featureType = featureType,
                                       verbose = verbose, use.gzip = use.gzip)
  }
}


write.junction.bed.file.GTFCONVERT <- function(file, trackLine = NULL, 
                                               output.format = c("GTF","GFF3"),
                                               chrom, chromStart, chromEnd, strand, featureName, featureScore, featureRGB, featureType,
                                               verbose = TRUE, use.gzip = use.gzip){
  
  output.format <- match.arg(output.format)
  in.data <- data.frame(chrom = chrom, chromStart = chromStart, chromEnd = chromEnd, strand = strand, 
                        featureName=featureName, featureScore = featureScore, featureRGB = featureRGB, featureType = featureType)
  
  

    if(use.gzip){
      gzf <- gzfile(paste0(file,""),"w")
    } else {
      gzf <- file(paste0(file,""),"w")
    }
  if(! is.null(trackLine)){
    write.table(trackLine,gzf, quote=FALSE,col.names=FALSE,row.names=FALSE)
  }
  
  if(nrow(in.data) > 0){
    out.GTF <- do.call(rbind.data.frame, lapply( 1:nrow(in.data), function(i){
      group <- if(output.format == "GTF"){
        paste0("gene_id ",featureName[i],"; transcript_id ",featureName[i])
      } else if(output.format == "GFF3"){
        paste0("ID=",featureName[i],";Name=",featureName[i])
      }
      if(featureType[i] == "exonic_part"){
        data.frame(chrom = chrom[i], source = "JunctionSeq", feature = "exon", 
                   start = chromStart[i]+1, end = chromEnd[i], score = featureScore[i], strand = strand[i],
                   frame = ".", group = group)
      } else {
        data.frame(chrom = rep(chrom[i],2), source = rep("JunctionSeq",2), feature = rep("exon",2), 
                   start = c(chromStart[i], chromEnd[i]+1), end = c(chromStart[i],chromEnd[i]+1), score = rep(featureScore[i],2), strand = rep(strand[i],2),
                   frame = rep(".",2), group = rep(group,2))
      }
    }))
    write.table(out.GTF, gzf, quote=FALSE, col.names=FALSE,row.names=FALSE,sep='\t')
  }
  
  close(gzf)
}

write.simple.table.gz <- function(d, file, use.gzip = TRUE, ...){
   if(use.gzip){
      gzf <- gzfile(paste0(file,".gz"),"w")
      write.table(d, gzf, ...)
      close(gzf)
   } else {
      write.table(d, file, ...)
   }
}



write.table.gz <- function(write.data, file, use.gzip = TRUE, sep = "	", quote=FALSE, row.names = FALSE, ...){
   if(use.gzip){
      gzf <- gzfile(paste0(file,".gz"),"w")
      write.fmt.table(write.data, gzf, row.names = row.names, quote = quote, sep = sep, ...)
      close(gzf)
   } else {
      write.fmt.table(write.data, file, row.names = row.names, quote = quote, sep = sep, ...)
   }
}

write.fmt.table <- function(d, file, row.names = FALSE, quote = FALSE, sep = sep, ...){
  if(row.names == FALSE){
    out.data <- d
  } else {
    if(row.names == TRUE){
      row.names.name <- "ROWNAME"
    } else {
      row.names.name <- row.names
    }
    
    out.data <- cbind.data.frame( row.names(d), d )
    
    names(out.data) <- c(row.names.name, names(d))
  }
  
  write.table(out.data, file = file, row.names = FALSE, quote = quote, sep = sep, ...)
}
