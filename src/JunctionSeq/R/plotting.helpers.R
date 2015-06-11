
########################################################
# UTILITY PLOTTING FUNCTIONS:
########################################################

plotTranscriptsOnly <- function(anno.data, geneID, debug.mode = TRUE, par.cex = 1, anno.cex.text = 1, ...)
{  
   #message("plotTranscriptsOnly cex = ",cex);

   draw.legend <- TRUE;

   rt.allExon <- which(anno.data$gene_id==geneID & anno.data$featureType == "exonic_part");
   rango <- 1:length(rt.allExon);
   

   
   #y.axis.title <- "";
   main.title <- paste0("Known transcripts for gene ",geneID);
   
   exoncol <- rep("#CCCCCC",length(rt.allExon))
   exonlty <- rep(1,length(exoncol));
   sig.feature <- rep(FALSE,length(exoncol));

   
   ################## DETERMINE THE LAYOUT OF THE PLOT DEPENDING OF THE OPTIONS THE USER PROVIDES ###########
      sub.allExon <- data.frame(start=anno.data$start[rt.allExon], end=anno.data$end[rt.allExon], chr=anno.data$chrom[rt.allExon], strand=anno.data$strand[rt.allExon])
      intervals<-(0:nrow(sub.allExon))/nrow(sub.allExon)


      
      rel.calc.min <- min(sub.allExon$start)
      rel.calc.max <- max(sub.allExon$end)

      rel <- (data.frame(sub.allExon$start, sub.allExon$end)) - rel.calc.min;
      rel <- rel/(rel.calc.max - rel.calc.min);

         transcripts <- sapply(sapply(anno.data$transcripts[rt.allExon],toString), function(x){strsplit(x, "+",fixed=TRUE)})
         trans <- Reduce(union, transcripts)
         if(length(trans) > 42){
            warning("This gene contains more than 42 transcripts annotated, only the first 42 will be plotted\n")
            trans <- trans[1:42];
         }
         mat <- 1:(min(length(trans), 42)) ## max support from transcripts is 45, which seems to be the max for the layout supported by graphics
         hei <- c(1, 1.5, rep(1.5, min(length(trans), 42)))

      layout(matrix(mat), heights=hei)
      
      p.values.labels <- rep("",length(rt.allExon));
   
      #########PLOT THE GENE MODEL:
      par(mar=c(0, 4, 0, 2))
      plot.new()
      par(cex = par.cex);
      # lines linking exons / splices to their column.
      segments(apply((rbind(rel[rango,2], rel[rango, 1])), 2, median), 0, apply(rbind(intervals[rango], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2)), 2, median), 1, col=exoncol, lty = exonlty, ...)
      
      par(mar=c(1.5, 4, 0, 2))
      drawGene.noSplices(rel.calc.min, rel.calc.max, tr.allExon=sub.allExon, par.cex = par.cex, anno.cex.text = anno.cex.text,...)
            for(i in 1:min(length(trans), 42)){
               logicexons <- sapply(transcripts, function(x){length(which(x==trans[i]))})
               tr <- data.frame(start = anno.data$start[rt.allExon][logicexons==1], end = anno.data$end[rt.allExon][logicexons==1], featureType = anno.data$featureType[rt.allExon][logicexons==1]);
               tr <- tr[tr$featureType == "exonic_part",]
               drawTranscript(rel.calc.min, rel.calc.max, tr=tr, rango=1:(length(anno.data$start[rt.allExon][logicexons==1])), exoncol=NULL, names=c(), trName=trans[i], par.cex = par.cex, anno.cex.text = anno.cex.text,sub.sig = data.frame(start=c(),end=c()))
            }
      
      axis1.main <- make.evenly.spaced.seq(rel.calc.min, rel.calc.max, 10);
      axis1.minor <- make.evenly.spaced.seq.minor(rel.calc.min, rel.calc.max, 10);
      
      axis(1,at=axis1.minor,labels=rep("",length(axis1.minor)),pos=0,lwd.ticks=0.2,padj=-0.7,tcl=-0.25, cex.axis = anno.cex.text, ...);
      axis(1,at=axis1.main,labels=axis1.main,pos=0,lwd.ticks=0.2,padj=-0.7,tcl=-1, cex.axis = anno.cex.text, ...);
}

########################################################
#FUNCTION TO MAKE THE AXIS OF THE VST VALUES
########################################################

makevstaxis <- function(min, ylimn, ecs, use.vst, use.log,plot.type, par.cex = 1, anno.cex.text = 1, draw.axis = 2, pos = 0, yAxisLabels.inExponentialForm = TRUE, ...)
{
   #message("makevstaxis cex = ",cex);

   if(use.vst){
   
  #  ylimn <- c(0.0156,39.89)

     range_scaled <- ylimn[2] - ylimn[1];

     decades_unscaled <- c(0,1);
     decades_scaled <- vst(decades_unscaled,ecs);
     #while(decades_scaled[1] > ylimn[1]){
     #  decades_unscaled <- c(decades_unscaled[1] / 10, decades_unscaled);
     #  decades_scaled <- vst(decades_unscaled,ecs);
     #}
     while(decades_scaled[length(decades_scaled)] < ylimn[2]){
       decades_unscaled <- c(decades_unscaled, decades_unscaled[length(decades_unscaled)] * 10);
       decades_scaled <- vst(decades_unscaled,ecs);
     }
     #first.decade.threshold <- (range_scaled * 0.025) + ylimn[1];
     #while(decades_scaled[1] < first.decade.threshold){
     #  decades_scaled <- decades_scaled[2:length(decades_scaled)];
     #  decades_unscaled <- decades_unscaled[2:length(decades_unscaled)];
     #}

     if(yAxisLabels.inExponentialForm & length(decades_unscaled) > 2){
       decades_labels <- c(expression(0),expression(1),
         sapply(decades_unscaled[-c(1,2)], function(D){
           exponent <- log10(D);
           substitute(10 ^ X, list(X = exponent));
         })
       );
     } else {
       decades_labels <- c(decades_unscaled);
     }
     
     decades_scaled <- vst(decades_unscaled,ecs);

     ticks_unscaled <- c();
     first.tick.threshold <- (range_scaled * 0.05) + ylimn[1];
     
     for(i in 1:length(decades_unscaled)){
       if(decades_scaled[i] >  first.tick.threshold & decades_unscaled[i] > 0) {
       ticks_unscaled <- c(ticks_unscaled,
                           decades_unscaled[i]*2,
                           decades_unscaled[i]*3,
                           decades_unscaled[i]*4,
                           decades_unscaled[i]*5,
                           decades_unscaled[i]*6,
                           decades_unscaled[i]*7,
                           decades_unscaled[i]*8,
                           decades_unscaled[i]*9
                           );
       }
     }
     ticks_scaled <- vst(ticks_unscaled,ecs);
   
     axis( draw.axis, at=ticks_scaled, labels=FALSE, las=2, pos=pos,tcl=-0.25, cex.axis = anno.cex.text, ...)
     axis( draw.axis, at=decades_scaled, labels=decades_labels, las=2, pos=pos,tcl=-0.5, cex.axis = anno.cex.text, ...)

#   minlog10 <- floor(log10(ylimn[1]));
#   maxlog10 <- ceiling(log10(ylimn[2]));
#   decades_unscaled <- 10^(minlog10:maxlog10);
#   min <- 1/ncol(count)
#   minlog10 <- floor( log10( 1/ncol(counts(ecs)) ) )
#   maxlog10 <- ceiling( log10( ylimn[2] ) )
#   ticks <- 10^seq( minlog10, maxlog10 )
#   decade_lengths <- ( vst(ticks, ecs)[ 2 : length(ticks) ] - vst(ticks, ecs)[ 1 : (length(ticks)-1) ] ) /
#      ( vst( ylimn[2], ecs) - vst( ylimn[1], ecs) )
#   ticks <- c( 0, ticks[ min( which( decade_lengths > .1 ) ) : length(ticks) ] )
#   ticks <- ticks[ min( which( decade_lengths > .1 ) ) : length(ticks) ]
#   axis( 2, at=vst(c(0, ticks), ecs), labels=c("",ticks), las=2, pos=0, ...)
#   for( i in minlog10 : (maxlog10-1) ) {
#      decade_length <- ( vst( 10^(i+1), ecs) - vst( 10^i, ecs) ) / ( vst( ylimn[2], ecs) - vst( ylimn[1], ecs) )
#      if( decade_length > .1 ) {
#         axis( 2, at = vst( 1:9 * 10^i, ecs), labels = FALSE, las=2, tcl=-.25, pos=0, ...)
#      }
#      if( decade_length > .4 & decade_length <= .6) {
#         axis( 2, at = vst( c(2,3,5) * 10^i, ecs ), labels = ( c(2,3,5) * 10^i ), las=2, tcl=-.25, pos=0, ...)
#      } else if( decade_length > .6 ) {
#         axis( 2, at = vst( c(1.5,2:9) * 10^i, ecs ), labels = ( c(1.5,2:9) * 10^i ), las=2, tcl=-.25, pos=0, ...)
#      }
#   }
#   axis( 2, at=vst(0, ecs), labels=FALSE, las=2, pos=0, ...)
  } else if(use.log) {
    if(plot.type == "rawCounts"){
      ticks <- rep(c(2,3,4,5,6,7,8,9),ceiling(ylimn[2])) * rep(10^(0:ceiling(ylimn[2]-1)),each=8)
      logticks <- log10(ticks);
      axis(draw.axis,at=c(INTERNAL.NINF.VALUE,0),labels=c(0,""), las=2, tcl=-0.5, pos=pos, cex.axis = anno.cex.text, ...)
      qorts.axis.break(draw.axis,breakpos=(INTERNAL.NINF.VALUE / 2),pos=pos, cex = anno.cex.text, yw = INTERNAL.NINF.VALUE * 0.3,...)
      
      #axis( 2, at=ticks_scaled, labels=FALSE, las=2, pos=0,tcl=0.25, cex.axis = anno.cex.text, ...)
      #axis( 2, at=decades_scaled, labels=decades_unscaled, las=2, pos=0,tcl=0.5, cex.axis = anno.cex.text, ...)
      if(yAxisLabels.inExponentialForm & ceiling(ylimn[2]) > 2){
        decade.labels <- c(expression(1),expression(10),sapply(2:ceiling(ylimn[2]), function(X){
          substitute(10 ^ x, list(x = X));
        }))
      } else {
        decade.labels <- (10^(0:ceiling(ylimn[2])))
      }
      
      axis(draw.axis,   at=logticks ,labels=FALSE, las=2, pos=pos, tcl=-0.25, cex.axis = anno.cex.text, ...)
      axis(draw.axis,   at=0:ceiling(ylimn[2]),labels=decade.labels, las=2, tcl=-0.5, pos=pos, cex.axis = anno.cex.text, ...)
    } else {
      ticks <- rep(c(2,3,4,5,6,7,8,9),ceiling(ylimn[2])) * rep(10^(0:ceiling(ylimn[2]-1)),each=8)
      logticks <- log10(ticks);
      axis(draw.axis,   at=c(INTERNAL.NINF.VALUE,0),labels=c("0",""),      las=2, tcl=-0.5,  pos=pos, cex.axis = anno.cex.text, ...)
      qorts.axis.break(draw.axis,breakpos=(INTERNAL.NINF.VALUE / 2),pos=pos, cex = anno.cex.text, yw = INTERNAL.NINF.VALUE * 0.3,...)
      #axis(2,   at = ((1:9)/10)*INTERNAL.NINF.VALUE, labels=FALSE, las=2, tcl=0.25, pos=0, cex.axis = anno.cex.text, ...);
      
      if(yAxisLabels.inExponentialForm & ceiling(ylimn[2]) > 2){
        decade.labels <- c(expression(1),expression(10),sapply(2:ceiling(ylimn[2]), function(X){
          substitute(10 ^ x, list(x = X));
        }))
      } else {
        decade.labels <- (10^(0:ceiling(ylimn[2])))
      }
      
      axis(draw.axis,   at=logticks ,labels=FALSE, las=2, pos=pos, tcl=-0.25, cex.axis = anno.cex.text, ...)
      axis(draw.axis,   at=0:ceiling(ylimn[2]),labels= decade.labels, las=2, tcl=-0.5, pos=pos, cex.axis = anno.cex.text, ...)
    }
  } else {
    axis(draw.axis, las=2, tcl=0.5, pos=pos, cex.axis = anno.cex.text, ...);
  }
}

makeGeneLevelAxis <- function(ylimn, geneColumn.TOP, use.vst, use.log, plot.type, par.cex = 1, anno.cex.text = 1, draw.axis = 2, pos = 1.04, geneLevel.rescaleFactor, yAxisLabels.inExponentialForm = TRUE, ...){
  if(use.vst){
    stop("Gene-Level data with vst-transformations is not currently supported!");
  } else if(use.log) {
    if(plot.type == "rawCounts"){
      ticks <- rep(c(2,3,4,5,6,7,8,9),ceiling(ylimn[2])) * rep(10^(0:ceiling(ylimn[2]-1)),each=8)
      logticks <- log10(ticks);
      axis(draw.axis,at=c(INTERNAL.NINF.VALUE,0),labels=c(0,""), las=2, tcl=-0.5, pos=pos, cex.axis = anno.cex.text, ...)
      qorts.axis.break(draw.axis,breakpos=(INTERNAL.NINF.VALUE / 2),pos=pos, cex = anno.cex.text, yw = INTERNAL.NINF.VALUE * 0.3,...)
      
      if(yAxisLabels.inExponentialForm & ceiling(ylimn[2]) > 2){
        decade.labels <- c(expression(1),expression(10),sapply(2:ceiling(ylimn[2]), function(X){
          substitute(10 ^ x, list(x = X));
        }))
      } else {
        decade.labels <- (10^(0:ceiling(ylimn[2])))
      }
      
      #axis( 2, at=ticks_scaled, labels=FALSE, las=2, pos=0,tcl=0.25, cex.axis = anno.cex.text, ...)
      #axis( 2, at=decades_scaled, labels=decades_unscaled, las=2, pos=0,tcl=0.5, cex.axis = anno.cex.text, ...)
      logticks <- logticks / geneLevel.rescaleFactor;
      decades <- 0:ceiling(ylimn[2]) / geneLevel.rescaleFactor;
      
      keep.ticks <- logticks <= geneColumn.TOP;
      keep.decades <- decades <= geneColumn.TOP;
      
      axis(draw.axis,   at = logticks[keep.ticks],labels=FALSE, las=2, pos=pos, tcl=-0.25, cex.axis = anno.cex.text, ...)
      axis(draw.axis,   at = decades[keep.decades], labels=decade.labels[keep.decades], las=2, tcl=-0.5, pos=pos, cex.axis = anno.cex.text, ...)
    } else {
      ticks <- rep(c(2,3,4,5,6,7,8,9),ceiling(ylimn[2])) * rep(10^(0:ceiling(ylimn[2]-1)),each=8)
      logticks <- log10(ticks);
      axis(draw.axis,   at=c(INTERNAL.NINF.VALUE,0),labels=c("0",""),      las=2, tcl=-0.5,  pos=pos, cex.axis = anno.cex.text, ...)
      qorts.axis.break(draw.axis,breakpos=(INTERNAL.NINF.VALUE / 2),pos=pos, cex = anno.cex.text, yw = INTERNAL.NINF.VALUE * 0.3,...)
      #axis(2,   at = ((1:9)/10)*INTERNAL.NINF.VALUE, labels=FALSE, las=2, tcl=0.25, pos=0, cex.axis = anno.cex.text, ...);
      
      if(yAxisLabels.inExponentialForm & ceiling(ylimn[2]) > 2){
        decade.labels <- c(expression(1),expression(10),sapply(2:ceiling(ylimn[2]), function(X){
          substitute(10 ^ x, list(x = X));
        }))
      } else {
        decade.labels <- (10^(0:ceiling(ylimn[2])))
      }
      logticks <- logticks / geneLevel.rescaleFactor;
      decades <- 0:ceiling(ylimn[2]) / geneLevel.rescaleFactor;
      
      keep.ticks <- logticks <= geneColumn.TOP;
      keep.decades <- decades <= geneColumn.TOP;
      
      axis(draw.axis,   at = logticks[keep.ticks],labels=FALSE, las=2, pos=pos, tcl=-0.25, cex.axis = anno.cex.text, ...)
      axis(draw.axis,   at = decades[keep.decades], labels=decade.labels[keep.decades], las=2, tcl=-0.5, pos=pos, cex.axis = anno.cex.text, ...)
    }
  } else {
    axis(draw.axis, las=2, tcl=0.5, pos=pos, cex.axis = anno.cex.text, ...);
  }
}

#######################
#FUNCTION TO DRAW THE EXPRESSION PLOTS:
#######################
drawPlot <- function(matr, ylimn,ecs, intervals, rango, fitExpToVar, numexons, textAxis, rt, color.count, 
                     colorlines,countbinIDs,use.vst, use.log,plot.type,main.title,draw.legend,color.key,
                     condition.names,p.values=NULL,draw.p.values=FALSE,plot.lwd = 1,axes.lwd = 1, 
                     anno.lwd =1, par.cex = 1, anno.cex.text = 1, anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                     fit.countbin.names = TRUE, fit.pvals = TRUE, stagger.pvals = TRUE,debug.mode = FALSE,
                     plot.gene.level.expression = FALSE, geneCount = NULL, color.geneCount = NULL, 
                     yAxisLabels.inExponentialForm = TRUE,
                     autofit.legend = TRUE, italicize.label = NULL, condition.legend.text = condition.legend.text,
                     rel = NULL, annolink.col = NULL, exonlty = NULL,
                     graph.margins = c(2,3,3,2),
                     ...)
{
   plot.junction.ids.vertically <- NULL;
   plot.junction.ids.vertically.threshold <- Inf;
   
   if(debug.mode) message("> Step 6.1");
   
   #message("drawPlot cex = ",cex);
   par(mar=graph.margins, cex = par.cex);
   plot.new();
   par(cex = par.cex);
   
   #message(paste("ylimn:",ylimn));
   plot.window(xlim=c(0, 1), ylim=ylimn);
   
   if(debug.mode) message("> Step 6.2");
   makevstaxis(1/ncol(matr), ylimn, ecs,use.vst=use.vst, use.log=use.log,plot.type=plot.type, lwd = axes.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, ...)
   title(main=main.title, cex.main = anno.cex.main, ...);
   if(debug.mode) message("> Step 6.3");
   
   if(plot.type == "logRawCounts")  segments(0,INTERNAL.NINF.VALUE,2,INTERNAL.NINF.VALUE,lty="dotted",col="black", lwd = plot.lwd,...); 
   if((! use.vst) & use.log ) segments(0,INTERNAL.NINF.VALUE,2,INTERNAL.NINF.VALUE,lty="dotted",col="black", lwd = plot.lwd,...); 
   segments(0,0,2,0,lty="dotted",col="black", lwd = plot.lwd,...);
   junctionColumns.RIGHT <- (intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2))[length(rango)];
   if(debug.mode) message("> Step 6.4");
   if(plot.gene.level.expression){
     geneColumn.RIGHT <- par("usr")[2];
     
     if(autofit.legend & draw.legend){
       legend.rect <- legend(geneColumn.RIGHT,  par("usr")[4],fill=color.key[names(condition.legend.text)],legend=condition.legend.text[condition.names], cex = anno.cex.text, bg = "white", box.lwd = axes.lwd, xjust = 1, yjust = 1, xpd = NA, plot=FALSE,...)$rect;
       #note: rect$w = width, rect$h = height, rect$left = left coord, rect$top = top coord.
       devlimits <- device.limits();
       junctionColumns.RIGHT.RESCALED <- devlimits[2] - legend.rect$w;
       if(junctionColumns.RIGHT.RESCALED < junctionColumns.RIGHT){
         if(debug.mode) message("Shifting plots to autofit the legend box.");
         intervals <- intervals * (junctionColumns.RIGHT.RESCALED / junctionColumns.RIGHT);
         junctionColumns.RIGHT <- junctionColumns.RIGHT.RESCALED;
       }
       geneColumn.TOP.PCT <- 0.98;
       geneColumn.TOP <- (legend.rect$top - legend.rect$h) * geneColumn.TOP.PCT;
     } else {
       geneColumn.TOP.PCT <- 0.9;
       geneColumn.TOP   <- par("usr")[4] * geneColumn.TOP.PCT;
     }
     geneColumn.LEFT  <- ((geneColumn.RIGHT - junctionColumns.RIGHT) * 0.25) + junctionColumns.RIGHT;
     geneColumn.MID   <- ((geneColumn.RIGHT - geneColumn.LEFT) *0.5 + geneColumn.LEFT);
     #rect(1, par("usr")[3], 1.01, par("usr")[4], xpd=NA, lwd = axes.lwd,xpd=NA, col="white",border="white",...);
     #end.plotting <- par("usr")[2];
     geneColumn.YMAX <- geneColumn.TOP * 0.96;
     
     rect(junctionColumns.RIGHT,    par("usr")[3], geneColumn.LEFT,    par("usr")[4], xpd=NA, lwd = axes.lwd,xpd=NA, border="white",col="white",...);
     rect(0,    par("usr")[3], junctionColumns.RIGHT,    par("usr")[4], xpd=NA, lwd = axes.lwd,xpd=NA,...);
     rect(geneColumn.LEFT, par("usr")[3], geneColumn.RIGHT, geneColumn.TOP, xpd=NA, lwd = axes.lwd,xpd=NA,...);
     if(draw.legend){
         #legend(geneColumn.LEFT,  par("usr")[4],fill=color.key,legend=condition.names, cex = anno.cex.text, bg = "white", box.lwd = axes.lwd, xjust = 0, yjust = 1, xpd = NA,...);
         #legend(geneColumn.RIGHT,  par("usr")[4],fill=color.key,legend=condition.names, cex = anno.cex.text, bg = "white", box.lwd = axes.lwd, xjust = 1, yjust = 1, xpd = NA,...);
         if(autofit.legend){
           #legend(geneColumn.RIGHT,  par("usr")[4],fill=color.key,legend=condition.names, cex = anno.cex.text, box.lwd = axes.lwd, xjust = 1, yjust = 1, xpd = NA, box.col = "transparent", bg = "transparent",...);
           legend(junctionColumns.RIGHT,  par("usr")[4],fill=color.key[names(condition.legend.text)],legend=condition.legend.text[condition.names], cex = anno.cex.text, box.lwd = axes.lwd, xjust = 0, yjust = 1, xpd = NA, box.col = "transparent", bg = "transparent",...);
         } else {
           legend(geneColumn.RIGHT,  par("usr")[4],fill=color.key[names(condition.legend.text)],legend=condition.legend.text[condition.names], cex = anno.cex.text, bg = "white", box.lwd = axes.lwd, xjust = 1, yjust = 1, xpd = NA, ...);
         }
     }
   } else {
     if(draw.legend){
        legend("topright",fill=color.key[names(condition.legend.text)],legend=condition.legend.text[condition.names], cex = anno.cex.text, bg = "transparent", box.lwd = axes.lwd, xpd = NA,...);
     }
   }
   

   
   if(debug.mode) message("> Step 6.5");
   
   #intervals<-(0:nrow(matr))/nrow(matr)
   middle <- apply(cbind(intervals[rango], (intervals[rango+1]-((intervals[rango+1])-intervals[rango])*0.2)), 1, median)
   matr <- rbind(matr, NA)
   j <- 1:ncol(matr)
   
   #abline(h = 0,lty="dotted",col="black", lwd = plot.lwd,...);
   segments(intervals[rango],matr[rango,j], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), matr[rango,j], col=color.count, lwd = plot.lwd,...)  #### line with the y level
   segments(intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), matr[rango,j], intervals[rango+1], matr[rango+1,j], col=color.count, lty="dotted", lwd =plot.lwd,...)  #### line joining the y levels
   abline(v=middle[rango], lty="dotted", col=colorlines, lwd = plot.lwd)
   mtext(textAxis, side=2, adj=0.5, line=1.5, outer=FALSE, cex = anno.cex.text * par("cex"),...)
   
   if(debug.mode) message("> Step 6.6");
   
   cex.countbinIDs <- anno.cex.axis;
   srt.countbinIDs <- 0;

   if(is.null(plot.junction.ids.vertically)){
     if(length(rt) > plot.junction.ids.vertically.threshold){
       plot.junction.ids.vertically <- TRUE;
     } else {
       plot.junction.ids.vertically <- FALSE;
     }
   }
   if(debug.mode) message("> Step 6.7");
   srt.countbinIDs <- 0;
   cex.countbinIDs <- anno.cex.axis;
   adj.countbinIDs <- c(0.5,0.5);
   tcl.countbinIDs <- -0.4
   if(fit.countbin.names){
     bin.width <- abs(intervals[1] - intervals[2]);
     cex.countbinIDs <- fit.character.vector(countbinIDs, min.width = 0.6 * bin.width, max.width = 0.9 * bin.width, max.width.per.char = 0.4 * bin.width)
     cex.countbinIDs <- min(cex.countbinIDs, anno.cex.axis);
     #message("cex.countbinIDs = ", cex.countbinIDs);
     if(cex.countbinIDs < anno.cex.axis / 4){
       srt.countbinIDs <- 90;
       #cex.countbinIDs <- anno.cex.axis / 2;
       cex.countbinIDs <- cex.countbinIDs * 3;
       tcl.countbinIDs <- -0.1;
     } else if(cex.countbinIDs < anno.cex.axis / 2){
       srt.countbinIDs <- 45;
       #cex.countbinIDs <- anno.cex.axis / 2;
       cex.countbinIDs <- cex.countbinIDs * 2;
       tcl.countbinIDs <- -0.25;
       #adj.countbinIDs[1] <- 1;
     }
   }
   if(is.null(italicize.label)){
     axis.font <- NULL;
   } else {
     axis.font <- ifelse(italicize.label,3,2);
   }
   axis.label.floor <- JS.axis(1, at=middle[1:length(rt)], labels=countbinIDs, tcl = tcl.countbinIDs,  lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, srt = srt.countbinIDs, font = axis.font, adj = adj.countbinIDs,...)# ,

   #if(plot.junction.ids.vertically){
     
     #JS.axis(1, at=middle[1:length(rt)], labels=countbinIDs, tcl = -0.5,  lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 3, mgp = c(3,0.55,0),...)# ,
     
     #if(is.null(italicize.label)){
     #  axis(1, at=middle[1:length(rt)], labels=countbinIDs, tcl = -0.25, lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 3,  line = 0, mgp = c(3,0.3,0),padj=0.5,hadj=1, ...)  #
     #} else {
     #  axis(1, at=middle[1:length(rt)][italicize.label], labels=countbinIDs[italicize.label], font = 3, tcl = -0.25, lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 3,  line = 0, mgp = c(3,0.3,0),padj=0.5,hadj=1, ...)  #
     #  axis(1, at=middle[1:length(rt)][! italicize.label], labels=countbinIDs[! italicize.label], font = 2, tcl = -0.25, lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 3,  line = 0, mgp = c(3,0.3,0),padj=0.5,hadj=1, ...)  #
     #}
   #} else {

     
     #if(is.null(italicize.label)){
     #  axis(1, at=middle[1:length(rt)], labels=countbinIDs, tcl = -0.5,  lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 0, mgp = c(3,0.55,0),padj=0.5, hadj=0.5,...)# ,
     #} else {
     #  axis(1, at=middle[1:length(rt)][italicize.label], labels=countbinIDs[italicize.label], font = 3, tcl = -0.5,  lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 0, mgp = c(3,0.55,0),padj=0.5, hadj=0.5,...)# ,
     #  axis(1, at=middle[1:length(rt)][! italicize.label], labels=countbinIDs[! italicize.label], font = 2, tcl = -0.5,  lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 0, mgp = c(3,0.55,0),padj=0.5, hadj=0.5,...)# ,
     #}
   #}
   if(debug.mode) message("> Step 6.8");
   
   #par usr: (x1, x2, y1, y2)
   #rect(xleft, ybottom, xright, ytop
   
   if(debug.mode) message("> Step 6.9");
   if(any(p.values != "")){
     if(draw.p.values){
        bin.width <- abs(intervals[1] - intervals[2]);
        if(fit.pvals){  
          cex.pvalues <- fit.character.vector(p.values, min.width = 1 * bin.width, max.width = 2 * bin.width, max.width.per.char = 1 * bin.width)
          cex.pvalues <- min(cex.pvalues, anno.cex.text);
        } else {
          cex.pvalues <- anno.cex.text;
        }
        pval.height <- max(strheight(p.values, cex = cex.pvalues));
        if(debug.mode) message("> Step 6.9b");
        if(stagger.pvals){
          pval.y <- rep(c(ylimn[2], ylimn[2] - pval.height, ylimn[2] - (2 * pval.height)), length(middle[rango]))[rango]
        } else {
          pval.y <- rep(ylimn[2],rango)
        }
        
        if(strwidth(p.values[1], cex = cex.pvalues) > bin.width  * 0.9){
          text(intervals[1], pval.y[1], labels = p.values[1], col = colorlines[1], cex = cex.pvalues, adj = 0, ...);
        } else {
          text(middle[1], pval.y[1], labels = p.values[1], col = colorlines[1], cex = cex.pvalues,  ...);
        }
        if(strwidth(p.values[length(p.values)], cex = cex.pvalues) > bin.width * 0.9){
          text(junctionColumns.RIGHT, pval.y[length(middle)], labels = p.values[length(middle)], col = colorlines[length(middle)], cex = cex.pvalues, adj = 1, ...);
        } else {
          text(middle[length(middle)], pval.y[length(middle)], labels = p.values[length(middle)], col = colorlines[length(middle)], cex = cex.pvalues, ...);
        }
        text(middle[2:(length(middle)-1)],pval.y[2:(length(middle)-1)],labels=p.values[2:(length(middle)-1)],col=colorlines[2:(length(middle)-1)], cex = cex.pvalues,...);
     }
   }
   if(debug.mode) message("> Step 6.10 (anno.cex.text = ",anno.cex.text,")");
   
   if(plot.gene.level.expression){
     #makevstaxis(1/ncol(matr), c(0,max(geneCount)), ecs,use.vst=use.vst, use.log=use.log,plot.type=plot.type, lwd = axes.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, draw.axis=4, pos = par("usr")[2], ...)
     geneLevel.max <- max(geneCount, na.rm=TRUE);
     #matr.max <- max(matr, na.rm=TRUE);
     
     geneLevel.rescaleFactor <- (geneLevel.max) / (geneColumn.YMAX);
     
     geneLevelRescaled <- ifelse(geneCount < 0, geneCount, geneCount / geneLevel.rescaleFactor);
     points(rep(geneColumn.MID,length(geneLevelRescaled)), geneLevelRescaled, col = color.geneCount, pch = 4, ...);
     
     segments(geneColumn.LEFT, geneLevelRescaled, geneColumn.RIGHT, geneLevelRescaled, col = color.geneCount, lwd = plot.lwd, ... );
     segments(junctionColumns.RIGHT,matr[rango[length(rango)],j],geneColumn.LEFT, geneLevelRescaled,col = color.geneCount,lwd = plot.lwd, lty = "dotted",...);
     
     makeGeneLevelAxis( c(0,geneColumn.TOP * geneLevel.rescaleFactor), geneColumn.TOP, use.vst=use.vst, use.log=use.log, plot.type=plot.type, lwd = axes.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, draw.axis=4, pos = geneColumn.RIGHT, geneLevel.rescaleFactor = geneLevel.rescaleFactor, mgp = c(3,0.5,0), ...);

     # srt.countbinIDs <- 0;
     # cex.countbinIDs <- anno.cex.axis;
     #adj.countbinIDs <- c(0.5,0.5);
     #tcl.countbinIDs <- -0.4
     JS.axis(1, at=geneColumn.MID, labels="GENE", tcl =  tcl.countbinIDs,  lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = min(anno.cex.axis, cex.countbinIDs * 1.5), srt = 0, font = 2, adj = adj.countbinIDs,...)# ,

     #axis(1, at=geneColumn.MID, labels="GENE", tcl = -0.5,  lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 0, mgp = c(3,0.55,0),padj=0.5, hadj=0.5,...)#
   } else {
     rect(0, par("usr")[3], par("usr")[2], par("usr")[4], xpd=NA, lwd = axes.lwd,xpd=NA,...);
   }
   

   if(debug.mode) message("> Step 6.11");
   #if(draw.box){
   #   box(lwd = plot.lwd, ...);
   #}
   #Removed from this section and replaced with a separate frame.
   #  It worked fine here as part of this plot, but 
   #segments(
   #        apply((rbind(rel[rango,2], rel[rango, 1])), 2, median), 
   #        device.limits()[3],
   #        middle, 
   #        axis.label.floor, 
   #        col = annolink.col, lty = exonlty, lwd = plot.lwd, cex = anno.cex.text,cex.axis=anno.cex.axis, cex.main=anno.cex.main, xpd=NA, ...) #col=colorlinesB,...)

   return(intervals);
}


#########################
#FUNCTION TO DRAW THE GENE MODELS AND TRANSCRIPTS:
#########################
drawGene <- function(minx, maxx, tr, tr.allExon, tr.allJunction, rango, rescale.iv = NULL, exoncol=NULL,allExon.exonCol=NULL, names, 
                     trName, exonlty, plot.lwd = 1, anno.lwd = 1, show.strand.arrows = 0, par.cex = 1, anno.cex.text = 1, 
                     arrows.length = 0.5,draw.untestable.annotation = TRUE,
                     draw.start.end.sites = TRUE, startSites = NULL, endSites = NULL,
                     ...)
{
   #message("drawGene cex = ",cex);
   #par(cex = par.cex);
   xpd <- NA
   plot.new()
   plot.window(xlim=c(minx, maxx), ylim=c(0, 1))
   ymax <- par("usr")[4];
   ymin <- par("usr")[3];
   
   rect.floor <- 0.05;
   rect.ceil  <- 0.75;
   startSite.marker.width <- 0.025;
   if(draw.start.end.sites){
     rect.floor <- 0.15;
     startSites.floor <- 0.05;
   }
   rect.mid   <- ((rect.ceil - rect.floor)/2) + rect.floor

   if(! is.null(rescale.iv)){
     tr.allExon$start <- (rescale.coords(tr.allExon$start,rescale.iv) * (maxx-minx)) + minx;
     tr.allExon$end   <- (rescale.coords(tr.allExon$end,rescale.iv) * (maxx-minx)) + minx;
     tr$start <- (rescale.coords(tr$start,rescale.iv) * (maxx-minx)) + minx;
     tr$end   <- (rescale.coords(tr$end,rescale.iv) * (maxx-minx)) + minx;
     
     #if(draw.untestable.annotation){
     if(nrow(tr.allJunction) > 0){
       tr.allJunction$start <- (rescale.coords(tr.allJunction$start,rescale.iv) * (maxx-minx)) + minx;
       tr.allJunction$end   <- (rescale.coords(tr.allJunction$end,rescale.iv) * (maxx-minx)) + minx;
     }
     
     if(draw.start.end.sites){
       startSites <- (rescale.coords(startSites,rescale.iv) * (maxx-minx)) + minx;
       endSites <- (rescale.coords(endSites,rescale.iv) * (maxx-minx)) + minx;
     }
   }
   
   #print("TR:");
   #print(tr);
   #print("TR.ALLEXON:");
   #print(tr.allExon);
   
   
   rect(tr.allExon$start, rect.floor, tr.allExon$end, rect.ceil, lty = 1, lwd = anno.lwd / 2,col=tr.allExon$fillColor,...)
   lines(c(minx,maxx),c(rect.mid,rect.mid),lty= 1, lwd = anno.lwd, ...);
   
   if(draw.start.end.sites){
     segments(startSites, startSites.floor,startSites,rect.floor, lwd = plot.lwd / 2,xpd=NA,...);
     segments(endSites,   startSites.floor,endSites,rect.floor, lwd = plot.lwd / 2,xpd=NA,...);
     
     segments(startSites, startSites.floor,startSites + (maxx-minx) * startSite.marker.width, startSites.floor, lwd = plot.lwd / 2,xpd=xpd,...);
     segments(endSites,   startSites.floor,endSites   - (maxx-minx) * startSite.marker.width, startSites.floor, lwd = plot.lwd / 2,xpd=xpd,...);
   }
   
   if(show.strand.arrows > 0){
      #if(is.null(strand)){
      #   message("Warning: Cannot draw strand arrows, strandedness is NULL! Skipping.");
      #} else {
         strand <- tr.allExon$strand[1];
         
         #even.spacing <- ((1:(show.strand.arrows)) / (show.strand.arrows+1)) * (maxx - minx) + minx;
         if(strand == "+"){
            #arrow.char <- ">";
            even.spacing <- ((0:(show.strand.arrows - 1)) / (show.strand.arrows)) * (maxx - minx) + minx
            arrows( even.spacing[-(show.strand.arrows)] ,rep(rect.mid,show.strand.arrows - 1), even.spacing[-1]  , rep(rect.mid,show.strand.arrows - 1), lwd = anno.lwd, length = arrows.length);
         } else if(strand == "-"){
            #arrow.char <- "<";
            even.spacing <- ((1:(show.strand.arrows)) / (show.strand.arrows)) * (maxx - minx) + minx
            arrows( even.spacing[-1] ,rep(rect.mid,show.strand.arrows - 1), even.spacing[-(show.strand.arrows)]  , rep(rect.mid,show.strand.arrows - 1), lwd = anno.lwd, length = arrows.length);
         }
         #points(even.spacing, rep(0.4,show.strand.arrows), pch = arrow.char, cex = anno.cex.text, ...);
      #}
   }

#   if(draw.untestable.annotation){
#      tr.splice <- tr.allJunction[, ];
#      start.splice <- tr.splice$start
#      end.splice <- tr.splice$end
#      middle.splice <- apply(rbind(start.splice, end.splice), 2, median)
#      colors.splice <- exoncol[! tr$is.exon]
#      segments(start.splice, 0.75, middle.splice, ymax, col=colors.splice, lty = exonlty[! tr$is.exon], lwd = plot.lwd,xpd=FALSE,...);
#      segments(middle.splice, ymax, end.splice, 0.75, col=colors.splice, lty = exonlty[! tr$is.exon], lwd = plot.lwd,xpd=FALSE,...);
#   }

   #if(any( ! tr$is.exon)){
   if(nrow(tr.allJunction) > 0){
      tr.splice <- tr.allJunction;
      #if(! draw.untestable.annotation){
      #  tr.splice <- tr.splice[(! tr.splice$is.testable) & (tr.splice$is.untestable),];
      #  message("############ REMOVING UNTESTED JCT!")
      #}
      #tr.splice <- tr
      start.splice <- tr.splice$start
      end.splice <- tr.splice$end
      middle.splice <- apply(rbind(start.splice, end.splice), 2, median)
      #colors.splice <- exoncol[! tr$is.exon]
      segments(start.splice, rect.ceil, middle.splice, ymax, col=tr.splice$lineColor, lty = tr.splice$lty, lwd = plot.lwd, xpd=xpd,...);
      segments(middle.splice, ymax, end.splice, rect.ceil, col=tr.splice$lineColor,   lty = tr.splice$lty, lwd = plot.lwd, xpd=xpd,...);
    }
   #}

   if( any(tr$is.exon) ){
      tr.ex <- tr[tr$is.exon, , drop = FALSE]
      middle.ex <- apply(rbind(tr.ex$start,tr.ex$end),2,median);
      colors.ex <- exoncol[tr$is.exon]
      segments(middle.ex, rect.ceil, middle.ex, ymax, col=colors.ex, lty = exonlty[tr$is.exon], lwd = plot.lwd, xpd = xpd,...);
   }

}
drawGene.noSplices <- function(minx, maxx, tr.allExon, exon.names=NULL, anno.lwd = 1, par.cex = 1, anno.cex.text = 1, ...)
{
   #message("drawGene.noSplices cex = ",cex);
   plot.new()
   par(cex = par.cex);
   plot.window(xlim=c(minx, maxx), ylim=c(0, 1))

   rect(tr.allExon[,1], 0, tr.allExon[,2], 0.75, col="lightgrey", lty = 1, lwd=anno.lwd,...)
   
   if(! is.null(exon.names)){
     text(apply(tr.allExon,1,mean),0.5,labels=exon.names,srt=90, cex = anno.cex.text);
   }
}

#drawTranscript(rel.calc.min, rel.calc.max, tr=tr, rango=1:(length(anno.data$start[rt.allExon][logicexons==1])), exoncol=NULL, names=c(), trName=trans[i], cex=0.8,sub.sig = sub.sig)external.margins = c(1,4,4,2),

drawTranscript <- function(minx, maxx, ymin, tr, tr.allJunction, rango, rescale.iv = NULL, names, trName, sub.sig, anno.lwd = 1, par.cex = 1, anno.cex.text = 1, anno.cex.TX.ID = anno.cex.text * 0.5,  ...)
{
   exon.height = 0.8;
   if(! is.null(trName)){
     id.ht <- strheight(trName, cex = anno.cex.TX.ID);
     if(id.ht < 0.75){
       exon.height <- 0.94 - id.ht;
     } else {
       exon.height <- 0.5;
     }
     text(minx,ymin+0.95,labels=trName,srt=0,xpd=NA,cex= anno.cex.TX.ID, adj=c(0,1),...);
     ##text(minx,0.5,labels=trName,srt=0,xpd=TRUE,cex= anno.cex.TX.ID, adj=c(0,0),...);
   }
   splice.top <- exon.height * 0.9;
   splice.mid <- splice.top / 2;
   
   #exoncol
   #message("drawTranscript cex = ",cex);
   #plot.new()
   #par(cex = par.cex, mar=c(0, external.margins[2], anno.cex.TX.ID, external.margins[4]));
   ##par(cex = par.cex, mar=c(0, external.margins[2], 0, external.margins[4]));
   #plot.window(xlim=c(minx, maxx), ylim=c(0, 1),...)
   
   #Remove adjacent exons:
   tr.raw <- tr;
   #print(tr.raw);
   
   is.adjacent <- tr$start[-1] == tr$end[-length(tr$end)];
   tr.mergedExons <- data.frame(start = tr.raw$start[1], end = tr.raw$end[1]);
   exon.breaks <- c();
   if(nrow(tr.raw) > 1){
     for(i in 2:nrow(tr.raw)){
       if( tr.raw$start[i] == tr.mergedExons$end[nrow(tr.mergedExons)] ){
         tr.mergedExons$end[nrow(tr.mergedExons)] <- tr.raw$end[i];
         exon.breaks <- c(exon.breaks, tr.raw$start[i]);
       } else {
         tr.mergedExons <- rbind(tr.mergedExons, c(tr.raw$start[i], tr.raw$end[i] ));
       }
     }
   } else {
     tr.mergedExons <- tr.raw;
   }
   
   if(nrow(tr.mergedExons) > 1){
     tr.splices <- data.frame(start = tr.mergedExons$end[-nrow(tr.mergedExons)], end = tr.mergedExons$start[-1]);
     tr.splices$lineColor <- rep("#DDDDDD",nrow(tr.splices));
     for(i in 1:nrow(tr.splices)){
        #color.sig <- color.sig | (is.non.na & sub.sig$start[i] == t.splices[,1] & sub.sig$end[i] == t.splices[,2]);
        matched.jct <- which(tr.allJunction$start == tr.splices$start[i] & tr.allJunction$end == tr.splices$end[i] );
        if(length(matched.jct) != 1){
          message("WARNNG: unmatched junction (len=",length(matched.jct),") (Jct #",i,")!\n","[",tr.splices$start[i],",",tr.splices$end[i],"]");
          print("####tr.splices:");
          print(tr.splices);
          print("####tr.allJunctions:");
          print(tr.allJunction);
          print("####tr.mergedExons:");
          print(tr.mergedExons);
          print("####tr.raw:");
          print(tr.raw);
          print("####exon.breaks:");
          print(exon.breaks);
        } else {
          tr.splices$lineColor[i] <- tr.allJunction$lineColor[matched.jct];
        }
     }
   }
   
   if(! is.null(rescale.iv)){
     tr.mergedExons$start <- (rescale.coords(tr.mergedExons$start,rescale.iv) * (maxx-minx)) + minx;
     tr.mergedExons$end   <- (rescale.coords(tr.mergedExons$end,  rescale.iv) * (maxx-minx)) + minx;
     tr.raw$start <- (rescale.coords(tr.raw$start,rescale.iv) * (maxx-minx)) + minx;
     tr.raw$end   <- (rescale.coords(tr.raw$end,  rescale.iv) * (maxx-minx)) + minx;     
     if(nrow(sub.sig) > 0){
       sub.sig$start <- (rescale.coords(sub.sig$start,rescale.iv) * (maxx-minx)) + minx;
       sub.sig$end   <- (rescale.coords(sub.sig$end,  rescale.iv) * (maxx-minx)) + minx;
     }
     if(nrow(tr.mergedExons) > 1){
       tr.splices$start <- (rescale.coords(tr.splices$start,rescale.iv) * (maxx-minx)) + minx;
       tr.splices$end   <- (rescale.coords(tr.splices$end,rescale.iv) * (maxx-minx)) + minx;
     }
     if(length(exon.breaks) > 0){
       exon.breaks <- (rescale.coords(exon.breaks,rescale.iv) * (maxx-minx)) + minx;
     }
   }
   
   #print("TR.RAW:");
   #print(tr.raw);
   #print("TR.MERGEDEXONS:");
   #print(tr.mergedExons);
   #}
   #splice.color <- ifelse(color.sig,"#F219ED","black");
   
   if(nrow(tr.mergedExons) > 1){
     zr <- apply(rbind(tr.splices$start, tr.splices$end), 2, median)
     segments(tr.splices$start, ymin+splice.mid, zr, ymin+splice.top,col=tr.splices$lineColor, lwd = anno.lwd,...)
     segments(zr, ymin+splice.top,tr.splices$end, ymin+splice.mid,col=tr.splices$lineColor, lwd = anno.lwd,...)
   }
   
   TX.EXON.BORDER.COLOR <- "black";
   rect(tr.raw$start, ymin, tr.raw$end, ymin+exon.height, col=tr$fillColor, border = "transparent", lwd = anno.lwd / 2, ...);
   rect(tr.mergedExons$start, ymin, tr.mergedExons$end, ymin+exon.height, col="transparent", border = TX.EXON.BORDER.COLOR, lwd = anno.lwd / 2, ...);
   
   if(length(exon.breaks) > 0){
     segments(exon.breaks,ymin,exon.breaks,ymin+exon.height, col = "black", lwd = anno.lwd / 2, lty = 3, ...)
   }
   #rect(tr[rango,1], 0.1, tr[rango,2], 0.5, col=exoncol, lwd = anno.lwd / 2)
   #zr <- apply(rbind(tr[rango, 2], tr[rango+1, 1]), 2, median)
   #segments(tr[rango,2], 0.3, zr, 0.45,col=splice.color, lwd = anno.lwd,...)
   #segments(zr, 0.45, tr[rango+1,1], 0.3,col=splice.color, lwd = anno.lwd,...)
   

}
#############################################


generate.interval.scale <- function(gene.features, exon.fraction){
  if(is.na(exon.fraction) | exon.fraction <= 0 | exon.fraction >= 1){ #if rescaling is off, do not rescale!
    iv <- data.frame(start = MIN, end = MAX, span = SPAN,
                       rescale.start = 0, rescale.end = 1);
    return(iv);
  }

  MAX <- max(c(gene.features$start, gene.features$end));
  MIN <- min(c(gene.features$start, gene.features$end));
  SPAN <- MAX - MIN;
  
  exon.iv.raw <- data.frame(start = gene.features$start[gene.features$is.exon], 
                            end = gene.features$end[gene.features$is.exon]);
  
  exon.iv <- exon.iv.raw[1,];
  if(nrow(exon.iv.raw) > 1){
    for(i in 2:nrow(exon.iv.raw)){
      if(exon.iv$end[nrow(exon.iv)] == exon.iv.raw$start[i]){
        exon.iv$end[nrow(exon.iv)] = exon.iv.raw$end[i];
      } else {
        exon.iv <- rbind(exon.iv, exon.iv.raw[i,]);
      }
    }
  }
  
  intr.iv <- data.frame(start = exon.iv$end[-nrow(exon.iv)], end = exon.iv$start[-1]);
  if(max(exon.iv$end) < MAX){
    intr.iv <- rbind.data.frame(intr.iv,data.frame(start = max(exon.iv$end), end = MAX));
  }
  if(min(exon.iv$start) > MIN){
    intr.iv <- rbind.data.frame(data.frame(start = MIN, end = min(exon.iv$start)),intr.iv);
  }
  exon.iv$span <- exon.iv$end - exon.iv$start;
  intr.iv$span <- intr.iv$end - intr.iv$start;
  
  if(nrow(intr.iv) > 0){
    exon.iv$normSpan <- exon.fraction * exon.iv$span / sum(exon.iv$span);
    intr.iv$normSpan <- (1 - exon.fraction) * intr.iv$span / sum(intr.iv$span);
    idx <- order(c(1:nrow(exon.iv),1:nrow(intr.iv)))
    iv <- rbind(exon.iv, intr.iv)[idx,, drop = FALSE];
    cumulativeSum <- cumsum(iv$normSpan);
    iv$rescale.start <- c(0,cumulativeSum[-length(cumulativeSum)]);
    iv$rescale.end <- cumulativeSum;
    iv$simple.normSpan <- iv$span / sum(iv$span);
    iv$simple.rescale.start <- (iv$start - MIN) / SPAN;
    iv$simple.rescale.end   <- (iv$end - MIN) / SPAN;

    return(iv);
  } else { #Special case: the entire gene is exonic. Use simple linear scale.
    iv <- data.frame(start = MIN, end = MAX, span = SPAN,
                       rescale.start = 0, rescale.end = 1);
    return(iv);
  }
}

rescale.coords <- function(x, rescale.iv){
  sapply(x, FUN=function(y){
    idy <- which(y >= rescale.iv$start & y <= rescale.iv$end)[1];
    pct <- (y - rescale.iv$start[idy]) / rescale.iv$span[idy]
    rescale.iv$rescale.start[idy] + rescale.iv$normSpan[idy] * pct;
  });
}

#############################################

#Expansion of axis.break from the plotrix package
# adds the ability to modify additional graphical parameters, and fixes a bug where
# this function globally overrides certain graphical parameters via par().

qorts.axis.break <- function (axis = 1, breakpos = NULL, pos = NA, bgcol = "white", 
    breakcol = "black", style = "slash", cex = 1, xw = NULL, yw = NULL,...) 
{
    brw <- strwidth("W",cex = cex, units="figure");
    brh <- strheight("W",cex = cex, units="figure");
    figxy <- par("usr")
    xaxl <- par("xlog")
    yaxl <- par("ylog")
    if(is.null(xw)) xw <- (figxy[2] - figxy[1]) * brw
    if(is.null(yw)) yw <- (figxy[4] - figxy[3]) * brh
    if (!is.na(pos)) 
        figxy <- rep(pos, 4)
    if (is.null(breakpos)) 
        breakpos <- ifelse(axis%%2, figxy[1] + xw * 2, figxy[3] + 
            yw * 2)
    if (xaxl && (axis == 1 || axis == 3)) 
        breakpos <- log10(breakpos)
    if (yaxl && (axis == 2 || axis == 4)) 
        breakpos <- log10(breakpos)
    switch(axis, 
        br <- c(breakpos - xw/2, figxy[3] - yw/2, breakpos + xw/2, figxy[3] + yw/2), 
        br <- c(figxy[1] - xw/2, breakpos - yw/2, figxy[1] + xw/2, breakpos + yw/2), 
        br <- c(breakpos - xw/2, figxy[4] - yw/2, breakpos + xw/2, figxy[4] + yw/2), 
        br <- c(figxy[2] - xw/2, breakpos - yw/2, figxy[2] + xw/2, breakpos + yw/2), 
        stop("Improper axis specification."))
    #old.xpd <- par("xpd")
    #par(xpd = TRUE)
    if (xaxl) 
        br[c(1, 3)] <- 10^br[c(1, 3)]
    if (yaxl) 
        br[c(2, 4)] <- 10^br[c(2, 4)]
    if (style == "gap") {
        if (xaxl) {
            figxy[1] <- 10^figxy[1]
            figxy[2] <- 10^figxy[2]
        }
        if (yaxl) {
            figxy[3] <- 10^figxy[3]
            figxy[4] <- 10^figxy[4]
        }
        if (axis == 1 || axis == 3) {
            rect(breakpos, figxy[3], breakpos + xw, figxy[4], 
                col = bgcol, border = bgcol)
            xbegin <- c(breakpos, breakpos + xw)
            ybegin <- c(figxy[3], figxy[3])
            xend <- c(breakpos, breakpos + xw)
            yend <- c(figxy[4], figxy[4])
            if (xaxl) {
                xbegin <- 10^xbegin
                xend <- 10^xend
            }
        }
        else {
            rect(figxy[1], breakpos, figxy[2], breakpos + yw, 
                col = bgcol, border = bgcol)
            xbegin <- c(figxy[1], figxy[1])
            ybegin <- c(breakpos, breakpos + yw)
            xend <- c(figxy[2], figxy[2])
            yend <- c(breakpos, breakpos + yw)
            if (xaxl) {
                xbegin <- 10^xbegin
                xend <- 10^xend
            }
        }
        #par(xpd = TRUE)
    }
    else {
        #rect(br[1], br[2], br[3], br[4], col = bgcol, border = bgcol)
        if (style == "slash") {
            if (axis == 1 || axis == 3) {
                xbegin <- c(breakpos - xw, breakpos)
                xend <- c(breakpos, breakpos + xw)
                ybegin <- c(br[2], br[2])
                yend <- c(br[4], br[4])
                if (xaxl) {
                  xbegin <- 10^xbegin
                  xend <- 10^xend
                }
            }
            else {
                xbegin <- c(br[1], br[1])
                xend <- c(br[3], br[3])
                ybegin <- c(breakpos - yw, breakpos)
                yend <- c(breakpos, breakpos + yw)
                if (yaxl) {
                  ybegin <- 10^ybegin
                  yend <- 10^yend
                }
            }
        }
        else {
            if (axis == 1 || axis == 3) {
                xbegin <- c(breakpos - xw/2, breakpos - xw/4, 
                  breakpos + xw/4)
                xend <- c(breakpos - xw/4, breakpos + xw/4, breakpos + 
                  xw/2)
                ybegin <- c(ifelse(yaxl, 10^figxy[3 + (axis == 
                  3)], figxy[3 + (axis == 3)]), br[4], br[2])
                yend <- c(br[4], br[2], ifelse(yaxl, 10^figxy[3 + 
                  (axis == 3)], figxy[3 + (axis == 3)]))
                if (xaxl) {
                  xbegin <- 10^xbegin
                  xend <- 10^xend
                }
            }
            else {
                xbegin <- c(ifelse(xaxl, 10^figxy[1 + (axis == 
                  4)], figxy[1 + (axis == 4)]), br[1], br[3])
                xend <- c(br[1], br[3], ifelse(xaxl, 10^figxy[1 + 
                  (axis == 4)], figxy[1 + (axis == 4)]))
                ybegin <- c(breakpos - yw/2, breakpos - yw/4, 
                  breakpos + yw/4)
                yend <- c(breakpos - yw/4, breakpos + yw/4, breakpos + 
                  yw/2)
                if (yaxl) {
                  ybegin <- 10^ybegin
                  yend <- 10^yend
                }
            }
        }
    }
    segments(xbegin, ybegin, xend, yend, col = breakcol, lty = 1, xpd = NA, ...)
    
    #par(xpd = FALSE)
}

#axis(1, at=geneColumn.MID, labels="GENE", tcl = -0.5,  lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 0, mgp = c(3,0.55,0),padj=0.5, hadj=0.5,...)#
#INCOMPLETE:
JS.axis.old <- function(side, at, labels = at, tcl = -0.5, lwd = 1, lwd.ticks = 1, cex.axis = 1, srt = 0, mgp = c(3,1,0), line = NA, pos = NA, adj = c(0.5,1.1), font = 1, ...){
  #axis(side, at, labels = FALSE, tcl = tcl, lwd = lwd, lwd.ticks = lwd.ticks, cex.axis = cex.axis, mgp = mgp, line = line, pos = pos, ...);
  
  if(side == 1){
    axis.height <- par("usr")[3];
    axis.tick.floor <- axis.height + (tcl * (2 * par("cxy")[2]));
    segments(at, axis.tick.floor,at, axis.height, lwd = lwd.ticks, cex = cex.axis, xpd = NA, ...);
    
    axis.label.height <- axis.height - (mgp[2] * (2 * par("cxy")[2]));
    
    text(x = at, y = axis.label.height, labels = labels, cex = cex.axis, adj = adj, srt = srt, xpd = NA, font = font,...);
    axis.floor <- axis.label.height - max(strheight(labels, cex = cex.axis))
    
    #message("Axis range: [",axis.label.height,",",axis.floor,"]");
    
    return(axis.floor);
  }
}
JS.axis <- function(side, at, labels = at, tcl = -0.5, lwd = 1, lwd.ticks = 1, cex.axis = 1, srt = 0, line = NA, pos = NA, adj = c(0.5,1.1), font = 1, ...){
  #axis(side, at, labels = FALSE, tcl = tcl, lwd = lwd, lwd.ticks = lwd.ticks, cex.axis = cex.axis, mgp = mgp, line = line, pos = pos, ...);
  
  if(side == 1){
    axis.height <- par("usr")[3];
    axis.tick.floor <- axis.height + (tcl * (2 * par("cxy")[2]));
    segments(at, axis.tick.floor,at, axis.height, lwd = lwd.ticks, cex = cex.axis, xpd = NA, ...);
    devlim <- device.limits();

    axis.label.height <- abs((axis.tick.floor - devlim[3])/2) + devlim[3];
    
    text(x = at, y = axis.label.height, labels = labels, cex = cex.axis, adj = adj, srt = srt, xpd = NA, font = font,...);
    axis.floor <- axis.label.height - max(strheight(labels, cex = cex.axis))
    
    #message("Axis range: [",axis.label.height,",",axis.floor,"]");
    
    return(axis.floor);
  }
}


autofit.strings <- function(xlim, y, at, labels, start.cex, buffer.char = "m", min.horiz.cex = 1, ...){
  
}



