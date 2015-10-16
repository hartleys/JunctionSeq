
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
drawPlot <- function(matr, ylimn,ecs, intervals, rango, fitExpToVar, numexons, textAxis, geneLevelAxisTitle = NULL, rt, color.count, 
                     colorlines,countbinIDs,use.vst, use.log,plot.type,main.title,draw.legend,color.key,
                     condition.names,p.values=NULL,draw.p.values=FALSE,plot.lwd = 1,axes.lwd = 1, 
                     anno.lwd =1, par.cex = 1, anno.cex.text = 1, anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                     fit.countbin.names = TRUE, fit.pvals = TRUE, stagger.pvals = 2,debug.mode = FALSE,
                     plot.gene.level.expression = FALSE, geneCount = NULL, color.geneCount = NULL, 
                     yAxisLabels.inExponentialForm = TRUE,
                     autofit.legend = TRUE, italicize.label = NULL, condition.legend.text = condition.legend.text,
                     annolink.col = NULL, exonlty = NULL,
                     graph.margins = c(2,3,3,2),
                     plotWindowXmax = 1.04,
                     fit.labels = TRUE,
                     above.pvals = FALSE,
                     ...)
{
   plot.junction.ids.vertically <- NULL;
   plot.junction.ids.vertically.threshold <- Inf;
   
   if(debug.mode) message("> Step 6.1");
   
   #message("drawPlot cex = ",cex);
   par(mar=graph.margins, cex = par.cex);
   plot.new();
   par(cex = par.cex);
   
     #str.ht <- strheight("X",units="figure",
     #plot.area.yaxs <- 1 - (par("plt")[1] + par("plt")[4])
     #spare.yaxs.margin <- 
     #TODO: add autofit. For now just hard-code it:
     spare.yaxs.margin <- (abs(ylimn[2] - ylimn[1])) * 0.08;
     spare.yaxs.margin <- (abs(ylimn[2] - ylimn[1]))*0.04;
   
   #message(paste("ylimn:",ylimn));
   plot.window(xlim=c(0, plotWindowXmax), ylim=c(ylimn[1], ylimn[2] + spare.yaxs.margin), xaxs = "i", yaxs="i");
   #plot.window(xlim=c(0, 1), ylim=ylimn, xaxs = "r");
   
   if(debug.mode) message("> Step 6.2");
   makevstaxis(1/ncol(matr), ylimn, ecs,use.vst=use.vst, use.log=use.log,plot.type=plot.type, lwd = axes.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, ...)
   title(main=main.title, cex.main = anno.cex.main, ...);
   if(debug.mode) message("> Step 6.3");
   
   if(plot.type == "logRawCounts")  segments(0,INTERNAL.NINF.VALUE,par("usr")[2]+1,INTERNAL.NINF.VALUE,lty="dotted",col="black", lwd = plot.lwd,...); 
   if((! use.vst) & use.log ) segments(0,INTERNAL.NINF.VALUE,par("usr")[2]+1,INTERNAL.NINF.VALUE,lty="dotted",col="black", lwd = plot.lwd,...); 
   segments(0,0,par("usr")[2]+1,0,lty="dotted",col="black", lwd = plot.lwd,...);
   junctionColumns.RIGHT <- (intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2))[length(rango)];
   if(debug.mode) message("> Step 6.4");

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
     
     
     
     
   if(plot.gene.level.expression){
     rect(junctionColumns.RIGHT,    par("usr")[3], geneColumn.LEFT,    par("usr")[4], xpd=NA, lwd = axes.lwd,xpd=NA, border="white",col="white",...);
     rect(0,    par("usr")[3], junctionColumns.RIGHT,    par("usr")[4], xpd=NA, lwd = axes.lwd,xpd=NA,...);
     rect(geneColumn.LEFT, par("usr")[3], geneColumn.RIGHT, geneColumn.TOP, xpd=NA, lwd = axes.lwd,xpd=NA,...);
   } else {
     rect(junctionColumns.RIGHT,    par("usr")[3], par("usr")[2],    par("usr")[4], xpd=NA, lwd = axes.lwd,xpd=NA, border="white",col="white",...);
     rect(0,    par("usr")[3], junctionColumns.RIGHT,    par("usr")[4], xpd=NA, lwd = axes.lwd,xpd=NA,...);
   }
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
   #else {
   #  if(draw.legend){
   #     legend("topright",fill=color.key[names(condition.legend.text)],legend=condition.legend.text[condition.names], cex = anno.cex.text, bg = "transparent", box.lwd = axes.lwd, xpd = NA,...);
   #  }
   #}
   

   
   if(debug.mode) message("> Step 6.5");
   
   #intervals<-(0:nrow(matr))/nrow(matr)
   middle <- apply(cbind(intervals[rango], (intervals[rango+1]-((intervals[rango+1])-intervals[rango])*0.2)), 1, median)
   matr <- rbind(matr, NA)
   j <- 1:ncol(matr)
   
   #abline(h = 0,lty="dotted",col="black", lwd = plot.lwd,...);
   segments(intervals[rango],matr[rango,j], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), matr[rango,j], col=color.count, lwd = plot.lwd,...)  #### line with the y level
   segments(intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), matr[rango,j], intervals[rango+1], matr[rango+1,j], col=color.count, lty="dotted", lwd =plot.lwd,...)  #### line joining the y levels
   abline(v=middle[rango], lty="dotted", col=colorlines, lwd = plot.lwd)
   
   cex.mainPlotAxis <- anno.cex.text;
   if(fit.labels){
     cex.mainPlotAxis <- shrink.character.vector.VERT(textAxis, cex.mainPlotAxis, abs(ylimn[2] - ylimn[1]) * 0.95);
   }
   text(device.limits()[1], ylimn[1] + abs(ylimn[2] - ylimn[1])*0.5, textAxis,           cex = cex.mainPlotAxis, xpd = NA, adj = c(0.5,1.1), srt = 90 , ...);
   if(plot.gene.level.expression){
     cex.geneLevelAxis <- anno.cex.text;
     if(fit.labels){
       cex.geneLevelAxis <- shrink.character.vector.VERT(geneLevelAxisTitle, cex.geneLevelAxis, abs(geneColumn.TOP - ylimn[1]) *0.95);
     }
     text(device.limits()[2], abs(geneColumn.TOP - ylimn[1]) *0.5 + ylimn[1], geneLevelAxisTitle, cex = cex.geneLevelAxis, xpd = NA, adj = c(0.5,1.1), srt = 270 , ...);
   }
   #title(ylab = textAxis, cex.lab = anno.cex.text, mgp = c(2.5,1,0), ...);
   #mtext(textAxis, side=2, adj=0.5, line=2.5, outer=FALSE, cex = anno.cex.text * par("cex"),...)
   
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
           cex.pvalues <- fit.character.vector(p.values, min.width = 1 * bin.width, max.width = 1.5 * bin.width, max.width.per.char = 0.5 * bin.width)
           cex.pvalues <- min(cex.pvalues, anno.cex.text);
         } else {
           cex.pvalues <- anno.cex.text;
         }
         pval.height <- max(strheight(p.values, cex = cex.pvalues));
         
         if(above.pvals){
           pval.ceiling <- par("usr")[4] + pval.height * stagger.pvals;
           pval.yadj <- 0;
         } else {
           pval.ceiling <- par("usr")[4] - pval.height*0.1;
           pval.yadj <- 1;
         }
         
        pval.y <- rep(pval.ceiling,length(rango))
        if(stagger.pvals > 1){
          tryCatch({
            if(debug.mode) message("> Step 6.9b");
            
            #CURRENTLY ONLY IMPLEMENTING 2-layer stagger! 3-layer stagger is DEPRECIATED.
            
            pval.drop <- rep(-pval.height,length(p.values));
              
              pw <- which(p.values != "");
              pwo <- pw;
              #print(pw);
              while(length(pw) > 0){
                if(any( c(pw[1] - 1, pw[1] + 1) %in% pwo ) && pw %% 2 == 0){
                  pval.drop[pw[1]] <- pval.drop[pw[1]] + pval.height;
                 # message("(",pw[1],")");
                } else{
                  #do nothing.
                }
                pw <- pw[-1];
              }
            pval.y <- pval.ceiling + pval.drop;
            #if(stagger.pvals == 3){
            #  pval.y <- rep(c(pval.ceiling, pval.ceiling - pval.height, pval.ceiling - (2 * pval.height)), length(middle[rango]))[rango]
            #} else if(stagger.pvals == 2){
            #  pval.y <- rep(c(pval.ceiling, pval.ceiling - pval.height), length(middle[rango]))[rango]
            #} else {
            #  pval.y <- rep(pval.ceiling,length(rango))
            #}
          }, error = function(e){
            message("WARNING: staggering p-values failed. Falling back to simpler method. This warning should never appear. If you see this, reporting this to the developer would be appreciated.");
          });
        }
        
        if(strwidth(p.values[1], cex = cex.pvalues) > bin.width  * 0.9){
          text(intervals[1], pval.y[1], labels = p.values[1], col = colorlines[1], cex = cex.pvalues, adj = c(0,pval.yadj), ...);
        } else {
          text(middle[1], pval.y[1], labels = p.values[1], col = colorlines[1], cex = cex.pvalues, adj = c(0.5,pval.yadj),  ...);
        }
        if(strwidth(p.values[length(p.values)], cex = cex.pvalues) > bin.width * 0.9){
          text(junctionColumns.RIGHT, pval.y[length(middle)], labels = p.values[length(middle)], col = colorlines[length(middle)], cex = cex.pvalues, adj = c(1,pval.yadj), xpd = NA, ...);
        } else {
          text(middle[length(middle)], pval.y[length(middle)], labels = p.values[length(middle)], col = colorlines[length(middle)], cex = cex.pvalues, adj = c(0.5,pval.yadj), xpd = NA, ...);
        }
        text(middle[2:(length(middle)-1)],pval.y[2:(length(middle)-1)],labels=p.values[2:(length(middle)-1)],col=colorlines[2:(length(middle)-1)], cex = cex.pvalues, adj = c(0.5,pval.yadj), xpd = NA,...);
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
     cex.gene.tick <- min(anno.cex.axis, cex.countbinIDs * 1.5);
     srt.gene.tick <- 0;
     if(strwidth("GENE",cex = cex.gene.tick)/2 > abs(geneColumn.MID - junctionColumns.RIGHT) * 0.9 ){
       cex.gene.tick <- cex.countbinIDs;
       srt.gene.tick <- srt.countbinIDs;
     }
     JS.axis(1, at=geneColumn.MID, labels="GENE", tcl =  tcl.countbinIDs,  lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.gene.tick, srt = srt.gene.tick, font = 2, adj = adj.countbinIDs,...)# ,
     #axis(1, at=geneColumn.MID, labels="GENE", tcl = -0.5,  lwd = axes.lwd, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 0, mgp = c(3,0.55,0),padj=0.5, hadj=0.5,...)#
   } else {
     #rect(0, par("usr")[3], par("usr")[2], par("usr")[4], xpd=NA, lwd = axes.lwd,xpd=NA,...);
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
                     trName, exonlty, plot.lwd = 1, anno.lwd = 1, show.strand.arrows = 0, geneStrand = ".", par.cex = 1, anno.cex.text = 1, 
                     draw.untestable.annotation = TRUE,
                     draw.start.end.sites = TRUE, startSites = NULL, endSites = NULL,
                     cex.arrows, chrom.label = "", label.chromosome = TRUE,
                     draw.curved.SJ = TRUE, draw.nested.SJ = TRUE,
                     draw.exon.lines = TRUE,
                     splice.junction.drawing.style = "hyperbola",
                     merge.exon.parts = TRUE,
                     plot.untestable.results = FALSE,
                     exon.height = 0.75,
                     INTERNAL.VARS = INTERNAL.VARS,
                     flip.splicing = FALSE, gapped.flip = FALSE,
                     ...)
{
   #message("drawGene cex = ",cex);
   #par(cex = par.cex);
   xpd <- NA
   plot.new()
   #plot.window(xlim=c(minx, maxx + ((maxx-minx)* 0.04)), ylim=c(0,1), xaxs = "i");
   plot.window(xlim=c(minx, maxx), ylim=c(0,1), xaxs = "i");
   #plot.window(xlim=c(minx, maxx), ylim=c(0, 1))
   ymax <- par("usr")[4];
   ymin <- par("usr")[3];
   startSiteMarks.lwd <- plot.lwd;
   arrows.lwd <- plot.lwd;
   centerLine.lwd <- plot.lwd;
   geneRect.lwd <- plot.lwd;
   
   start.end.sites.height <- 0.2;

   if(! flip.splicing){
     rect.floor <- 0.1 * (exon.height);
     rect.ceil  <- exon.height;
     splice.ceil <- (0.975 * (1-exon.height)) + exon.height;
     splice.floor <- rect.ceil;
     
     startSite.marker.width <- strwidth("M",cex = anno.cex.text);
     startSite.angle.width <- startSite.marker.width * 1/3;
     if(draw.start.end.sites){
       rect.floor <- start.end.sites.height * (exon.height);
       startSites.floor <- 0.05 * exon.height;
       endSites.floor <- (rect.floor - startSites.floor) * (1/3) + startSites.floor;
     }
     rect.mid   <- ((rect.ceil - rect.floor)/2) + rect.floor
   } else {
     rect.floor <- (1-exon.height);
     rect.ceil  <- (0.975 * exon.height) + rect.floor;
     
     splice.ceil <- rect.floor
     splice.floor <- 0.05 * (1-exon.height);
     
     startSite.angle.width <-  strwidth("M",cex = anno.cex.text)  * 1/3;
     startSite.marker.width <- strwidth("M",cex = anno.cex.text)  * 1/3;
     
     if(draw.start.end.sites){
       startSites.floor <- rect.floor - (0.1 * exon.height);
       #endSites.floor <- (rect.floor - startSites.floor) * (1/3) + startSites.floor;
       endSites.floor <- startSites.floor
       
       if(gapped.flip){
       #proposed alt:
         rect.floor <- (1-exon.height) + start.end.sites.height * exon.height;
         startSites.floor <- 0.05 * exon.height + (1-exon.height);
         #endSites.floor <- (rect.floor - startSites.floor) * (1/3) + startSites.floor;
         endSites.floor <- startSites.floor
       }
     }
   }
   rect.mid   <- ((rect.ceil - rect.floor)/2) + rect.floor
   
   #message("rescale.iv:"); print(rescale.iv);
   #message("RAW:"); print(tr.allExon);
   if(! is.null(rescale.iv)){
     tr.allExon$start <- (rescale.coords(tr.allExon$start,rescale.iv) * (maxx-minx)) + minx;
     tr.allExon$end   <- (rescale.coords(tr.allExon$end,rescale.iv) * (maxx-minx)) + minx;
     tr$start <- (rescale.coords(tr$start,rescale.iv) * (maxx-minx)) + minx;
     tr$end   <- (rescale.coords(tr$end,rescale.iv) * (maxx-minx)) + minx;
     
     #print(class(tr.allExon$start));
     
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
   #message("RESCALE:"); print(tr.allExon);
   
   
   if(merge.exon.parts & nrow(tr.allExon) > 1 & any(tr.allExon$start %in% tr.allExon$end) ){
     #merge adjacent exonic parts:
     tr.mergedExon <- data.frame(start = tr.allExon$start[1], end = tr.allExon$end[1]);
     tr.breakLine <- c();
     
     for(i in 2:nrow(tr.allExon)){
       j <- nrow(tr.mergedExon);
       if(tr.mergedExon$end[j] == tr.allExon$start[i]){
         tr.mergedExon$end[j] <- tr.allExon$end[i];
         tr.breakLine <- c(tr.breakLine, tr.allExon$start[i]);
       } else {
         tr.mergedExon[j+1,] <- c(tr.allExon$start[i], tr.allExon$end[i]);
       }
     }
     rect(tr.allExon$start, rect.floor, tr.allExon$end, rect.ceil, lty = 1, lwd = geneRect.lwd ,col=tr.allExon$fillColor, border= "transparent", xpd = xpd,...);
     rect(tr.mergedExon$start, rect.floor, tr.mergedExon$end, rect.ceil, lty = 1, lwd = geneRect.lwd, col = "transparent", border = "black", xpd = xpd, ...);
     segments(tr.breakLine, rect.floor, tr.breakLine, rect.ceil, lty = 3, lwd = geneRect.lwd, col = "gray33", ...);
   } else {
     rect(tr.allExon$start, rect.floor, tr.allExon$end, rect.ceil, lty = 1, lwd = geneRect.lwd ,col=tr.allExon$fillColor, xpd = xpd,...)
   }
   lines(c(minx,maxx),c(rect.mid,rect.mid),lty= 1, lwd = centerLine.lwd, xpd = xpd, ...);
   
   if(draw.start.end.sites){
     segments(startSites + startSite.angle.width, startSites.floor,startSites,rect.floor, lwd = startSiteMarks.lwd,xpd=NA,...);
     segments(endSites   - startSite.angle.width,   endSites.floor,  endSites,  rect.floor, lwd = startSiteMarks.lwd,xpd=NA,...);
     
     segments(startSites + startSite.angle.width, startSites.floor,  startSites + startSite.marker.width, startSites.floor, lwd = startSiteMarks.lwd,xpd=xpd,...);
     segments(endSites   - startSite.angle.width,   endSites.floor,  endSites   - startSite.marker.width,   endSites.floor, lwd = startSiteMarks.lwd,xpd=xpd,...);
   }
   
   if(show.strand.arrows > 0){
      #if(is.null(strand)){
      #   message("Warning: Cannot draw strand arrows, strandedness is NULL! Skipping.");
      #} else {
         strand <- tr.allExon$strand[1];
         
         if(cex.arrows == "auto"){
             arrow.height <- strheight(">", cex = 1);
             if(arrow.height > abs(rect.ceil - rect.floor) * 0.75 ){
               cex.arrows <- (abs(rect.ceil - rect.floor) * 0.75) / arrow.height
             } else {
               cex.arrows <- 1;
             }
         }
         if(show.strand.arrows == 1){
           if(strand == "+" | geneStrand == "+"){
             arrow.X <- par("usr")[2] + par("cxy")[1]*1.5
             lines(c(par("usr")[2],arrow.X),c(rect.mid, rect.mid), lwd = arrows.lwd, xpd = NA, ...);
             JS.arrowChars( arrow.X , rect.mid, "right", arrow.cex = cex.arrows, lwd = arrows.lwd, xpd = NA, ...);
           } else if(strand == "-" | geneStrand == "-"){
             arrow.X <- par("usr")[1] - par("cxy")[1]*1.5
             lines(c(par("usr")[1],arrow.X),c(rect.mid, rect.mid), lwd = arrows.lwd, xpd = NA, ...);
             JS.arrowChars( arrow.X , rect.mid, "left", arrow.cex = cex.arrows, lwd = arrows.lwd, xpd = NA, ...);
           }
           
         } else {
           arrow.X <-((1:(show.strand.arrows - 1)) / (show.strand.arrows)) * (maxx - minx) + minx
           if(strand == "+" | geneStrand == "+"){
              JS.arrowChars( arrow.X ,rep(rect.mid,length(arrow.X)), "right", arrow.cex = cex.arrows, lwd = arrows.lwd, ...);
           } else if(strand == "-" | geneStrand == "-"){
              JS.arrowChars( arrow.X ,rep(rect.mid,length(arrow.X) - 1), "left", arrow.cex = cex.arrows, lwd = arrows.lwd, ...);
           }
         }
         
         
         #even.spacing <- ((1:(show.strand.arrows)) / (show.strand.arrows+1)) * (maxx - minx) + minx;

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
      if(draw.nested.SJ){
        if(is.null(INTERNAL.VARS[["SJ.ceiling"]])){
          tr.splice$span <- tr.splice$end - tr.splice$start;
          #print(tr.splice);
          nh <- get.nested.heights(as.data.frame(tr.splice), splice.floor, splice.ceil);
          SJ.ceiling <- nh[["SJ.ceiling"]]
        } else {
          SJ.ceiling <- INTERNAL.VARS[["SJ.ceiling"]];
          SJ.ceiling <- SJ.ceiling * abs(splice.ceil - splice.floor) + splice.floor;
          SJ.ceiling <- SJ.ceiling[match(names(SJ.ceiling), tr.splice$featureID)];
        }
      } else {
        SJ.ceiling <- splice.ceil;
      }
      if(any(tr.splice$is.plotted)){
        segments(middle.splice[tr.splice$is.plotted], SJ.ceiling[tr.splice$is.plotted], middle.splice[tr.splice$is.plotted], ymax, col=tr.splice$lineColor[tr.splice$is.plotted], lty = tr.splice$lty[tr.splice$is.plotted], lwd = plot.lwd, xpd=xpd,...);
      }
      
      if(draw.curved.SJ){
        if(flip.splicing){
           drawHalfLoops("up",start.splice, end.splice, splice.ceil - SJ.ceiling + splice.floor, splice.ceil, col = tr.splice$lineColor, lwd = plot.lwd, xpd = xpd, lty = tr.splice$lty, style = splice.junction.drawing.style, ...);
        } else {
           drawHalfLoops("down",start.splice, end.splice, splice.floor,SJ.ceiling, col = tr.splice$lineColor, lwd = plot.lwd, xpd = xpd, lty = tr.splice$lty, style = splice.junction.drawing.style, ...);
        }
      } else {
        if(flip.splicing){
          segments(start.splice, splice.ceil - SJ.ceiling + splice.floor, middle.splice, splice.ceil, col=tr.splice$lineColor, lty = tr.splice$lty, lwd = plot.lwd, xpd=xpd,...);
          segments(middle.splice, splice.ceil, end.splice, splice.ceil - SJ.ceiling + splice.floor, col=tr.splice$lineColor,   lty = tr.splice$lty, lwd = plot.lwd, xpd=xpd,...);
        } else {
          segments(start.splice, splice.floor, middle.splice, SJ.ceiling, col=tr.splice$lineColor, lty = tr.splice$lty, lwd = plot.lwd, xpd=xpd,...);
          segments(middle.splice, SJ.ceiling, end.splice, splice.floor, col=tr.splice$lineColor,   lty = tr.splice$lty, lwd = plot.lwd, xpd=xpd,...);
        }

      }
    }
   #}

   if( any(tr$is.exon) & draw.exon.lines){
      tr.ex <- tr[tr$is.exon, , drop = FALSE]
      middle.ex <- apply(rbind(tr.ex$start,tr.ex$end),2,median);
      colors.ex <- exoncol[tr$is.exon]
      segments(middle.ex, rect.ceil, middle.ex, ymax, col=colors.ex, lty = exonlty[tr$is.exon], lwd = plot.lwd, xpd = xpd,...);
   }
   

}


JS.arrowChars <- function(xcenter, ycenter, direction, arrow.cex, lwd, ...){
  arrow.width = strwidth(">", cex= arrow.cex);
  arrow.height = strheight(">", cex = arrow.cex);
  if(direction == "left"){
    segments(xcenter,ycenter,xcenter + arrow.width,ycenter + arrow.height / 2, lwd = lwd, lend = 0, ...);
    segments(xcenter,ycenter,xcenter + arrow.width,ycenter - arrow.height / 2, lwd = lwd, lend = 0, ...);
  } else if(direction == "right"){
    segments(xcenter - arrow.width,ycenter + arrow.height / 2, xcenter, ycenter, lwd = lwd, lend = 0, ...);
    segments(xcenter - arrow.width,ycenter - arrow.height / 2, xcenter, ycenter, lwd = lwd, lend = 0, ...);
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

drawTranscript <- function(minx, maxx, ymin, tr, tr.allJunction, rango, rescale.iv = NULL, names, trName, trStrand = ".", sub.sig, anno.lwd = 1, par.cex = 1, anno.cex.text = 1, anno.cex.TX.ID = anno.cex.text * 0.5, cex.arrows = 1, draw.strand = FALSE,  ...)
{
   if(cex.arrows == "auto") cex.arrows <- 1;
   
   exon.height = 0.8;
   if(! is.null(trName)){
     id.ht <- strheight(trName, cex = anno.cex.TX.ID);
     if(id.ht < 0.75){
       exon.height <- 0.75 - id.ht;
     } else {
       exon.height <- 0.25;
     }
     
     ##text(minx,0.5,labels=trName,srt=0,xpd=TRUE,cex= anno.cex.TX.ID, adj=c(0,0),...);
   }
   splice.top <- exon.height * 0.9;
   splice.mid <- splice.top / 2;
   #Change? Splice junctions are now plotted at straight lines in TX annotation?
   #splice.top <- splice.mid;
   
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
     segments(tr.splices$start, ymin+splice.mid, zr, ymin+splice.top,col=tr.splices$lineColor, lwd = anno.lwd, xpd = NA,...)
     segments(zr, ymin+splice.top,tr.splices$end, ymin+splice.mid,col=tr.splices$lineColor, lwd = anno.lwd, xpd = NA,...)
   }
   
   TX.EXON.BORDER.COLOR <- "black";
   rect(tr.raw$start, ymin, tr.raw$end, ymin+exon.height, col=tr$fillColor, border = "transparent", lwd = anno.lwd / 2, xpd = NA, ...);
   rect(tr.mergedExons$start, ymin, tr.mergedExons$end, ymin+exon.height, col="transparent", border = TX.EXON.BORDER.COLOR, lwd = anno.lwd / 2, xpd = NA, ...);
   
   if(draw.strand){
     arrow.width <- convertHeightToWidth(exon.height) / 2;
     #strwidth(">", cex= cex.arrows) / 2;
     #message("trStrand = \"",trStrand,"\"");
     if(is.null(trStrand)){
       #Do Nothing
     } else if(trStrand == "+"){
       #message("+!");
       lastEnd <- tr.raw$end[length(tr.raw$end)];

       polygon(c(lastEnd,lastEnd, lastEnd + arrow.width), 
               #c(ymin + exon.height/3,ymin + exon.height * 2/3, exon.height/2 + ymin), 
               c(ymin,ymin + exon.height, exon.height/2 + ymin), 
               col=tr$fillColor[length(tr$fillColor)], border = TX.EXON.BORDER.COLOR, xpd = NA, lwd = anno.lwd / 2)
       #polygon(c(lastEnd,lastEnd, lastEnd + arrow.width), 
       #        c(ymin,ymin + exon.height, exon.height/2 + ymin),
       #        #c(ymin + exon.height/4,ymin + exon.height * 3/4, exon.height/2 + ymin), 
       #        col=tr$fillColor[length(tr$fillColor)], border = tr$fillColor[length(tr$fillColor)], xpd = NA, lwd = anno.lwd / 2)
       #lines(c(lastEnd,lastEnd + arrow.width), c(ymin, exon.height/2 + ymin),               col = TX.EXON.BORDER.COLOR, xpd = NA, lwd = anno.lwd / 2)
       #lines(c(lastEnd,lastEnd + arrow.width), c(ymin + exon.height, exon.height/2 + ymin), col = TX.EXON.BORDER.COLOR, xpd = NA, lwd = anno.lwd / 2)
     } else if (trStrand == "-"){
       #message("-!");
       firstStart <- tr.raw$start[1];
       polygon(c(firstStart,firstStart, firstStart - arrow.width),
               c(ymin,ymin + exon.height, exon.height/2 + ymin), 
               col=tr$fillColor[1], border = TX.EXON.BORDER.COLOR, xpd = NA, lwd = anno.lwd / 2)
       #polygon(c(firstStart,firstStart, firstStart - arrow.width),
       #        c(ymin,ymin + exon.height, exon.height/2 + ymin),
       #        #c(ymin + exon.height/4,ymin + exon.height * 3/4, exon.height/2 + ymin),
       #        col=tr$fillColor[1], border = tr$fillColor[1], xpd = NA, lwd = anno.lwd / 2)
       #lines(c(firstStart,firstStart - arrow.width), c(ymin, exon.height/2 + ymin),               col = TX.EXON.BORDER.COLOR, xpd = NA, lwd = anno.lwd / 2)
       #lines(c(firstStart,firstStart - arrow.width), c(ymin + exon.height, exon.height/2 + ymin), col = TX.EXON.BORDER.COLOR, xpd = NA, lwd = anno.lwd / 2)
     }
   }
   
   if(length(exon.breaks) > 0){
     segments(exon.breaks,ymin,exon.breaks,ymin+exon.height, col = "black", lwd = anno.lwd / 2, lty = 3, ...)
   }
   #rect(tr[rango,1], 0.1, tr[rango,2], 0.5, col=exoncol, lwd = anno.lwd / 2)
   #zr <- apply(rbind(tr[rango, 2], tr[rango+1, 1]), 2, median)
   #segments(tr[rango,2], 0.3, zr, 0.45,col=splice.color, lwd = anno.lwd,...)
   #segments(zr, 0.45, tr[rango+1,1], 0.3,col=splice.color, lwd = anno.lwd,...)
   
   id.wd <- strwidth(trName, cex = anno.cex.TX.ID);
   id.leftX <- tr.raw$start[1];
   if(id.leftX + id.wd > maxx){
     id.leftX <- maxx - id.wd;
   }
   
   text(id.leftX,ymin + 0.8,labels=trName,srt=0,xpd=NA,cex= anno.cex.TX.ID, adj=c(0,1),...);
}
#############################################


generate.interval.scale <- function(gene.features, exon.fraction, rescaleFunction = c("sqrt","log","linear","34root"), debug.mode = FALSE){
  rescaleFunction <- match.arg(rescaleFunction);
  if(rescaleFunction == "34root"){
    resFunc <- function(x){ x ^ 0.75 }
  } else if(rescaleFunction == "log"){
    resFunc <- function(x){ ifelse(x == 0, 0, log(x)) }
  } else if(rescaleFunction == "sqrt"){
    resFunc <- function(x){ sqrt(x) }
  } else {
    resFunc <- function(x){x}
  }
  
  MAX <- max(c(gene.features$start, gene.features$end));
  MIN <- min(c(gene.features$start, gene.features$end));
  SPAN <- MAX - MIN;
  
  if(is.na(exon.fraction) | exon.fraction <= 0 | exon.fraction >= 1){ #if rescaling is off, do not rescale!
    iv <- data.frame(start = MIN, end = MAX, span = SPAN,
                       rescale.start = 0, rescale.end = 1, normSpan = 1);
    return(iv);
  }
  
  exon.iv <- data.frame(start = gene.features$start[gene.features$is.exon], 
                            end = gene.features$end[gene.features$is.exon]);
  exon.iv$span <- exon.iv$end - exon.iv$start;
  #print(exon.iv);
  intr.iv <- data.frame(start = exon.iv$end[-nrow(exon.iv)], end = exon.iv$start[-1]);
  #print(intr.iv);
  if(max(exon.iv$end) < MAX){
    intr.iv <- rbind.data.frame(intr.iv,data.frame(start = max(exon.iv$end), end = MAX));
  }
  if(min(exon.iv$start) > MIN){
    intr.iv <- rbind.data.frame(data.frame(start = MIN, end = min(exon.iv$start)),intr.iv);
  }
  exon.iv$span <- exon.iv$end - exon.iv$start;
  intr.iv$span <- intr.iv$end - intr.iv$start;
  
  #intr.iv <- intr.iv[intr.iv$span > 0,,drop=FALSE];
  #print("intr.iv:");
  #print(intr.iv);
  
  if((nrow(intr.iv) > 0) && sum(intr.iv$span) > 0 ){
    exon.iv$type <- "E";
    intr.iv$type <- "I";
    #if(logRescale){
    #  exon.iv$normSpan <- exon.fraction * log(exon.iv$span) / sum(log(exon.iv$span));
    #  intr.iv$normSpan <- (1 - exon.fraction) * log(intr.iv$span) / sum(log(intr.iv$span));
    #} else {
    #  exon.iv$normSpan <- exon.fraction * exon.iv$span / sum(exon.iv$span);
    #  intr.iv$normSpan <- (1 - exon.fraction) * intr.iv$span / sum(intr.iv$span);
    #}
    #resFunc
    exon.iv$normSpan <- exon.fraction * resFunc(exon.iv$span) / sum(resFunc(exon.iv$span));
    intr.iv$normSpan <- (1 - exon.fraction) * resFunc(intr.iv$span) / sum(resFunc(intr.iv$span));
    
    #alternate intron and exon:
    
    iv <- rbind(exon.iv, intr.iv);
    iv <- iv[order(iv$start, iv$end),,drop=FALSE];
    cumulativeSum <- cumsum(iv$normSpan);
    
    iv$rescale.start <- c(0,cumulativeSum[-length(cumulativeSum)]);
    iv$rescale.end <- cumulativeSum;
    iv$normSpan <- iv$rescale.end - iv$rescale.start;
    iv$simple.normSpan <- iv$span / sum(iv$span);
    iv$simple.rescale.start <- (iv$start - MIN) / SPAN;
    iv$simple.rescale.end   <- (iv$end - MIN) / SPAN;
    return(iv);
  } else { #Catch special case: the entire gene is exonic:
    exon.iv$type <- "E";
    exon.iv$normSpan <- 1 * resFunc(exon.iv$span) / sum(resFunc(exon.iv$span));
    
    iv <- exon.iv;
    cumulativeSum <- cumsum(iv$normSpan);
    iv$rescale.start <- c(0,cumulativeSum[-length(cumulativeSum)]);
    iv$rescale.end <- cumulativeSum;
    iv$normSpan <- iv$rescale.end - iv$rescale.start;
    iv$simple.normSpan <- iv$span / sum(iv$span);
    iv$simple.rescale.start <- (iv$start - MIN) / SPAN;
    iv$simple.rescale.end   <- (iv$end - MIN) / SPAN;
    return(iv);
  }
}

 # #Remove adjacent exons? REMOVED.
 # #exon.iv <- exon.iv.raw[1,];
 # #if(nrow(exon.iv.raw) > 1){
 # #  for(i in 2:nrow(exon.iv.raw)){
 # #    if(exon.iv$end[nrow(exon.iv)] == exon.iv.raw$start[i]){
 # #      exon.iv$end[nrow(exon.iv)] = exon.iv.raw$end[i];
 # #    } else {
 # #      exon.iv <- rbind(exon.iv, exon.iv.raw[i,]);
 # #    }
 # #  }
 # #}
 # exon.iv <- exon.iv.raw;
 # 
 # intr.iv <- data.frame(start = exon.iv$end[-nrow(exon.iv)], end = exon.iv$start[-1]);
 # if(max(exon.iv$end) < MAX){
 #   intr.iv <- rbind.data.frame(intr.iv,data.frame(start = max(exon.iv$end), end = MAX));
 # }
 # if(min(exon.iv$start) > MIN){
 #   intr.iv <- rbind.data.frame(data.frame(start = MIN, end = min(exon.iv$start)),intr.iv);
 # }
 # exon.iv$span <- exon.iv$end - exon.iv$start;
 # intr.iv$span <- intr.iv$end - intr.iv$start;
 # 
 # #intr.iv <- intr.iv[intr.iv$span > 0,,drop=FALSE];
 # #print("intr.iv:");
 # #print(intr.iv);
 # 
 # if((nrow(intr.iv) > 0) && sum(intr.iv$span) > 0 ){
 #   exon.iv$type <- "E";
 #   intr.iv$type <- "I";
 #   #if(logRescale){
 #   #  exon.iv$normSpan <- exon.fraction * log(exon.iv$span) / sum(log(exon.iv$span));
 #   #  intr.iv$normSpan <- (1 - exon.fraction) * log(intr.iv$span) / sum(log(intr.iv$span));
 #   #} else {
 #   #  exon.iv$normSpan <- exon.fraction * exon.iv$span / sum(exon.iv$span);
 #   #  intr.iv$normSpan <- (1 - exon.fraction) * intr.iv$span / sum(intr.iv$span);
 #   #}
 #   #resFunc
 #   exon.iv$normSpan <- exon.fraction * resFunc(exon.iv$span) / sum(resFunc(exon.iv$span));
 #   intr.iv$normSpan <- (1 - exon.fraction) * resFunc(intr.iv$span) / sum(resFunc(intr.iv$span));
 #   
 #   #alternate intron and exon:
 #   idx <- order(c(1:nrow(exon.iv),1:nrow(intr.iv)))
 #   iv <- rbind(exon.iv, intr.iv)[idx,, drop = FALSE];
 #   cumulativeSum <- cumsum(iv$normSpan);
 #   
 #   iv$rescale.start <- c(0,cumulativeSum[-length(cumulativeSum)]);
 #   iv$rescale.end <- cumulativeSum;
 #   iv$normSpan <- iv$rescale.end - iv$rescale.start;
 #   iv$simple.normSpan <- iv$span / sum(iv$span);
 #   iv$simple.rescale.start <- (iv$start - MIN) / SPAN;
 #   iv$simple.rescale.end   <- (iv$end - MIN) / SPAN;
 #   return(iv);
 # } else { #Catch special case: the entire gene is exonic:
 #   exon.iv$type <- "E";
 #   exon.iv$normSpan <- 1 * resFunc(exon.iv$span) / sum(resFunc(exon.iv$span));
 #   
 #   iv <- exon.iv;
 #   cumulativeSum <- cumsum(iv$normSpan);
 #   iv$rescale.start <- c(0,cumulativeSum[-length(cumulativeSum)]);
 #   iv$rescale.end <- cumulativeSum;
 #   iv$normSpan <- iv$rescale.end - iv$rescale.start;
 #   iv$simple.normSpan <- iv$span / sum(iv$span);
 #   iv$simple.rescale.start <- (iv$start - MIN) / SPAN;
 #   iv$simple.rescale.end   <- (iv$end - MIN) / SPAN;
 #   return(iv);
 # }
#}


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

drawHalfLoops <- function(dir = c("down","up","left","right"), x0,x1,y0,y1,sampleCt = 100, col = "black", lwd = 1, lty = 1, style = c("hyperbola","ellipse","triangular","line"), ...){
  dir <- match.arg(dir);
  style <- match.arg(style);
  
  len <- max(sapply(list(x0,x1,y0,y1),length));
  if(length(x0) == 1) x0 <- rep(x0,len);
  if(length(x1) == 1) x1 <- rep(x1,len);
  if(length(y0) == 1) y0 <- rep(y0,len);
  if(length(y1) == 1) y1 <- rep(y1,len);
  if(length(col) == 1) col <- rep(col, len);
  if(length(lwd) == 1) lwd <- rep(lwd,len);
  if(length(lty) == 1) lwd <- rep(lty,len);
  for(i in 1:len){
    drawHalfLoop(dir = dir, x0 = x0[i], x1 = x1[i], y0 = y0[i], y1 = y1[i], sampleCt = sampleCt, col = col[i], lwd = lwd[i], lty = lty[i], style = style, ...);
  }
}

#This handy little function draws a bracket or "half-loop" in a number of different styles. 
# \item{dir}{ The "open" direction of the half-loop. Not all directions are currently implemented for all styles. }
# \item{x0,x1,y0,y1}{ The x and y limits of the loop.}
# \item{sampleCt}{ The number of points to draw to form the shape. Higher sampleCt values will produce smoother curves. Not applicable to all styles.  }
# \item{style}{ the half-loop style. hyperbola produces a pair of mirrored hyperbolas with the ends straightened out.
#               ellipse produces an ellipse over the region.
#               triangular produces a isosceles triangle over the region.
#               line produces a simple line segment along the "open" edge.
#               Some styles will be affected by the styleFactor parameter.}
# \item{styleFactor}{
#                Currently only implemented for the hyperbola style. The default is c(8,8,0.5). The first number is the maximum X value, the second value is the maximum Y value.
#                The third value determines the proportion of the vertical space that is straightened out.
#                This function then calculates a curve for Y = 1 / X over the range c(1 / styleFactor[2], styleFactor[1]). The curve is then mirrored and refitted alongside straight
#                lines to fit in the assigned region.
#}
# \item{...}{ graphical parameters to be passed to the lines function, such as lty, lwd, and col.}
# Note this function is a bit quick-and-dirty, and is definitely NOT optimized for performance.
#
# Future functionality (WIP): "bracket" style.
# 

drawHalfLoop <- function(dir = c("down","up","left","right"), x0,x1,y0,y1, sampleCt = NULL, style = c("hyperbola","ellipse","triangular","line"), styleFactor = NULL, ...){
  dir <- match.arg(dir);
  style <- match.arg(style);
  #message("style is: \"",style);
  
  if(style == "ellipse"){
    if(is.null(sampleCt)){
      sampleCt <- 100;
    }
    if(dir == "down"){
      theta <- (0:sampleCt / sampleCt) * pi
      A <- abs(x1-x0)/2;
      B <- abs(y1-y0);
      X <- A * cos(theta) + abs(x1-x0)/2 + x0;
      Y <- B * sin(theta) + y0;
      lines(X,Y, ...)
    } else if(dir == "up"){
      theta <- ((0:sampleCt) / sampleCt) * pi + pi;
      A <- abs(x1-x0)/2;
      B <- abs(y1-y0);
      X <- A * cos(theta) + abs(x1-x0)/2 + x0;
      Y <- B * sin(theta) + y1;
      lines(X,Y, ...)
    } else {
      stop(paste0("Half-loop direction: \"",dir,"\" not supported for style \"",style,"\"."));
    }
  } else if(style == "hyperbola"){
    if(is.null(styleFactor)){
      styleFactor <- c(8,8, 0.7);
    }
    if(is.null(sampleCt)){
      sampleCt <- 100;
    }
    if(dir == "down"){
      maxRawX <- styleFactor[1];
      maxRawY <- styleFactor[2];
      curvePct <- styleFactor[3];
      minRawX <- 1 / maxRawY;
      X.raw <- ((1:sampleCt) / sampleCt) * (maxRawX - minRawX) + minRawX;
      #X.raw <- ((1:sampleCt) / sampleCt) * styleFactor;
      Y.raw <- 1 / X.raw;
      Y.raw <- c(max(Y.raw) * (1 / curvePct), Y.raw);
      Y.raw <- max(Y.raw) - Y.raw;
      X <- c(X.raw, X.raw + max(X.raw) - min(X.raw));
      X <- c(min(X), X, max(X));
      X <- X / abs(min(X) - max(X));
      X <- X - min(X);
      X <- X * abs(x1-x0) + x0;
      Y <- c(Y.raw, rev(Y.raw));
      Y <- Y / abs(min(Y) - max(Y));
      Y <- Y - min(Y);
      Y <- Y * abs(y1-y0) + y0;
      X[2] <- X[1];
      X[length(X)-1] <- X[length(X)];
      #advlines(x, y, col, lty, lwd, secondary = FALSE, secondary.col = col, secondary.alpha = 100, secondary.lty = 1, secondary.lwd = lwd / 2, ...)
      advlines(X,Y, secondary = TRUE, ...);
    } else if(dir == "up"){
      maxRawX <- styleFactor[1];
      maxRawY <- styleFactor[2];
      curvePct <- styleFactor[3];
      minRawX <- 1 / maxRawY;
      X.raw <- ((1:sampleCt) / sampleCt) * (maxRawX - minRawX) + minRawX;
      #X.raw <- ((1:sampleCt) / sampleCt) * styleFactor;
      Y.raw <- 1 / X.raw;
      Y.raw <- c(max(Y.raw) * (1 / curvePct), Y.raw);
      Y.raw <- max(Y.raw) - Y.raw;
      X <- c(X.raw, X.raw + max(X.raw) - min(X.raw));
      X <- c(min(X), X, max(X));
      X <- X / abs(min(X) - max(X));
      X <- X - min(X);
      X <- X * abs(x1-x0) + x0;
      Y <- c(Y.raw, rev(Y.raw));
      Y <- Y / abs(min(Y) - max(Y));
      Y <- Y - min(Y);
      Y <- Y * abs(y1-y0) + y0;
      X[2] <- X[1];
      X[length(X)-1] <- X[length(X)];
      Y <- y1 - Y + y0;
      lines(X,Y, ...);
    } else{
      stop(paste0("Half-loop direction: \"",dir,"\" not supported for style \"",style,"\"."));
    }
  } else if(style == "triangular"){
    if(dir == "down"){
      xM <- abs(x1-x0)/2 + x0;
      lines(c(x0,xM,x1), c(y0,y1,y0), ...);
    } else if(dir == "up"){
      xM <- abs(x1-x0)/2 + x0;
      lines(c(x0,xM,x1), c(y1,y0,y1), ...);
    } else if(dir == "left"){
      yM <- abs(y1-y0)/2 + y0;
      lines(c(x0,x1,x0), c(y0,yM,y1), ...);
    } else if(dir == "right"){
      yM <- abs(y1-y0)/2 + y0;
      lines(c(x1,x0,x1), c(y0,yM,y1), ...);
    } else {
      stop(paste0("Half-loop direction: \"",dir,"\" not supported for style \"",style,"\"."));
    }
    
    #stop(paste0("Half-loop style: \"",style,"\" not (yet) implemented."));
  } else if(style == "line"){
    if(dir == "down"){
      lines(c(x0,x1),c(y0,y0),...)
    } else if(dir == "up"){
      lines(c(x0,x1),c(y1,y1),...)
    } else if(dir == "left"){
      lines(c(x0,x0),c(y0,y1),...)
    } else if(dir == "right"){
      lines(c(x1,x1),c(y0,y1),...)
    } else {
      stop(paste0("Half-loop direction: \"",dir,"\" not supported for style \"",style,"\"."));
    }
  } else {
    stop(paste0("Half-loop style: \"",style,"\" not supported."));
  }
}


get.connection.cramp <- function(connection.lines.bottom, connection.lines.top, intron.break.threshold){
   if(is.na(intron.break.threshold)){
     FALSE;
   } else {
     #This only applies when there's lots of features:
     if(length(connection.lines.bottom) < 30){
       FALSE;
     } else {
       return(any(abs(connection.lines.top - connection.lines.bottom) > intron.break.threshold))
     }
   }
}

#add.intron.break <- function(rescale.iv, 


#Nests splices inside one another:
# Modified Version: calculates nesting when junctions overlap.
get.nested.heights <- function(tr.splice, ymin, ymax, verbose = TRUE, debug.mode = FALSE){
  tryCatch({
    if(verbose) message("Starting nested heights...");
    get.nested.heights.advanced(tr.splice,ymin,ymax,verbose,debug.mode);
  }, error = function(e){
    message("Construction of nested splices failed. Falling back to simpler nesting method...");
    message("---");
    print(e);
    message("---");
    tryCatch({
      get.nested.heights.FALLBACK(tr.splice,ymin,ymax,verbose,debug.mode);
    }, error = function(e){
      message("Fallback nested splices failed. Falling back to un-nested splices.");
      SJ.ceiling = rep(ymax,nrow(tr.splice))
      names(SJ.ceiling) <- tr.splice$featureID;
      return( list( SJ.ceiling = SJ.ceiling, maxDepth = 0 ));
    });
  });
}

#Nests splices inside one another:
# Any splice site is guaranteed to have a lower ceiling than any splice site that contains it.
# The height of a splice site's ceiling depends on the number of splice sites above and below it.
# Splice sites that mutually overlap but are not subsets of one another do not affect one another's height.

get.nested.heights.FALLBACK <- function(tr.splice, ymin, ymax, verbose = TRUE, debug.mode = FALSE){
        num.over <- sapply(1:nrow(tr.splice), function(i){
          sum(tr.splice$start[i] >= tr.splice$start & tr.splice$end[i] <= tr.splice$end ) - 1;
        })
        num.under <- sapply(1:nrow(tr.splice), function(i){
          sum(tr.splice$start[i] <= tr.splice$start & tr.splice$end[i] >= tr.splice$end ) - 1;
        });
        SJ.ceiling <- (num.under + 1) / (num.over + num.under + 1)
        SJ.ceiling <- SJ.ceiling * abs(ymax-ymin) + ymin;
        names(SJ.ceiling) <- tr.splice$featureID;
        #message("(",ymin,",",ymax,")"," = [",paste0(SJ.ceiling,collapse=","),"]");
        return(list(  SJ.ceiling = SJ.ceiling, maxDepth = max(num.over)  ));
}



get.nested.heights.advanced <- function(tr.splice, ymin, ymax, verbose = TRUE, debug.mode = FALSE){
  tr <- data.frame(ID = 1:nrow(tr.splice), start = tr.splice$start, end = tr.splice$end, span = tr.splice$end - tr.splice$start);
  
  if(debug.mode) print(tr);
  
  nhh <- get.nested.heights.helper(tr, currDepth = 1, debug.mode = debug.mode);
  if(debug.mode) print(nhh);
  
  tr <- get.nested.heights.depths(nhh, tr);
  
  num.under <- tr$depthBelow - tr$depth;
  num.over <- tr$depth;
  
  SJ.ceiling <- (num.under + 1) / (num.over + num.under + 1)
  SJ.ceiling <- SJ.ceiling * abs(ymax-ymin) + ymin;
  names(SJ.ceiling) <- tr.splice$featureID;
  #message("(",ymin,",",ymax,")"," = [",paste0(SJ.ceiling,collapse=","),"]");
  #print(cbind(tr,SJ.ceiling = SJ.ceiling));
  return(list(  SJ.ceiling = SJ.ceiling, maxDepth = max(num.over)  ));
}

repString <- function(s,n){
  paste0(rep(s,n),collapse="");
}

###WIP!!!!


get.nested.heights.helper <- function(tr, currDepth = 1, debug.mode = FALSE){
  #if(currDepth > 5) stop();

  if(debug.mode) message(repString("|",currDepth),"^^^^");
  if(debug.mode) message(repString("|",currDepth),"[", paste0(tr$ID,collapse=","),"]");
  intersections <- lapply(1:nrow(tr), function(i){
    tr$start[i] <= tr$end & tr$end[i] >= tr$start
  });
  ixMatrix <- as.matrix(do.call(rbind, intersections));
  notDiag <- as.matrix(do.call(rbind, lapply(1:nrow(tr), function(i){
    out <- rep(TRUE, nrow(tr));
    out[i] <- FALSE;
    out;
  })));
  ixMatrix <- ixMatrix & notDiag;
  numIX <- rowSums( ixMatrix )
  
  if(debug.mode) message(repString("|",currDepth),"numIX=",paste0(numIX,collapse=","));
  
  is.ungrouped <- rep(TRUE,nrow(tr));
  is.alone <- sapply(intersections, function(I){ sum(I) == 1 });
  intersection.groups <- as.list((tr$ID)[is.alone]);
  is.in.group <- lapply(intersection.groups, function(ig){ tr$ID %in% ig });
  
  is.ungrouped[is.alone] <- FALSE;
  
  if(all(is.alone)){ 
    if(debug.mode) message(repString("|",currDepth),"vvvv");
    return(tr$ID);
  }

  for(i in 1:nrow(tr)){
    if(is.ungrouped[i]){
      g <- length(is.in.group) + 1;
      ix <- intersections[[i]];
      is.in.group[[g]] <- ix;
      unchecked <- rep(TRUE, nrow(tr));
      unchecked[i] <- FALSE;
      while( any(unchecked & is.in.group[[g]]) ){
        j <- min(which(unchecked & is.in.group[[g]]));
        is.in.group[[g]] <- is.in.group[[g]] | intersections[[j]];
        unchecked[j] <- FALSE;
      }
      is.ungrouped[is.in.group[[g]]] <- FALSE;
      intersection.groups[[length(intersection.groups)+1]] <- tr$ID[is.in.group[[g]]];
      #message("new group: (",i,"): ",paste0(which(is.in.group[[g]]), collapse=","));
    }
  }
  
  if(debug.mode) message(repString("|",currDepth),"GRP: [",paste0(lapply(is.in.group, function(ig){ paste0(tr$ID[ig], collapse=",")}), collapse = "] [" ) ,"]");
  
  out <- lapply(is.in.group, function(ig){
    if(sum(ig) == 1){
      tr$ID[ig];
    } else {
      #order by # intersections, then span:
      if(debug.mode) message(repString("|",currDepth)," For GRP: [",paste0(tr$ID[ig],collapse=","),"]");
      maxIX <- max(numIX[ig]);
      #if(debug.mode) message(repString("|",currDepth)," maxIX = ", maxIX);
      if( sum(numIX == maxIX & ig) == 1){
        ex <- which(numIX == maxIX  & ig);
        #if(debug.mode) message(repString("|",currDepth)," (A) ex = tr[",ex,"] = ",tr$ID[ex]);
      } else {
        ex <- which( numIX == maxIX  & ig );
        ex <- ex[which.max(tr$span[ex])];
        #if(debug.mode) message(repString("|",currDepth)," (B) ex = tr[",ex,"] = ",tr$ID[ex]);
      }
      ig[ex] <- FALSE;
      if(debug.mode) message(repString("|",currDepth)," extract: ",tr$ID[ex]);
      if(sum(ig) == 1){
        list(parent = tr$ID[ex], children = tr$ID[ig]);
      } else {
        list(parent = tr$ID[ex], children = get.nested.heights.helper(tr[ig,,drop=FALSE], currDepth + 1, debug.mode));
      }
    }
  })
  if(debug.mode) message(repString("|",currDepth),"vvvv");
  if(length(out) == 1){
    out[[1]];
  } else {
    out;
  }
}

get.nested.heights.depths <- function(nhh, tr){
  depth <- rep(-1, nrow(tr));
  depthBelow <- rep(-1, nrow(tr));
  parents <- list();

  for(i in 1:nrow(tr)){
    currDepth <- 0;
    currList <- nhh;
    currParents <- list();
    
    while(is.list(currList)){
      if((! is.null(names(currList))) && names(currList)[1] == "parent"){
        if(i == currList$parent){
          currList <- NULL;
          depth[i] <- currDepth;
        } else {
          currParents <- c(currParents, currList$parent);
          currDepth <- currDepth + 1;
          currList <- currList$children;
        }
      } else {
        contains.i <- which(sapply(currList, function(L){ any(unlist(L) == tr$ID[i]) }))
        currList <- currList[[contains.i]];
      }
    }
    depth[i] <- currDepth;
    parents[[i]] <- currParents;
  }
  parentCt <- sapply(parents, function(p){length(p)});
  for(i in 1:nrow(tr)){
    containsi <- sapply(parents,function(p){ any(p == i) });
    if(any(containsi)){
      depthBelow[i] <- max(c(parentCt[containsi],depth[i]));
    } else {
      depthBelow[i] <- max(0,depth[i]);
    }
  }
  tr$depth <- depth;
  tr$depthBelow <- depthBelow;
  tr$parentCt <- parentCt;
  tr$lvl <- ifelse(depth == 0 & depthBelow == 0, 1, depth / depthBelow);
  tr$parents <- sapply(parents,function(p){ paste0(p,collapse=",")})
  return(tr);
}