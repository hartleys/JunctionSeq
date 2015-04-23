
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

makevstaxis <- function(min, ylimn, ecs, use.vst, use.log,plot.type, par.cex = 1, anno.cex.text = 1, ...)
{
   #message("makevstaxis cex = ",cex);

   if(use.vst){
   
  #  ylimn <- c(0.0156,39.89)

     range_scaled <- ylimn[2] - ylimn[1];

     decades_unscaled <- c(1);
     decades_scaled <- vst(decades_unscaled,ecs);
     while(decades_scaled[1] > ylimn[1]){
       decades_unscaled <- c(decades_unscaled[1] / 10, decades_unscaled);
       decades_scaled <- vst(decades_unscaled,ecs);
     }
     while(decades_scaled[length(decades_scaled)] < ylimn[2]){
       decades_unscaled <- c(decades_unscaled, decades_unscaled[length(decades_unscaled)] * 10);
       decades_scaled <- vst(decades_unscaled,ecs);
     }
     first.decade.threshold <- (range_scaled * 0.025) + ylimn[1];
     while(decades_scaled[1] < first.decade.threshold){
       decades_scaled <- decades_scaled[2:length(decades_scaled)];
       decades_unscaled <- decades_unscaled[2:length(decades_unscaled)];
     }

     decades_unscaled <- c(0,decades_unscaled);
     decades_scaled <- vst(decades_unscaled,ecs);

     ticks_unscaled <- c();
     first.tick.threshold <- (range_scaled * 0.05) + ylimn[1];
     
     for(i in 1:length(decades_unscaled)){
       if(decades_scaled[i] >  first.tick.threshold) {
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
   
     axis( 2, at=ticks_scaled, labels=FALSE, las=2, pos=0,tcl=0.25, cex.axis = anno.cex.text, ...)
     axis( 2, at=decades_scaled, labels=decades_unscaled, las=2, pos=0,tcl=0.5, cex.axis = anno.cex.text, ...)

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
      axis(2,at=c(INTERNAL.NINF.VALUE,0),labels=c(0,""), las=2, tcl=-.25, pos=0, cex.axis = anno.cex.text, ...)
      qorts.axis.break(2,breakpos=(INTERNAL.NINF.VALUE / 2),pos=0,...)
      
      #axis( 2, at=ticks_scaled, labels=FALSE, las=2, pos=0,tcl=0.25, cex.axis = anno.cex.text, ...)
      #axis( 2, at=decades_scaled, labels=decades_unscaled, las=2, pos=0,tcl=0.5, cex.axis = anno.cex.text, ...)
      
      axis(2,   at=logticks ,labels=FALSE, las=2, pos=0, tcl=0.25, cex.axis = anno.cex.text, ...)
      axis(2,   at=0:ceiling(ylimn[2]),labels=(10^(0:ceiling(ylimn[2]))), las=2, tcl=0.5, pos=0, cex.axis = anno.cex.text, ...)
    } else {
      ticks <- rep(c(2,3,4,5,6,7,8,9),ceiling(ylimn[2])) * rep(10^(0:ceiling(ylimn[2]-1)),each=8)
      logticks <- log10(ticks);
      axis(2,   at=c(INTERNAL.NINF.VALUE,0),labels=c("0",""),      las=2, tcl=0.5,  pos=0, cex.axis = anno.cex.text, ...)
      qorts.axis.break(2,breakpos=(INTERNAL.NINF.VALUE / 2),pos=0,...)
      #axis(2,   at = ((1:9)/10)*INTERNAL.NINF.VALUE, labels=FALSE, las=2, tcl=0.25, pos=0, cex.axis = anno.cex.text, ...);
      
      axis(2,   at=logticks ,labels=FALSE, las=2, pos=0, tcl=0.25, cex.axis = anno.cex.text, ...)
      axis(2,   at=0:ceiling(ylimn[2]),labels=(10^(0:ceiling(ylimn[2]))), las=2, tcl=0.5, pos=0, cex.axis = anno.cex.text, ...)
    }
  } else {
    axis(2, las=2, tcl=0.5, pos=0, cex.axis = anno.cex.text, ...);
  }
}

#######################
#FUNCTION TO DRAW THE EXPRESSION PLOTS:
#######################
drawPlot <- function(matr, ylimn,ecs, intervals, rango, fitExpToVar, numexons, textAxis, rt, color.count, 
                     colorlines,countbinIDs,use.vst, use.log,plot.type,main.title,draw.legend,color.key,
                     condition.names,p.values=NULL,draw.p.values=FALSE,plot.lwd = 1,axes.lwd = 1, 
                     anno.lwd =1, par.cex = 1, anno.cex.text = 1, anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                     fit.countbin.names = TRUE, fit.pvals = TRUE, stagger.pvals = TRUE,
                     ...)
{
   #message("drawPlot cex = ",cex);
   plot.new()
   par(cex = par.cex);
   
   #message(paste("ylimn:",ylimn));
   plot.window(xlim=c(0, 1), ylim=ylimn)
   makevstaxis(1/ncol(matr), ylimn, ecs,use.vst=use.vst, use.log=use.log,plot.type=plot.type, lwd = axes.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, ...)
   title(main=main.title, cex = anno.cex.text, ...);

   intervals<-(0:nrow(matr))/nrow(matr)
   middle <- apply(cbind(intervals[rango], (intervals[rango+1]-((intervals[rango+1])-intervals[rango])*0.2)), 1, median)
   matr <- rbind(matr, NA)
   j <- 1:ncol(matr)
   segments(0,0,2,0,lty="dotted",col="black", lwd = plot.lwd,...);
   #abline(h = 0,lty="dotted",col="black", lwd = plot.lwd,...);
   segments(intervals[rango],matr[rango,j], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), matr[rango,j], col=color.count, lwd = plot.lwd,...)  #### line with the y level
   segments(intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2), matr[rango,j], intervals[rango+1], matr[rango+1,j], col=color.count, lty="dotted", lwd =plot.lwd,...)  #### line joining the y levels
   abline(v=middle[rango], lty="dotted", col=colorlines, lwd = plot.lwd)
   mtext(textAxis, side=2, adj=0.5, line=1.5, outer=FALSE, cex = anno.cex.text * par("cex"),...)
   
   cex.countbinIDs <- anno.cex.axis;
   if(fit.countbin.names){
     bin.width <- abs(intervals[1] - intervals[2]);
     cex.countbinIDs <- fit.character.vector(countbinIDs, min.width = 0.6 * bin.width, max.width = 0.8 * bin.width, max.width.per.char = 0.3 * bin.width)
     cex.countbinIDs <- min(cex.countbinIDs, anno.cex.axis);
   }
   
   if(length(rt) > 15){
     axis(1, at=middle[1:length(rt)], labels=countbinIDs, tcl = 0, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 3, padj=0.5,hadj=1, line = 0, mgp = c(3,0,0),...)
   } else {
     axis(1, at=middle[1:length(rt)], labels=countbinIDs, tcl = 0, lwd.ticks = axes.lwd, cex.axis = cex.countbinIDs, las = 0, padj=0.5,hadj=0.5,...)
   }
   
   
   #par usr: (x1, x2, y1, y2)
   #rect(xleft, ybottom, xright, ytop
   rect(0, par("usr")[3], par("usr")[2], par("usr")[4], xpd=FALSE, lwd = axes.lwd,...);
   
   if(draw.p.values){
      if(fit.pvals){
        bin.width <- abs(intervals[1] - intervals[2]);
        cex.pvalues <- fit.character.vector(p.values, min.width = 1 * bin.width, max.width = 2 * bin.width, max.width.per.char = 1 * bin.width)
        cex.pvalues <- min(cex.pvalues, anno.cex.text);
      } else {
        cex.pvalues <- anno.cex.text;
      }
      pval.height <- max(strheight(p.values, cex = cex.pvalues));
      
      if(stagger.pvals){
        pval.y <- rep(c(ylimn[2], ylimn[2] - pval.height, ylimn[2] - (2 * pval.height)), length(middle[rango]))[rango]
      } else {
        pval.y <- rep(ylimn[2],rango)
      }
      
      text(middle[rango],pval.y,labels=p.values,col=colorlines, cex = cex.pvalues);
   }

   if(draw.legend){
     legend("topright",fill=color.key,legend=condition.names, cex = anno.cex.text, bg = "transparent");
   }
   #if(draw.box){
   #   box(lwd = plot.lwd, ...);
   #}
}

#########################
#FUNCTION TO DRAW THE GENE MODELS AND TRANSCRIPTS:
#########################
drawGene <- function(minx, maxx, tr, tr.allExon, rango, rescale.iv = NULL, exoncol=NULL,allExon.exonCol=NULL, names, 
                     trName, exonlty, anno.lwd = 1, show.strand.arrows = 0, par.cex = 1, anno.cex.text = 1, arrows.length = 0.5,
                     ...)
{
   #message("drawGene cex = ",cex);
   #par(cex = par.cex);
   plot.new()
   plot.window(xlim=c(minx, maxx), ylim=c(0, 1))
   ymax <- par("usr")[4];
   ymin <- par("usr")[3];
   
   rect.floor <- 0.05;
   
   if(! is.null(rescale.iv)){
     tr.allExon$start <- (rescale.coords(tr.allExon$start,rescale.iv) * (maxx-minx)) + minx;
     tr.allExon$end   <- (rescale.coords(tr.allExon$end,rescale.iv) * (maxx-minx)) + minx;
     tr$start <- (rescale.coords(tr$start,rescale.iv) * (maxx-minx)) + minx;
     tr$end   <- (rescale.coords(tr$end,rescale.iv) * (maxx-minx)) + minx;
   }
   
   rect(tr.allExon$start, rect.floor, tr.allExon$end, 0.75, lty = 1, lwd = anno.lwd,col=allExon.exonCol,...)
   lines(c(minx,maxx),c(0.4,0.4),lty= 1, lwd = anno.lwd);
   
   if(show.strand.arrows != 0){
      #if(is.null(strand)){
      #   message("Warning: Cannot draw strand arrows, strandedness is NULL! Skipping.");
      #} else {
         strand <- tr.allExon$strand[1];
         
         if(strand == "+"){
            even.spacing <- ((0:(show.strand.arrows - 1)) / (show.strand.arrows)) * (maxx - minx) + minx
            arrows( even.spacing[-(show.strand.arrows)] ,rep(0.4,show.strand.arrows - 1), even.spacing[-1]  , rep(0.4,show.strand.arrows - 1), lwd = anno.lwd, length = arrows.length);
         } else if(strand == "-"){
            even.spacing <- ((1:(show.strand.arrows)) / (show.strand.arrows)) * (maxx - minx) + minx
            arrows( even.spacing[-1] ,rep(0.4,show.strand.arrows - 1), even.spacing[-(show.strand.arrows)]  , rep(0.4,show.strand.arrows - 1), lwd = anno.lwd, length = arrows.length);
         }
      #}
   }

   if( any(tr$is.exon) ){
      tr.ex <- tr[tr$is.exon, ]
      middle.ex <- apply(rbind(tr.ex$start,tr.ex$end),2,median);
      colors.ex <- exoncol[tr$is.exon]
      segments(middle.ex, 0.75, middle.ex, 1, col=colors.ex, lty = exonlty[tr$is.exon], lwd = anno.lwd,...);
   }
   if(any( ! tr$is.exon)){
      tr.splice <- tr[! tr$is.exon, ];
      #tr.splice <- tr
      start.splice <- tr.splice$start
      end.splice <- tr.splice$end
      middle.splice <- apply(rbind(start.splice, end.splice), 2, median)
      colors.splice <- exoncol[! tr$is.exon]
      segments(start.splice, 0.75, middle.splice, ymax, col=colors.splice, lty = exonlty[! tr$is.exon], lwd = anno.lwd,xpd=FALSE,...);
      segments(middle.splice, ymax, end.splice, 0.75, col=colors.splice, lty = exonlty[! tr$is.exon], lwd = anno.lwd,xpd=FALSE,...);
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

#drawTranscript(rel.calc.min, rel.calc.max, tr=tr, rango=1:(length(anno.data$start[rt.allExon][logicexons==1])), exoncol=NULL, names=c(), trName=trans[i], cex=0.8,sub.sig = sub.sig)

drawTranscript <- function(minx, maxx, tr, rango, rescale.iv = NULL, exoncol=NULL, names, trName, sub.sig, anno.lwd = 1, par.cex = 1, anno.cex.text = 1, ...)
{
   #message("drawTranscript cex = ",cex);
   plot.new()
   par(cex = par.cex);
   plot.window(xlim=c(minx, maxx), ylim=c(0, 1),...)
   
   if(! is.null(rescale.iv)){
     tr$start <- (rescale.coords(tr$start,rescale.iv) * (maxx-minx)) + minx;
     tr$end   <- (rescale.coords(tr$end,  rescale.iv) * (maxx-minx)) + minx;
   }
   
   t.splices <- cbind(tr[rango, 2], tr[rango+1, 1]);
   color.sig <- rep(FALSE,length(tr[rango,1]));
   is.non.na <- ! (is.na(t.splices[,1]) | is.na(t.splices[,2]));
   if(length(sub.sig$start) > 0){
   for(i in 1:(length(sub.sig$start))){
      curr.sig <- sub.sig[i,];
      color.sig <- color.sig | (is.non.na & sub.sig$start[i] == t.splices[,1] & sub.sig$end[i] == t.splices[,2]);
   }}
   splice.color <- ifelse(color.sig,"#F219ED","black");

   rect(tr[rango,1], 0, tr[rango,2], 1, col=exoncol)
   zr <- apply(rbind(tr[rango, 2], tr[rango+1, 1]), 2, median)

   segments(tr[rango,2], 0.5, zr, 0.65,col=splice.color, lwd = anno.lwd,...)
   segments(zr, 0.65, tr[rango+1,1], 0.5,col=splice.color, lwd = anno.lwd,...)
   
   if(! is.null(trName)){
     #text(minx - ((maxx - minx)*0.02),0.5,labels=trName,srt=45,xpd=FALSE,cex=0.5);
     text(minx - ((maxx - minx)*0.02),0.5,labels=trName,srt=45,xpd=TRUE,cex=0.5 * anno.cex.text,...);
   }
}

#############################################


generate.interval.scale <- function(gene.features, exon.fraction){
  MAX <- max(c(gene.features$start, gene.features$end));
  MIN <- min(c(gene.features$start, gene.features$end));
  SPAN <- MAX - MIN;
  
  exon.iv <- data.frame(start = gene.features$start[gene.features$is.exon], 
                          end = gene.features$end[gene.features$is.exon]);
  for(i in 2:nrow(exon.iv)){
    if(i > nrow(exon.iv)) break;
    if(exon.iv$end[i-1] == exon.iv$start[i]){
      exon.iv$end[i-1] = exon.iv$end[i];
      exon.iv <- exon.iv[-i,];
    }
  }
  intr.iv <- data.frame(start = exon.iv$end[-nrow(exon.iv)], end = exon.iv$start[-1]);
  if(max(exon.iv$end) < MAX){
    intr.iv <- rbind(intr.iv,c(start = max(exon.iv$end), end = MAX));
  }
  if(min(exon.iv$start) > MIN){
    intr.iv <- rbind(c(start = MIN, end = min(exon.iv$start)),intr.iv);
  }
  exon.iv$span <- exon.iv$end - exon.iv$start;
  intr.iv$span <- intr.iv$end - intr.iv$start;
  if(is.na(exon.fraction) | exon.fraction <= 0 | exon.fraction >= 1){
    idx <- order(c(1:nrow(exon.iv),1:nrow(intr.iv)))
    iv <- rbind(exon.iv, intr.iv)[idx,];
    iv$normSpan <- iv$span / sum(iv$span);
    iv$rescale.start <- (iv$start - MIN) / SPAN;
    iv$rescale.end   <- (iv$end - MIN) / SPAN;
    return(iv);
  } else {
    exon.iv$normSpan <- exon.fraction * exon.iv$span / sum(exon.iv$span);
    intr.iv$normSpan <- (1 - exon.fraction) * intr.iv$span / sum(intr.iv$span);
    idx <- order(c(1:nrow(exon.iv),1:nrow(intr.iv)))
    iv <- rbind(exon.iv, intr.iv)[idx,];
    cumulativeSum <- cumsum(iv$normSpan);
    iv$rescale.start <- c(0,cumulativeSum[-length(cumulativeSum)]);
    iv$rescale.end <- cumulativeSum;
    iv$simple.normSpan <- iv$span / sum(iv$span);
    iv$simple.rescale.start <- (iv$start - MIN) / SPAN;
    iv$simple.rescale.end   <- (iv$end - MIN) / SPAN;

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
    breakcol = "black", style = "slash", brw = 0.02,...) 
{
    figxy <- par("usr")
    xaxl <- par("xlog")
    yaxl <- par("ylog")
    xw <- (figxy[2] - figxy[1]) * brw
    yw <- (figxy[4] - figxy[3]) * brw
    if (!is.na(pos)) 
        figxy <- rep(pos, 4)
    if (is.null(breakpos)) 
        breakpos <- ifelse(axis%%2, figxy[1] + xw * 2, figxy[3] + 
            yw * 2)
    if (xaxl && (axis == 1 || axis == 3)) 
        breakpos <- log10(breakpos)
    if (yaxl && (axis == 2 || axis == 4)) 
        breakpos <- log10(breakpos)
    switch(axis, br <- c(breakpos - xw/2, figxy[3] - yw/2, breakpos + 
        xw/2, figxy[3] + yw/2), br <- c(figxy[1] - xw/2, breakpos - 
        yw/2, figxy[1] + xw/2, breakpos + yw/2), br <- c(breakpos - 
        xw/2, figxy[4] - yw/2, breakpos + xw/2, figxy[4] + yw/2), 
        br <- c(figxy[2] - xw/2, breakpos - yw/2, figxy[2] + 
            xw/2, breakpos + yw/2), stop("Improper axis specification."))
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
        rect(br[1], br[2], br[3], br[4], col = bgcol, border = bgcol)
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
    segments(xbegin, ybegin, xend, yend, col = breakcol, lty = 1, xpd = FALSE, ...)
    #par(xpd = FALSE)
}
