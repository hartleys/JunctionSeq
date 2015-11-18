
simpleReportMem <- function(){
     if(requireNamespace("pryr", quietly=TRUE)){
       message("Mem used:");
       print(pryr::mem_used());
     }
}

advlines <- function(x, y, col = "black", lty = 1, lwd = 1, secondary = FALSE, secondary.col = col, secondary.alpha = 100, secondary.lty = 1, secondary.lwd = lwd / 4, ...){
  if(secondary){
    lines(x, y,col= color2transparentVector(secondary.col, secondary.alpha), lty = secondary.lty, lwd = secondary.lwd, ...);
  }
  lines(x, y,col=col, lty = lty, lwd = lwd,...);
}

REPORTMEM.VERBOSE <- FALSE;

reportMem <- function(jscs){
    if(requireNamespace("pryr", quietly=TRUE)){
      message("     Total mem_used(): ",pryr::mem_used());
    }
     
    if(REPORTMEM.VERBOSE){
     #message("     object.size(jscs): ",format(object.size(jscs),units="auto"));
      message("     Memory Usage:");
      message("        jscs: ", format(object.size(jscs), units="auto"));
      message("        fData(jscs): ",format(object.size(fData(jscs)),"auto"));
      message("        pData(jscs): ",format(object.size(pData(jscs)),"auto"));
      
      message("        jscs@fittedMu: ",format(object.size(jscs@fittedMu),"auto"));
      message("        counts(jscs): ", format(object.size(counts(jscs)), units="auto"));
      message("        jscs@geneCountData: ", format(object.size(jscs@geneCountData), units="auto"));
      message("        jscs@countVectors: ", format(object.size(jscs@countVectors), units="auto"));      
      message("        jscs@flatGffData: ", format(object.size(jscs@flatGffData), units="auto"));
      message("        jscs@DESeqDataSet : ", format(object.size(jscs@DESeqDataSet), units="auto"));
    }
}

make.progress.report.fcn <- function(maxVal, numReports, reportStringPrefix){
  reportIndices <- pretty(c(1,maxVal), numReports);
  return(
    function(i){
      if(any(i == reportIndices)){
        message(paste0(reportStringPrefix,i," of ",maxVal,"","(",date(),")"));
      }
    }
  );
}

plotting.limits <- function(){
   usr <- par("usr");
   x.log <- par("xlog");
   y.log <- par("ylog");
   x.out <- c(usr[1],usr[2]);
   if(x.log) x.out <- 10 ^ x.out;
   y.out <- c(usr[3],usr[4]);
   if(y.log) y.out <- 10 ^ y.out;
   out <- c(x.out, y.out);
   return(out);
}

device.limits <- function(){
    usr <- par("usr");
    plt <- par("plt");
    x.log <- par("xlog");
    y.log <- par("ylog");

    x.plotfrac <- plt[2] - plt[1];
    x.range <- usr[2] - usr[1];
    x.adjust <- x.range / x.plotfrac;
    x.out <- c(usr[1] - (plt[1] * x.adjust), usr[2] + ((1 - plt[2]) * x.adjust));
    if(x.log) x.out <- 10 ^ x.out;

    y.plotfrac <- plt[4] - plt[3];
    y.range <- usr[4] - usr[3];
    y.adjust <- y.range / y.plotfrac;
    y.out <- c(usr[3] - (plt[3] * y.adjust), usr[4] + ((1 - plt[4]) * y.adjust));
    if(y.log) y.out <- 10 ^ y.out;

    out <- c(x.out,y.out);
    return(out);
}


#Helper function to remember what par("cxy") does:
get.character.dim <- function(){
  return(par("cxy"));
}

################################

timestamp <- function(){
   return(Sys.time());
}

timediff <- function(ts){
   return(as.numeric(Sys.time() - ts, units = "secs"));
}

getTimeAndDiff <- function(ts = NULL){
  if(is.null(ts)){
    return(paste0("[time: ",Sys.time(),"]"));
  } else {
    nts <- Sys.time();
    elapsed <- as.numeric(nts - ts, units = "secs");
    if(elapsed < 1){
      elapsed <- as.character(round(elapsed,digits=2));
    } else {
      elapsed <- floor(elapsed);
    }
    return(paste0("[time: ",nts,"],[elapsed: ",elapsed," secs]"));
  }
}

reportTimeAndDiff <- function(ts = NULL){
  message(getTimeAndDiff(ts));
}

################################################


fit.title <- function(title.text){
  plt <- par("plt");
  usr <- par("usr");
  
  x.dist <- abs(usr[2] - usr[1]);
  x.dist.frac <- min(plt[1], 1-plt[2]);
  x.frac <- (1 - x.dist.frac - x.dist.frac);
  extra.frac <- 1 + (1 - x.frac);
  out.dist <- x.dist * extra.frac;  
  
  default.width <- strwidth(title.text, cex = par("cex.main"))
  if(default.width > out.dist){
    return( fit.character.vector.helper(title.text, curr.cex = par("cex.main"), min.width = out.dist * 0.8, max.width = out.dist * 0.975, max.width.per.char = Inf) );
  } else {
    return(par("cex.main"));
  }
}

fit.vert <- function(title.text, default.cex = par("cex.ylab")){
  plt <- par("plt");
  usr <- par("usr");
  
  x.dist <- abs(usr[4] - usr[3]);
  x.frac <- abs(plt[4] - plt[3]);
  extra.frac <- 1 + (1 - x.frac);
  out.dist <- x.dist * extra.frac;  
  
  default.width <- strwidth(title.text, cex = default.cex)
  if(default.width > out.dist){
    return( fit.character.vector.helper(title.text, curr.cex = default.cex, min.width = out.dist * 0.8, max.width = out.dist, max.width.per.char = Inf) );
  } else {
    return(default.cex);
  }
}

#DEBUGMODE <- TRUE;

convertWidthToHeight <- function(x){
  width.per.inch <-  strwidth("X",units="user") / strwidth("X",units="inches");
  height.per.inch <- strheight("X",units="user") / strheight("X",units="inches");
  return(x / width.per.inch * height.per.inch);
}
convertHeightToWidth <- function(x){
  width.per.inch <-  strwidth("X",units="user") / strwidth("X",units="inches");
  height.per.inch <- strheight("X",units="user") / strheight("X",units="inches");
  return(x / height.per.inch * width.per.inch);
}


#Alternative version of strheight. Reports the height of a string when plotted with srt = 90. In other words the string width rotated by 90 degrees (in user coordinates).
strheightSRT90 <- function(s, cex = NULL, ...){
   return( strwidth(s, cex = cex, ...) *  strwidth("X",units="inches", cex = cex, ...) / strwidth("X",units="user", cex = cex, ...) * strheight("X",units="user", cex = cex, ...) / strheight("X",units="inches", cex = cex, ...) );
}

shrink.character.vector.VERT <- function(strs, curr.cex, max.height){
  curr.height <- max(strheightSRT90(strs, cex = curr.cex))
  while(curr.height > max.height){
    curr.cex <- curr.cex * 0.9;
    curr.height <- max(strheightSRT90(strs, cex = curr.cex));
    #message("Resizing. CEX = ",curr.cex, ", width = ",curr.width);
  }
  return(curr.cex);
}


shrink.character.vector <- function(strs, curr.cex, max.width){
  curr.width <- max(strwidth(strs, cex = curr.cex))
  while(curr.width > max.width){
    curr.cex <- curr.cex * 0.9;
    curr.width <- max(strwidth(strs, cex = curr.cex));
    #message("Resizing. CEX = ",curr.cex, ", width = ",curr.width);
  }
  return(curr.cex);
}

fit.character.vector <- function(strs, min.width = 0.6, max.width = 0.95, max.width.per.char = 0.15){
   curr.cex <- 1;
   return(fit.character.vector.helper(strs, curr.cex = curr.cex, min.width = min.width, max.width = max.width, max.width.per.char = max.width.per.char));
}


fit.character.vector.helper <- function(strs, curr.cex, min.width, max.width, max.width.per.char){
   strs <- strs[nchar(strs) > 0];
   
   curr.width <- max(strwidth(strs, cex = curr.cex));
   strs.nchar <- max(nchar(strs));
   
   curr.width.per.char <- curr.width / strs.nchar;
   
   desired.width <- ((max.width - min.width) / 2) + min.width;
   new.cex <- curr.cex * (desired.width / curr.width);
   
   new.width <- max(strwidth(strs, cex = new.cex));
   new.width.per.char <- new.width / strs.nchar;
   
   if(new.width.per.char > max.width.per.char){
      desired.width.per.char <- max.width.per.char;
      new.cex.perchar <- curr.cex * (desired.width.per.char / curr.width.per.char);
      return(new.cex.perchar);
   } else {
      return(new.cex);
   }
}



################################################


SUPPORTED_PLOTTING_DEVICE_LIST = c("png","x11","current","CairoPNG","svg","tiff","cairo_ps");

getPlottingDeviceFileExtension <- function(d = c("png","x11","current","CairoPNG","svg","tiff","cairo_ps") ){
  use.plotting.device <- match.arg(d);
  if(d == "png" | d == "CairoPNG"){
    return(".png");
  } else if(d == "x11" | d == "current"){
    return("");
  } else if(d == "svg"){
    return(".svg");
  } else if(d == "tiff"){
    return(".tiff");
  } else if(d == "cairo_ps"){
    return(".ps");
  } else {
    return("");
  }
}

strStartsWith <- function(s, prefix){
  substr(s,1,nchar(prefix)) == prefix;
}
f.na <- function(x){
  ifelse(is.na(x), FALSE,x);
}


overmerge.list <- function(list.old,list.new){
  list.out <- list.old;
  if(length(list.new) > 0){
    for(i in 1:length(list.new)){
      list.out[[names(list.new)[i]]] <- list.new[[i]];
    }
  }
  return(list.out);
}


getPlottingDeviceFunc <- function(use.plotting.device = c("png","x11","current","CairoPNG","svg","tiff","cairo_ps"), 
                                  base.plot.height, 
                                  base.plot.width, 
                                  base.plot.units = "px", 
                                  plotting.device.params = list()){

   use.plotting.device <- match.arg(use.plotting.device);
   if(use.plotting.device == "x11"){
       plotting.device.params <- overmerge.list(list(pointsize = 12), plotting.device.params);
       if(base.plot.units == "px"){ unitmod <- 150; 
       } else if(base.plot.units == "in"){ unitmod <- 1; 
       } else { stop("the x11 device only supports inches.") }
       
       plotdevfunc <- function(filename, heightMult, widthMult){
         plotting.device.params[["height"]] <- heightMult * base.plot.height / unitmod;
         plotting.device.params[["width"]] <- widthMult * base.plot.width / unitmod;
         do.call(x11,plotting.device.params);
       }
       closefunc <- function(){
         #do nothing
       }
   } else if(use.plotting.device == "current"){
       plotdevfunc <- function(filename, heightMult, widthMult){
         #do nothing
       }
       closefunc <- function(){
         #do nothing
       }
#   } else if(use.plotting.device == "RSvgDevice"){
#     package.found <- suppressMessages(suppressWarnings(require("RSvgDevice")));
#     warning("RSvgDevice is NOT CURRENTLY SUPPORTED. Errors will likely follow.");
#     if(package.found){
#        plotdevfunc <- function(filename, heightMult, widthMult){
#          plotting.device.params[["height"]] <- heightMult * base.plot.height;
#          plotting.device.params[["width"]] <- widthMult * base.plot.width;
#          plotting.device.params[["file"]] <- paste0(filename,".png");
#          
#          legal.params <- c("file", "height", "width", "bg", "fg", "onefile", "xmlHeader");
#          plotting.device.params <- plotting.device.params[ names(plotting.device.params) %in% legal.params ]
#          
#          do.call(RSvgDevice::devSVG,plotting.device.params);
#        }
#        closefunc <- function(){
#          dev.off();
#        }
#     } else {
#       stop("Package RSvgDevice not found! Install package RSvgDevice or use a different plotting device!");
#     }
   } else if(use.plotting.device == "CairoPNG"){
     plotting.device.params <- overmerge.list(list(pointsize = 18, res = 150), plotting.device.params);
     
     cairo.package.found <- requireNamespace("Cairo", quietly=TRUE)
     if(cairo.package.found){
       plotdevfunc <- function(filename, heightMult, widthMult){
         plotting.device.params[["height"]] <- heightMult * base.plot.height;
         plotting.device.params[["width"]] <- widthMult * base.plot.width;
         plotting.device.params[["units"]] <- base.plot.units;
         plotting.device.params[["filename"]] <- paste0(filename,".png");
         do.call(Cairo::CairoPNG,plotting.device.params);
       }
       closefunc <- function(){
         dev.off();
       }
     } else {
       stop("Package Cairo not found! Install package Cairo or use a different plotting device!");
     }
   } else if(use.plotting.device == "CairoSVG"){
     plotting.device.params <- overmerge.list(list(pointsize = 18), plotting.device.params);
     
     cairo.package.found <- requireNamespace("Cairo", quietly=TRUE)
     #warning("Note: R device CairoSVG has known issues that make labels unreadable on certain renderers.");
     if(cairo.package.found){
       plotdevfunc <- function(filename, heightMult, widthMult){
         if(base.plot.units == "px"){ unitmod <- 150; } else { unitmod <- 1; }
         plotting.device.params[["height"]] <- heightMult * base.plot.height / unitmod;
         plotting.device.params[["width"]] <- widthMult * base.plot.width / unitmod;
         #plotting.device.params[["units"]] <- base.plot.units;
         plotting.device.params[["file"]] <- paste0(filename,".svg");
         do.call(Cairo::CairoSVG,plotting.device.params);
       }
       closefunc <- function(){
         dev.off();
       }
     } else {
       stop("Package Cairo not found! Install package Cairo or use a different plotting device!");
     }
   } else if(use.plotting.device == "png"){
     plotting.device.params <- overmerge.list(list(pointsize = 18, res = 150), plotting.device.params);
     if(capabilities()[["png"]]){
       plotdevfunc <- function(filename, heightMult, widthMult){
         plotting.device.params[["height"]] <- heightMult * base.plot.height;
         plotting.device.params[["width"]] <- widthMult * base.plot.width;
         plotting.device.params[["units"]] <- base.plot.units;
         plotting.device.params[["filename"]] <- paste0(filename,".png");
         do.call(png,plotting.device.params);
       }
       closefunc <- function(){
         dev.off();
       }
     } else {
       stop("png functionality disabled on this installation of R. Reinstall/recompile R with png support, or use a different plotting device!");
     }
   } else if(use.plotting.device == "tiff"){
     plotting.device.params <- overmerge.list(list(pointsize = 18, res = 150, compression = "lzw"), plotting.device.params);
     if(capabilities()[["tiff"]]){
       plotdevfunc <- function(filename, heightMult, widthMult){
         plotting.device.params[["height"]] <- heightMult * base.plot.height;
         plotting.device.params[["width"]] <- widthMult * base.plot.width;
         plotting.device.params[["units"]] <- base.plot.units;
         plotting.device.params[["filename"]] <- paste0(filename,".tiff");
         do.call(tiff,plotting.device.params);
       }
       closefunc <- function(){
         dev.off();
       }
     } else {
       stop("tiff functionality disabled on this installation of R. Reinstall/recompile R with tiff support, or use a different plotting device!");
     }
   } else if(use.plotting.device == "svg"){
     #message("Note: the R svg device relies on certain external software packages that have known errors on certain versions of linux.");
     plotting.device.params <- overmerge.list(list(pointsize = 18), plotting.device.params);
     if(capabilities()[["cairo"]]){
       plotdevfunc <- function(filename, heightMult, widthMult){
         if(base.plot.units == "px"){ unitmod <- 150; } else { unitmod <- 1; }
         plotting.device.params[["height"]] <- heightMult * base.plot.height / unitmod;
         plotting.device.params[["width"]] <- widthMult * base.plot.width / unitmod;
         #plotting.device.params[["units"]] <- base.plot.units;
         plotting.device.params[["filename"]] <- paste0(filename,".svg");
         do.call(svg,plotting.device.params);
       }
       closefunc <- function(){
         dev.off();
       }
     } else {
       stop("cairo svg functionality disabled on this installation of R. Reinstall/recompile R with svg support, or use a different plotting device!");
     }
   } else if(use.plotting.device == "cairo_ps"){
     plotting.device.params <- overmerge.list(list(pointsize = 18), plotting.device.params);
     if(capabilities()[["cairo"]]){
       plotdevfunc <- function(filename, heightMult, widthMult){
         if(base.plot.units == "px"){ unitmod <- 150; } else { unitmod <- 1; }
         plotting.device.params[["height"]] <- heightMult * base.plot.height / unitmod;
         plotting.device.params[["width"]] <- widthMult * base.plot.width / unitmod;
         #plotting.device.params[["units"]] <- base.plot.units;
         plotting.device.params[["filename"]] <- paste0(filename,".ps");
         do.call(cairo_ps,plotting.device.params);
       }
       closefunc <- function(){
         dev.off();
       }
     } else {
       stop("cairo cairo_ps functionality disabled on this installation of R. Reinstall/recompile R with cairo_ps support, or use a different plotting device!");
     }
   } else {
     stop("Unrecognized plotting device name: ",use.plotting.device,"\n   Supported devices are: [",paste0(SUPPORTED_PLOTTING_DEVICE_LIST, collapse = ","),"]");
   }
   return(list(plotdevfunc,closefunc));
}

getMyApply <- function(nCores = 1, verbose = TRUE, allowWindowsMulticore = TRUE, testCapability = TRUE){
   if(! is.numeric(nCores)){
     if(is(nCores, "BiocParallelParam")){
       myApply <- function(X, FUN){ BiocParallel::bplapply( X, FUN, BPPARAM = nCores ) };
     } else {
       stop("Fatal Error: nCores must be either an integer or a BiocParallelParam object");
     }
   } else if(nCores > 1){
     #multicore.package.found <- suppressMessages(suppressWarnings(require("parallel")));
     #The older version was "enhanced" by BiocParallel. It is now mandatory as of JunctionSeq v0.3.72.
     #BiocParallel.package.found <- TRUE; #suppressMessages(suppressWarnings(require("BiocParallel")));
    
     if( Sys.info()[['sysname']] == 'Windows' ){
       message(">>> NOTE: Microsoft windows detected. As of BiocParallel v1.2.0 and R 3.1.1, simple multicore forking is not supported on windows. ");
       message("          JunctionSeq will fall-back to single-core operation if necessary.");
     } 
     
     #if(BiocParallel.package.found) {
       message("    [[Using package \"BiocParallel\" for parallelization. (Using ",nCores," cores)]]");
       if( Sys.info()[['sysname']] == 'Windows' ){
         message(">>> WARNING: attempting to use BiocParallel for multicore functionality. However: On windows machines some versions of BiocParallel appear to run very slowly and do not appear to actually use multiple cores.");
       }
       myApply <- function(X, FUN){ BiocParallel::bplapply( X, FUN, BPPARAM = BiocParallel::MulticoreParam(workers = nCores) ) };
     #} else if(multicore.package.found){
     #  if( Sys.info()[['sysname']] == 'Windows' ){
     #    message(">>> WARNING: Microsoft windows detected, and package BiocParallel not found, and nCores > 1:");
     #    message(">>>    Currently, multicore operation is not supported on windows without the BiocParallel package!");
     #    message(">>>    Falling back to single-core operation!");
     #    warning("Multicore lapply unavailable. Falling back to single-core operations.");
     #    myApply <- lapply;
     #  } else {
     #    message("    [[Using package \"parallel\" for parallelization. (Using ",nCores," cores)]]");
     #    myApply <- function(X, FUN){ parallel::mclapply( X, FUN, mc.cores=nCores ) };
     #  }
     #} else {
     #  message(">>> WARNING: Neither package BiocParallel nor package parallel found, but nCores > 1.");
     #  message(">>>          Paralell operations not supported without these packages!");
     #  message(">>>          Falling back to single-core operation!");
     #  warning("Multicore lapply unavailable. Falling back to single-core operations.");
     #  myApply <- lapply;
     #}
   } else {
     myApply <- lapply;
   }
   
   if(testCapability){
     myApply <- tryCatch({
         test.run <- myApply(1:10, FUN = function(x){2 * x});
         myApply;
       }, error = function(e){
         message(">>> WARNING: Attempted to run some form of multicore lapply, but it threw an error (likely due to an OS conflict).");
         message("    Error follows: ",e);
         message(">>> Falling back to single-core operations!");
         warning("Multicore lapply unavailable. Falling back to single-core operations.");
         lapply;
       }
     )
   }
   
   return(myApply);
}

##########################################################################
######### Front-end Utility Functions:
##########################################################################

add.gene.name.annotation <- function(geneIDs,geneNames,merged.results.data){
  out.gene.names <- sapply(merged.results.data$geneID,FUN=function(x){
    gene.index <- which(x == geneIDs);
    if(length(gene.index) == 0){
      return("-");
    } else if(length(gene.index) == 1){
      return(geneNames[gene.index])
    } else {
      return("ERROR!");
    }
  });
  
  merged.results.data$geneName <- out.gene.names;
  return(merged.results.data);
}

get.sig.genes <- function(merged.results.data, FDR.threshold=0.01){
  sig.features <- which(merged.results.data$padjust < 0.001);
  sig.genes    <- unique(as.character(merged.results.data$geneID[sig.features]))
}


make.evenly.spaced.seq <- function(start,end, approx.ct){
   len <- end - start;
   tickLen <- signif(len / approx.ct,digits=1);
   
   return(seq(floor(start / tickLen)*tickLen,ceiling(end / tickLen)*tickLen,by=tickLen));
}
make.evenly.spaced.seq.minor <- function(start,end, approx.ct, num.ticks.per.main.tick = 10){
   len <- end - start;
   tickLen <- signif(len / approx.ct,digits=1);
   
   return(seq(floor(start / tickLen)*tickLen,ceiling(end / tickLen)*tickLen,by= tickLen / num.ticks.per.main.tick));
}


##########################################################################
##########################################################################
##########################################################################
######### UTILITY FUNCTIONS:
##########################################################################
##########################################################################
##########################################################################

color2transparentVector <- function(c,t){
   sapply(c, FUN = color2transparent, t = t);
}
color2transparent <- function(c,t){
   r <- col2rgb(c,alpha=TRUE);
   return(rgb(r[1],r[2],r[3],t,maxColorValue = 255));
}

apply2 <- function( X, MARGIN, FUN, ... ) {
   if( length(MARGIN) > 0 ) 
      apply( X, MARGIN, FUN, ... ) 
   else 
      FUN( X, ... ) }
      
      
INTERNAL.NINF.VALUE <- -0.2;