
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


SUPPORTED_PLOTTING_DEVICE_LIST = c("CairoPNG","png","svg","tiff","cairo_ps");

strStartsWith <- function(s, prefix){
  substr(s,1,nchar(prefix)) == prefix;
}
f.na <- function(x){
  ifelse(is.na(x), F,x);
}

getPlottingDeviceFunc <- function(use.plotting.device, 
                                  base.plot.height, 
                                  base.plot.width, 
                                  base.plot.units = "px", 
                                  plotting.device.params = list()){
   
   if(use.plotting.device == "RSvgDevice"){
     package.found <- suppressMessages(suppressWarnings(require("RSvgDevice")));
     warning("RSvgDevice is NOT CURRENTLY SUPPORTED. Errors will likely follow.");
     if(package.found){
        plotdevfunc <- function(filename, heightMult, widthMult){
          plotting.device.params[["height"]] <- heightMult * base.plot.height;
          plotting.device.params[["width"]] <- widthMult * base.plot.width;
          plotting.device.params[["file"]] <- paste0(filename,".png");
          
          legal.params <- c("file", "height", "width", "bg", "fg", "onefile", "xmlHeader");
          plotting.device.params <- plotting.device.params[ names(plotting.device.params) %in% legal.params ]
          
          do.call(devSVG,plotting.device.params);
        }
        closefunc <- function(){
          dev.off();
        }
     } else {
       stop("Package RSvgDevice not found! Install package RSvgDevice or use a different plotting device!");
     }
   } else if(use.plotting.device == "CairoPNG"){
     cairo.package.found <- suppressMessages(suppressWarnings(require("Cairo")));
     if(cairo.package.found){
       plotdevfunc <- function(filename, heightMult, widthMult){
         plotting.device.params[["height"]] <- heightMult * base.plot.height;
         plotting.device.params[["width"]] <- widthMult * base.plot.width;
         plotting.device.params[["units"]] <- base.plot.units;
         plotting.device.params[["filename"]] <- paste0(filename,".png");
         do.call(CairoPNG,plotting.device.params);
       }
       closefunc <- function(){
         dev.off();
       }
     } else {
       stop("Package Cairo not found! Install package Cairo or use a different plotting device!");
     }
   } else if(use.plotting.device == "CairoSVG"){
     cairo.package.found <- suppressMessages(suppressWarnings(require("Cairo")));
     #warning("Note: R device CairoSVG has known issues that make labels unreadable on most renderers.");
     if(cairo.package.found){
       plotdevfunc <- function(filename, heightMult, widthMult){
         if(base.plot.units == "px"){ unitmod <- 150; } else { unitmod <- 1; }
         plotting.device.params[["height"]] <- heightMult * base.plot.height / unitmod;
         plotting.device.params[["width"]] <- widthMult * base.plot.width / unitmod;
         #plotting.device.params[["units"]] <- base.plot.units;
         plotting.device.params[["file"]] <- paste0(filename,".svg");
         do.call(CairoSVG,plotting.device.params);
       }
       closefunc <- function(){
         dev.off();
       }
     } else {
       stop("Package Cairo not found! Install package Cairo or use a different plotting device!");
     }
   } else if(use.plotting.device == "png"){
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
     #warning("Note: R device svg has known issues that make labels unreadable on most renderers.");
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

getMyApply <- function(nCores = 1, verbose = TRUE){
   if(nCores > 1){
     multicore.package.found <- suppressMessages(suppressWarnings(require("parallel")));
     BiocParallel.package.found <- suppressMessages(suppressWarnings(require("BiocParallel")));
     
     if(BiocParallel.package.found) {
       message("using package \"BiocParallel\" for parallelization. (Using ",nCores," cores)");
       myApply <- function(X, FUN){ BiocParallel::bplapply( X, FUN, BPPARAM = BiocParallel::MulticoreParam(workers = nCores) ) };
     } else if(multicore.package.found){
       message("using package \"parallel\" for parallelization. (Using ",nCores," cores)");
       myApply <- function(X, FUN){ parallel::mclapply( X, FUN, mc.cores=nCores ) };
     } else {
       warning("Package parallel not found! Set parameter nCores to 1. Falling back to single-core operation!");
       myApply <- lapply;
     }
   } else {
     myApply <- lapply;
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