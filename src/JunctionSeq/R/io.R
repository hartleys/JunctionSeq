

########################################################################################
### INPUT:



########################################################################################
### OUTPUT:

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
                                      group.RGB = NULL, 
                                      #nongroup.RGB = NULL, 
                                      use.score = FALSE, 
                                      FDR.threshold = 0.05, 
                                      count.digits = 1, 
                                      includeGeneID = FALSE,
                                      includeLocusID = TRUE,
                                      includeGroupID = TRUE){
  
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
  
  res.data <- res.data[keep.features,];
  estimates <- estimates[keep.features,];
  
  chrom <- res.data$chr;
  chromStart <- res.data$start;
  chromEnd <- res.data$end;
  strand <- res.data$strand;
  
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
  
  write.junction.bed.file(
    file = file,
    trackLine = trackLine,
    chrom = rep(chrom, each=numcond),
    chromStart = rep(chromStart,each=numcond),
    chromEnd = rep(chromEnd, each = numcond),
    strand = rep(strand, each = numcond),
    featureName = out.featureNames,
    featureScore = rep(score,each = numcond),
    featureRGB = rep(group.RGB, nrow(estimates))
  );
}

writeSigBedTrack <- function(file, 
                               jscs, 
                               trackLine = "track name='sigJct' description='Significant Splice Junction Loci' useScore=1 visibility=3",
                               only.sig = TRUE, only.testable = TRUE, 
                               sig.RGB = "255,0,0", nonsig.RGB = "0,0,0", 
                               use.score = TRUE, FDR.threshold = 0.05, 
                               pval.digits = 4, 
                               includeGeneID = FALSE,
                               includeLocusID = TRUE){ #fullFeatureID includeFeatureID
  res.data <- fData(jscs);
  if(only.testable){
    res.data <- res.data[res.data$testable,,drop=FALSE];
  }
  if(only.sig){
    res.data <- res.data[res.data$padjust < FDR.threshold,,drop=FALSE];
  }
  
  chrom <- res.data$chr;
  chromStart <- res.data$start;
  chromEnd <- res.data$end;
  strand <- res.data$strand;
  
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
  
  
  
  write.junction.bed.file(
    file = file,
    trackLine = trackLine,
    chrom = chrom,
    chromStart = chromStart,
    chromEnd = chromEnd,
    strand = strand,
    featureName = featureName,
    featureScore = score,
    featureRGB = featureRGB);
}

write.junction.bed.file <- function(file, trackLine = NULL, chrom, chromStart, chromEnd, strand, featureName, featureScore, featureRGB){
  out.bed <- data.frame(
  chrom = chrom,
  chromStart = chromStart - 1,
  chromEnd = chromEnd + 1,
  featureName = featureName,
  featureScore = featureScore,
  strand = strand,
  thickStart = chromStart,
  thickEnd = chromEnd,
  itemRGB = featureRGB,
  blockCount = rep(2,length(chrom)),
  blockSizes = rep("1,1",length(chrom)),
  blockStarts = paste0("0,",chromEnd - (chromStart - 1))
  );
  
  gzf <- gzfile(paste0(file,""),"w");
  if(! is.null(trackLine)){
    write.table(trackLine,gzf, quote=F,col.names=F,row.names=F);
  }
  
  write.table(out.bed, gzf, quote=F, col.names=F,row.names=F,sep='\t');
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
