
##########################################################################
######### Expanding, Extracting, and Annotating DEXSeq results
##########################################################################

#      designColumns = "character",
#      dispFitCoefs = "numeric",
#      formulas = "list",
#      annotationFile = "character",
#      geneCountData = "matrix",
#      countVectors = "matrix",
#      plottingEstimates = "list",
#      plottingEstimatesVST = "list",
#      modelFitForHypothesisTest = "list",
#      modelFitForEffectSize = "list",
#      modelCoefficientsSample    = "list", #Currently unused.
#      modelCoefficientsGene      = "list"  #Currently unused.
#writeCompleteResults(jscs,outfile.prefix="/cluster/ifs/users/mullikin/Klein/steve/projects/14-09-kleinPineal5/JunctionSeq/run2/pos/Ctrl_DvN/testOut/test");


#DEPRECIATED:
generateCompleteResults <- function(jscs, outfile.prefix = NULL, verbose = TRUE, nCores=1, gzip.output = TRUE
            ){
   message("WARNING: function generateCompleteResults is DEPRECIATED! Use writeCompleteResults instead!");
   warning("WARNING: function generateCompleteResults is DEPRECIATED! Use writeCompleteResults instead!");
   if(verbose){
       message("> STARTING generateCompleteResults (",date(),")");
       message(paste("> gcr: nCores:",nCores));
       message(paste("> gcr: gzip.output:",gzip.output));
   }
   keep.debug.model.data <- TRUE;
   
   #require(statmod)
   #require(DEXSeq)
   ##require(multicore)
   if(verbose) { message(paste0("> gcr: Starting generateCompleteResults: ",date(),".")) }
   if(verbose) { message(paste0("> gcr: generating.all.expression.estimates... ",date(),".")) }
   allExprEstOut <- generateAllExpressionEstimates(jscs,verbose=verbose,nCores=nCores);
   jscs <- allExprEstOut[[1]];
   expr.data <- allExprEstOut[[2]];
   if(verbose) { message(paste0("> gcr: generating.all.expression.estimates done.",date(),".")) }
   
   if(! is.null(outfile.prefix)) write.table.gz(expr.data,file=paste(outfile.prefix,"expression.data.txt",sep=''),use.gzip=gzip.output,row.names=F,col.names=T,quote=F);
   if(verbose) message(paste0("> gcr: Written Expression Estimates. ",date()));

   if(verbose) message(paste0("> gcr: merge 1: ",date()));
   if(! all( expr.data$featureID == row.names(featureData(jscs)@data) )){
      stop("FATAL ERROR: all( expr.data$featureID == row.names(featureData(jscs)@data) ) is not TRUE!");
   }
   fn <- names(expr.data)
   fd <- featureData(jscs)@data;
   featureChar <- substr(fd$countbinID,1,1);
   #fd$baseMean <- rowMeans(counts(jscs, normalized=TRUE))
   fd$featureType <- ifelse(featureChar == "E","exonic_part",ifelse(featureChar == "J", "splice_site", "novel_splice_site"))
   names(fd)[names(fd) == "chr"] <- "chrom";
   
   merged.data.1 <- cbind.data.frame(featureID = expr.data[,1], fd, expr.data[,4:length(fn)])
   
   if(verbose) message(paste0("> gcr: merge 1 complete ",date()));
   
   if(keep.debug.model.data){
     merged.data.1 <- cbind.data.frame(merged.data.1, jscs@countVectors);
   }

   if(! is.null(outfile.prefix)){
      write.table.gz(merged.data.1,file=paste(outfile.prefix,"complete.results.txt",sep=''),use.gzip = gzip.output,quote=F,row.names=F)
   } else {
      if(verbose) { message(paste0("> gcr: skipping write of merged data, because outfile.prefix is NULL. ",date(),".")) }
   }
   
   if(verbose) { message(paste0("> gcr: done. ",date())) }
   return(list(jscs,merged.data.1));
}
