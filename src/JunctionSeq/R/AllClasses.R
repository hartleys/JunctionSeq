setClass( "JunctionSeqCountSet",
   contains = "eSet",
   representation = representation(
      designColumns = "character",
      dispFitCoefs = "numeric",
      fittedMu = "matrix",
      dispFunctionType = "list",
      dispFunction = "function",
      dispFunctionJct  = "function",
      dispFunctionExon = "function",
      formulas = "list",
      annotationFile = "character",
      geneCountData = "matrix",
      countVectors = "matrix",
      altSizeFactors = "data.frame",
      plottingEstimates = "list",
      plottingEstimatesVST = "list",
      geneLevelPlottingEstimates = "list",
      modelFitForHypothesisTest = "list", #Currently unused.
      modelFitForEffectSize = "list", #Currently unused.
      flatGffData = "data.frame",
      flatGffGeneData = "list",
      analysisType = "character",
      DESeqDataSet = "DESeqDataSet",
      modelCoefficientsSample    = "list", #Currently unused.
      modelCoefficientsGene      = "list"  #Currently unused.
   ),
   prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), JunctionSeqCountSet = "0.0.5" ) ) )
)

#setMethod("show","JunctionSeqCountSet",
#   function(object){
#      cat("QoRTs_QC_Results object:\n");
#      
#      cat("Memory Usage:");
#      cat("object: ", format(object.size(object), units="Mb");
#      cat("fData(object): ",format(object.size(fData(object)),"Mb");
#      cat("pData(object): ",format(object.size(pData(object)),"Mb");
#      
#      cat("object@fittedMu: ",format(object.size(object@fittedMu),"Mb");
#      cat("counts(object): ", format(object.size(counts(object)), units="Mb");
#      cat("object@geneCountData: ", format(object.size(object@geneCountData), units="Mb");
#      cat("object@countVectors: ", format(object.size(object@countVectors), units="Mb");      
#      cat("object@flatGffData: ", format(object.size(object@flatGffData), units="Mb");
#
#      #Add more printing stuff here?
#   }
#);


#save.JunctionSeqCountSet <- function(jscs){
#  
#}

#load.JunctionSeqCountSet <- function(file){
#  
#}

makeDESeqDataSetFromJSCS <- function(jscs, test.formula1){
  countData <- jscs@countVectors;
  colData <- rbind.data.frame(
                              cbind.data.frame(data.frame(sample = rownames(pData(jscs))) , pData(jscs) ),
                              cbind.data.frame(data.frame(sample = rownames(pData(jscs))) , pData(jscs) )
  );
  colData$countbin <- factor(c(  rep("this",ncol(countData)/2)  ,   rep("others",ncol(countData)/2)   ), levels = c("this","others"));
  #print("colData:");
  #print(colData);
  for(i in 1:ncol(colData)){
    colData[[i]] <- factor(colData[[i]]);
  }
  colData <- DataFrame(colData);
  
  jscs@DESeqDataSet <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = test.formula1, ignoreRank = TRUE);
  return(jscs);
}

getModelFit <- function(jscs, featureID, geneID, countbinID, modelFitType = c("fitH0","fitH1","fitEffect")){
  if(is.null(featureID)){
    if(is.null(geneID) | is.null(countbinID)){
      stop("ERROR: getModelFit(): either featureID or geneID AND countbinID must be set!");
    } else {
      featureID <- paste0(geneID, ":", countbinID);
    }
  }
  
  out <- list();
  if(any(modelFitType == "fitH0")){
    if(is.null(jscs@modelFitForHypothesisTest)){
      stop("ERROR: getModelFit(): no modelFitForHypothesisTest found. Run testJunctionsForDJU().");
    }
    out[["fitH0"]] <- jscs@modelFitForHypothesisTest[[featureID]][["fitH0"]];
  }
  if(any(modelFitType == "fitH1")){
    if(is.null(jscs@modelFitForHypothesisTest)){
      stop("ERROR: getModelFit(): no modelFitForHypothesisTest found. Run testJunctionsForDJU().");
    }
    out[["fitH1"]] <- jscs@modelFitForHypothesisTest[[featureID]][["fitH1"]];
  }
  if(any(modelFitType == "fitEffect")){
    if(is.null(jscs@modelFitForEffectSize)){
      stop("ERROR: getModelFit(): no modelFitForHypothesisTest found. Run testJunctionsForDJU().");
    }
    out[["fitEffect"]] <- jscs@modelFitForEffectSize[[featureID]];
  }
  return(out);
}


newJunctionSeqCountSet <- function( countData, geneCountData, design, geneIDs, countbinIDs, featureIntervals=NULL, transcripts=NULL){

   countData <- as.matrix( countData )
   if( any( round( countData ) != countData ) ){
      stop( "The countData is not integer." )}
   mode( countData ) <- "integer"

   geneCountData <- as.matrix( geneCountData )
   if( any( round( geneCountData ) != geneCountData ) ){
      stop( "The geneCountData is not integer." )}
   mode( geneCountData ) <- "integer"
   
   if( is( design, "matrix" ) ){
      design <- as.data.frame( design )}

   rownames(countData) <- paste(geneIDs, countbinIDs, sep=":")
   if( any( duplicated(rownames(countData) ) ) ) {
      stop("The geneIDs or countbinIDs are not unique")
   }
   if(any(duplicated(rownames(geneCountData)))){
      stop("The geneIDs in the geneCountData are not unique");
   }
   if(any(duplicated(colnames(geneCountData)))){
      stop("The sample ID's in the geneCountData are not unique");
   }

   phenoData <- annotatedDataFrameFrom( countData, byrow=FALSE )
   featureData <- annotatedDataFrameFrom( countData, byrow=TRUE )

   phenoData$sizeFactor <- rep( NA_real_, ncol(countData) )
   varMetadata( phenoData )[ "sizeFactor", "labelDescription" ] <- "size factor (relative estimate of sequencing depth)"

   geneIDs <- as.factor( geneIDs )
   if( length(geneIDs) != nrow(countData) )
      stop( "geneIDs must be of the same length as the number of columns in countData")
   
   featureData$geneID <- geneIDs
   varMetadata( featureData )[ "geneID", "labelDescription" ] <- "ID of gene to which the exon belongs"

   countbinIDs <- as.character( countbinIDs )
   if( length(countbinIDs) != nrow(countData) )
      stop( "countbinIDs must be of the same length as the number of columns in countData")

   featureData$countbinID <- countbinIDs
   varMetadata( featureData )[ "countbinID", "labelDescription" ] <- "feature ID (unique only within a gene)"

   if( is.null(featureIntervals) ){
      featureIntervals <- data.frame(
         chr    = rep( NA_character_, nrow( countData ) ),
         start  = rep( NA_integer_,   nrow( countData ) ),
         end    = rep( NA_integer_,   nrow( countData ) ),
         strand = rep( NA_character_, nrow( countData ) ) ) }

   featureData$testable <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "testable", "labelDescription" ] <- "slot indicating if an feature should be considered in the test."
   featureData$status <- rep( "TBD", nrow( countData ) )
   
   featureData$allZero <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "allZero", "labelDescription" ] <- "slot indicating if the feature count is zero across all samples."

      
   featureData$status <- rep( "TBD", nrow( countData ) )
   varMetadata( featureData )[ "status", "labelDescription" ] <- "Feature status (either 'OK' or the reason that the feature is untestable)."
   
   featureData$baseMean <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "baseMean", "labelDescription" ] <- "Mean normalized counts across all samples."
   featureData$baseVar <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "baseVar", "labelDescription" ] <- "Simple variance of normalized counts across all samples."
   
   featureData$dispBeforeSharing <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "dispBeforeSharing", "labelDescription" ] <- "feature dispersion (Cox-Reid estimate)"

   featureData$dispFitted <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "dispFitted", "labelDescription" ] <- "Fitted mean-variance estimate.";

   featureData$dispersion <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "dispersion", "labelDescription" ] <- "Final dispersion estimate."

   featureData$pvalue <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "pvalue", "labelDescription" ] <- "p-value from testForDEU"

   featureData$padjust <- rep( NA_real_, nrow( countData ) )
   varMetadata( featureData )[ "padjust", "labelDescription" ] <- "BH adjusted p-value"
   

   featureIntervals <- as.data.frame( featureIntervals )

   # in case it was a GRanges object before, change the colname:
   if( "seqnames" %in% colnames(featureIntervals) ){
      colnames(featureIntervals)[ colnames(featureIntervals) == "seqnames" ] <- "chr"   }

   if( !all( c( "chr", "start", "end", "strand" ) %in% colnames(featureIntervals) ) ){
      stop( "featureIntervals must be a data frame with columns 'chr', 'start', 'end', and 'strand'." )}

   if(is.null(transcripts)){
      transcripts <- rep(NA_character_, nrow( countData ) )}

   if(!is(transcripts, "character")){
      stop("transcript information must be a character vector")}

   featureData$chr    <- factor( featureIntervals$chr )
   featureData$start  <- featureIntervals$start
   featureData$end    <- featureIntervals$end
   featureData$strand <- factor( featureIntervals$strand )
   featureData$transcripts <- transcripts
   varMetadata( featureData )[ "chr",    "labelDescription" ] <- "chromosome of feature"
   varMetadata( featureData )[ "start",  "labelDescription" ] <- "start of feature"
   varMetadata( featureData )[ "end",    "labelDescription" ] <- "end of feature"
   varMetadata( featureData )[ "strand", "labelDescription" ] <- "strand of feature"
   varMetadata( featureData )[ "transcripts", "labelDescription" ] <- "transcripts in which this feature is contained"

   featureData$baseMean <- rep( NA_real_, nrow( countData ) );
   varMetadata( featureData )[ "baseMean", "labelDescription" ] <- "The mean normalized counts across all samples.";
   

   if( is( design, "data.frame" ) || is( design, "AnnotatedDataFrame" ) ) {
      stopifnot( nrow( design ) == ncol( countData ) )
      stopifnot( all( unlist( lapply(design, class) ) == "factor" ) )
      design <- as( design, "AnnotatedDataFrame" )
      dimLabels(design) <- dimLabels(phenoData)
      rownames( pData(design) ) <- rownames( pData(phenoData) )
      phenoData <- combine( phenoData, design )
      rvft <- c( `_all` = NA_character_ )
      designColumns <- varLabels(design)
   } else {
      design <- factor( design, levels=unique(design))
      stopifnot( length( design ) == ncol( countData ) )
      phenoData$`condition` <- factor( design )
      varMetadata( phenoData )[ "condition", "labelDescription" ] <- "experimental condition, treatment or phenotype"
      designColumns <- "condition"
   }
   jscs <- new( "JunctionSeqCountSet",
      assayData = assayDataNew( "environment", counts=countData ),
      phenoData = phenoData,
      featureData = featureData,
      designColumns = designColumns,
      dispFitCoefs = c( NA_real_, NA_real_ ),
      geneCountData = geneCountData
      )
   jscs
}

setValidity( "JunctionSeqCountSet", function( object ) {

   if( !all( object@designColumns %in% names(pData(object)) ) )
      return( "Not all designColumns appear in phenoData." )

   if( ! "sizeFactor" %in% names(pData(object)) )
      return( "phenoData does not contain a 'sizeFactor' column.")
   if( ! is( pData(object)$`sizeFactor`, "numeric" ) )
      return( "The 'sizeFactor' column in phenoData is not numeric." )

   if( ! "geneID" %in% names(fData(object)) )
      return( "featureData does not contain a 'geneID' column.")
   if( ! is( fData(object)$geneID, "factor" ) )
      return( "The 'geneID' column in fData is not a factor." )

   if( ! "countbinID" %in% names(fData(object)) )
      return( "featureData does not contain an 'countbinID' column.")
   if( ! is( fData(object)$countbinID, "character" ) )
      return( "The 'countbinID' column in fData is not a character vector." )

   if( ! "chr"  %in% names(fData(object)) )
      return( "featureData does not contain a 'chr' column.")
   if( ! is( fData(object)$chr, "factor" ) )
      return( "The 'chr' column in fData is not a factor." )

   if( ! "start"  %in% names(fData(object)) )
      return( "featureData does not contain a 'start' column.")
   if( ! is( fData(object)$start, "integer" ) )
      return( "The 'start' column in fData is not integer." )

   if( ! all(featureNames(object) == paste(geneIDs(object), countbinIDs(object), sep=":") ) )
      return( "The featureNames do not match with the geneIDs:countbinIDs" )

   if( ! all(rownames( counts(object) ) == featureNames(object) ) )
      return( "The rownames of the count matrix do not coincide with the featureNames" )

   if( ! all(rownames( fData( object ) ) == featureNames( object ) ) )
      return( "The rownames of the featureData do not coincide with the featureNames" )

   

   if( ! "end"  %in% names(fData(object)) )
      return( "featureData does not contain a 'end' column.")
   if( ! is( fData(object)$end, "integer" ) )
      return( "The 'end' column in fData is not integer." )

   if( ! "strand"  %in% names(fData(object)) )
      return( "featureData does not contain a 'strand' column.")
   if( ! is( fData(object)$strand, "factor" ) )
      return( "The 'strand' column in fData is not a factor." )
   if( !is(fData(object)$dispersion, "numeric")){
      return( "The 'dispersion' is not numeric")}
   if( !is(fData(object)$dispFitted, "numeric")){
      return( "The 'dispFitted' is not numeric")}
   if( !is(fData(object)$dispBeforeSharing, "numeric")){
      return( "The 'dispBeforeSharing' column is not numeric")}
   if( !is(fData(object)$pvalue, "numeric")){
      return( "The 'pvalue' values are not numeric")}
   if( !is(fData(object)$padjust, "numeric")){
      return( "The 'padjust' values are not numeric")}
   if( !is.integer( assayData(object)[["counts"]] ) )
      return( "The count data is not in integer mode." )

   if( any( assayData(object)[["counts"]] < 0 ) )
      return( "The count data contains negative values." )

   if( length( object@dispFitCoefs ) != 2 )
      return( "dispFitCoefs is not a vector of length 2." )

   TRUE
} )


setMethod("counts", signature(object="JunctionSeqCountSet"),
  function( object, normalized=FALSE) {
    cds <- object
    if(!normalized){
      assayData(cds)[["counts"]]
    } else {
      if(any(is.na( sizeFactors(cds)))) {
         stop( "Please first calculate size factors or set normalized=FALSE")
      } else {
         t(t( assayData( cds )[["counts"]] ) / sizeFactors(cds) )
      }
   }
})
setReplaceMethod("counts", signature(object="JunctionSeqCountSet", value="matrix"),
  function( object, value ) {
   cds <- object
   assayData(cds)[[ "counts" ]] <- value
   validObject(cds)
   cds
})

setMethod("sizeFactors",  signature(object="JunctionSeqCountSet"),
  function(object) {
   cds <- object
   sf <- pData(cds)$sizeFactor
   names( sf ) <- colnames( counts(cds) )
   sf
})

setReplaceMethod("sizeFactors",  signature(object="JunctionSeqCountSet", value="numeric"),
  function(object, value ) {
   cds <- object
   pData(cds)$sizeFactor <- value
   validObject( cds )
   cds
})

setMethod("design", signature(object="JunctionSeqCountSet"),
   function( object, drop=TRUE, asAnnotatedDataFrame=FALSE ) {
      cds <- object
      if( asAnnotatedDataFrame )
         return( phenoData(cds)[, cds@designColumns ] )
      ans <- pData(cds)[, cds@designColumns, drop=FALSE ]
      if( ncol(ans) == 1 && drop ) {
         ans <- ans[,1]
         names(ans) <- colnames( counts(cds) ) }
      else
         rownames( ans ) <- colnames( counts(cds) )
      ans
})
setReplaceMethod("design", signature(object="JunctionSeqCountSet"),
      function( object, value ) {
      cds <- object
      ## Is it multivariate or just a vector?
      if( ncol(cbind(value)) > 1 )
         value <- as( value, "AnnotatedDataFrame" )
      else {
         value <- new( "AnnotatedDataFrame",
            data = data.frame( condition = value ) )
         varMetadata( value )[ "condition", "labelDescription" ] <-
            "experimental condition, treatment or phenotype" }

      rownames( pData(value) ) <- rownames( pData(cds) )
      dimLabels( value ) <- dimLabels( phenoData(cds) )
      phenoData(cds) <- combine(
         phenoData(cds)[ , !( colnames(pData(cds)) %in% cds@designColumns ), drop=FALSE ],
         value )
      cds@designColumns <- colnames( pData(value) )
      validObject(cds)
      cds
})

#setMethod("design",       signature(cds="JunctionSeqCountSet"), .design)
#setMethod("design<-",     signature(cds="JunctionSeqCountSet"), `.design<-`)
#setMethod("conditions",   signature(cds="JunctionSeqCountSet"), .design)
#setMethod("conditions<-", signature(cds="JunctionSeqCountSet"), `.design<-`)

geneIDs <- function( ecs ) {
   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
   g <- fData(ecs)$geneID
   names(g) <- rownames( counts(ecs) )
   g
}

`geneIDs<-` <- function( ecs, value ) {
   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
   fData(ecs)$geneID <- value
   validObject(ecs)
   ecs
}

countbinIDs <- function( ecs ) {
   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
   g <- fData(ecs)$countbinID
   names(g) <- rownames( counts(ecs) )
   g
}

`countbinIDs<-` <- function( ecs, value ) {
   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
   fData(ecs)$countbinID <- value
   validObject(ecs)
   ecs
}


subsetByGenes <- function( ecs, genes ) {
   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
   stopifnot( all( genes %in% levels(geneIDs(ecs)) ) )
   ecs2 <- ecs[ as.character(geneIDs(ecs)) %in% genes, ]
   ecs2
}

geneCountTable <- function( ecs ) {
   stopifnot( is( ecs, "JunctionSeqCountSet" ) )
   do.call( rbind,
      tapply( seq_len(nrow(ecs)), geneIDs(ecs), function(rows)
         colSums( counts(ecs)[rows,,drop=FALSE] ) ) )
}

DEUresultTable <- function(ecs)
{
   result <- data.frame(
      geneID=geneIDs(ecs),
      countbinID=countbinIDs(ecs),
      dispersion=featureData(ecs)$dispersion,
      pvalue=fData(ecs)$pvalue,
      padjust=fData(ecs)$padjust,
      baseMean=rowMeans(counts(ecs, normalized=TRUE)))

   extracol <- regexpr("log2fold", colnames(fData(ecs)))==1
   if(any(extracol)){
      w <- which(extracol)
      result <- data.frame(result, fData(ecs)[,w])
      colnames(result)[7:(6+length(w))] <- colnames(fData(ecs))[w]
   }
   result
}

