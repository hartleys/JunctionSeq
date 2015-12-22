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
      plottingEstimatesVST = "list", #DEPRECIATED! VST-xform is fast enough that it's better to calculate them as needed.
      geneLevelPlottingEstimates = "list",
      modelFitForHypothesisTest = "list", #USUALLY unused.
      modelFitForEffectSize = "list", #USUALLY unused.
      flatGffData = "data.frame",
      flatGffGeneData = "list",
      analysisType = "character",
      DESeqDataSet = "DESeqDataSet",
      modelCoefficientsSample    = "list", #USUALLY unused.
      modelCoefficientsGene      = "list"  #USUALLY unused.
   ),
   prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), JunctionSeqCountSet = "0.0.5" ) ) )
)

################################
##Accessors and Replacers:
##  These were too much of a pain to use. Just use the slots. They really shouldn't be used
##  directly by end-users anyways.
################################
##
##
#setGeneric("jscs.designColumns",   function(object,...)       standardGeneric("jscs.designColumns"))
#setGeneric("jscs.designColumns<-", function(object,...,value) standardGeneric("jscs.designColumns<-"))
#setMethod("jscs.designColumns", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@designColumns;
#})
#setReplaceMethod("jscs.designColumns", signature(object="JunctionSeqCountSet", value="character"),
#  function( object, value ) {
#   object@designColumns <- value;
#})
################################
#setGeneric("jscs.dispFitCoefs",   function(object,...)       standardGeneric("jscs.dispFitCoefs"))
#setGeneric("jscs.dispFitCoefs<-", function(object,...,value) standardGeneric("jscs.dispFitCoefs<-"))
#setMethod("jscs.dispFitCoefs", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@dispFitCoefs;
#})
#setReplaceMethod("jscs.dispFitCoefs", signature(object="JunctionSeqCountSet", value="numeric"),
#  function( object, value ) {
#   object@dispFitCoefs <- value;
#})
################################
#setGeneric("jscs.fittedMu",   function(object,...)       standardGeneric("jscs.fittedMu"))
#setGeneric("jscs.fittedMu<-", function(object,...,value) standardGeneric("jscs.fittedMu<-"))
#setMethod("jscs.fittedMu", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@fittedMu;
#})
#setReplaceMethod("jscs.fittedMu", signature(object="JunctionSeqCountSet", value="matrix"),
#  function( object, value ) {
#   object@fittedMu <- value;
#})
################################
#setGeneric("jscs.dispFunctionType",   function(object,...)       standardGeneric("jscs.dispFunctionType"))
#setGeneric("jscs.dispFunctionType<-", function(object,...,value) standardGeneric("jscs.dispFunctionType<-"))
#setMethod("jscs.dispFunctionType", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@dispFunctionType;
#})
#setReplaceMethod("jscs.dispFunctionType", signature(object="JunctionSeqCountSet", value="list"),
#  function( object, value ) {
#   object@dispFunctionType <- value;
#})
################################
#setGeneric("jscs.dispFunction",   function(object,...)       standardGeneric("jscs.dispFunction"))
#setGeneric("jscs.dispFunction<-", function(object,...,value) standardGeneric("jscs.dispFunction<-"))
#setMethod("jscs.dispFunction", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@dispFunction;
#})
#setReplaceMethod("jscs.dispFunction", signature(object="JunctionSeqCountSet", value="function"),
#  function( object, value ) {
#   object@dispFunction <- value;
#})
################################
#setGeneric("jscs.dispFunctionJct",   function(object,...)       standardGeneric("jscs.dispFunctionJct"))
#setGeneric("jscs.dispFunctionJct<-", function(object,...,value) standardGeneric("jscs.dispFunctionJct<-"))
#setMethod("jscs.dispFunctionJct", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@dispFunctionJct;
#})
#setReplaceMethod("jscs.dispFunctionJct", signature(object="JunctionSeqCountSet", value="function"),
#  function( object, value ) {
#   object@dispFunctionJct <- value;
#})
################################
#setGeneric("jscs.dispFunctionExon",   function(object,...)       standardGeneric("jscs.dispFunctionExon"))
#setGeneric("jscs.dispFunctionExon<-", function(object,...,value) standardGeneric("jscs.dispFunctionExon<-"))
#setMethod("jscs.dispFunctionExon", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@dispFunctionExon;
#})
#setReplaceMethod("jscs.dispFunctionExon", signature(object="JunctionSeqCountSet", value="function"),
#  function( object, value ) {
#   object@dispFunctionExon <- value;
#})
################################
#setGeneric("jscs.formulas",   function(object,...)       standardGeneric("jscs.formulas"))
#setGeneric("jscs.formulas<-", function(object,...,value) standardGeneric("jscs.formulas<-"))
#setMethod("jscs.formulas", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@formulas;
#})
#setReplaceMethod("jscs.formulas", signature(object="JunctionSeqCountSet", value="list"),
#  function( object, value ) {
#   object@formulas <- value;
#})
################################
#setGeneric("jscs.annotationFile",   function(object,...)       standardGeneric("jscs.annotationFile"))
#setGeneric("jscs.annotationFile<-", function(object,...,value) standardGeneric("jscs.annotationFile<-"))
#setMethod("jscs.annotationFile", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@annotationFile;
#})
#setReplaceMethod("jscs.annotationFile", signature(object="JunctionSeqCountSet", value="character"),
#  function( object, value ) {
#   object@annotationFile <- value;
#})
################################
#setGeneric("jscs.geneCountData",   function(object,...)       standardGeneric("jscs.geneCountData"))
#setGeneric("jscs.geneCountData<-", function(object,...,value) standardGeneric("jscs.geneCountData<-"))
#setMethod("jscs.geneCountData", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@geneCountData;
#})
#setReplaceMethod("jscs.geneCountData", signature(object="JunctionSeqCountSet", value="matrix"),
#  function( object, value ) {
#   object@geneCountData <- value;
#})
################################
#setGeneric("jscs.countVectors",   function(object,...)       standardGeneric("jscs.countVectors"))
#setGeneric("jscs.countVectors<-", function(object,...,value) standardGeneric("jscs.countVectors<-"))
#setMethod("jscs.countVectors", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@countVectors;
#})
#setReplaceMethod("jscs.countVectors", signature(object="JunctionSeqCountSet", value="matrix"),
#  function( object, value ) {
#   object@countVectors <- value;
#})
################################
#setGeneric("jscs.altSizeFactors",   function(object,...)       standardGeneric("jscs.altSizeFactors"))
#setGeneric("jscs.altSizeFactors<-", function(object,...,value) standardGeneric("jscs.altSizeFactors<-"))
#setMethod("jscs.altSizeFactors", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@altSizeFactors;
#})
#setReplaceMethod("jscs.altSizeFactors", signature(object="JunctionSeqCountSet", value="data.frame"),
#  function( object, value ) {
#   object@altSizeFactors <- value;
#})
################################
#setGeneric("jscs.plottingEstimates",   function(object,...)       standardGeneric("jscs.plottingEstimates"))
#setGeneric("jscs.plottingEstimates<-", function(object,...,value) standardGeneric("jscs.plottingEstimates<-"))
#setMethod("jscs.plottingEstimates", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@plottingEstimates;
#})
#setReplaceMethod("jscs.plottingEstimates", signature(object="JunctionSeqCountSet", value="list"),
#  function( object, value ) {
#   object@plottingEstimates <- value;
#})
################################
#setGeneric("jscs.geneLevelPlottingEstimates",   function(object,...)       standardGeneric("jscs.geneLevelPlottingEstimates"))
#setGeneric("jscs.geneLevelPlottingEstimates<-", function(object,...,value) standardGeneric("jscs.geneLevelPlottingEstimates<-"))
#setMethod("jscs.geneLevelPlottingEstimates", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@geneLevelPlottingEstimates;
#})
#setReplaceMethod("jscs.geneLevelPlottingEstimates", signature(object="JunctionSeqCountSet", value="list"),
#  function( object, value ) {
#   object@geneLevelPlottingEstimates <- value;
#})
################################
#setGeneric("jscs.modelFitForHypothesisTest",   function(object,...)       standardGeneric("jscs.modelFitForHypothesisTest"))
#setGeneric("jscs.modelFitForHypothesisTest<-", function(object,...,value) standardGeneric("jscs.modelFitForHypothesisTest<-"))
#setMethod("jscs.modelFitForHypothesisTest", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@modelFitForHypothesisTest;
#})
#setReplaceMethod("jscs.modelFitForHypothesisTest", signature(object="JunctionSeqCountSet", value="list"),
#  function( object, value ) {
#   object@modelFitForHypothesisTest <- value;
#})
################################
#setGeneric("jscs.modelFitForEffectSize",   function(object,...)       standardGeneric("jscs.modelFitForEffectSize"))
#setGeneric("jscs.modelFitForEffectSize<-", function(object,...,value) standardGeneric("jscs.modelFitForEffectSize<-"))
#setMethod("jscs.modelFitForEffectSize", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@modelFitForEffectSize;
#})
#setReplaceMethod("jscs.modelFitForEffectSize", signature(object="JunctionSeqCountSet", value="list"),
#  function( object, value ) {
#   object@modelFitForEffectSize <- value;
#})
################################
#setGeneric("jscs.flatGffData",   function(object,...)       standardGeneric("jscs.flatGffData"))
#setGeneric("jscs.flatGffData<-", function(object,...,value) standardGeneric("jscs.flatGffData<-"))
#setMethod("jscs.flatGffData", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@flatGffData;
#})
#setReplaceMethod("jscs.flatGffData", signature(object="JunctionSeqCountSet", value="data.frame"),
#  function( object, value ) {
#   object@flatGffData <- value;
#})
################################
#setGeneric("jscs.flatGffGeneData",   function(object,...)       standardGeneric("jscs.flatGffGeneData"))
#setGeneric("jscs.flatGffGeneData<-", function(object,...,value) standardGeneric("jscs.flatGffGeneData<-"))
#setMethod("jscs.flatGffGeneData", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@flatGffGeneData;
#})
#setReplaceMethod("jscs.flatGffGeneData", signature(object="JunctionSeqCountSet", value="list"),
#  function( object, value ) {
#   object@flatGffGeneData <- value;
#})
################################
#setGeneric("jscs.analysisType",   function(object,...)       standardGeneric("jscs.analysisType"))
#setGeneric("jscs.analysisType<-", function(object,...,value) standardGeneric("jscs.analysisType<-"))
#setMethod("jscs.analysisType", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@analysisType;
#})
#setReplaceMethod("jscs.analysisType", signature(object="JunctionSeqCountSet", value="character"),
#  function( object, value ) {
#   object@analysisType <- value;
#})
################################
#setGeneric("jscs.DESeqDataSet",   function(object,...)       standardGeneric("jscs.DESeqDataSet"))
#setGeneric("jscs.DESeqDataSet<-", function(object,...,value) standardGeneric("jscs.DESeqDataSet<-"))
#setMethod("jscs.DESeqDataSet", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@DESeqDataSet;
#})
#setReplaceMethod("jscs.DESeqDataSet", signature(object="JunctionSeqCountSet", value="DESeqDataSet"),
#  function( object, value ) {
#   object@DESeqDataSet <- value;
#})
################################
#setGeneric("jscs.modelCoefficientsSample",   function(object,...)       standardGeneric("jscs.modelCoefficientsSample"))
#setGeneric("jscs.modelCoefficientsSample<-", function(object,...,value) standardGeneric("jscs.modelCoefficientsSample<-"))
#setMethod("jscs.modelCoefficientsSample", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@modelCoefficientsSample;
#})
#setReplaceMethod("jscs.modelCoefficientsSample", signature(object="JunctionSeqCountSet", value="list"),
#  function( object, value ) {
#   object@modelCoefficientsSample <- value;
#})
################################
#setGeneric("jscs.modelCoefficientsGene",   function(object,...)       standardGeneric("jscs.modelCoefficientsGene"))
#setGeneric("jscs.modelCoefficientsGene<-", function(object,...,value) standardGeneric("jscs.modelCoefficientsGene<-"))
#setMethod("jscs.modelCoefficientsGene", signature(object="JunctionSeqCountSet"),
#  function( object) {
#    object@modelCoefficientsGene;
#})
#setReplaceMethod("jscs.modelCoefficientsGene", signature(object="JunctionSeqCountSet", value="list"),
#  function( object, value ) {
#   object@modelCoefficientsGene <- value;
#})
################################

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
  rownames(colData) <- colnames(countData);
  
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

