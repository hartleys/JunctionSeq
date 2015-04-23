



##########################################################################
##########################################################################
##########################################################################
######### Front-end Interface functions:
##########################################################################
##########################################################################
##########################################################################





##########################################################################
######### Running JunctionSeq Analyses:
##########################################################################

runJunctionSeqAnalyses <- function(sample.files, sample.names, condition, 
                                flat.gff.file = NULL, 
                                outfile.prefix = NULL,  saveState = FALSE,
                                use.splice.sites = TRUE, use.novel.splice.sites = TRUE, use.exons = FALSE, 
                                nCores = 1,
                                use.covars = NULL, 
                                test.formula0 = formula(~ sample + countbin), 
                                test.formula1 = formula(~ sample + countbin + condition : countbin),
                                effect.formula = formula(~ condition + countbin + condition : countbin),
                                gzip.output = TRUE,
                                test.aggregated.genes = FALSE,
                                use.alternate.method = TRUE,
                                verbose = TRUE){
  keep.debug.model.data <- TRUE;
  gtf.format <- TRUE;
  if(verbose){
    message("> STARTING runJunctionSeqAnalyses (",date(),")");
    message(paste("> rda: sample.files: ", paste0(sample.files,collapse=", ")));
    message(paste("> rda: sample.names: ",  paste0(sample.names,collapse=", ")));
    message(paste("> rda: condition: ",    paste0(condition,collapse=", ")));
    message(paste("> rda: use.splice.sites: ",use.splice.sites));
    message(paste("> rda: use.novel.splice.sites: ",use.novel.splice.sites));
    message(paste("> rda: use.exons: ",use.exons));
    message(paste("> rda: nCores: ",nCores));
    message(paste("> rda: use.covars: ",use.covars));
    message(paste("> rda: test.formula0: ",test.formula0));
    message(paste("> rda: test.formula1: ",test.formula1));
    message(paste("> rda: gzip.output: ",gzip.output));
    message(paste("> rda: use.alternate.method: ",use.alternate.method));
    message(paste("> rda: test.aggregated.genes: ",test.aggregated.genes ));
  }
  #require("DEXSeq");
  
  if(! is.factor(condition)){
     condition <- factor(condition, levels = sort(unique(condition)));
  }
  
  design <- data.frame(condition = condition);
  if(! is.null(use.covars)){
    message(paste0("using covars:"," ",date()));
    if(class(use.covars) != "data.frame"){
       stop(paste0("FATAL ERROR: use.covars must be a data.frame! Instead it appears to be: ",class(use.covars)));
    }
    for(i in 1:ncol(use.covars)){
      if(! is.factor(use.covars[[i]]) ){
        use.covars[[i]] <- factor(use.covars[[i]], levels = sort(unique(use.covars[[i]]))  );
      }
    }
    
    design <- data.frame(cbind(design,use.covars));
    for(i in 1:length(names(use.covars))){
      message(paste0("      covar: ",names(use.covars)[i]));
      message(paste0(c("      ",paste0(use.covars[,i],collapse=", "))));
    }
    names(design) <- c("condition",names(use.covars));
  }
  row.names(design) <- sample.names;

    if(verbose) { message(paste0("> rda: Reading Count files..."," ",date(),".")) }
  jscs = readJunctionSeqCounts(countfiles = as.character(sample.files),
                             samplenames = sample.names,
                             design = design,
                             flat.gff.file = flat.gff.file, 
                             verbose = verbose,
                             use.splice.sites = use.splice.sites,
                             use.novel.splice.sites = use.novel.splice.sites,
                             use.exons = use.exons
                             );
  #gct <- htsc[[1]];
  #jscs <- htsc[[2]];
  if(verbose) {  message(paste0("> rda: Count files read."," ",date(),".")) }
    
    if(verbose) { message(paste0("> rda: Estimating Size Factors..."," ",date(),".")) }
  jscs <- estimateSizeFactors(jscs)
    if(verbose) { message(paste0("> rda: Size Factors Done. Size Factors are:",".")) 
                  message(paste0("> rda: ",paste0(names(sizeFactors(jscs)),collapse=",")))
                  message(paste0("> rda: ",paste0(sizeFactors(jscs),collapse=",")) )
                  message(paste0("> rda: Estimating Dispersions..."," ",date(),"."))
                }
  #jscs <- estimateDispersions(jscs, nCores = nCores, formula=test.formula1);

  #if(use.gene.counts.variant){
    jscs <- estimateJunctionSeqDispersions(jscs, nCores = nCores, test.formula1=test.formula1, use.alternate.method  =  use.alternate.method, test.aggregated.genes = test.aggregated.genes, verbose = verbose);
  #} else {
  #  jscs <- estimateDispersions(jscs, nCores = nCores, formula=test.formula1);
  #}

  #jscs <- estimateDispersions(jscs, nCores = nCores);
    if(verbose) { message(paste0("> rda: Dispersions estimated."," ",date(),".")) }
    if(verbose) { message(paste0("> rda: Fitting Dispersion Fcn..."," ",date(),".")) }
  jscs <- fitDispersionFunction(jscs)
  
    if(verbose) { message(paste0("> rda: Dispersions Fcn Fitted."," ",date(),".")) }
    if(verbose) { message(paste0("> rda: Testing for DEU..."," ",date(),".")) }
  
  #if(use.gene.counts.variant){
    jscs <- testJunctionsForDJU(jscs, nCores = nCores, test.formula0 = test.formula0, test.formula1 = test.formula1, use.alternate.method  =  use.alternate.method, verbose = verbose); 
  #} else {
  #  jscs <- testForDEU(jscs, nCores = nCores, formula0 = test.formula0, formula1 = test.formula1); 
  #}
  #jscs <- testForDEU(jscs, nCores = nCores, formula0 = test.formula0, formula1 = test.formula1); 
  #jscs <- testForDEU(jscs, nCores = nCores); 
  
  
    if(verbose) { message(paste0("> rda: DEU tests complete."," ",date(),".")) }
    if(verbose) { message(paste0("> rda: Estimating effect sizes using effects models..."," ",date(),".")) }
  jscs <- estimateEffectSizes( jscs , effect.formula = effect.formula)
    if(verbose) { message(paste0("> rda: Effect Sizes estimated.",".")) }
    if(verbose) { message(paste0("> rda: Generating results table..."," ",date(),".")) }
  
  ###Removed! This is now depreciated!
  ###  if(verbose) { message(paste0("> rda: DEU tests complete."," ",date(),".")) }
  ###  if(verbose) { message(paste0("> rda: Estimating log2FoldChanges..."," ",date(),".")) }
  ###jscs <- estimatelog2FoldChanges( jscs , estimate.fc.using.genewide.model = estimate.fc.using.genewide.model)
  ###  if(verbose) { message(paste0("> rda: FC estimated.",".")) }
  ###  if(verbose) { message(paste0("> rda: Generating results table..."," ",date(),".")) }
  ###res <- DEUresultTable( jscs )
  ###  if(verbose) { message(paste0("> rda: Results table generated."," ",date(),".")) }

  sampleNames(jscs) = as.character(sample.names)

  if(saveState & ! is.null(outfile.prefix)){
    if(verbose) { message(paste0("> rda: Saving state..."," ",date(),".")) }
    save(jscs,file=paste(outfile.prefix,"jscs.RData",sep=''));
    #save(res,file=paste(outfile.prefix,"res.RData",sep='')); 
    if(verbose) { message(paste0("> rda: State saved."," ",date(),".")) }
  } else {
    if(verbose) { message(paste0("> rda: Skipping save state because saveState = ", saveState, " and/or is.null(outfile.prefix) = ", is.null(outfile.prefix),". ",date(),".")) }
  }
  #if(! is.null(outfile.prefix)){
  #  res.sort <- res[order(res$padjust),]
  #  write.table.gz(res.sort,file=paste(outfile.prefix,".results.txt",sep=''), use.gzip = gzip.output, row.names = "feature.ID");
  #} else {
  #  if(verbose) { message(paste0("> rda: Skipping save results because is.null(outfile.prefix) = ", is.null(outfile.prefix),". ",date(),".")) }
  #}
  #if(verbose) { message(paste0("> rda: Done. ",date(),".")) }
  #return(list(jscs,res));
  
  if(verbose) message("> FINISHED runJunctionSeqAnalyses (",date(),")");
  return(jscs);
}

##########################################################################
######### 
##########################################################################


estimateEffectSizes <- function(jscs, 
                        effect.formula = formula(~ condition + countbin + condition : countbin),
                        nCores=1, 
                        dispColumn="dispersion",
                        verbose = TRUE){
  runOnFeatures = seq_len(nrow(fData(jscs)));
  fitExpToVar <- "condition";
  varlist <- rownames(attr(terms(effect.formula), "factors"));
  covarlist <- varlist[ ! varlist %in% c(fitExpToVar, "countbin") ]
  
  myApply <- getMyApply(nCores);
  
  modelFrame <- constructModelFrame( jscs )
  modelFrame$countbin <- factor(modelFrame$countbin, levels = c("others","this"));
  mm <- rmDepCols( model.matrix( effect.formula, modelFrame ) );
  conditionLevels <- levels(modelFrame[[fitExpToVar]]);
  conditionCt <- length(conditionLevels);
                  
  mdl.out <- myApply(runOnFeatures, function(i){
    if( verbose & i %% 1000 == 0 ){
       message(paste0("-------> estimateEffectSizes: (Calculating effect size and predicted values for feature ",i," of ",nrow(jscs),")","(",date(),")"));
    }
    geneID <- fData(jscs)$geneID[i];
    countbinID <- fData(jscs)$countbinID[i];
    countVector <- jscs@countVectors[i,];
    
    if(! fData(jscs)$testable[i]){
      return(list(
        logFCs = rep(NA, conditionCt - 1),
        fitEffect = NA,
        exprEstimate = rep(NA,conditionCt),
        exprEstimateVST = rep(NA,conditionCt),
        otherExprEstimate = rep(NA,conditionCt),
        otherExprEstimateVST = rep(NA,conditionCt),
        relExprEstimate = rep(NA,conditionCt),
        relExprEstimateVST = rep(NA,conditionCt),
        normCounts = countVector * modelFrame$sizeFactor,
        normCountsVST = vst(countVector * modelFrame$sizeFactor, jscs)
      ));
    } else {
      disp <- fData(jscs)[i, dispColumn];
      fit <- try( {
        glmnb.fit( mm,  countVector, dispersion = disp, offset = log( modelFrame$sizeFactor ) );
      });
      if( any(inherits( fit, "try-error" ) )) {
        warning( sprintf("glmnb.fit failed for %s:%s\n", as.character( geneID ), countbinID) );
        return(list(fitConverged = FALSE, fit = NULL));

        return(list(
          logFCs = rep(NA, conditionCt - 1),
          fitEffect = NA,
          exprEstimate = rep(NA,conditionCt),
          exprEstimateVST = rep(NA,conditionCt),
          otherExprEstimate = rep(NA,conditionCt),
          otherExprEstimateVST = rep(NA,conditionCt),
          relExprEstimate = rep(NA,conditionCt),
          relExprEstimateVST = rep(NA,conditionCt),
          normCounts = countVector * modelFrame$sizeFactor,
          normCountsVST = vst(countVector * modelFrame$sizeFactor, jscs)
        ));
      }
      
      coefs <- arrangeCoefs( effect.formula, modelFrame, mm, fit = fit, insertValues = TRUE );
      
      predictedEstimates <- getPredictedEstimates(coefs = coefs, 
                                                  forVarName = fitExpToVar, 
                                                  forVarValues = conditionLevels, 
                                                  selectVarName = "countbin", 
                                                  selectVarValue = "this", 
                                                  averageVarNames = covarlist);
      
      predictedEstimatesGene <- getPredictedEstimates(coefs = coefs, 
                                                  forVarName = fitExpToVar, 
                                                  forVarValues = conditionLevels, 
                                                  selectVarName = "countbin", 
                                                  selectVarValue = "others", 
                                                  averageVarNames = covarlist);
      
      #meanCondition <- mean(sapply(predictedEstimates, function(pe){
      #    pe[[fitExpToVar]];
      #}));
      
      countBinThisBeta <- predictedEstimates[[1]][["countbin"]];
      countBinOthersBeta <- predictedEstimatesGene[[1]][["countbin"]];
      interceptBeta <- predictedEstimates[[1]][["(Intercept).(Intercept)"]];
      
      relativeEstimates <- predictedEstimates;
      for(j in seq_len(length(predictedEstimates))){
        relativeEstimates[[j]] <- predictedEstimates[[j]] - predictedEstimatesGene[[j]];
      }
      relativeLogExprDiff <- sapply(relativeEstimates, sum);
      relativeLogExprEstimate <- sapply(relativeEstimates, sum) - countBinThisBeta + countBinOthersBeta;
      #logFCs <- relativeLogExprEstimate[-1] - relativeLogExprEstimate[1];
      logFCs <- relativeLogExprDiff[-1] - relativeLogExprDiff[1];
      
      exprEstimate <- exp(sapply(predictedEstimates, sum));
      otherExprEstimate <- exp(sapply(predictedEstimatesGene, sum));
      #relExprEstimate <- exp(relativeLogExprEstimate + interceptBeta + countBinThisBeta);
      relExprEstimate <- exp((relativeLogExprDiff - mean(relativeLogExprDiff)) + mean(log(exprEstimate)));
      
      return(list(
        logFCs = logFCs,
        fitEffect = fit,
        exprEstimate = exprEstimate,
        exprEstimateVST = vst(exprEstimate, jscs),
        otherExprEstimate = otherExprEstimate,
        otherExprEstimateVST = vst(otherExprEstimate, jscs),
        relExprEstimate = relExprEstimate,
        relExprEstimateVST = vst(relExprEstimate, jscs),
        normCounts = countVector * modelFrame$sizeFactor,
        normCountsVST = vst(countVector * modelFrame$sizeFactor, jscs)
      ));
    }
  });
  
  logFCs <- do.call( rbind.data.frame , lapply(mdl.out, "[[", "logFCs") );
  #Change the base to base-2:
  logFCs <- apply(logFCs, c(1,2), function(x){x / log(2)});
  
  
  colnames(logFCs) <- paste0("log2FC(",conditionLevels[-1],"/",conditionLevels[1],")");
  for(colname in colnames(logFCs)){
    fData(jscs)[[colname]] <- logFCs[,colname];
  }
  
  modelFitForEffectSize <- lapply(mdl.out, "[[", "fitEffect");
  names(modelFitForEffectSize) <- featureNames( jscs );
  jscs@modelFitForEffectSize <- modelFitForEffectSize;
  
  exprEstimate         <- extractPlottingEstimates_helper(jscs, mdl.out, "exprEstimate",         paste0("expr","_",conditionLevels));
  exprEstimateVST      <- extractPlottingEstimates_helper(jscs, mdl.out, "exprEstimateVST",      paste0("exprVST","_",conditionLevels));
  otherExprEstimate    <- extractPlottingEstimates_helper(jscs, mdl.out, "otherExprEstimate",    paste0("geneExpr","_",conditionLevels));
  otherExprEstimateVST <- extractPlottingEstimates_helper(jscs, mdl.out, "otherExprEstimateVST", paste0("geneExprVST","_",conditionLevels));
  
  relExprEstimate      <- extractPlottingEstimates_helper(jscs, mdl.out, "relExprEstimate",      paste0("relExpr","_",conditionLevels));
  relExprEstimateVST   <- extractPlottingEstimates_helper(jscs, mdl.out, "relExprEstimateVST",   paste0("relExprVST","_",conditionLevels));
  
  sampleNames <- colnames(counts(jscs));
  normCounts           <- extractPlottingEstimates_helper(jscs, mdl.out, "normCounts",   c(paste0("normCount_",sampleNames), paste0("normGeneCount_",sampleNames)));
  normCountsVST        <- extractPlottingEstimates_helper(jscs, mdl.out, "normCountsVST",   c(paste0("normCountsVST_",sampleNames), paste0("normGeneCountsVST_",sampleNames)));
  
  jscs@plottingEstimates    <- list(exprEstimate = exprEstimate,
                                     geneExprEstimate = otherExprEstimate,
                                     relExprEstimate = relExprEstimate, 
                                     normCounts = normCounts);
  jscs@plottingEstimatesVST <- list(exprEstimateVST = exprEstimateVST, 
                                     geneExprEstimateVST = otherExprEstimateVST, 
                                     relExprEstimateVST = relExprEstimateVST, 
                                     normCountsVST = normCountsVST);
  return(jscs);
}

extractPlottingEstimates_helper <- function(jscs, mdl.out, estName, columnNames){
  x <- do.call( rbind.data.frame , lapply(mdl.out, "[[", estName) );
  rownames(x) <- featureNames( jscs );
  colnames(x) <- columnNames;
  return(x);
}

getPredictedEstimates <- function(coefs, forVarName, forVarValues, selectVarName, selectVarValue, averageVarNames){
  predictedEstimates <- lapply(forVarValues, function(pvn){
    pvToAdd <- sapply(coefs, function(coef){
      if("(Intercept)" %in% names(dimnames(coef))){
        return(coef[1]);
      }
      coef <- applyByDimname(coef, dimname = selectVarName, FUN = function(cf){
        cf[[selectVarValue]];
      });
      coef <- applyByDimname(coef, dimname = forVarName, FUN = function(cf){
        cf[[pvn]];
      });
      for(covar in averageVarNames){
        coef <- applyByDimname(coef, dimname = covar, FUN = function(cf){
          mean(cf);
        });
      }
      #message("class(coef) = ", class(coef));
      #message("length(coef) = ", length(coef));
      as.vector(coef);
    });
    pvToAdd;
  });
  names(predictedEstimates) <- forVarValues;
  return(predictedEstimates);
}



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

writeCompleteResults <- function(jscs, outfile.prefix, 
                            gzip.output = TRUE, FDR.threshold = 0.05,
                            save.allGenes = TRUE, save.sigGenes = TRUE, save.fit = TRUE, save.VST = TRUE,
                            save.bedTracks = TRUE,
                            save.jscs = FALSE,
                            verbose = TRUE
                            ){
   if(verbose){
       message("> STARTING writeCompleteResults (",date(),")");
       message(paste("> wcr: outfile.prefix:",outfile.prefix));
       message(paste("> wcr: FDR.threshold:",FDR.threshold));
       message(paste("> wcr: save.allGenes:",save.allGenes));
       message(paste("> wcr: save.sigGenes:",save.sigGenes));
       message(paste("> wcr: save.fit:",save.fit));
       message(paste("> wcr: save.VST:",save.VST));
   }
   keep.debug.model.data <- TRUE;
   
   countVectors <- jscs@countVectors;
   sampleNames <- sampleNames(jscs@phenoData);
   colnames(countVectors) <- c(paste0("rawCounts_",sampleNames), paste0("rawGeneCounts_",sampleNames));
   
   expression.data <- cbind.data.frame(featureID = rownames(fData(jscs)), countVectors, do.call(cbind.data.frame,jscs@plottingEstimates)  );
   colnames(expression.data) <- c("featureID", colnames(countVectors), unlist( lapply( jscs@plottingEstimates, colnames ) ) )
   expression.data.vst <- cbind.data.frame(featureID = rownames(fData(jscs)), do.call(cbind.data.frame,jscs@plottingEstimatesVST)  );
   colnames(expression.data.vst) <- c("featureID", unlist( lapply( jscs@plottingEstimatesVST, colnames ) ) )
   
   if(save.allGenes & verbose) message("> wcr: Writing results for ",length(unique(as.character(fData(jscs)$geneID)))," genes.");
   if(save.allGenes & verbose) message("> wcr:     Found ",nrow(fData(jscs))," counting bins belonging to these genes.");
   if(save.allGenes) write.simple.table.gz(expression.data,     file=paste0(outfile.prefix,"allGenes.expression.data.txt"),      use.gzip=gzip.output,row.names=F,col.names=T,quote=F, sep = '\t');
   if(save.allGenes & save.VST) write.simple.table.gz(expression.data.vst, file=paste0(outfile.prefix,"allGenes.expression.data.VST.txt"),  use.gzip=gzip.output,row.names=F,col.names=T,quote=F, sep = '\t');
   
   if(save.allGenes) write.table.gz(fData(jscs),         file=paste0(outfile.prefix,"allGenes.results.txt"), use.gzip=gzip.output,row.names="featureID",col.names=T,quote=F);
   
   modelFitForHypothesisTest <- jscs@modelFitForHypothesisTest;
   modelFitForEffectSize <- jscs@modelFitForEffectSize;
   if(save.fit) save( modelFitForHypothesisTest, file = paste0(outfile.prefix, "modelFitForHypothesisTest.RData") );
   if(save.fit) save( modelFitForEffectSize,     file = paste0(outfile.prefix, "modelFitForEffectSize.RData") );
   
   if(save.sigGenes){
     sig.features <- which(fData(jscs)$padjust < FDR.threshold);
     gene.list <- unique(as.character(fData(jscs)$geneID[sig.features]));
     if(verbose) message("> wcr: Writing results for ", length(gene.list), " genes with 1 or more significant junctions (at adjusted-p-value threshold ", FDR.threshold,")");
     sig.rows <- which( as.character(fData(jscs)$geneID) %in% gene.list );
     if(verbose) message("> wcr:     Found ", length(sig.rows), " counting bins belonging to those genes.");
     
     
                  write.simple.table.gz(expression.data[sig.rows,, drop=FALSE],     file=paste0(outfile.prefix,"sigGenes.expression.data.txt"),      use.gzip=gzip.output,row.names=F,col.names=T,quote=F, sep = '\t');
     if(save.VST) write.simple.table.gz(expression.data.vst[sig.rows,, drop=FALSE], file=paste0(outfile.prefix,"sigGenes.expression.data.VST.txt"),  use.gzip=gzip.output,row.names=F,col.names=T,quote=F, sep = '\t');
                  write.table.gz(fData(jscs)[sig.rows,, drop=FALSE],         file=paste0(outfile.prefix,"sigGenes.results.txt"), use.gzip=gzip.output,row.names="featureID",col.names=T,quote=F);
   }
   
   if(save.bedTracks){
     if(save.allGenes){
       writeExprBedTrack(paste0(outfile.prefix,"allGenes.coverage.bed.gz"), 
                         jscs = jscs, 
                         only.with.sig.gene = FALSE,
                         trackLine = paste0("track name='JctExprAll' description='Splice Junction Coverage Estimates, by group' itemRgb='On' visibility=3"));
     }
     if(save.sigGenes){
       writeExprBedTrack(paste0(outfile.prefix,"sigGenes.coverage.bed.gz"), 
                         jscs = jscs, 
                         trackLine = paste0("track name='JctExprSig' description='Splice Junction Coverage Estimates, by group' itemRgb='On' visibility=3"));
       writeSigBedTrack(paste0(outfile.prefix,"sigGenes.pvalues.bed.gz"), 
                        jscs = jscs,
                        trackLine = paste0("track name='JctPvals' description='Significant Splice Junctions' useScore=1 visibility=3")
                        );
     }
   }
   
   if(verbose) message("> DONE writeCompleteResults (",date(),")");
   
   if(save.jscs) save( jscs, file = paste0(outfile.prefix, "jscs.RData") );
}

generateCompleteResults <- function(jscs, outfile.prefix = NULL, verbose = TRUE, nCores=1, gzip.output = TRUE
						#use.splice.sites = TRUE, use.novel.splice.sites = TRUE, use.exons = FALSE,  #,keep.debug.model.data = FALSE
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
   #fd$meanBase <- rowMeans(counts(jscs, normalized=TRUE))
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

##########################################################################
######### Generating Plots From Results:
##########################################################################

buildAllPlots <- function(jscs,
                          outfile.prefix = "./",
                          flat.gff.data = NULL, flat.gff.file = NULL, 
                          gene.list = NULL, FDR.threshold = 0.01, 
                          use.plotting.device = "png", 
                          use.vst=FALSE,use.log = TRUE, truncateBelowOne = TRUE, 
                          exon.rescale.factor = 0.3,
                          ma.plot=FALSE, variance.plot=FALSE,
                          with.TX=TRUE,without.TX=TRUE,
                          expr.plot=TRUE,normCounts.plot=TRUE,
                          rExpr.plot=TRUE,rawCounts.plot=TRUE,
                          colorRed.FDR.threshold = FDR.threshold, 
                          color=NULL,
                          plot.exon.results = FALSE, 
                          plot.splice.results = TRUE, 
                          plot.novel.splice.results = plot.splice.results,
                          plot.untestable.results = FALSE,
                          plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                          par.cex = 1, anno.cex.text = 2, anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                          drawCoordinates = TRUE,
                          arrows.length = 0.25,
                          base.plot.height = 22.222, base.plot.width = 22.222, 
                          base.plot.units = "in", 
                          plotting.device.params = list(pointsize = 18, res = 72), 
                          number.plots = TRUE,
                          openPlottingDeviceFunc = NULL, closePlottingDeviceFunc = NULL,
                          verbose=TRUE,
                          ...){
  condition <- jscs@phenoData$condition;
  
  if(is.null(openPlottingDeviceFunc) & is.null(closePlottingDeviceFunc)){
    devFunctions <- getPlottingDeviceFunc(use.plotting.device = use.plotting.device,
                                            base.plot.height = base.plot.height,
                                            base.plot.width = base.plot.width,
                                            base.plot.units = base.plot.units,
                                            plotting.device.params = plotting.device.params);
    openPlottingDeviceFunc <- devFunctions[[1]];
    closePlottingDeviceFunc <- devFunctions[[2]];
  }
  
  
  gtf.format <- TRUE;
  
  if(is.null(gene.list)){
    sig.features <- which(fData(jscs)$padjust < FDR.threshold);
    gene.list <- unique(as.character(fData(jscs)$geneID[sig.features]));
    if(verbose) message("> buildAllPlots: Found ", length(gene.list), " significant genes to plot, at adjusted-p-value threshold ", FDR.threshold);
  } else {
    if(verbose) message("> buildAllPlots: Found ", length(gene.list), " genes to plot.");
  }
  
  if(length(gene.list) > 0){
    if(is.null(flat.gff.data)){
      if(is.null(flat.gff.file)){
        stop("FATAL ERROR: buildAllPlots: either flat.gff.file or flat.gff.data must be specified!");
      }
      flat.gff.data <- readAnnotationData(flat.gff.file);
    }

    #numcond <- length(levels(condition));

    if(verbose) message("> buildAllPlots: Starting plotting...");

    if(variance.plot){
      openPlottingDeviceFunc(paste(outfile.prefix,"dispersion-plot","",sep=""),heightMult=1,widthMult=1);
      plotDispEsts( jscs, cex = anno.cex.text * 0.45 );
      dev.off();
    }

    if(ma.plot){
      #which.is.fc <- which(substr(names(merged.results.data),1,8) == "log2fold")
      ##print("which is fc:");
      ##print(substr(names(merged.results.data),1,8));
      #ma.frame <- data.frame(baseMean=merged.results.data$meanBase, log2FoldChange = merged.results.data[,which.is.fc], padj = merged.results.data$padjust);

      fc.cols <- which(strStartsWith(names(fData(jscs)), "log2FC("));

      for(fc.name in names(fData(jscs))[fc.cols]){
        if(verbose) message("> buildAllPlots: Generating MA-Plot (",fc.name,")");

        fc.title <- sub("/","vs",fc.name, fixed = TRUE);
        openPlottingDeviceFunc(paste(outfile.prefix,"ma-plot-",fc.title,sep=""),heightMult=1,widthMult=1);
        plotMA( jscs, FDR.threshold=colorRed.FDR.threshold, fc.name = fc.name, par.cex = 2 * par.cex, text.cex = anno.cex.text / 2, points.cex = anno.cex.text / 2, ... );
        dev.off();
      }
    }

    FDR <- colorRed.FDR.threshold

    total.gene.ct <- length(gene.list);
    number.width <- floor(log10(total.gene.ct)) + 1;
    if( number.plots ){
      geneNum.strings <- paste0(formatC(1:total.gene.ct, width = number.width, format = 'd', flag = '0'), "-");
    } else {
      geneNum.strings <- rep("",total.gene.ct);
    }

    #if(verbose){
    #   message(paste0("total gene ct: ",total.gene.ct,", column width: ",number.width));
    #}

    geneNum <- 1;
    for(geneID in gene.list){
      message(paste("> buildAllPlots: starting geneID:",geneID));
      geneNum.string <- geneNum.strings[geneNum];

      buildAllPlotsForGene(geneID = geneID, jscs = jscs,
                            outfile.prefix = paste0(outfile.prefix,geneNum.string),
                            flat.gff.data = flat.gff.data, 
                            use.plotting.device = use.plotting.device,
                            use.vst=use.vst, use.log = use.log, truncateBelowOne = truncateBelowOne,
                            exon.rescale.factor = exon.rescale.factor,
                            with.TX=with.TX,without.TX=without.TX,
                            expr.plot=expr.plot,normCounts.plot=normCounts.plot,
                            rExpr.plot=rExpr.plot,rawCounts.plot=rawCounts.plot,
                            colorRed.FDR.threshold = colorRed.FDR.threshold, 
                            color= color,
                            plot.exon.results = plot.exon.results, 
                            plot.splice.results = plot.splice.results,
                            plot.novel.splice.results = plot.novel.splice.results,
                            plot.untestable.results = plot.untestable.results,
                            plot.lwd = plot.lwd, axes.lwd = axes.lwd, anno.lwd = anno.lwd, 
                            drawCoordinates = drawCoordinates,arrows.length = arrows.length,
                            par.cex = par.cex,
                            anno.cex.text = anno.cex.text,
                            anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,
                            base.plot.height = base.plot.height, base.plot.width = base.plot.width,
                            base.plot.units = base.plot.units,
                            plotting.device.params = plotting.device.params,
                            verbose=verbose,
                            ...);

      geneNum <- geneNum + 1;
    }
    if(verbose) message("> buildAllPlots: Plotting complete.");
  }
  if(verbose) message("> buildAllPlots: Plotting and data writing complete.");
}

buildAllPlotsForGene <- function(geneID,jscs,
                          outfile.prefix = "./", 
                          flat.gff.data = NULL, flat.gff.file = NULL,
                          use.plotting.device = "png", 
                          use.vst=FALSE, use.log = TRUE, truncateBelowOne = TRUE,
                          exon.rescale.factor = 0.3,
                          with.TX=TRUE,without.TX=TRUE,
                          expr.plot=TRUE,normCounts.plot=TRUE,
                          rExpr.plot=TRUE,rawCounts.plot=TRUE,
                          colorRed.FDR.threshold = 0.01, 
                          color=NULL,
                          plot.exon.results = FALSE, 
                          plot.splice.results = TRUE, 
                          plot.novel.splice.results = plot.splice.results,
                          plot.untestable.results = FALSE,
                          plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                          par.cex = 1, anno.cex.text = 2,
                          anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                          drawCoordinates = TRUE,
                          arrows.length = 0.25,
                          base.plot.height = 22.222, base.plot.width = 22.222, 
                          base.plot.units = "in", 
                          plotting.device.params = list(pointsize = 18, res = 72),
                          openPlottingDeviceFunc = NULL, closePlottingDeviceFunc = NULL,
                          verbose=TRUE,  ...){
    message(paste("starting buildAllPlotsForGene() for geneID:",geneID));
    condition <- jscs@phenoData$condition
    
    #if("HIDDENPARAM_FIX_SVG_RHEL5" %in% names(plotting.device.params)){
    #   closeDeviceFunc <- function(filename){
    #     dev.off();
    #     #Hack (doesn't work), attempting to workaround a flaw in one particular linux distribution where svg files are rendered incorrectly.
    #     #svgfixfile <- system.file("extdata/FIX_SVG.pl", package="JunctionSeq", mustWork=TRUE);
    #     #svgfixcommand <- paste("perl ",svgfixfile," ", filename, " > ", filename);
    #     #system(svgfixcommand);
    #   };
    #   plotting.device.params <- plotting.device.params[- which(names(plotting.device.params) == "HIDDENPARAM_FIX_SVG_RHEL5")];
    #} else {
    #   closeDeviceFunc <- function(filename){
    #     dev.off();
    #   }
    #}
    #if(is.null(openPlottingDeviceFunc)){
    #    openPlottingDeviceFunc <- getPlottingDeviceFunc(use.plotting.device = use.plotting.device,
    #                                          base.plot.height = base.plot.height,
    #                                          base.plot.width = base.plot.width,
    #                                          base.plot.units = base.plot.units,
    #                                          plotting.device.params = plotting.device.params);
    #}
    
    if(is.null(openPlottingDeviceFunc) & is.null(closePlottingDeviceFunc)){
      devFunctions <- getPlottingDeviceFunc(use.plotting.device = use.plotting.device,
                                              base.plot.height = base.plot.height,
                                              base.plot.width = base.plot.width,
                                              base.plot.units = base.plot.units,
                                              plotting.device.params = plotting.device.params);
      openPlottingDeviceFunc <- devFunctions[[1]];
      closePlottingDeviceFunc <- devFunctions[[2]];
    }

    
    FDR <- colorRed.FDR.threshold;
    if(is.null(flat.gff.data)){
       if(is.null(flat.gff.file)){
          stop("FATAL ERROR: buildAllPlotsForGene: either flat.gff.file or flat.gff.data must be specified!");
       }
       flat.gff.data <- readAnnotationData(flat.gff.file);
    }
    if(expr.plot){
      plot.type <- "expr"
      if(with.TX){
        outfile <- paste(outfile.prefix,geneID,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1.5,widthMult=1);
        plotJunctionSeqResultsForGene(geneID, jscs, flat.gff.data, colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd , par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.splice.results = plot.splice.results , plot.novel.splice.results = plot.novel.splice.results , plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix,geneID,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=1);
        plotJunctionSeqResultsForGene(geneID, jscs, flat.gff.data, colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.splice.results = plot.splice.results , plot.novel.splice.results = plot.novel.splice.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length,...)
        closePlottingDeviceFunc();  
      }
    }
    if(normCounts.plot){
      plot.type <- "normCounts"
      if(with.TX){
        outfile <- paste(outfile.prefix,geneID,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1.5,widthMult=1);
        plotJunctionSeqResultsForGene(geneID, jscs, flat.gff.data, colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.splice.results = plot.splice.results , plot.novel.splice.results = plot.novel.splice.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix,geneID,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=1);
        plotJunctionSeqResultsForGene(geneID, jscs, flat.gff.data, colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.splice.results = plot.splice.results , plot.novel.splice.results = plot.novel.splice.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length,...)
        closePlottingDeviceFunc();  
      }
    }
    if(rExpr.plot){
      plot.type <- "rExpr"
      if(with.TX){
        outfile <- paste(outfile.prefix,geneID,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1.5,widthMult=1);
        plotJunctionSeqResultsForGene(geneID, jscs, flat.gff.data, colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.splice.results = plot.splice.results , plot.novel.splice.results = plot.novel.splice.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix,geneID,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=1);
        plotJunctionSeqResultsForGene(geneID, jscs, flat.gff.data, colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.splice.results = plot.splice.results , plot.novel.splice.results = plot.novel.splice.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length,...)
        closePlottingDeviceFunc();  
      }
    }
    if(rawCounts.plot){
      plot.type <- "rawCounts"
      if(with.TX){
        outfile <- paste(outfile.prefix,geneID,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1.5,widthMult=1);
        plotJunctionSeqResultsForGene(geneID, jscs, flat.gff.data, colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=FALSE,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.splice.results = plot.splice.results , plot.novel.splice.results = plot.novel.splice.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix,geneID,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=1);
        plotJunctionSeqResultsForGene(geneID, jscs, flat.gff.data, colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=FALSE,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.splice.results = plot.splice.results , plot.novel.splice.results = plot.novel.splice.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length,...)
        closePlottingDeviceFunc();  
      }
    }
}

#read.gtf.annotation <- function(gtf.file){
#   gff.data <- read.table(gtf.file,header=FALSE,sep='	');
#   names(gff.data) <- c("chrom","src","feature","start","end","score","strand","frame","attrString");
#}


#readAnnotationData <- function(tabbed.anno.file){
#  return(read.table( anno.file, header=TRUE, stringsAsFactors=FALSE));
  #note file columns are:
  #featureName	featureType	chrom	start	end	strand	gene_id	part_number	transcripts
#}

#read.HTSeqCounts.old <- function( countfiles, design, flattenedfile=NULL )
#{
#   lf <- lapply( countfiles, function(x)
#      read.table( x, header=FALSE,stringsAsFactors=FALSE ) )
#   if( !all( sapply( lf[-1], function(x) all( x$V1 == lf[1]$V1 ) ) ) )
#      stop( "Count files have differing gene ID column." )
#   dcounts <- sapply( lf, `[[`, "V2" )
#   rownames(dcounts) <- lf[[1]][,1]
#   rownames(dcounts)
#   dcounts <- dcounts[ substr(rownames(dcounts),1,1)!="_", ]
#   rownames(dcounts) <- sub(":", ":E", rownames(dcounts))
#   colnames(dcounts) <- countfiles
#   splitted <- strsplit(rownames(dcounts), ":")
#   exons <- sapply(splitted, "[[", 2)
#   genesrle <- sapply( splitted, "[[", 1)
#   if(!is.null(flattenedfile)){
#      aggregates<-read.delim(flattenedfile, stringsAsFactors=FALSE, header=FALSE)
#      colnames(aggregates)<-c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
#      aggregates<-aggregates[which(aggregates$class =="exonic_part"),]
#      aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
#      aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
#      transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
#      transcripts <- gsub("\\+", ";", transcripts)
#      featureIDs <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", aggregates$attr)
#      exoninfo<-data.frame(chr=aggregates$chr, start=aggregates$start, end=aggregates$end, strand=aggregates$strand)
#      rownames( exoninfo ) <- paste( aggregates$gene_id, featureIDs, sep=":E" )
#      names(transcripts) <- rownames(exoninfo)
#      if (!all( rownames(dcounts) %in% rownames(exoninfo) )){
#         stop("Count files do not correspond to the flattened annotation file")
#      }
#      matching <- match(rownames(dcounts), rownames(exoninfo))
#      jscs<-newJunctionSeqCountSet(countData=dcounts, design=design, geneIDs=genesrle, featureIDs=exons, exonIntervals=exoninfo[matching,], transcripts=transcripts[matching])
#      jscs@annotationFile <- flattenedfile
#      pData(jscs)$countfiles <- countfiles
#      return(jscs)
#   }else{
#      return(newJunctionSeqCountSet(countData=dcounts, design=design, geneIDs=genesrle, featureIDs=exons))
#   }
#}


##########################################################################
##########################################################################
##########################################################################
######### Top-Level Coefficient calculation functions
##########################################################################
##########################################################################
##########################################################################


generateAllExpressionEstimates.v2 <- function(jscs,nCores = 1,fitExpToVar="condition",verbose=TRUE){

   myApply <- getMyApply(nCores);

   levelCt <- length(levels(jscs@phenoData$condition))
   sampleCt <- length(sampleNames(jscs@phenoData))
   
   out <- getAllData2(jscs, runOnFeatures = seq_len(nrow(fData(jscs))),
                            fitExpToVar=fitExpToVar, 
                            formula1 = formula(paste0("~ ",fitExpToVar," + countbin + ",fitExpToVar," : countbin")),
                            nCores=nCores, dispColumn="dispersion",verbose = verbose);
   
   out <- cbind.data.frame(featureID = rownames(fData(jscs)),
                           geneID = fData(jscs)$geneID,
                           countbinID = fData(jscs)$countbinID,
                           out);
   
   out;
}



generateAllExpressionEstimates <- function(jscs,nCores = 1,fitExpToVar="condition",verbose=TRUE){
   #require(statmod)
   #require(DEXSeq)
   
   #require(multicore)
   
   myApply <- getMyApply(nCores);

   levelCt <- length(levels(jscs@phenoData$condition))
   sampleCt <- length(sampleNames(jscs@phenoData))
   
   out <- rbind( c(0.0,0.0,0.0,rep(0.0,levelCt),rep(0.0,levelCt),rep(0.0,sampleCt),rep(0.0,levelCt),rep(0.0,levelCt),rep(0.0,sampleCt)))
   #out <- rbind( c(0,0,0,rep(0,levelCt),rep(0,levelCt),rep(0,sampleCt),rep(0,levelCt),rep(0,levelCt),rep(0,sampleCt)))

   out <- data.frame(out)[0,]
      
   #names.out <- c("featureID","geneID","countbinID",
   #                paste("exprVST",sort(levels(jscs@phenoData$condition)),sep='_'),paste("expr",sort(levels(jscs@phenoData$condition)),sep='_'),
   #                paste("rExprVST",sort(levels(jscs@phenoData$condition)),sep='_'),paste("rExpr",sort(levels(jscs@phenoData$condition)),sep='_'),
   #                paste("normCountVST",sampleNames(jscs@phenoData),sep='_'),paste("normCount",sampleNames(jscs@phenoData),sep='_'),paste("rawCount",sampleNames(jscs@phenoData),sep='_')
   #                );
   names.out <- c("featureID","geneID","countbinID",
                   paste("exprVST",(levels(jscs@phenoData$condition)),sep='_'),
                   paste("expr",(levels(jscs@phenoData$condition)),sep='_'),
                   paste("rExprVST",(levels(jscs@phenoData$condition)),sep='_'),
                   paste("rExpr",(levels(jscs@phenoData$condition)),sep='_'),
                   paste("normCountVST",sampleNames(jscs@phenoData),sep='_'),
                   paste("normCount",sampleNames(jscs@phenoData),sep='_'),
                   paste("rawCount",sampleNames(jscs@phenoData),sep='_')
                   );
   #message(dim(out))
   gene.list <- unique(jscs@featureData$geneID);
   blank.line.length <- length(levels(jscs@phenoData$condition)) * 4;
   #rep(NA,length(levels(jscs@phenoData$condition)) * 4);
   
   if(verbose){
      message(paste("Starting expression estimate generation on",length(gene.list),"top-level features (eg genes)."));
   }
   get.one.genes.data <- function(geneIndex){
      geneID <- gene.list[geneIndex];
      es <- fitAndArrangeCoefs( jscs, geneID, frm=as.formula(paste("count ~", fitExpToVar,  "* countbin")) )
      if(! is.null(es)){
          #If the model converged properly, then output the results:
          curr <- data.frame(getAllData(jscs,es,geneID))
          if(geneIndex - 1 %% 10000 == 0 & verbose){
            message("colnames:");
            message(paste0(colnames(curr),collapse=" "));
          }
          names(curr) <- names.out;
          return(  list(TRUE , curr)  );
      } else {
          #If the model does not converge, return NA for all fitted values.
          if(verbose){
            message(paste("glm fit failed for gene", geneID, " index: ",geneIndex));
          }
          curr <- featureData(jscs)@data[featureData(jscs)@data$geneID == geneID,c("geneID","countbinID"), drop=FALSE]
          rn <- row.names(curr);
          curr <- cbind.data.frame(rn,curr);
          names(curr)[1] <- "featureID";
          exonCt <- dim(curr)[1]
          NA.buffer <- matrix(rep(NA,blank.line.length * exonCt), nrow=exonCt, ncol=blank.line.length );
          #blankList <- lapply(1:exonCt, function(x){ blank.line });
          norCounts <- cbind.data.frame(get.norcounts.data(jscs,geneID),get.norcounts.data(jscs,geneID,vst.xform=FALSE), get.rawcounts.data(jscs,geneID))          
          curr <- cbind.data.frame(curr,NA.buffer,norCounts);
          names(curr) <- names.out;
          row.names(curr) <- rn;
          return( list(FALSE, curr) );
      }
   }
   
   pair.buffer <- myApply(as.list(1:length(gene.list)), FUN=get.one.genes.data);
   
   genes.data.list <- lapply(pair.buffer, "[[", 2);
   did.gene.converge <- sapply(pair.buffer, "[[", 1);
   
   if(verbose){
     message(paste("generated expression data for: ",length(genes.data.list),"top-level features (ie genes). ",sum(did.gene.converge)," converged. ",sum(! did.gene.converge)," failed."));
     #message("head(genes.data.list):");
     #print(head(genes.data.list));
   }

   out <- do.call(rbind.data.frame,genes.data.list);
   names(out) <- names.out;
   
   #############################
   #NEW CODE: save data to jscs!
   message("Saving data to JunctionSeqCountSet. (",date(),")");
   conditionLevels <- levels(jscs@phenoData$condition);
   names(conditionLevels) <- conditionLevels;
   sampleNames <- colnames(counts(jscs));
   
   exprEstimate          <- out[,  strStartsWith(names(out), "expr_") ];
   exprEstimateVST       <- out[,  strStartsWith(names(out), "exprVST_") ];
   otherExprEstimate     <- do.call(cbind.data.frame, lapply( conditionLevels, function(x){ rep(NA,nrow(jscs))} )  );
   otherExprEstimateVST  <- do.call(cbind.data.frame, lapply( conditionLevels, function(x){ rep(NA,nrow(jscs))} )  );
   relExprEstimate       <- out[,  strStartsWith(names(out), "rExpr_") ];
   relExprEstimateVST    <- out[,  strStartsWith(names(out), "rExprVST_") ];
   
   normCounts      <- cbind( out[,  strStartsWith(names(out),"normCount_") ],     do.call(cbind.data.frame, lapply( 1:sampleCt, function(x){ rep(NA,nrow(jscs))} )  )   );
   normCountsVST    <- cbind( out[,  strStartsWith(names(out),"normCountVST") ],   do.call(cbind.data.frame, lapply( 1:sampleCt, function(x){ rep(NA,nrow(jscs))} )  )   );
   
   colnames(otherExprEstimate) <- paste0("geneExpr","_",conditionLevels);
   colnames(otherExprEstimateVST) <- paste0("geneExprVST","_",conditionLevels);
   colnames(relExprEstimate) <- paste0("relExpr","_",conditionLevels);
   colnames(relExprEstimateVST) <- paste0("relExprVST","_",conditionLevels);
   
   colnames(normCounts) <- c(paste0("normCount_",sampleNames), paste0("normGeneCount_",sampleNames));
   colnames(normCountsVST) <- c(paste0("normCountsVST_",sampleNames), paste0("normGeneCountsVST_",sampleNames));
   
   jscs@plottingEstimates    <- list(exprEstimate = exprEstimate,
                                     geneExprEstimate = otherExprEstimate,
                                     relExprEstimate = relExprEstimate, 
                                     normCounts = normCounts);
   jscs@plottingEstimatesVST <- list(exprEstimateVST = exprEstimateVST, 
                                     geneExprEstimateVST = otherExprEstimateVST, 
                                     relExprEstimateVST = relExprEstimateVST, 
                                     normCountsVST = normCountsVST);
   message("Done saving data to JunctionSeqCountSet. (",date(),")");
   #############################
   
   if(verbose){
     message(paste("Collapsed expression data into a single data frame. generateAllExpressionEstimates() complete."));
   }
   #out <- out[! is.na(out$geneID) ,];
   
   #for(i in 4:length(names(out))){
   #  out[,i] <- as.numeric(as.character(out[,i]));
   #}
   
   list(jscs, out);
}

##########################################################################
##########################################################################
##########################################################################
######### Plotting Functions:
##########################################################################
##########################################################################
##########################################################################

plotDispEsts <- function( jscs, ymin, linecol="#ff000080",
  xlab = "mean of normalized counts", ylab = "dispersion",
  log = "xy", cex = 0.45, ... )
{
  px = rowMeans( counts( jscs, normalized=TRUE ) )
  sel = (px>0)
  px = px[sel]

  py = fData(jscs)$dispBeforeSharing[sel]
  if(missing(ymin))
      ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)

  plot(px, pmax(py, ymin), xlab=xlab, ylab=ylab,
    log=log, pch=ifelse(py<ymin, 6, 16), cex=cex, ... )
  xg = 10^seq( -.5, 5, length.out=100 )
  fun = function(x) { jscs@dispFitCoefs[1] + jscs@dispFitCoefs[2] / x }
  lines( xg, fun(xg), col=linecol, lwd=4)
}

plotMA <- function(jscs, 
                           FDR.threshold = 0.05, 
                           fc.name = NULL,
                           fc.thresh = 1,
                           use.pch = 19,
                           smooth.nbin = 256, 
                           ylim = c( 1 / 1000,1000),
                           use.smoothScatter = TRUE,
                           label.counts = TRUE,
                           label.axes = c(TRUE,TRUE,FALSE,FALSE), 
                           show.labels = TRUE,
                           par.cex = 1, points.cex = 1, text.cex = 1,
                           lines.cex = 8,
                           anno.lwd = 2,
                           miniTicks = TRUE, ...){

  sig.thresh <- FDR.threshold;
  color.orange <- "red";
  fdata <- fData(jscs);
  
  if(is.null(fc.name)){
    fc.names <- names(fdata)[which(strStartsWith(names(fdata), "log2FC("))];
    if(length(fc.names) != 1){
      stop("Error: you must specify the name of the column containing the desired fold change variable.");
    }
    fc.name <- fc.names[1];
  }
  if( strStartsWith(fc.name, "log2FC(") ){
    lvlsString <- substr(fc.name, 8, nchar(fc.name) - 1);
    lvls <- strsplit(lvlsString,"/",fixed=TRUE)[[1]];
  } else {
    lvls <- NULL;
    label.counts <- FALSE;
  }
  
  X <- fData(jscs)$meanBase;
  Y <- fData(jscs)[[fc.name]];
  pvals <- fData(jscs)[["padjust"]];
  
  xlim <- c(log10(1),log10( max(X,na.rm=TRUE) ))
  
  if(is.null(ylim)){
    ymax <- 2 ^ max(abs(Y),na.rm=TRUE);
    ylim <- c(1 / ymax, ymax);
    
    ymax.exp <- 2 ^ ymax;
  }
  ylim.exp <- log2( ylim );
  
  markerlines.lty <- "dashed";
  
  point.color <- sapply(pvals,FUN=function(p){
    if(is.na(p)){
      "grey"
    } else if(p < sig.thresh){
      "red"
    } else {
      "black"
    }
  });
  point.pch <- sapply(Y,function(fc){
    if(is.na(fc)){
      46
    } else if(fc > ylim.exp[2]){
      94
    } else if(fc < ylim.exp[1]){
      118
    } else {
      use.pch
    }
  })
  limited.log2fc <- sapply(Y,function(fc){
    if(is.na(fc)){
      fc
    } else if(fc > ylim.exp[2]){
      ylim.exp[2]
    } else if(fc < ylim.exp[1]){
      ylim.exp[1]
    } else {
      fc
    }
  })

  point.color[point.color == "red"] <- ifelse(f.na((2 ^ abs(limited.log2fc)) < fc.thresh)[point.color == "red"], color.orange,"red");
  is.sig <- f.na(pvals < sig.thresh);
  is.over <- f.na(Y > ylim.exp[2] | Y < ylim.exp[1]);

  par(cex = par.cex);
  if(use.smoothScatter){
    smoothScatter(log10(X),Y,nbin = smooth.nbin, nrpoints = 0,
                 xlim = xlim,ylim=ylim.exp,xlab="",ylab="",axes=F, ...);
    points(log10(X[is.over & (! is.sig)]),limited.log2fc[is.over & (! is.sig)], col= point.color[is.over & (! is.sig)],pch= point.pch[is.over & (! is.sig)],cex= points.cex, ...);
  } else {
    plot(log10(X[! is.sig]),limited.log2fc[! is.sig], col= point.color[! is.sig], pch= point.pch[! is.sig],cex= points.cex,
                 xlim = xlim,ylim=ylim.exp,xlab="",ylab="",axes=F, ...);
  }
  
  points(log10(X[is.sig]),limited.log2fc[is.sig], col= point.color[is.sig],pch= point.pch[is.sig],cex= points.cex, ...);
  abline(h=0,col=rgb(255,0,0,100,maxColorValue=255), lwd = anno.lwd, cex = lines.cex)

  box(lwd = anno.lwd);

  if(miniTicks){
    log10.label <- c(-3,-2,-1,0,1,2,3);
    log10.expression.label <- c(expression(10^"-3"), expression(10^"-2"), expression(10^"-1"), 1, expression(10), expression(10^2), expression(10^3));
  } else {
    log10.label <- c(-3,-2,-1,-log10(4),0,log10(4),1,2,3);
    log10.expression.label <- c(expression(10^"-3"), expression(10^"-2"), expression(10^"-1"),expression(1/4), 1, 4, expression(10), expression(10^2), expression(10^3));
  }
  #log10.expression.label <- c(expression(10^-3), expression(10^-2), expression(10^-1),"\u00BC", 1, 4, expression(10), expression(10^2), expression(10^3));
  log10.at <- log10.label * log2(10);

  miniTicks.temp <- c(2,3,4,5,6,7,8,9)
  log10.miniTicks <- log10(c(miniTicks.temp, miniTicks.temp * 10, miniTicks.temp * 100, 1 / miniTicks.temp, 1 / (miniTicks.temp * 10), 1 / (miniTicks.temp * 100) ));
  miniTicks.pts <- log10.miniTicks * log2(10);

  #axis(1);
  if(miniTicks){
     x.at <- log10(c(1,10,100,1000,10000,100000));
     x.label <- c(1, expression(10), expression(10^2), expression(10^3), expression(10^4), expression(10^5));
  } else {
     x.at <- log10(c(1,5,10,100,1000,10000,100000));
     x.label <- c(1,5, expression(10), expression(10^2), expression(10^3), expression(10^4), expression(10^5));
  }

  miniTicks.x <- log10(c(miniTicks.temp, miniTicks.temp * 10, miniTicks.temp * 100, miniTicks.temp * 1000, miniTicks.temp * 10000, miniTicks.temp * 100000, miniTicks.temp * 1000000));

  #axis(2);
  #axis(2,at=log10.at,label=log10.label);

  if(label.axes[1]){
    axis(1,at=x.at,labels=x.label, lwd = anno.lwd, lwd.ticks = anno.lwd,  cex.axis = text.cex);
  } else {
    axis(1,at=x.at,labels=FALSE, lwd = anno.lwd, lwd.ticks = anno.lwd, cex.axis = text.cex);  
  }
  if(miniTicks ) axis(1,at=miniTicks.x,labels=FALSE, lwd = anno.lwd / 2, lwd.ticks = anno.lwd / 2, tcl = -0.25, cex.axis = text.cex);
  if(label.axes[2]){
    axis(2,at=log10.at,labels=log10.expression.label,las=1, lwd = anno.lwd, lwd.ticks = anno.lwd, cex.axis = text.cex);
    if(miniTicks ) axis(2,at = miniTicks.pts, labels = F, lwd = anno.lwd/2, lwd.ticks = anno.lwd / 2, tcl = -0.25)
  } else {
    axis(2,at=log10.at,labels=FALSE,las=1, lwd = anno.lwd, lwd.ticks = anno.lwd, tcl = -0.25, cex.axis = text.cex);
    if(miniTicks ) axis(2,at = miniTicks.pts, labels = F, lwd = anno.lwd/2, lwd.ticks = anno.lwd / 2, tcl = -0.2)
  }
  if(label.axes[4]){
    axis(4,at=log10.at,labels=log10.expression.label,las=1, lwd = anno.lwd, lwd.ticks = anno.lwd, cex.axis = text.cex);
    if(miniTicks ) axis(4,at = miniTicks.pts, labels = F, lwd = anno.lwd/2, lwd.ticks = anno.lwd / 2, tcl = -0.25)
  } else {
    axis(4,at=log10.at,labels=FALSE,las=1, lwd = anno.lwd, lwd.ticks = anno.lwd, tcl = -0.25, cex.axis = text.cex);
    if(miniTicks ) axis(4,at = miniTicks.pts, labels = F, lwd = anno.lwd/2, lwd.ticks = anno.lwd / 2, tcl = -0.2) 
  }

  box(lwd=anno.lwd,cex=lines.cex);
  
  if(label.counts){
    if(is.na(fc.thresh)){
      sig.count.plus <- sum(f.na(Y > 0 & pvals < sig.thresh));
      sig.count.minus <- sum(f.na(Y < 0 & pvals < sig.thresh));
    } else {
      sig.count.plus <- sum(f.na(Y > log2(fc.thresh) & pvals < sig.thresh));
      sig.count.minus <- sum(f.na(Y < -log2(fc.thresh) & pvals < sig.thresh));
      
      if(fc.thresh != 1){
        abline(h=log2(fc.thresh),col="gray",lty=markerlines.lty, lwd = anno.lwd, cex = lines.cex);
        abline(h=-log2(fc.thresh),col="gray",lty=markerlines.lty, lwd = anno.lwd, cex = lines.cex);
      }
    }
    legend("topright",   legend=c(paste0(sig.count.plus, " higher in ",lvls[1]," (FC > ",fc.thresh,")")),bty="n",cex=text.cex);
    legend("bottomright",legend=c(paste0(sig.count.minus," higher in ",lvls[2]," (FC < ",1 / fc.thresh,")")),bty="n",cex=text.cex);
  }
  
  if(show.labels){
    title(main = "MA Plot",
          xlab = "Mean Normalized Coverage",
          ylab = "Fold Change",
          cex.main = 1.2 * text.cex,
          cex.lab = 1 * text.cex
          );
  }
}


plotMA.OLDVER <- function(jscs, FDR.threshold, fc.name = NULL, ylim = c(-5,5), ...){
    fdata <- fData(jscs);
    
    if(is.null(fc.name)){
      fc.names <- names(fdata)[which(strStartsWith(names(fdata), "log2FC("))];
      if(length(fc.names) != 1){
        stop("Error: you must specify the name of the column containing the desired fold change variable.");
      }
      fc.name <- fc.names[1];
    }

    #which.is.fc <- which(substr(names(merged.results.data),1,8) == "log2fold")
    #fc.name <- names(merged.results.data)[fc.col];
    ma.frame <- data.frame(baseMean=fdata$meanBase);
    ma.frame[[fc.name]] <- fdata[[fc.name]];
    ma.frame[["padj"]] <- fdata$padjust;
    
    #, log2FoldChange = merged.results.data[,which.is.fc], padj = merged.results.data$padjust);
    #CairoPNG(paste(outfile.prefix,"ma.plot",".png",sep=""),height=2500,width=1600,pointsize = 18);
    plotMA_HELPER( ma.frame, FDR=FDR.threshold, ylim=ylim, ylab = fc.name, ...)
    #dev.off();
}


plotMA_HELPER <- function(ma.frame, ylim,FDR=0.05,
  col = ifelse(ma.frame$padj>=FDR, "gray32", "red3"),
  linecol = "#ff000080",
  xlab = "mean of normalized counts", ylab = expression(log[2]~fold~change),
  log = "x", cex=0.45, ...)
{
  #if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x))))
  #  stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")

  ma.frame = subset(ma.frame, ma.frame$baseMean != 0)
  py = ma.frame[,2];
  if(missing(ylim)){
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  }
  plot(ma.frame$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  abline(h=0, lwd=4, col=linecol);
}

readJunctionSeqCounts <- function(countfiles = NULL, countdata = NULL,
                                  samplenames,  design, 
                                  flat.gff.file=NULL,
                                  use.splice.sites = TRUE, use.novel.splice.sites = TRUE, use.exons = FALSE,
                                  verbose = TRUE)
{
   if(verbose) {
      message("-> STARTING readJunctionSeqCounts (",date(),")");
      #message("---> RJSC: countfiles: ",paste0(countfiles, collapse=","));
      #message("---> RJSC: samplenames: ",paste0(samplenames, collapse=","));
      #message("---> RJSC: flat.gff.file: ",flat.gff.file);
      #message("---> RJSC: use.splice.sites:",use.splice.sites);
      #message("---> RJSC: use.novel.splice.sites:",use.novel.splice.sites);
      #message("---> RJSC: use.exons:",use.exons);
   }
   
   if((is.null(countfiles) & is.null(countdata))){
     stop("Fatal error: Either countfiles OR countdata must be set! Both are null!");
   }
   if(  (!is.null(countfiles)) & (!is.null(countdata))   ){
     stop("Fatal error: Either countfiles OR countdata must be set! Both are non-null!");
   }
   
   stopifnot( class(design) == "data.frame" );
   
   for(i in 1:ncol(design)){
     if( ! is.factor(design[[i]])){
        stop("ERROR: design must be a data.frame composed entirely of factors!");
     }
   }
   
   if(! is.null(countfiles)){
     lf <- lapply( countfiles, function(x)
        read.table( x, header=FALSE,stringsAsFactors=FALSE ) )
   } else {
     lf <- countdata;
   }
   
   if( !all( sapply( lf[-1], function(x) all( x$V1 == lf[1]$V1 ) ) ) )
      stop( "Count files have differing gene ID column." )
   if(verbose) message("---> File read complete.");  

   dcounts <- sapply( lf, `[[`, "V2" )
   rownames(dcounts) <- lf[[1]][,1]
   #rownames(dcounts)
   dcounts <- dcounts[ substr(rownames(dcounts),1,1)!="_", ]

   bin.type <- sapply( rownames(dcounts), 
                     function(x){
                        substr(strsplit(x, ":",fixed=TRUE)[[1]][2],1,1)
                     });
   
   if(verbose) message(paste0("---> Extracted counts. Found ",dim(dcounts)[1]," features so far."));  
   
   geneCountTable <- dcounts[bin.type == "A",, drop=FALSE];
   rownames(geneCountTable) <- sapply(strsplit(rownames(geneCountTable), ":"),"[[",1);
   colnames(geneCountTable) <- as.character(samplenames);
   use.bins <- bin.type != "A"
   
   if(verbose) message(paste0("---> Extracted gene-level counts. Found: ",dim(geneCountTable)[1], " aggregate genes."));
   if(verbose) message(paste0("---> Removed aggregate gene features. Found: ",sum(use.bins), " features to be included so far."))
   
   if(! use.splice.sites & ! use.novel.splice.sites & ! use.exons){
      stop("FATAL ERROR: At least one of: use.splice.sites, use.novel.splice.sites, or use.exons must be set to TRUE. Otherwise you've got no data to test!");
   }
   if(! use.splice.sites){
      use.bins <- use.bins & bin.type != "J" & bin.type != "N";
      if(verbose) message(paste0("---> Removed splice junction features. Found: ",sum(use.bins), " features to be included so far."))
      
   } else if(! use.novel.splice.sites){
      #dcounts <- dcounts[bin.type != "N",]
      use.bins <- use.bins & bin.type != "N";
      if(verbose) message(paste0("---> Removed novel splice junction features. Found: ",sum(use.bins), " features to be included so far."))
   }
   if(! use.exons){
      use.bins <- use.bins & bin.type != "E";
      if(verbose) message(paste0("---> Removed exon features. Found: ",sum(use.bins), " features to be included so far."))
   }
   if(verbose) message(paste0("---> Final feature count: ",sum(use.bins), " features to be included in the analysis."))
   dcounts <- dcounts[use.bins,, drop=FALSE];
   bin.type <- bin.type[use.bins];
   
   if(verbose) message("---> Extracted feature counts.");  
   
   #rownames(dcounts) <- sub(":", ":E", rownames(dcounts))
   colnames(dcounts) <- as.character(samplenames);
   splitted <- strsplit(rownames(dcounts), ":")
   exons <- sapply(splitted, "[[", 2)
   genesrle <- sapply( splitted, "[[", 1)
   
   if(verbose) message("---> counts complete.");  
   
   if(! is.null(flat.gff.file)){
      if(verbose) message("-----> reading annotation...");  
      anno.data <- readAnnotationData(flat.gff.file);
      if(verbose) message("-----> formatting annotation...");  
      #read.table(flat.gff.file,stringsAsFactors=FALSE,header=TRUE);
      #featureName     featureType     chrom   start   end     strand  gene_id part_number     transcripts
      
      exoninfo <- data.frame(chr = anno.data$chrom, start = anno.data$start, end = anno.data$end, strand = anno.data$strand);
      if(verbose) message("-----> initial generation...");  
      rownames(exoninfo) <- anno.data$featureName;
      
      transcripts <- anno.data$transcripts;
      transcripts <- gsub("\\+", ";", transcripts)
      names(transcripts) <- rownames(exoninfo)
      
      matching <- match(rownames(dcounts), rownames(exoninfo))
      if(any(is.na(matching))){
         stop("FATAL ERROR! Annotation file appears to be missing information! Are you sure you are using the correct flattened annotation file, created by prepare_annotation_with_splices.py?");
      }
      if(verbose) message("-----> creating jscs...");  
      jscs <- newJunctionSeqCountSet(countData=dcounts, geneCountData = geneCountTable, design=design, geneIDs=genesrle, countbinIDs=exons, featureIntervals=exoninfo[matching,], transcripts=transcripts[matching])
      
      featureChar <- substr(fData(jscs)$countbinID,1,1);
      fData(jscs)[["featureType"]] <- ifelse(featureChar == "E","exonic_part",ifelse(featureChar == "J", "splice_site", "novel_splice_site"))
      
      jscs@annotationFile <- flat.gff.file
      pData(jscs)$countfiles <- countfiles
      
      if(verbose) message("-----> generating count vectors... (",date(),")");
      jscs@countVectors <- getAllJunctionSeqCountVectors(jscs);
      if(verbose) message("-----> count vectors generated (",date(),")");
      
      if(verbose) message("-> FINISHED readJunctionSeqCounts (",date(),")");  
      return(jscs);
      #return(list(geneCountTable, jscs))

   }else{
      #return(list(geneCountTable, newJunctionSeqCountSet(countData=dcounts, design=design, geneIDs=genesrle, countbinIDs=exons)))
      if(verbose) message("-> FINISHED readJunctionSeqCounts (",date(),")"); 
      message("Warning: flat gff annotation not set! While technically optional, running without the annotation data may make interpretation of the data difficult.");
      warning("flat gff annotation not set! While technically optional, running without the annotation data may make interpretation of the data difficult.");
      jscs <- newJunctionSeqCountSet(countData=dcounts, design=design, geneIDs=genesrle, countbinIDs=exons);
      
      if(verbose) message("-----> generating count vectors... (",date(),")");
      jscs@countVectors <- getAllJunctionSeqCountVectors(jscs);
      if(verbose) message("-----> count vectors generated (",date(),")");
      
      return(jscs);
   }
}


#setMethod("estimateDispersions", signature(object="JunctionSeqCountSet"),
estimateJunctionSeqDispersions <- function( jscs, 
                                            test.formula1 = formula(~ sample + countbin + condition : countbin),
                                            minCount=10, nCores=1, test.aggregated.genes = FALSE, 
                                            use.alternate.method  = TRUE,
                                            verbose = TRUE){
   if(verbose) message("---> STARTING estimateJunctionSeqDispersions: (",date(),")");
   stopifnot( inherits( jscs, "JunctionSeqCountSet" ) )
   if( all( is.na( sizeFactors( jscs )) ) ){
     stop("Please calculate size factors before estimating dispersions\n")
   }
   fData(jscs)$meanBase <- rowMeans(counts(jscs, normalized=TRUE));
   
   testable <- rowSums(counts(jscs)) >= minCount;
   if(verbose) message("-----> ejsd: ",sum(testable)," counting bins dropped due to low coverage (total count < ",minCount,")");
   
   is.aggregate <- grepl("+", geneIDs(jscs), fixed=TRUE)
   if(verbose) message("-----> ejsd: ",sum(is.aggregate)," counting bins are from ",length(unique(grep("+", geneIDs(jscs), fixed=TRUE)))," multigene aggregates (ie. overlapping genes).");
   
   if(! test.aggregated.genes){
     if(verbose) {
       message("-----> ejsd: ",sum(is.aggregate & testable)," counting bins dropped because they were from multigene aggregates.");
       message("             (To include these aggregate genes, set test.aggregated.genes = TRUE)");
     }
     testable <- testable & (! is.aggregate);
   }
   
   onlyOneJunctionDropCt <- 0;
   for( r in split( seq_len(nrow(jscs)), geneIDs(jscs) ) ) {
     if( sum( testable[r] ) == 1 ){
       testable[r] <- FALSE;
       onlyOneJunctionDropCt <- onlyOneJunctionDropCt + 1;
     }
   }
   message("-----> ejsd: ",sum(onlyOneJunctionDropCt)," counting bins dropped because they belong to a gene with only 1 testable bin.");
   
   fData(jscs)$testable <- testable;
   
   myApply <- getMyApply(nCores);

   jscs@formulas[["formulaDispersion"]] <- deparse(test.formula1)

   rows <- seq_len(nrow(jscs))
   modelFrame <- constructModelFrame( jscs )
   mm <- rmDepCols( model.matrix( test.formula1, modelFrame ) )
   disps <- myApply( rows, function(i) {
      if( verbose & i %% 1000 == 0 ){
         message(paste0("-------> estimateJunctionSeqDispersions: (Calculating dispersion for feature ",i," of ",nrow(jscs),")","(",date(),")"));
      }
         
      if( fData(jscs)$testable[i] ) {
         #a <- try( estimateFeatureDispersion( jscs, geneIDs(jscs)[i], countbinIDs(jscs)[i], modelFrame, mm , use.alternate.method  = use.alternate.method ), silent=TRUE)
         a <- try( estimateFeatureDispersionFromRow( jscs, i, modelFrame, mm , use.alternate.method  = use.alternate.method ), silent=TRUE)
         if( inherits( a, "try-error" ) ) {
            warning( paste0( sprintf("Unable to estimate dispersions for %s:%s", as.character( geneIDs(jscs)[i] ), countbinIDs(jscs)[i]),
                         "\nReason: ", a
                    ))
            NA 
         } else{
            a 
         }
      }else{
          NA 
      }
   })
   names(disps) <- featureNames(jscs)
   fData(jscs)[names(disps), "dispBeforeSharing"] <- unlist(disps)
   fData(jscs)$testable[which( is.na( fData(jscs)$dispBeforeSharing ) )] <- FALSE
   if(verbose) message("---> FINISHED estimateJunctionSeqDispersions (",date(),")");
   jscs
} #)



testJunctionsForDJU <- function( jscs,
                                test.formula0 = formula(~ sample + countbin), 
                                test.formula1 = formula(~ sample + countbin + condition : countbin),
                                dispColumn="dispersion", nCores=1 , use.alternate.method = TRUE, #keep.debug.model.data = FALSE,
                                verbose = TRUE){
  keep.debug.model.data <- TRUE;
  stopifnot( inherits( jscs, "JunctionSeqCountSet" ) )
   if( all( is.na( sizeFactors( jscs )))) {
     stop("Please calculate size factors before estimating dispersions\n")
   } 
   if( all( is.na( fData(jscs)[,dispColumn] ) ) ){
     stop("Please estimate dispersions before calling this function\n")
   }
   myApply <- getMyApply(nCores);
   
   jscs@formulas[["test.formula0"]] <- deparse(test.formula0)
   jscs@formulas[["test.formula1"]] <- deparse(test.formula1)
   rows <- seq_len(nrow(jscs));
   #names(rows) <- featureNames( jscs );
   modelFrame <- constructModelFrame( jscs )
   mm0 <- rmDepCols( model.matrix( test.formula0, modelFrame ) )
   mm1 <- rmDepCols( model.matrix( test.formula1, modelFrame ) )
   
#########################################
   fitExpToVar <- "condition"
   keepCoefs <- which(attr(mm1,"assign") == length( attr(terms(test.formula1),"term.labels"))  )
   #keepCoefs <- which( substr(colnames(mm1),1,13) == "countbinthis:");
   keepCoefNames <- paste0("HtestCoef(",colnames(mm1)[keepCoefs],")");
   conditionLevels <- levels(modelFrame[[fitExpToVar]]);
#########################################
   
   mdl.out <- myApply( rows,
     function(i) {
       if( verbose & i %% 1000 == 0 ){
         message(paste0("-------> testJunctionsForDJU: (testing for DJU on feature ",i," of ",nrow(jscs),")","(",date(),")"));
       }
          
       if( fData(jscs)$testable[i] ) {
          #a <- try( 
          #            testFeatureForDJU( jscs, gct, geneIDs(jscs)[i], countbinIDs(jscs)[i], modelFrame, mm0, mm1, fData(jscs)[i, dispColumn] , use.alternate.method = use.alternate.method) 
          #         )
          #if( any(inherits( a, "try-error" ) )) {
          #   warning( sprintf("Unable to calculate p-values for %s:%s\n", as.character( geneIDs(jscs)[i] ), countbinIDs(jscs)[i]) )
          #   #list(coefficient = NA,pval = NA, disp = NA, countVector = NA);
          #   list(NA,NA);
          #} else {
          #   a ;
          #}
          testFeatureForDJU.fromRow(test.formula1, jscs, i, modelFrame, mm0, mm1, fData(jscs)[i, dispColumn] , keepCoefs = keepCoefs, use.alternate.method = use.alternate.method);
          #testFeatureForDJU(test.formula1, jscs, geneIDs(jscs)[i], countbinIDs(jscs)[i], modelFrame, mm0, mm1, fData(jscs)[i, dispColumn] , keepCoefs = keepCoefs, use.alternate.method = use.alternate.method) 
       } else {
         return(list(coefficient = rep(NA,length(keepCoefs)), 
                     logFC = rep(NA, length(conditionLevels) - 1),
                     pval = NA, 
                     disp = NA, 
                     countVector = jscs@countVectors[i,], 
                     fit = list(fitH0 = NA, fitH1 = NA) 
                ));
       }
     })
    
    #list(coefficient = coefficient, logFC = logFC,pval = pval, disp = disp, countVector = countVector)
    pvals <- sapply(mdl.out, "[[", "pval");
    coefficient <- do.call(rbind.data.frame, lapply(mdl.out,"[[", "coefficient"));
    logFC <- do.call(rbind.data.frame, lapply(mdl.out,"[[", "logFC"));
    
    modelFits <- lapply(mdl.out, "[[", "pval");
    names(modelFits) <- featureNames( jscs );
    jscs@modelFitForHypothesisTest <- modelFits;
    
    #if(keep.debug.model.data){
    #   countVectors <- t(as.matrix(sapply(mdl.out,"[[",4)));
    #   rownames(countVectors) <- featureNames(jscs);
    #   colnames(countVectors) <- paste0("modelVectorCT_",as.character(modelFrame$sample),"_",as.character(modelFrame$countbin));
    #   jscs@countVectors  <- countVectors;
    #}
       
    #names( pvals ) <- featureNames( jscs );
    
    message("dim(coefficient) = ",paste0(dim(coefficient),collapse=","));
    message("dim(logFC) = ",paste0(dim(logFC),collapse=","));
    #message("dim(pvals) = ",paste0(dim(pvals),collapse=","));
    
    #if(length(keepCoefs) > 1){
    #rownames(coefficient) <- featureNames( jscs );
    colnames(coefficient) <- keepCoefNames
    #rownames(logFC) <- featureNames( jscs );
    conditionLevels <- levels(pData(jscs)[["condition"]]);
    colnames(logFC) <- paste0("logFC(",conditionLevels[-1],"/",conditionLevels[1],")")
    #colnames(pvals) <- "pvalue";
    
    #fData(jscs) <- cbind.data.frame(fData(jscs), coefficient, logFC);
    for(CN in colnames(coefficient)){
      fData(jscs)[[CN]] <- coefficient[[CN]];
    }
    
    fData(jscs)$pvalue <- pvals;
    fData(jscs)$padjust <- p.adjust( fData(jscs)$pvalue, method="BH" )

    #   fData(jscs) <- cbind.data.frame(fData(jscs), coefficient);
    #} else {
    #   names( coefficient ) <- featureNames( jscs );
    #   fData(jscs)[[ keepCoefNames ]] <- coefficient;
    #}
    
    #fData(jscs)[names(pvals), "pvalue"] <- pvals;
    
    jscs;
}

readAnnotationData <- function(flat.gff.file){
     aggregates<-read.delim(flat.gff.file, stringsAsFactors=FALSE, header=FALSE)
     colnames(aggregates)<-c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
     aggregates<-aggregates[which(aggregates$class != "aggregate_gene"),]
     aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
     gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
     transcripts <- gsub(".*tx_set\\s(\\S+).*", "\\1", aggregates$attr)
     part_number <- gsub(".*num\\s(\\S+).*", "\\1", aggregates$attr)
     
     featureCode <- ifelse(aggregates$class == "exonic_part","E", ifelse(aggregates$class == "splice_site","J", ifelse(aggregates$class == "novel_splice_site","N","?")));
     
     featureName <- paste0(gene_id,":",featureCode,part_number);
     
     out <- data.frame(featureName = featureName,
                       featureType = aggregates$class,
                       chrom = aggregates$chr,
                       start = aggregates$start,
                       end = aggregates$end,
                       strand = aggregates$strand,
                       gene_id = gene_id,
                       part_number = part_number,
                       transcripts = transcripts
                       );
                       
     out$start <- as.integer(out$start - 1);
     return(out);

}


#extractSigResults <- function(jscs, outfile.prefix = NULL, verbose = FALSE,nCores=1, sig.threshold = 0.05, use.adjusted.pvals = TRUE,
#						use.splice.sites = TRUE, use.novel.splice.sites = TRUE, use.exons = FALSE, gzip.output = TRUE){
#   #require(statmod)
#   #require(DEXSeq)
#   ##require(multicore)
#   if(verbose) { message(paste0("> gcr: Starting generateCompleteResults: ",date(),".")) }
#   if(verbose) { message(paste0("> gcr: generating.all.expression.estimates... ",date(),".")) }
#   
#   feature.is.sig <- fData(jscs)$pvalue < sig.threshold;
#   if(use.adjusted.pvals){
#      feature.is.sig <- fData(jscs)$padjust < sig.threshold;
#   }
#   feature.is.sig <- ifelse(is.na(feature.is.sig),FALSE,feature.is.sig);
#   
#   sig.gene.list <- unique(fData(jscs)$geneID[ feature.is.sig ]);
#   
#   expr.data <- generateSigExpressionEstimates(jscs,verbose=verbose,nCores=nCores, sig.gene.list = sig.gene.list);
#   if(verbose) { message(paste0("> gcr: generating.all.expression.estimates done.",date(),".")) }
#   
#   #if(! is.null(outfile.prefix)) write.table.gz(expr.data,file=paste(outfile.prefix,".expression.data.txt",sep=''),use.gzip=gzip.output,row.names=F,col.names=T,quote=F);
#   if(verbose) message(paste0("> gcr: Generated Expression Estimates. ",date()));
#
#   if(verbose) message(paste0("> gcr: merge 1: ",date()));
#   #if(! all( expr.data$featureID == row.names(featureData(jscs)@data) )){
#   #   stop("FATAL ERROR: all( expr.data$featureID == row.names(featureData(jscs)@data) ) is not TRUE!");
#   #}
#   expr.match <- match(expr.data$featureID, row.names(featureData(jscs)@data));
#   
#   fn <- names(expr.data)
#   fd <- featureData(jscs)@data[expr.match,];
#   featureChar <- substr(fd$countbinID,1,1);
#   fd$meanBase <- rowMeans(counts(jscs, normalized=TRUE))[expr.match];
#   fd$featureType <- ifelse(featureChar == "E","exonic_part",ifelse(featureChar == "J", "splice_site", "novel_splice_site"))
#   names(fd)[names(fd) == "chr"] <- "chrom";
#   
#   merged.data.1 <- cbind.data.frame(featureID = expr.data[,1], fd, expr.data[,4:length(fn)])
#   
#   if(verbose) message(paste0("> gcr: merge 1 complete",date()));
#
#   if(! is.null(outfile.prefix)){
#      write.table.gz(merged.data.1,file=paste(outfile.prefix,".merged.out.txt",sep=''),use.gzip = gzip.output,quote=F,row.names=F)
#   } else {
#      if(verbose) { message(paste0("> gcr: skipping write of merged data, because outfile.prefix is NULL. ",date(),".")) }
#   }
#   
#   return(merged.data.1);
#}
#ENSRNOG00000002095


plotJunctionSeqResultsForGene <- function(geneID, jscs, flat.gff.data, 
                              colorRed.FDR.threshold=0.05,
                              plot.type = "expr", 
                              displayTranscripts = FALSE,
                              color = NULL, 
                              use.vst = FALSE, use.log = FALSE,truncateBelowOne = TRUE,
                              exon.rescale.factor = 0.3,
                              label.p.vals = TRUE, 
                              plot.lwd = 3,axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                              par.cex = 1, anno.cex.text = 1,
                              anno.cex.axis=anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                              fit.countbin.names = TRUE,
                              debug.mode = FALSE,
                              plot.exon.results = FALSE, plot.splice.results = TRUE, plot.novel.splice.results = plot.splice.results, 
                              plot.untestable.results = FALSE,
                              show.strand.arrows = 20,
                              sort.features = TRUE,
                              drawCoordinates = TRUE,
                              arrows.length = 0.25,
                              title.main=NULL, title.ylab=NULL,
                              verbose=TRUE,
                              ...)
{
   FDR <- colorRed.FDR.threshold;
   merged.data <- fData(jscs);
   condition <- jscs@phenoData$condition;
   
   if(verbose){
      message("> plotJunctionSeqResultsForGene(): ", geneID, ", plot.type: ", plot.type);
   }

   rt <- merged.data$geneID == geneID;
   if(! plot.exon.results){
     rt <- rt & merged.data$featureType != "exonic_part" ;
   }
   if(! plot.splice.results){
     rt <- rt & merged.data$featureType != "splice_site" ;
   }
   if(! plot.novel.splice.results){
     rt <- rt & merged.data$featureType != "novel_splice_site" ;
   }
   if(! plot.untestable.results){
     rt <- rt & merged.data$testable;
   }
   
   rt <- which(rt);
   
   if(sort.features){
     rt <- rt[order(merged.data$start[rt],merged.data$end[rt])];
   }
   
   #rt.jscs <- which(featureData(jscs)$geneID == geneID)
   rango <- 1:length(rt)
   draw.legend <- TRUE;
   #condition.names <- sort(levels(condition));
   condition.names <- levels(condition);
   sample.names <- sampleNames(jscs@phenoData);
   colorcode.title <- TRUE;
   numcond <- length(condition.names);

if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached step 2.");

   rt.allExon <- which(flat.gff.data$gene_id==geneID & flat.gff.data$featureType == "exonic_part");
   rango.allExon <- 1:length(rt.allExon);

   #x.lim.coord <- c(,);

   # print("rt.allExon:");
   # print(rt.allExon);
   
   ####### DETERMINE COLORS, IF THE USER DOES NOT PROVIDE ONE PER SAMPLE THE COUNT WILL OBTAIN THEM CORRESPONDING TO THEIR DESIGN ####
   ##### determine colors if not provided by user ######
   if(is.null(color)){
      #color<-rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond)), maxColorValue=255, alpha=175)
      default.color.list <- c("red","blue","orange","green3","purple","cyan1", "magenta","yellow3","tan4")

      if(numcond > length(default.color.list)){
         message("Too many condition values, the default color selection may not look good! Set your own colors by setting the \"color\" parameter.");
         color<-rgb(colorRamp(c("#D7191C", "#FFFFBF", "#2B83BA"))(seq(0, 1, length.out=numcond)), maxColorValue=255, alpha=175)
      } else {
         color <- color2transparentVector(default.color.list[1:numcond], t = 175);
      }
   }
if(debug.mode & verbose) message(">    pJSRforGene(): ","color = ", paste0(color, collapse=","));
if(debug.mode & verbose) message(">    pJSRforGene(): ","length(color) = ", length(color));
if(debug.mode & verbose) message(">    pJSRforGene(): ","length(condition.names) = ", length(condition.names));
if(debug.mode & verbose) message(">    pJSRforGene(): ","condition.names = ", paste0(condition.names,collapse=","));
   
   names(color) <- condition.names;
   #print(color);
   
   y.axis.title <- "";
   main.title <- "";
   
if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached step 3.");

   #For reference:
   #jscs@plottingEstimates    <- list(exprEstimate = exprEstimate,
   #                                  geneExprEstimate = otherExprEstimate,
   #                                  relExprEstimate = relExprEstimate, 
   #                                  normCounts = normCounts);
   #jscs@plottingEstimatesVST <- list(exprEstimateVST = exprEstimateVST, 
   #                                   geneExprEstimateVST = otherExprEstimateVST, 
   #                                  relExprEstimateVST = relExprEstimateVST, 
   #                                  normCountsVST = normCountsVST);

   if(truncateBelowOne){
     convertY <- function(y){
       ifelse(y < 0, ifelse(is.infinite(y), INTERNAL.NINF.VALUE, (1-exp(y)) * INTERNAL.NINF.VALUE), y);
     }
   } else {
     convertY <- function(y){
       y
     }
   }

   #count <- merged.data[rt,c("expr_1","expr_2")]
   if(plot.type == "rExpr"){
     #count <- if(use.vst) merged.data[rt,paste("rExprVST_",condition.names,sep='')]
     #         else merged.data[rt,paste("rExpr_",condition.names,sep='')] 
     count <- if(use.vst) { 
                jscs@plottingEstimatesVST[["relExprEstimateVST"]][rt,, drop=FALSE];
              } else  if(use.log) {
                apply(log10(jscs@plottingEstimates[["relExprEstimate"]][rt,, drop=FALSE]),c(1,2),FUN=convertY);
              } else {
                jscs@plottingEstimates[["relExprEstimate"]][rt,, drop=FALSE];
              }
     
     color.count <- rep(color[condition.names],each=nrow(count))
     y.axis.title <- "Normalized Coverage, excluding gene-wise expression";
     main.title <- paste0("Normalized Coverage, Relative to Whole-Gene Expression (",geneID,")");
   } else if(plot.type == "normCounts"){
     #count <- if(use.vst) merged.data[rt,paste("normCountVST_",sample.names,sep='')]
     #         else merged.data[rt,paste("normCount_",sample.names,sep='')]
     count <- if(use.vst){ 
                jscs@plottingEstimatesVST[["normCountsVST"]][rt,1:length(sample.names), drop=FALSE] 
              } else if(use.log) {
                apply(log10(jscs@plottingEstimates[["normCounts"]][rt,1:length(sample.names), drop=FALSE]),c(1,2),FUN=convertY)
              } else {
                jscs@plottingEstimates[["normCounts"]][rt,1:length(sample.names), drop=FALSE]
              }
                    
     color.count <- rep(color[as.character(condition)], each=nrow(count))
     y.axis.title <- "Normalized Counts";
     main.title <- paste0("Normalized Counts (",geneID,")");
   } else if(plot.type == "rawCounts"){
     count <- if(use.vst){ 
                #jscs@plottingEstimatesVST[["normCountsVST"]][rt,1:length(sample.names), drop=FALSE];
                message("VST-transformed raw counts not supported!");
              } else if(use.log) {
                apply(log10(jscs@countVectors[rt,1:length(sample.names), drop=FALSE]),c(1,2),FUN=function(y){max(INTERNAL.NINF.VALUE,y)})
              } else {
                jscs@countVectors[rt,1:length(sample.names), drop=FALSE]
              }
     color.count <- rep(color[as.character(condition)], each=nrow(count))
     y.axis.title <- "log10 Raw Counts";
     main.title <- paste0("Raw Counts (",geneID,")");
   } else if(plot.type == "expr"){
     count <- if(use.vst){ 
                jscs@plottingEstimatesVST[["exprEstimateVST"]][rt,, drop=FALSE];
              } else if(use.log) {
                apply(log10(jscs@plottingEstimates[["exprEstimate"]][rt,, drop=FALSE]),c(1,2),FUN=convertY);
              } else {
                jscs@plottingEstimates[["exprEstimate"]][rt,, drop=FALSE];
              }
     color.count <- rep(color[condition.names],each=nrow(count))
     y.axis.title <- "Normalized Coverage";
     main.title <- paste0("Normalized Coverage (",geneID,")");
     #print(count);
   } else {
     stop(paste0("FATAL ERROR: Unknown plot type! plot.type = \"",plot.type,"\""));
   }
   if(! is.null(title.ylab)) y.axis.title <- title.ylab;
   if(! is.null(title.main)) main.title <- title.main;
   
   count <- as.matrix(count);
   
   #if(use.vst) main.title <- paste(main.title," (VST-transformed)",sep='');

if(debug.mode & verbose) message(">    pJSRforGene(): ", "Reached step 4.");

   #if(length(levels(condition)) == 2){
   #   if(colorcode.title){
   #     paste(main.title,"\n",sep='');
   #     for(i in 2:length(names(color))){
   #        
   #     }
   #     
   #     
   #   } else {
   #     main.title <- paste(main.title,"\n",condition.names[1]," vs ",condition.names[2],"",sep='')
   #   }   
   #}
   
   intervals<-(0:nrow(count))/nrow(count)
   #numcond<-length(unique(design(jscs, drop=FALSE)[[fitExpToVar]]))
   numexons<-nrow(count)
   each <- merged.data$padjust[rt]
   #exoncol<-ifelse(each<=FDR, "#8B0000", "dark green")
   #exoncol[is.na(exoncol)]<-"black"
   #colorlines <- ifelse(each<=FDR, "#FF000060", "lightgrey")
   exoncol<-ifelse(each<=FDR, "#F219ED", "#CCCCCC")
   exonlty <- rep(1,length(exoncol));
   exonlty[as.character(merged.data$featureType[rt]) == "novel_splice_site"] <- 2 
   exoncol[is.na(exoncol)]<-color2transparent("#CCCCCC",25)
   
   sig.feature <- which(ifelse(is.na(each),FALSE,ifelse(each <= FDR, TRUE,FALSE)));
   is.sig.feature <- ifelse(is.na(each),FALSE,ifelse(each <= FDR, TRUE,FALSE));
   #print(paste("Index of sig splices: ", sig.feature));
   #print(merged.data[rt,][sig.feature,]);

if(debug.mode & verbose) message(">    pJSRforGene(): ", "Reached step 5.");

   colorlines <- ifelse(each<=FDR, "#F219ED60", "#B3B3B360")   # vertical dashed lines
   colorlines[is.na(colorlines)] <- "#B3B3B360"
   colorlinesB <- ifelse(each<=FDR, "#9E109B", "#666666")  # slanted solid lines
   colorlinesB[is.na(colorlinesB)] <- "#666666"
   
   #exrt <- which(merged.data$geneID==geneID & merged.data$featureType=="exonic_part");
   
   ################## DETERMINE THE LAYOUT OF THE PLOT DEPENDING OF THE OPTIONS THE USER PROVIDES ###########
      sub <- data.frame(start=merged.data$start[rt], end=merged.data$end[rt], chr=merged.data$chr[rt], 
                      strand=merged.data$strand[rt], is.exon = (merged.data$featureType[rt] == "exonic_part"), featureType = merged.data$featureType[rt],
                      stringsAsFactors=F);
      sub.allExon <- data.frame(start=flat.gff.data$start[rt.allExon], 
                                end=flat.gff.data$end[rt.allExon], 
                                chr=flat.gff.data$chrom[rt.allExon], 
                                strand=flat.gff.data$strand[rt.allExon], 
                                is.exon = (flat.gff.data$featureType[rt.allExon] == "exonic_part"), 
                                featureID = as.character(flat.gff.data$featureName[rt.allExon]), stringsAsFactors=F);

   sig.exon.names <- sub$featureID[is.sig.feature & sub$is.exon];
   allExon.isSig <- sub.allExon$featureID %in% sig.exon.names;
   allExon.exonCol <- ifelse(allExon.isSig, "#F219ED", "#CCCCCC");

   #TEMPORARY DEBUGGING MESSAGES:
if(debug.mode & verbose){
     message(">    Debugging Info:");
     message("        Exons:");
     message("           ", paste0(names(sub.allExon),collapse='\t'));
     for(i in 1:length(sub.allExon$start)){
       message("           ", paste0(sub.allExon[i,],collapse='\t'));
     }
     message(">       dim(count) = ", paste0(dim(count), collapse=","));
     message(">       rt: ", paste0(rt, collapse=","));
     message(">       rango: ", paste0(rango, collapse=","));
     message(">       ncol(count): ", ncol(count));
}
   #

      sub.sig <- sub[sig.feature,, drop=FALSE];
       #print(sub.sig);
       #print("sub.allExon");
       #print(sub.allExon);
      rel.calc.min <- min(sub$start,sub.allExon$start)
      rel.calc.max <- max(sub$end, sub.allExon$end)
      
      if(is.na(exon.rescale.factor) | exon.rescale.factor <= 0 | exon.rescale.factor >= 1){
        rel <- (data.frame(sub$start, sub$end))-rel.calc.min;
        rel <- rel/(rel.calc.max - rel.calc.min);
        rescale.iv <- NULL;
      } else {
        #message("attempting rescale...");
        rescale.iv <- generate.interval.scale(
          data.frame(
            start = c(sub.allExon$start, sub$start),
            end = c(sub.allExon$end, sub$end),
            is.exon = c(sub.allExon$is.exon, sub$is.exon)
          ),exon.rescale.factor);
        rel <- data.frame(start = rescale.coords(sub$start,rescale.iv), 
                          end   = rescale.coords(sub$end,  rescale.iv));
        #message("rescale.iv:");
        #print(rescale.iv);
        #message("rel:");
        #print(rel);
      }

      
      
      
      if(displayTranscripts==TRUE){
         transcripts <- sapply(sapply(flat.gff.data$transcripts[rt.allExon],toString), function(x){strsplit(x, "+",fixed=TRUE)})
         trans <- Reduce(union, transcripts)
         if(length(trans) > 42){
            warning("This gene contains more than 42 transcripts annotated, only the first 42 will be plotted\n")
            trans <- trans[1:42];
         }
         mat <- 1:(3+min(length(trans), 42)) ## max support from transcripts is 45, which seems to be the max for the layout supported by graphics
         hei<-c(8, 1, 1.5, rep(1.5, min(length(trans), 42)))
         
      }else{
         mat<-1:3
         hei<-c(5, 1, 1.5)
         
      }
      hei <- c(hei, .2)
      mat <- c(mat, length(mat)+1)
      layout(matrix(mat), heights=hei)
      par(mar=c(2, 4, 4, 2), cex = par.cex)

if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached step 6.");

      ylimn <- c(min(min(count,na.rm=TRUE),0), max(count, na.rm=TRUE));
      #if(plot.type == "rawCounts" & use.log) ylimn[1] <- INTERNAL.NINF.VALUE;
      if((! use.vst) & use.log ) ylimn[1] <- INTERNAL.NINF.VALUE;
      
      p.values.labels <- ifelse(each<=FDR, format(each,digits=3), "");

      drawPlot(matr=count, ylimn,jscs, 
               intervals, rango, textAxis=y.axis.title, 
               rt=rt, color.count=color.count, 
               colorlines=colorlines, 
               countbinIDs = merged.data$countbinID[rt],
               use.vst=use.vst, use.log = use.log,plot.type=plot.type,
               main.title=main.title,draw.legend=draw.legend,
               color.key=color,condition.names=condition.names,
               p.values=p.values.labels,draw.p.values=TRUE, 
               plot.lwd=plot.lwd, axes.lwd = axes.lwd, 
               anno.lwd = anno.lwd, par.cex = par.cex, 
               anno.cex.text = anno.cex.text,
               anno.cex.axis=anno.cex.axis, 
               anno.cex.main=anno.cex.main,
               fit.countbin.names = fit.countbin.names,
               ...);
      if(plot.type == "logRawCounts")  segments(0,INTERNAL.NINF.VALUE,2,INTERNAL.NINF.VALUE,lty="dotted",col="black", lwd = plot.lwd,...); 
      if((! use.vst) & use.log ) segments(0,INTERNAL.NINF.VALUE,2,INTERNAL.NINF.VALUE,lty="dotted",col="black", lwd = plot.lwd,...); 
      
if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached step 7.");

      #box()
      #########PLOT THE GENE MODEL:
      par(mar=c(0, 4, 0, 2), cex = par.cex)
      
      
      plot.new()
      #par();
      #plot.new(cex = anno.cex.text)
      # lines linking exons / splices to their column.
      segments(
                    apply((
                            rbind(rel[rango,2], rel[rango, 1])
                           ), 2, median), 
                    0, #par("usr")[3], (old version: lines are connected.)
                    apply(
                            rbind(intervals[rango], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2)), 
                            2, 
                            median), 
                    1, col=exoncol, lty = exonlty, lwd = anno.lwd, cex = anno.cex.text,cex.axis=anno.cex.axis, cex.main=anno.cex.main, xpd=FALSE, ...) #col=colorlinesB,...)
      #axes.lwd = axes.lwd, anno.lwd = anno.lwd,
      par(mar=c(0, 4, 0, 2), cex = par.cex);
      
      drawGene(rel.calc.min, rel.calc.max, tr=sub, tr.allExon=sub.allExon, rango, rescale.iv = rescale.iv, exoncol=exoncol,allExon.exonCol=allExon.exonCol, names, trName="Gene model", anno.cex.text=anno.cex.text, par.cex = par.cex, exonlty = exonlty, anno.lwd=anno.lwd, show.strand.arrows = show.strand.arrows, cex.axis=anno.cex.axis, cex.main=anno.cex.main, arrows.length = arrows.length,...)
      
         if(displayTranscripts){
            for(i in 1:min(length(trans), 42)){
               logicexons <- sapply(transcripts, function(x){length(which(x==trans[i]))})
               tr <- data.frame(start = flat.gff.data$start[rt.allExon][logicexons==1], end = flat.gff.data$end[rt.allExon][logicexons==1], featureType = flat.gff.data$featureType[rt.allExon][logicexons==1], stringsAsFactors = F);
               tr <- tr[tr$featureType == "exonic_part",, drop=FALSE]
               
               curr.exoncol <- ifelse(allExon.isSig[logicexons==1],"#F219ED", "#CCCCCC");
               
               drawTranscript(rel.calc.min, rel.calc.max, tr=tr, rango=1:(length(flat.gff.data$start[rt.allExon][logicexons==1])), rescale.iv=rescale.iv, exoncol=curr.exoncol, names=c(), trName=trans[i], par.cex = par.cex, anno.cex.text = anno.cex.text,sub.sig = sub.sig,anno.lwd=anno.lwd, cex.axis=anno.cex.axis, cex.main=anno.cex.main,...)
            }
         }

   if(drawCoordinates){
     if(! is.null(rescale.iv)){
       pretty.x <- pretty(c(rel.calc.min,rel.calc.max),n=5);
       pretty.x <- pretty.x[pretty.x > rel.calc.min & pretty.x < rel.calc.max];
       rescaled.pretty.x <- rescale.coords(pretty.x,rescale.iv);
       #if(any(is.na(rescaled.pretty.x))){
       #  print(pretty.x);
       #  print(rescaled.pretty.x);
       #  print(rescale.iv);
       #}
       if(min(rescaled.pretty.x) > 0.05){
         pretty.x <- c(rel.calc.min, pretty.x);
         rescaled.pretty.x <- c(0,rescaled.pretty.x);
       }
       if(max(rescaled.pretty.x) < 0.95){
         pretty.x <- c(pretty.x, rel.calc.max);
         rescaled.pretty.x <- c(rescaled.pretty.x,1);
       }
       rescaled.pretty.x <- rescaled.pretty.x * (rel.calc.max-rel.calc.min) + rel.calc.min;
     } else {
       pretty.x <- pretty(c(rel.calc.min,rel.calc.max),n=5);
     }
     #axis(2,   at=logticks ,labels=FALSE, las=2, pos=0, tcl=0.25, cex.axis = anno.cex.text, ...)
     axis(1, at = rescaled.pretty.x,labels=pretty.x,cex.axis=anno.cex.text, ...);
   }

if(debug.mode & verbose) message("> plotJunctionSeqResultsForGene(): "," Done.");
      
      #Depreciated?
      #axis1.main <- make.evenly.spaced.seq(rel.calc.min, rel.calc.max, 10);
      #axis1.minor <- make.evenly.spaced.seq.minor(rel.calc.min, rel.calc.max, 10);
      
      #axis(1,at=axis1.minor,labels=rep("",length(axis1.minor)),pos=0,lwd.ticks=0.2,padj=-0.7,tcl=-0.25, lwd = axes.lwd, cex = anno.cex.text,cex.axis=cex.axis, cex.main=cex.main,...);
      #axis(1,at=axis1.main,labels=axis1.main,pos=0,lwd.ticks=0.2,padj=-0.7,tcl=-1, lwd = axes.lwd, cex = anno.cex.text,cex.axis=cex.axis, cex.main=cex.main,...);
}


