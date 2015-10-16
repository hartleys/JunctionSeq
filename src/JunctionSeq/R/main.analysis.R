
##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


estimateSizeFactors <- function( jscs , method.sizeFactors = c("byGenes","byCountbins"), replicateDEXSeqBehavior.useRawBaseMean = FALSE, calcAltSF = TRUE, verbose = FALSE){
   #Updated version, uses gene-level counts rather than sub-feature counts:
   #  This is preferred, since using the sub-feature counts in the normalization
   #  effectively over-weights genes with numerous annotated sub-features.
   #In practice, the results are almost always functionally-identical in most datasets.
   stopifnot( is( jscs, "JunctionSeqCountSet") );
   
   method.sizeFactors <- match.arg(method.sizeFactors);
   counts.byGenes <- jscs@geneCountData;
   counts.byCountbins <- counts(jscs);
   
   attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.sizeFactors = method.sizeFactors);
   
   sizeFactors.byGenes <- estimateSizeFactorsForMatrix(counts.byGenes);
   sizeFactors.byCountbins <- estimateSizeFactorsForMatrix(counts.byCountbins);
   
   if(method.sizeFactors == "byGenes"){
     sizeFactors(jscs) <- sizeFactors.byGenes;
   } else {
     sizeFactors(jscs) <- sizeFactors.byCountbins;
   }
   
   sizeFactors(jscs@DESeqDataSet) <- rep(sizeFactors(jscs),2);
   
   fData(jscs)$baseMean <- rowMeans(counts(jscs, normalized= ! replicateDEXSeqBehavior.useRawBaseMean));
   fData(jscs)$baseVar <- rowVars(counts(jscs, normalized=! replicateDEXSeqBehavior.useRawBaseMean));
   
   if(calcAltSF){
       altSF <- data.frame(sizeFactors.byGenes = sizeFactors.byGenes, sizeFactors.byCountbins = sizeFactors.byCountbins);
       rownames(altSF) <- pData(jscs)$sample.ID;
       jscs@altSizeFactors <- altSF;
   } else {
       jscs@altSizeFactors <- data.frame();
   }
   
   return(jscs);
   #Old version, uses feature counts.
   #stopifnot( is( jscs, "JunctionSeqCountSet") );
   #geomeans <- exp( rowMeans( log( counts(jscs) ) ) );
   #sizeFactors(jscs) <- apply( counts(jscs), 2, function(cnts){
   #   median( ( cnts / geomeans )[ geomeans>0 ] )
   #});
   #return(jscs);                                                      
}



##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################

fitDispersionFunction <- function(jscs , 
                                  method.GLM = c(c("advanced","DESeq2-style"), c("simpleML","DEXSeq-v1.8.0-style")),
                                  method.dispFit = c("parametric", "local", "mean"), 
                                  method.dispFinal = c("shrink","max","fitted","noShare"),
                                  fitDispersionsForExonsAndJunctionsSeparately = TRUE, 
                                  #advancedMode = TRUE, 
                                  verbose = TRUE){
   method.GLM <- match.arg(method.GLM);
   method.dispFit <- match.arg(method.dispFit);
   method.dispFinal <- match.arg(method.dispFinal);
   
   attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.GLM = method.GLM);
   attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.dispFit = method.dispFit);
   attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.dispFinal = method.dispFinal);
   
     if(method.dispFinal != "noShare"){
       jscs <- fitDispersionFunction_advancedMode(jscs, fitType = method.dispFit, verbose = verbose, finalDispersionMethod = method.dispFinal, fitDispersionsForExonsAndJunctionsSeparately = fitDispersionsForExonsAndJunctionsSeparately);
     } else {
       fData(jscs)$dispersion <- fData(jscs)$dispBeforeSharing;
     }
     mcols(jscs@DESeqDataSet)$dispersion <- fData(jscs)$dispersion;
     return(jscs)
   

}



##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


estimateEffectSizes <- function(jscs, 
                        method.expressionEstimation = c("feature-vs-gene","feature-vs-otherFeatures"),
                        effect.formula = formula(~ condition + countbin + condition : countbin),
                        geneLevel.formula = formula(~ condition), calculate.geneLevel.expression = TRUE,
                        keep.estimation.fit = FALSE,
                        nCores=1, 
                        dispColumn="dispersion",
                        verbose = TRUE){
  stopifnot(is(jscs, "JunctionSeqCountSet"))
  runOnFeatures = seq_len(nrow(fData(jscs)));
  names(runOnFeatures) <- rownames(fData(jscs));
  method.expressionEstimation <- match.arg(method.expressionEstimation);
  #method.GLM <- match.arg(method.GLM);
  fitExpToVar <- "condition";
  varlist <- rownames(attr(terms(effect.formula), "factors"));
  covarlist <- varlist[ ! varlist %in% c(fitExpToVar, "countbin") ]
    
  jscs@formulas[["effect.formula"]] <- effect.formula;
  attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.expressionEstimation = method.expressionEstimation);
  
  myApply <- getMyApply(nCores);
  
  modelFrame <- constructModelFrame( jscs )
  modelFrame$countbin <- factor(modelFrame$countbin, levels = c("others","this"));
  mm <- rmDepCols( model.matrix( effect.formula, modelFrame ) );
  conditionLevels <- levels(modelFrame[[fitExpToVar]]);
  conditionCt <- length(conditionLevels);
  
  #To implement later: progress function!
  #prog.fcn <- make.progress.report.fcn(length(runOnFeatures), 20, "-------> estimateEffectSizes: (Calculating effect size and predicted values for feature ");
  #if(method.GLM %in% c("simpleML","DEXSeq-v1.8.0-style")){
  if(method.expressionEstimation == "feature-vs-gene"){
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
            logFCvst = rep(NA, conditionCt - 1),
            fitEffect = NA,
            exprEstimate = rep(NA,conditionCt),
            otherExprEstimate = rep(NA,conditionCt),
            relExprEstimate = rep(NA,conditionCt)
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
              logFCvst = rep(NA, conditionCt - 1),
              fitEffect = NA,
              exprEstimate = rep(NA,conditionCt),
              otherExprEstimate = rep(NA,conditionCt),
              relExprEstimate = rep(NA,conditionCt)
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

          countBinThisBeta <- predictedEstimates[[1]][["countbin"]];
          countBinOthersBeta <- predictedEstimatesGene[[1]][["countbin"]];
          interceptBeta <- predictedEstimates[[1]][["(Intercept).(Intercept)"]];

          relativeEstimates <- predictedEstimates;
          for(j in seq_len(length(predictedEstimates))){
            relativeEstimates[[j]] <- predictedEstimates[[j]] - predictedEstimatesGene[[j]];
          }
          relativeLogExprDiff <- sapply(relativeEstimates, sum);
          relativeLogExprEstimate <- sapply(relativeEstimates, sum) - countBinThisBeta + countBinOthersBeta;
          logFCs <- relativeLogExprDiff[-1] - relativeLogExprDiff[1];
          logFCvst <- log(vst(exp(relativeLogExprDiff[-1]), jscs) / 
                          vst(exp(relativeLogExprDiff[1]),jscs)
                          );

          exprEstimate <- exp(sapply(predictedEstimates, sum));
          otherExprEstimate <- exp(sapply(predictedEstimatesGene, sum));
          relExprEstimate <- exp((relativeLogExprDiff - mean(relativeLogExprDiff)) + mean(log(exprEstimate)));
          if(! keep.estimation.fit) fit <- "FIT_NOT_SAVED";
          return(list(
            logFCs = logFCs,
            logFCvst = logFCvst,
            fitEffect = fit,
            exprEstimate = exprEstimate,
            otherExprEstimate = otherExprEstimate,
            relExprEstimate = relExprEstimate
          ));
        }
      });
  } else {  
       mdl.out <- estimateEffectSizes.BigModel(jscs, 
                            conditionCt = conditionCt, conditionLevels = conditionLevels,
                            effect.formula = effect.formula,
                            nCores=nCores, 
                            dispColumn=dispColumn,
                            keep.estimation.fit = keep.estimation.fit,
                            verbose = verbose);
  }
  
  logFCs <- do.call( rbind.data.frame , lapply(mdl.out, "[[", "logFCs") );
  logFCvst <- do.call( rbind.data.frame , lapply(mdl.out, "[[", "logFCvst") );
  #Change the base to base-2:
  logFCs <- apply(logFCs, c(1,2), function(x){x / log(2)});
  logFCvst <- apply(logFCvst, c(1,2), function(x){x / log(2)});

  colnames(logFCs) <-   paste0("log2FC(",   conditionLevels[-1],"/",conditionLevels[1],")");
  colnames(logFCvst) <- paste0("log2FCvst(",conditionLevels[-1],"/",conditionLevels[1],")");
  for(colname in colnames(logFCs)){
    fData(jscs)[[colname]] <- logFCs[,colname];
  }
  for(colname in colnames(logFCvst)){
    fData(jscs)[[colname]] <- logFCvst[,colname];
  }

  modelFitForEffectSize <- lapply(mdl.out, "[[", "fitEffect");
  #names(modelFitForEffectSize) <- featureNames( jscs );
  jscs@modelFitForEffectSize <- modelFitForEffectSize;

  exprEstimate         <- extractPlottingEstimates_helper(jscs, mdl.out, "exprEstimate",         paste0("expr","_",conditionLevels));
  otherExprEstimate    <- extractPlottingEstimates_helper(jscs, mdl.out, "otherExprEstimate",    paste0("geneExpr","_",conditionLevels));
  relExprEstimate      <- extractPlottingEstimates_helper(jscs, mdl.out, "relExprEstimate",      paste0("relExpr","_",conditionLevels));
  sampleNames <- colnames(counts(jscs));
  normCounts <- as.matrix( t( t(jscs@countVectors) / rep(sizeFactors(jscs),2) ) )
  rownames(normCounts) <- featureNames( jscs );
  colnames(normCounts) <-  c(paste0("normCount_",sampleNames), paste0("normGeneCount_",sampleNames));
  jscs@plottingEstimates    <- list(exprEstimate = exprEstimate,
                                    geneExprEstimate = otherExprEstimate,
                                    relExprEstimate = relExprEstimate, 
                                    normCounts = normCounts);
  if(calculate.geneLevel.expression){
    if(verbose) message(paste0("-------> estimateEffectSizes: Estimating gene-level expression."));
    
    #if(method.expressionEstimation != "feature-vs-gene"){
    #  if(verbose) message(paste0("---------> estimateEffectSizes: Re-setting gene-level counts to equal the sum of all bins. (",date(),")"));
    #  gcd <- do.call(rbind.data.frame, lapply(1:nrow(jscs@geneCountData), function(i){
    #    
    #  }));
    #  jscs@geneCountData <- 
    #  if(verbose) message(paste0("---------> (done) (",date(),")"));
    #}
    
    modelFrame <- cbind(
               sample = sampleNames(jscs),
               design(jscs, drop=FALSE),
               sizeFactor = sizeFactors(jscs)
               )
    mm <- rmDepCols( model.matrix( geneLevel.formula, modelFrame ) );
    
    geneLevelEstModeled <- do.call(rbind.data.frame,myApply(seq_len(nrow(jscs@geneCountData)), function(i){
      if( verbose & i %% 100 == 0 ){
         message(paste0("-------> estimateEffectSizes: (Calculating gene-level effect size and predicted values for gene ",i," of ",nrow(jscs@geneCountData),")","(",date(),")"));
      }
      geneID <- rownames(jscs@geneCountData)[i];
      countVector <- jscs@geneCountData[i,];
      #Old version, only works for parametric fits
      #fitted.dispersion <- jscs@dispFitCoefs[1] + jscs@dispFitCoefs[2] / mean(countVector / modelFrame$sizeFactor);
      fitted.dispersion <- jscs@dispFunction( mean( countVector / modelFrame$sizeFactor ) ) ;
      fit <- try( {
        glmnb.fit( mm,  countVector, dispersion = fitted.dispersion, offset = log( modelFrame$sizeFactor ) );
      });
      if( any(inherits( fit, "try-error" ) )) {
        warning( sprintf("glmnb.fit failed for %s:\n", as.character( geneID )) );
        return(rep(NA,length(conditionLevels)));
      } else {
        coefs <- arrangeCoefs( geneLevel.formula, modelFrame, mm, fit = fit, insertValues = TRUE );
        predictedEstimatesGene <- getPredictedEstimatesGeneLevel(coefs = coefs, 
                                                  forVarName = fitExpToVar, 
                                                  forVarValues = conditionLevels,
                                                  averageVarNames = covarlist);
        exp(sapply(predictedEstimatesGene, sum));
      }
    }));
    colnames(geneLevelEstModeled) <- paste0("geneLevelEst_",conditionLevels);
    rownames(geneLevelEstModeled) <- rownames(jscs@geneCountData);
    
    jscs@geneLevelPlottingEstimates <- list(geneLevelEstModeled = geneLevelEstModeled);
    
    #jscs@geneLevelPlottingEstimates <- list(geneLevelEstSimple = geneLevelEstSimple, geneLevelEstAveraged = geneLevelEstAveraged, geneLevelEstModeled = geneLevelEstModeled);

    ##conditionLevels <- levels(modelFrame[[fitExpToVar]]);
    ##conditionCt <- length(conditionLevels);
    ##
    ##colDataNames <- rownames(attr(terms(geneLevel.formula), "factors"));
    ##colData <- pData(jscs)[, colnames(pData(jscs)) %in% colDataNames];
    ##dds <- DESeqDataSetFromMatrix(countData = jscs@geneCountData,
    ##                              colData = colData,
    ##                              design = geneLevel.formula);
    ##dds <- DESeq2::DESeq(dds);
    ##ddsRes <- DESeq2::results(dds, cooksCutoff = FALSE);
    ##
    ##geneLevel.exprEstimate <- 
  }
  en <- colnames(jscs@plottingEstimates[["exprEstimate"]])
  fData(jscs)[,en] <- jscs@plottingEstimates[["exprEstimate"]][,en];
  
  tryCatch({
      if(verbose) message("-------> estimateEffectSizes: Starting gene-wise p-adjust. (",date(),")");
      genewise.sig <- JS.perGeneQValue(pvals = fData(jscs)$pvalue, wTest = fData(jscs)$testable, fData(jscs)$geneID);
      fData(jscs)$geneWisePadj <- sapply(as.character(fData(jscs)$geneID), function(g){ if(any(g == names(genewise.sig))){ genewise.sig[[g]]; } else { NA;} });
      if(verbose) message("-------> estimateEffectSizes: Finished gene-wise p-adjust. (",date(),")");
  }, error = function(e){
    message("WARNING: Failed gene-level padjust.");
    print(e);
  })
  
  return(jscs);
}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


estimateEffectSizes.BigModel <- function(jscs, conditionCt, conditionLevels,
                        effect.formula = formula(~ condition + countbin + condition : countbin),
                        nCores=1, 
                        keep.estimation.fit = FALSE,
                        dispColumn="dispersion",
                        verbose = TRUE
                        ) {
   
   myApply <- getMyApply(nCores);
   gene.list <- unique(jscs@featureData$geneID);
   gene.idx <- seq_len(length(gene.list));
   names(gene.idx) <- gene.list;
   
   if(verbose){
      message(paste("Starting expression estimate generation on",length(gene.list),"top-level features (eg genes)."));
   }
   get.one.genes.data <- function(geneIndex){
      if( verbose & geneIndex %% 1000 == 0 ){
         message(paste0("-------> estimateEffectSizes.BigModel: (Calculating effect size and predicted values for gene ",geneIndex," of ",length(gene.list),")","(",date(),")"));
      }
      geneID <- gene.list[geneIndex];
      rowIdx <- which( fData(jscs)$geneID == geneID );
      mdl.data <- fitAndArrangeCoefs( jscs, geneID, frm = effect.formula , set.na.dispersions = 1e-8, tol = 0.1)

      if(! is.null(mdl.data)) {
          fit <- mdl.data[["fit"]];
          coefs <- mdl.data[["coefs"]];
          #If the model converged properly, then output the results:
          exprEstimate <-  exp(as.matrix( t( getEffectsForPlotting(coefs, averageOutExpression=FALSE, groupingVar="condition") ) ));
          rExprEstimate <- as.matrix( t( getEffectsForPlotting(coefs, averageOutExpression=TRUE,  groupingVar="condition") ) );
          otherExprEstimate <- matrix(NA, ncol = ncol(exprEstimate), nrow = nrow(exprEstimate));
          
          logFCs <-   rExprEstimate[,-1,drop=FALSE] - rExprEstimate[,1,drop=FALSE];
          logFCvst <- log(vst(exp(rExprEstimate[,-1,drop=FALSE]), jscs) / 
                          vst(exp(rExprEstimate[, 1,drop=FALSE]),  jscs)
                          );
          rExprEstimate <- exp(rExprEstimate);
          
          colnames(logFCs) <-   paste0("log2FC(",   conditionLevels[-1],"/",conditionLevels[1],")");
          colnames(logFCvst) <- paste0("log2FCvst(",conditionLevels[-1],"/",conditionLevels[1],")");
          
          if(! keep.estimation.fit) fit <- "FIT_NOT_SAVED";
          
          return( list(logFCs = as.matrix(logFCs),
                       logFCvst = as.matrix(logFCvst),
                       fitEffect = fit,
                       exprEstimate = exprEstimate,
                       otherExprEstimate = otherExprEstimate,
                       relExprEstimate = rExprEstimate
          ));
      } else {
          #If the model does not converge, return NA for all fitted values.
          if(verbose){
            message(paste("glm fit failed for gene", geneID, " index: ",geneIndex));
          }
          logFCs <- matrix(NA, ncol = conditionCt - 1, nrow = length(rowIdx))
          logFCvst <- matrix(NA, ncol = conditionCt - 1, nrow = length(rowIdx))
          colnames(logFCs) <-   paste0("log2FC(",   conditionLevels[-1],"/",conditionLevels[1],")");
          colnames(logFCvst) <- paste0("log2FCvst(",conditionLevels[-1],"/",conditionLevels[1],")");
          
          return( list(
            logFCs = logFCs,
            logFCvst = logFCvst,
            fitEffect = NA,
            exprEstimate = matrix(NA, ncol = conditionCt, nrow = length(rowIdx)),
            otherExprEstimate = matrix(NA, ncol = conditionCt, nrow = length(rowIdx)),
            relExprEstimate = matrix(NA, ncol = conditionCt, nrow = length(rowIdx))
          ));
      }
   }
   mdl.out <- myApply(gene.idx, get.one.genes.data);
   
   failed.fit <- is.na(lapply(mdl.out,"[[","fitEffect"));
   has.testable.feature <- sapply(gene.idx, function(geneIndex){
     geneID <- gene.list[geneIndex];
     rowIdx <- which( fData(jscs)$geneID == geneID );
     any(fData(jscs)$testable[rowIdx]);
   });
   if(verbose){
     message(paste("-------> estimateEffectSizes.BigModel: glm fit failed for ",sum(failed.fit), " out of ",length(failed.fit)," genes."));
     message(paste("-------> estimateEffectSizes.BigModel: glm fit failed for ",sum(failed.fit & has.testable.feature), " out of ",sum(has.testable.feature)," genes with 1 or more testable features."));
   }
   
   #DEBUGGING:
   #logFCs <- do.call( rbind.data.frame , lapply(mdl.out, "[[", "logFCs") );
   #message("dim(logFCs): ", paste0(dim(logFCs),collapse=","));
   #message("sum(nrow(exprEstimate)): ", sum(sapply(sapply(mdl.out,"[[","exprEstimate"), nrow
   #                                         )     )
   #);
   #write.table(logFCs,"test.txt");
   
   return(mdl.out);
}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


extractPlottingEstimates_helper <- function(jscs, mdl.out, estName, columnNames){
  x <- do.call( rbind , lapply(mdl.out, "[[", estName) );
  rownames(x) <- featureNames( jscs );
  colnames(x) <- columnNames;
  return(x);
}

getPredictedEstimatesGeneLevel <- function(coefs, forVarName, forVarValues, averageVarNames){
  predictedEstimates <- lapply(forVarValues, function(pvn){
    pvToAdd <- sapply(coefs, function(coef){
      if("(Intercept)" %in% names(dimnames(coef))){
        return(coef[1]);
      }
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

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


#setMethod("estimateDispersions", signature(object="JunctionSeqCountSet"),
estimateJunctionSeqDispersions <- function( jscs, 
                                            method.GLM = c(c("advanced","DESeq2-style"), c("simpleML","DEXSeq-v1.8.0-style")),
                                            test.formula1 = formula(~ sample + countbin + condition : countbin),
                                            meanCountTestableThreshold="auto", 
                                            nCores=1, 
                                            use.multigene.aggregates = FALSE, 
                                            verbose = TRUE){
   if(verbose) message("---> STARTING estimateJunctionSeqDispersions: (",date(),")");
   stopifnot( inherits( jscs, "JunctionSeqCountSet" ) )
   if( all( is.na( sizeFactors( jscs )) ) ){
     stop("Please calculate size factors before estimating dispersions\n")
   }
   method.GLM <- match.arg(method.GLM);
   fData(jscs)$status <- rep("OK",nrow(fData(jscs)));
   #fData(jscs)$baseMean <- rowMeans(  );
   
   attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.GLM = method.GLM);
   
   testable <- ! fData(jscs)$allZero;
   fData(jscs)$status[! testable] <- "ALL_ZERO";
   if(meanCountTestableThreshold != "auto"){
     fData(jscs)$status[testable & fData(jscs)$baseMean <= meanCountTestableThreshold] <- "LOW_COUNT";
     testable <- testable & fData(jscs)$baseMean > meanCountTestableThreshold;
     if(verbose) message("-----> ejsd: ",length(testable) - sum(testable)," counting bins dropped due to low coverage (mean normalized count < ",meanCountTestableThreshold,")");
   }
   

   #is.aggregate <- grepl("+", geneIDs(jscs), fixed=TRUE)
   #if(verbose) message("-----> ejsd: ",sum(is.aggregate)," counting bins are from ",length(unique(grep("+", geneIDs(jscs), fixed=TRUE)))," multigene aggregates (ie. overlapping genes).");
   #if(! use.multigene.aggregates){
   #  if(verbose) {
   #    message("-----> ejsd: ",sum(is.aggregate & testable)," counting bins dropped because they were from multigene aggregates.");
   #    message("             (To include these aggregate genes, set test.aggregated.genes = TRUE)");
   #  }
   #  fData(jscs)$status[ testable & is.aggregate ] <- "IS_MULTIGENE_AGGREGATE";
   #  testable <- testable & (! is.aggregate);
   #}

   #onlyOneJunctionDropCt <- 0;
   #for( r in split( seq_len(nrow(jscs)), geneIDs(jscs) ) ) {
   #  if( sum( testable[r] ) == 1 ){
   #    fData(jscs)$status[ r ] <- "ONLY_ONE_BIN_TESTABLE";
   #    testable[r] <- FALSE;
   #    onlyOneJunctionDropCt <- onlyOneJunctionDropCt + 1;
   #  }
   #}
   #message("-----> ejsd: ",sum(onlyOneJunctionDropCt)," counting bins dropped because they belong to a gene with only 1 testable bin.");
   fData(jscs)$testable <- testable;
   if(verbose) message("-----> ejsd: ",sum(testable)," counting bins are marked 'testable'. across ",length(unique(fData(jscs)$geneID[fData(jscs)$testable]))," genes.");
   if(verbose) message("             (",sum(fData(jscs)$featureType[fData(jscs)$testable] == "exonic_part")," exonic regions, ",
                                        sum(fData(jscs)$featureType[fData(jscs)$testable] == "splice_site")," known junctions, ",
                                        sum(fData(jscs)$featureType[fData(jscs)$testable] == "novel_splice_site")," novel junctions)");
   
   jscs@formulas[["formulaDispersion"]] <- deparse(test.formula1)
   
   if(method.GLM %in% c("simpleML","DEXSeq-v1.8.0-style")){
     myApply <- getMyApply(nCores);
     rows <- seq_len(nrow(jscs))
     modelFrame <- constructModelFrame( jscs )
     mm <- rmDepCols( model.matrix( test.formula1, modelFrame ) )
     disps <- myApply( rows, function(i) {
        if( verbose & i %% 1000 == 0 ){
           message(paste0("-------> estimateJunctionSeqDispersions: (Calculating dispersion for feature ",i," of ",nrow(jscs),")","(",date(),")"));
        }

        if( fData(jscs)$testable[i] ) {
           a <- try( estimateFeatureDispersionFromRow( jscs, i, modelFrame, mm ), silent=TRUE)
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

     message("-----> ejsd: Calculating mu...");
     mu <- myApply( rows, function(i) {
        if( verbose & i %% 1000 == 0 ){
           message(paste0("-------> ejsd: (Calculating fitted values for feature ",i," of ",nrow(jscs),")","(",date(),")"));
        }

        if( fData(jscs)$testable[i] ) {
           disp <- fData(jscs)$dispBeforeSharing[i];
           count <- jscs@countVectors[i,];
           fit <- try( glmnb.fit( mm, count, dispersion = disp, offset = log( modelFrame$sizeFactor ) ), silent=TRUE)
           if( inherits( fit, "try-error" ) ) {
              warning( paste0( sprintf("Unable to estimate mu for %s:%s", as.character( geneIDs(jscs)[i] ), countbinIDs(jscs)[i]),
                           "\nReason: ", fit
                      ))
              rep(NA, nrow(mm));
           } else {
              fit[["fitted.values"]];
           }
        } else{
            rep(NA, nrow(mm)); 
        }
     })
     mu <- do.call(rbind, mu);
     colnames(mu) <- colnames(jscs@countVectors);
     rownames(mu) <- rownames(jscs@countVectors);
     jscs@fittedMu <- as.matrix(mu);
     
     #mcols(jscs@DESeqDataSet)$dispGeneEst <- fData(jscs)$dispBeforeSharing;
     #mcols(jscs@DESeqDataSet)$baseMean <- rowMeans(counts(jscs, normalized=TRUE));
     #mcols(jscs@DESeqDataSet)$baseVar <- rowVars(counts(jscs, normalized=TRUE));
     mcols(jscs@DESeqDataSet)$dispGeneEst <- fData(jscs)$dispBeforeSharing;
     mcols(jscs@DESeqDataSet)$baseMean <- fData(jscs)$baseMean;
     mcols(jscs@DESeqDataSet)$baseVar <- fData(jscs)$baseVar;
   } else {
     dds <- jscs@DESeqDataSet
     dds <- estimateUnsharedDispersions( dds, formula=test.formula1, BPPARAM=MulticoreParam(workers=nCores), quiet= !verbose);
     fData(jscs)$dispBeforeSharing <- mcols(dds)$dispGeneEst;
     mcols(dds)$baseMean <- fData(jscs)$baseMean;
     mcols(dds)$baseVar <- fData(jscs)$baseVar;
     mcols(dds)$allZero <- fData(jscs)$allZero;
     
     jscs@DESeqDataSet <- dds;
     jscs@fittedMu <- assays(jscs@DESeqDataSet)[["mu"]]
   }
   if(verbose) message("-----> ejsd: Dispersion estimation failed for ",sum(is.na( fData(jscs)$dispBeforeSharing[fData(jscs)$testable] )) ," out of ",sum(fData(jscs)$testable)," 'testable' counting bins. Setting these features to be 'untestable'");
   fData(jscs)$status[is.na( fData(jscs)$dispBeforeSharing ) & fData(jscs)$testable] <- "DISPERSION_EST_FAILED";
   fData(jscs)$testable[which( is.na( fData(jscs)$dispBeforeSharing ) )] <- FALSE;
   
   if(verbose) message("---> FINISHED estimateJunctionSeqDispersions (",date(),")");
   jscs
} #)

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


testForDiffUsage <- function( jscs,
                                test.formula0 = formula(~ sample + countbin), 
                                test.formula1 = formula(~ sample + countbin + condition : countbin),
                                method.GLM = c(c("advanced","DESeq2-style"), c("simpleML","DEXSeq-v1.8.0-style")),
                                dispColumn="dispersion", nCores=1 , 
                                keep.hypothesisTest.fit = FALSE,
                                meanCountTestableThreshold = "auto",
                                optimizeFilteringForAlpha = 0.01,
                                method.cooksFilter = TRUE, cooksCutoff,
                                pAdjustMethod = "BH",
                                verbose = TRUE){
  method.GLM <- match.arg(method.GLM);
  attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.GLM = method.GLM);
  
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
   
   
   if(verbose & INTERNALDEBUGMODE) simpleReportMem();
   if(method.GLM %in% c("advanced","DESeq2-style")){
     BPPARAM <- MulticoreParam(workers = nCores);
     reducedModelMatrix <- mm0;
     fullModelMatrix <- mm1;
     object <- jscs@DESeqDataSet
     mcols(object)$allZero <- ! fData(jscs)$testable;
     
     if(verbose & INTERNALDEBUGMODE) simpleReportMem();
     splitParts <- sort(
         rep(seq_len(BPPARAM$workers), 
         length.out=nrow(object) ) )
     splitObject <- split( object, splitParts )
     
     if(verbose & INTERNALDEBUGMODE) simpleReportMem();
     message(paste0("-------> testJunctionsForDiffUsage: Starting hypothesis test iteration. ","(",date(),")"));
     splitObject <- bplapply( splitObject,
       function(x){
         x <- nbinomLRT( x, reduced = reducedModelMatrix, full=fullModelMatrix )
     }, BPPARAM=BPPARAM )
     if(verbose & INTERNALDEBUGMODE) simpleReportMem();
     message(paste0("-------> testJunctionsForDiffUsage: Finished hypothesis test iteration. ","(",date(),")"));
     
     mergeObject <- do.call(rbind, splitObject)
     #if(any(rownames(object) != rownames(mergeObject))){
     #  message(paste0("-------> testJunctionsForDiffUsage: Reordering output..."));
     if(verbose & INTERNALDEBUGMODE) simpleReportMem();
     matchedNames <- match( rownames(object), rownames(mergeObject)) 
     if(verbose & INTERNALDEBUGMODE) simpleReportMem();
     #} else {
     #  message(paste0("-------> testJunctionsForDiffUsage: Output is already ordered."));
     #  matchedNames <- 1:nrow(object);
     #}
     mcols(object) <- mcols( mergeObject )[matchedNames,]
     assays(object) <- assays(mergeObject[matchedNames,])
     extraAttributes <- setdiff( names( attributes(splitObject[[1]]) ),  names( attributes(object) ) )
     if(verbose & INTERNALDEBUGMODE) simpleReportMem();
     
     for( atr in extraAttributes ){
       attr( object, atr ) <- attr( splitObject[[1]], atr )
     }
     message(paste0("-------> testJunctionsForDiffUsage: Finished compiling hypothesis test results. ","(",date(),")"));
     jscs@DESeqDataSet <- object;
     
     fData(jscs)$pvalue <- mcols(jscs@DESeqDataSet)$LRTPvalue;
     fData(jscs)$testable <- fData(jscs)$testable & (! is.na(fData(jscs)$pvalue));
     fData(jscs)$padjust_noFilter <- rep(NA,nrow(fData(jscs)));
     fData(jscs)$padjust_noFilter[fData(jscs)$testable] <- p.adjust(fData(jscs)$pvalue[fData(jscs)$testable] , method = pAdjustMethod);
     #fData(jscs)$padjust[fData(jscs)$testable] <- p.adjust(fData(jscs)$pvalue[fData(jscs)$testable] , method = "fdr");
     if(! all(is.na(mcols(jscs@DESeqDataSet)$maxCooks))){
       fData(jscs)$maxCooks <- mcols(jscs@DESeqDataSet)$maxCooks;
     } else {
       if(verbose) message("---> tJfDU(): No non-NA maxCooks values. Ignoring cooks.");
     }
     
     

        modelFrame <- constructModelFrame( jscs )
        mm <- rmDepCols( model.matrix( test.formula1, modelFrame ) )
        cfp <- calc.filtered.adjusted.p(
                        testable = fData(jscs)$testable, 
                        status = fData(jscs)$status,
                        filter = fData(jscs)$baseMean, 
                        pvalue = fData(jscs)$pvalue, 
                        maxCooks = fData(jscs)$maxCooks,
                        alpha = optimizeFilteringForAlpha,
                        dispModelMatrix = mm,
                        cooksFilter = method.cooksFilter,
                        cooksCutoff = cooksCutoff,
                        pAdjustMethod = pAdjustMethod,
                        independentFiltering = (meanCountTestableThreshold == "auto"),
                        verbose = verbose);
        fData(jscs)$padjust <- cfp[["res"]]$padjust;
        fData(jscs)$status <- cfp[["res"]]$status;
        fData(jscs)$testable <- cfp[["res"]]$testable;
        attr(jscs, "filterThreshold") <- cfp[["filterThreshold"]];
        attr(jscs, "filterNumRej") <- cfp[["filterNumRej"]];

     if(verbose & INTERNALDEBUGMODE) simpleReportMem();

     return(jscs);
   } else if(method.GLM %in% c("simpleML","DEXSeq-v1.8.0-style")){
     mdl.out <- myApply( rows,
       function(i) {
         if( verbose & i %% 1000 == 0 ){
           message(paste0("-------> testJunctionsForDiffUsage: (testing for DJU on feature ",i," of ",nrow(jscs),")","(",date(),")"));
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
            out <- testFeatureForDJU.fromRow(test.formula1, jscs, i, modelFrame, mm0, mm1, fData(jscs)[i, dispColumn] , keepCoefs = keepCoefs);
            out[["fit"]] <- "FIT_NOT_SAVED";
            out;
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
      fData(jscs)$simple_padjust <- p.adjust( fData(jscs)$pvalue, method=pAdjustMethod )
      
      #if(meanCountTestableThreshold == "auto"){

      #}
      #   fData(jscs) <- cbind.data.frame(fData(jscs), coefficient);
      #} else {
      #   names( coefficient ) <- featureNames( jscs );
      #   fData(jscs)[[ keepCoefNames ]] <- coefficient;
      #}

      #fData(jscs)[names(pvals), "pvalue"] <- pvals;
      
      return(jscs);
    }
}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


get.filtered.padjust <- function(jscs, optimizeFilteringForAlpha = 0.01, 
                                 method.cooksFilter = TRUE,
                                 cooksCutoff,
                                 pAdjustMethod = "BH",
                                 meanCountTestableThreshold = "auto",
                                 verbose = verbose){
        modelFrame <- constructModelFrame( jscs )
        mm <- rmDepCols( model.matrix( formula(jscs@formulas[["formulaDispersion"]]), modelFrame ) );
        
        attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.cooksFilter = method.cooksFilter);
        
        testable <- fData(jscs)$testable;
        baseMean <- fData(jscs)$baseMean;
        if(meanCountTestableThreshold != "auto"){
          if(verbose) message(">     (Explicit baseMean filter. Removing ",sum(testable & (! f.na(baseMean > meanCountTestableThreshold)))," features with baseMean <= ",meanCountTestableThreshold,")");
          testable <- testable & f.na(baseMean > meanCountTestableThreshold);
        }
        
        cfp <- calc.filtered.adjusted.p(
                        testable = testable, 
                        status = fData(jscs)$status,
                        filter = baseMean, 
                        pvalue = fData(jscs)$pvalue, 
                        maxCooks = fData(jscs)$maxCooks,
                        alpha = optimizeFilteringForAlpha,
                        dispModelMatrix = mm,
                        cooksFilter = method.cooksFilter,
                        cooksCutoff = cooksCutoff,
                        pAdjustMethod = pAdjustMethod,
                        independentFiltering = (meanCountTestableThreshold == "auto"),
                        verbose = verbose);
        return(cfp);
        
}
