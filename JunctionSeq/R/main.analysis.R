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

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################

estimateSizeFactors <- function(...){
  estimateJunctionSeqSizeFactors(...)
}

estimateJunctionSeqSizeFactors <- function( jscs , method.sizeFactors = c("byGenes","byCountbins"), replicateDEXSeqBehavior.useRawBaseMean = FALSE, calcAltSF = TRUE, verbose = FALSE){
   #Updated version, uses gene-level counts rather than sub-feature counts:
   #  This is preferred, since using the sub-feature counts in the normalization
   #  effectively over-weights genes with numerous annotated sub-features.
   #In practice, the results are almost always functionally-identical in most datasets.
   stopifnot( is( jscs, "JunctionSeqCountSet") )
   

   method.sizeFactors <- match.arg(method.sizeFactors)
   counts.byGenes <- jscs@geneCountData
   counts.byCountbins <- counts(jscs)
   
   attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.sizeFactors = method.sizeFactors)
   attr(jscs,"CallStack") <- c(attr(jscs,"CallStack"), list(deparse(match.call())))
   
   sizeFactors.byGenes <- estimateSizeFactorsForMatrix(counts.byGenes)
   sizeFactors.byCountbins <- estimateSizeFactorsForMatrix(counts.byCountbins)
   
   if(method.sizeFactors == "byGenes"){
     sizeFactors(jscs) <- sizeFactors.byGenes
   } else {
     sizeFactors(jscs) <- sizeFactors.byCountbins
   }
   
   sizeFactors(jscs@DESeqDataSet) <- rep(sizeFactors(jscs),2)
   
   if(replicateDEXSeqBehavior.useRawBaseMean){
     warning("Option replicateDEXSeqBehavior.useRawBaseMean replicates a bug in an old version of DEXSeq. It is only intended for testing purposes. NOT FOR GENERAL USE!")
   }
   
   fData(jscs)$baseMean <- rowMeans(counts(jscs, normalized= ! replicateDEXSeqBehavior.useRawBaseMean))
   fData(jscs)$baseVar <- rowVars(counts(jscs, normalized=! replicateDEXSeqBehavior.useRawBaseMean))
   
   if(calcAltSF){
       altSF <- data.frame(sizeFactors.byGenes = sizeFactors.byGenes, sizeFactors.byCountbins = sizeFactors.byCountbins)
       rownames(altSF) <- pData(jscs)$sample.ID
       jscs@altSizeFactors <- altSF
   } else {
       jscs@altSizeFactors <- data.frame()
   }
   
   return(jscs);                         
}



##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################
fitDispersionFunction <- function(...){ #DEPRECATED!
   fitJunctionSeqDispersionFunction(...)
}

fitJunctionSeqDispersionFunction <- function(jscs , 
                                  method.GLM = c(c("advanced","DESeq2-style"), c("simpleML","DEXSeq-v1.8.0-style")),
                                  method.dispFit = c("parametric", "local", "mean"), 
                                  method.dispFinal = c("shrink","max","fitted","noShare"),
                                  fitDispersionsForExonsAndJunctionsSeparately = TRUE, 
                                  #advancedMode = TRUE, 
                                  verbose = TRUE){
   method.GLM <- match.arg(method.GLM)
   method.dispFit <- match.arg(method.dispFit)
   method.dispFinal <- match.arg(method.dispFinal)
   attr(jscs,"CallStack") <- c(attr(jscs,"CallStack"), list(deparse(match.call())))
   
   attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.GLM.DispFit = method.GLM)
   attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.dispFit = method.dispFit)
   attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.dispFinal = method.dispFinal)
   
     if(method.dispFinal != "noShare"){
       jscs <- fitDispersionFunction_advancedMode(jscs, fitType = method.dispFit, verbose = verbose, finalDispersionMethod = method.dispFinal, fitDispersionsForExonsAndJunctionsSeparately = fitDispersionsForExonsAndJunctionsSeparately)
     } else {
       fData(jscs)$dispersion <- fData(jscs)$dispBeforeSharing
     }
     mcols(jscs@DESeqDataSet)$dispersion <- fData(jscs)$dispersion
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
   attr(jscs,"CallStack") <- c(attr(jscs,"CallStack"), list(deparse(match.call())))
  
  runOnFeatures = seq_len(nrow(fData(jscs)))
  names(runOnFeatures) <- rownames(fData(jscs))
  method.expressionEstimation <- match.arg(method.expressionEstimation)
  fitExpToVar <- "condition"
  varlist <- rownames(attr(terms(effect.formula), "factors"))
  covarlist <- varlist[ ! varlist %in% c(fitExpToVar, "countbin") ]
    
  jscs@formulas[["effect.formula"]] <- effect.formula
  attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.expressionEstimation = method.expressionEstimation)
  
  myApply <- getMyApply(nCores)
  
  modelFrame <- constructModelFrame( jscs )
  modelFrame$countbin <- factor(modelFrame$countbin, levels = c("others","this"))
  mm <- rmDepCols( model.matrix( effect.formula, modelFrame ) )
  conditionLevels <- levels(modelFrame[[fitExpToVar]])
  conditionCt <- length(conditionLevels)
  
  #To implement later: progress function!
  #prog.fcn <- make.progress.report.fcn(length(runOnFeatures), 20, "-------> estimateEffectSizes: (Calculating effect size and predicted values for feature ")
  #if(method.GLM %in% c("simpleML","DEXSeq-v1.8.0-style")){
  if(method.expressionEstimation == "feature-vs-gene"){
      mdl.out <- myApply(runOnFeatures, function(i){
        if( verbose && i %% 1000 == 0 ){
           message(paste0("-------> estimateEffectSizes: (Calculating effect size and predicted values for feature ",i," of ",nrow(jscs),")","(",date(),")"))
        }
        geneID <- fData(jscs)$geneID[i]
        countbinID <- fData(jscs)$countbinID[i]
        countVector <- jscs@countVectors[i,]

        if(! fData(jscs)$testable[i]){
          return(list(
            logFCs = rep(NA, conditionCt - 1),
            logFCvst = rep(NA, conditionCt - 1),
            fitEffect = NA,
            exprEstimate = rep(NA,conditionCt),
            otherExprEstimate = rep(NA,conditionCt),
            relExprEstimate = rep(NA,conditionCt)
          ))
        } else {
          disp <- fData(jscs)[i, dispColumn]
          fit <- try( {
            glmnb.fit( mm,  countVector, dispersion = disp, offset = log( modelFrame$sizeFactor ) )
          })
          if( any(inherits( fit, "try-error" ) )) {
            warning( sprintf("glmnb.fit failed for %s:%s\n", as.character( geneID ), countbinID) )
            return(list(fitConverged = FALSE, fit = NULL))

            return(list(
              logFCs = rep(NA, conditionCt - 1),
              logFCvst = rep(NA, conditionCt - 1),
              fitEffect = NA,
              exprEstimate = rep(NA,conditionCt),
              otherExprEstimate = rep(NA,conditionCt),
              relExprEstimate = rep(NA,conditionCt)
            ))
          }

          coefs <- arrangeCoefs( effect.formula, modelFrame, mm, fit = fit, insertValues = TRUE )

          predictedEstimates <- getPredictedEstimates(coefs = coefs, 
                                                      forVarName = fitExpToVar, 
                                                      forVarValues = conditionLevels, 
                                                      selectVarName = "countbin", 
                                                      selectVarValue = "this", 
                                                      averageVarNames = covarlist)

          predictedEstimatesGene <- getPredictedEstimates(coefs = coefs, 
                                                      forVarName = fitExpToVar, 
                                                      forVarValues = conditionLevels, 
                                                      selectVarName = "countbin", 
                                                      selectVarValue = "others", 
                                                      averageVarNames = covarlist)

          countBinThisBeta <- predictedEstimates[[1]][["countbin"]]
          countBinOthersBeta <- predictedEstimatesGene[[1]][["countbin"]]
          interceptBeta <- predictedEstimates[[1]][["(Intercept).(Intercept)"]]

          relativeEstimates <- predictedEstimates
          for(j in seq_len(length(predictedEstimates))){
            relativeEstimates[[j]] <- predictedEstimates[[j]] - predictedEstimatesGene[[j]]
          }
          relativeLogExprDiff <- sapply(relativeEstimates, sum)
          relativeLogExprEstimate <- sapply(relativeEstimates, sum) - countBinThisBeta + countBinOthersBeta
          logFCs <- relativeLogExprDiff[-1] - relativeLogExprDiff[1]
          logFCvst <- log(vst(exp(relativeLogExprDiff[-1]), jscs) / 
                          vst(exp(relativeLogExprDiff[1]),jscs)
                          )

          exprEstimate <- exp(sapply(predictedEstimates, sum))
          otherExprEstimate <- exp(sapply(predictedEstimatesGene, sum))
          relExprEstimate <- exp((relativeLogExprDiff - mean(relativeLogExprDiff)) + mean(log(exprEstimate)))
          if(! keep.estimation.fit) fit <- "FIT_NOT_SAVED"
          return(list(
            logFCs = logFCs,
            logFCvst = logFCvst,
            fitEffect = fit,
            exprEstimate = exprEstimate,
            otherExprEstimate = otherExprEstimate,
            relExprEstimate = relExprEstimate
          ))
        }
      })
  } else {  
       mdl.out <- estimateEffectSizes.BigModel(jscs, 
                            conditionCt = conditionCt, conditionLevels = conditionLevels,
                            effect.formula = effect.formula,
                            nCores=nCores, 
                            dispColumn=dispColumn,
                            keep.estimation.fit = keep.estimation.fit,
                            verbose = verbose)
  }
  #if(verbose) message("   Note: dim(mdl.out) = [",paste0(dim(mdl.out),collapse=","),"]");
  if(verbose) message("   Note: length(mdl.out) = ",length(mdl.out));
  
  logFCs <- do.call( rbind.data.frame , lapply(mdl.out, "[[", "logFCs") )
  logFCvst <- do.call( rbind.data.frame , lapply(mdl.out, "[[", "logFCvst") )
  #Change the base to base-2:
  logFCs <- apply(logFCs, c(1,2), function(x){x / log(2)})
  logFCvst <- apply(logFCvst, c(1,2), function(x){x / log(2)})

  colnames(logFCs) <-   paste0("log2FC(",   conditionLevels[-1],"/",conditionLevels[1],")")
  colnames(logFCvst) <- paste0("log2FCvst(",conditionLevels[-1],"/",conditionLevels[1],")")
  for(colname in colnames(logFCs)){
    fData(jscs)[[colname]] <- logFCs[,colname]
  }
  for(colname in colnames(logFCvst)){
    fData(jscs)[[colname]] <- logFCvst[,colname]
  }

  modelFitForEffectSize <- lapply(mdl.out, "[[", "fitEffect")
  jscs@modelFitForEffectSize <- modelFitForEffectSize

  exprEstimate         <- extractPlottingEstimates_helper(jscs, mdl.out, "exprEstimate",         paste0("expr","_",conditionLevels))
  otherExprEstimate    <- extractPlottingEstimates_helper(jscs, mdl.out, "otherExprEstimate",    paste0("geneExpr","_",conditionLevels))
  relExprEstimate      <- extractPlottingEstimates_helper(jscs, mdl.out, "relExprEstimate",      paste0("relExpr","_",conditionLevels))
  sampleNames <- colnames(counts(jscs))
  normCounts <- as.matrix( t( t(jscs@countVectors) / rep(sizeFactors(jscs),2) ) )
  rownames(normCounts) <- featureNames( jscs )
  colnames(normCounts) <-  c(paste0("normCount_",sampleNames), paste0("normGeneCount_",sampleNames))
  jscs@plottingEstimates    <- list(exprEstimate = exprEstimate,
                                    geneExprEstimate = otherExprEstimate,
                                    relExprEstimate = relExprEstimate, 
                                    normCounts = normCounts)
  if(calculate.geneLevel.expression){
    if(verbose) message(paste0("-------> estimateEffectSizes: Estimating gene-level expression."))
    modelFrame <- cbind(
               sample = sampleNames(jscs),
               design(jscs, drop=FALSE),
               sizeFactor = sizeFactors(jscs)
               )
    mm <- rmDepCols( model.matrix( geneLevel.formula, modelFrame ) )
    
    geneLevelEstModeled <- do.call(rbind.data.frame,myApply(seq_len(nrow(jscs@geneCountData)), function(i){
      if( verbose && i %% 100 == 0 ){
         message(paste0("-------> estimateEffectSizes: (Calculating gene-level effect size and predicted values for gene ",i," of ",nrow(jscs@geneCountData),")","(",date(),")"))
      }
      geneID <- rownames(jscs@geneCountData)[i]
      countVector <- jscs@geneCountData[i,]
      #Old version, only works for parametric fits
      #fitted.dispersion <- jscs@dispFitCoefs[1] + jscs@dispFitCoefs[2] / mean(countVector / modelFrame$sizeFactor)
      fitted.dispersion <- jscs@dispFunction( mean( countVector / modelFrame$sizeFactor ) ) 
      fit <- try( {
        glmnb.fit( mm,  countVector, dispersion = fitted.dispersion, offset = log( modelFrame$sizeFactor ) )
      })
      if( any(inherits( fit, "try-error" ) )) {
        warning( sprintf("glmnb.fit failed for %s:\n", as.character( geneID )) )
        return(rep(NA,length(conditionLevels)))
      } else {
        coefs <- arrangeCoefs( geneLevel.formula, modelFrame, mm, fit = fit, insertValues = TRUE )
        predictedEstimatesGene <- getPredictedEstimatesGeneLevel(coefs = coefs, 
                                                  forVarName = fitExpToVar, 
                                                  forVarValues = conditionLevels,
                                                  averageVarNames = covarlist)
        exp(sapply(predictedEstimatesGene, sum))
      }
    }))
    colnames(geneLevelEstModeled) <- paste0("geneLevelEst_",conditionLevels)
    rownames(geneLevelEstModeled) <- rownames(jscs@geneCountData)
    
    jscs@geneLevelPlottingEstimates <- list(geneLevelEstModeled = geneLevelEstModeled)
    
  }
  en <- colnames(jscs@plottingEstimates[["exprEstimate"]])
  fData(jscs)[,en] <- jscs@plottingEstimates[["exprEstimate"]][,en]
  
  tryCatch({
      if(verbose) message("-------> estimateEffectSizes: Starting gene-wise p-adjust. (",date(),")")
      genewise.sig <- JS.perGeneQValue(pvals = fData(jscs)$pvalue, wTest = fData(jscs)$testable, fData(jscs)$geneID)
      fData(jscs)$geneWisePadj <- sapply(as.character(fData(jscs)$geneID), function(g){ if(any(g == names(genewise.sig))){ genewise.sig[[g]]; } else { NA;} })
      if(verbose) message("-------> estimateEffectSizes: Finished gene-wise p-adjust. (",date(),")")
  }, error = function(e){
    message("WARNING: Failed gene-level padjust.")
    print(e)
  })
  
  return(jscs)
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
   attr(jscs,"CallStack") <- c(attr(jscs,"CallStack"), list(deparse(match.call())))
   myApply <- getMyApply(nCores)
   gene.list <- unique(jscs@featureData$geneID)
   gene.idx <- seq_len(length(gene.list))
   names(gene.idx) <- gene.list
   
   if(verbose){
      message(paste("Starting expression estimate generation on",length(gene.list),"top-level features (eg genes)."))
   }
   get.one.genes.data <- function(geneIndex){
      if( verbose && geneIndex %% 1000 == 0 ){
         message(paste0("-------> estimateEffectSizes.BigModel: (Calculating effect size and predicted values for gene ",geneIndex," of ",length(gene.list),")","(",date(),")"))
      }
      geneID <- gene.list[geneIndex]
      rowIdx <- which( fData(jscs)$geneID == geneID )
      mdl.data <- fitAndArrangeCoefs( jscs, geneID, frm = effect.formula , set.na.dispersions = 1e-8, tol = 0.1)

      if(! is.null(mdl.data)) {
          fit <- mdl.data[["fit"]]
          coefs <- mdl.data[["coefs"]]
          #If the model converged properly, then output the results:
          exprEstimate <-  exp(as.matrix( t( getEffectsForPlotting(coefs, averageOutExpression=FALSE, groupingVar="condition") ) ))
          rExprEstimate <- as.matrix( t( getEffectsForPlotting(coefs, averageOutExpression=TRUE,  groupingVar="condition") ) )
          otherExprEstimate <- matrix(NA, ncol = ncol(exprEstimate), nrow = nrow(exprEstimate))
          
          logFCs <-   rExprEstimate[,-1,drop=FALSE] - rExprEstimate[,1,drop=FALSE]
          logFCvst <- log(vst(exp(rExprEstimate[,-1,drop=FALSE]), jscs) / 
                          vst(exp(rExprEstimate[, 1,drop=FALSE]),  jscs)
                          )
          rExprEstimate <- exp(rExprEstimate)
          
          colnames(logFCs) <-   paste0("log2FC(",   conditionLevels[-1],"/",conditionLevels[1],")")
          colnames(logFCvst) <- paste0("log2FCvst(",conditionLevels[-1],"/",conditionLevels[1],")")
          
          if(! keep.estimation.fit) fit <- "FIT_NOT_SAVED"
          
          return( list(logFCs = as.matrix(logFCs),
                       logFCvst = as.matrix(logFCvst),
                       fitEffect = fit,
                       exprEstimate = exprEstimate,
                       otherExprEstimate = otherExprEstimate,
                       relExprEstimate = rExprEstimate
          ))
      } else {
          #If the model does not converge, return NA for all fitted values.
          if(verbose){
            message(paste("glm fit failed for gene", geneID, " index: ",geneIndex))
          }
          logFCs <- matrix(NA, ncol = conditionCt - 1, nrow = length(rowIdx))
          logFCvst <- matrix(NA, ncol = conditionCt - 1, nrow = length(rowIdx))
          colnames(logFCs) <-   paste0("log2FC(",   conditionLevels[-1],"/",conditionLevels[1],")")
          colnames(logFCvst) <- paste0("log2FCvst(",conditionLevels[-1],"/",conditionLevels[1],")")
          
          return( list(
            logFCs = logFCs,
            logFCvst = logFCvst,
            fitEffect = NA,
            exprEstimate = matrix(NA, ncol = conditionCt, nrow = length(rowIdx)),
            otherExprEstimate = matrix(NA, ncol = conditionCt, nrow = length(rowIdx)),
            relExprEstimate = matrix(NA, ncol = conditionCt, nrow = length(rowIdx))
          ))
      }
   }
   mdl.out <- myApply(gene.idx, get.one.genes.data)
   
   failed.fit <- is.na(lapply(mdl.out,"[[","fitEffect"))
   has.testable.feature <- sapply(gene.idx, function(geneIndex){
     geneID <- gene.list[geneIndex]
     rowIdx <- which( fData(jscs)$geneID == geneID )
     any(fData(jscs)$testable[rowIdx])
   })
   if(verbose){
     message(paste("-------> estimateEffectSizes.BigModel: glm fit failed for ",sum(failed.fit), " out of ",length(failed.fit)," genes."))
     message(paste("-------> estimateEffectSizes.BigModel: glm fit failed for ",sum(failed.fit & has.testable.feature), " out of ",sum(has.testable.feature)," genes with 1 or more testable features."))
   }
   
   
   return(mdl.out)
}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


extractPlottingEstimates_helper <- function(jscs, mdl.out, estName, columnNames){
  x <- do.call( rbind , lapply(mdl.out, "[[", estName) )
  rownames(x) <- featureNames( jscs )
  colnames(x) <- columnNames
  return(x)
}

getPredictedEstimatesGeneLevel <- function(coefs, forVarName, forVarValues, averageVarNames){
  predictedEstimates <- lapply(forVarValues, function(pvn){
    pvToAdd <- sapply(coefs, function(coef){
      if("(Intercept)" %in% names(dimnames(coef))){
        return(coef[1])
      }
      coef <- applyByDimname(coef, dimname = forVarName, FUN = function(cf){
        cf[[pvn]]
      })
      for(covar in averageVarNames){
        coef <- applyByDimname(coef, dimname = covar, FUN = function(cf){
          mean(cf)
        })
      }
      as.vector(coef)
    })
    pvToAdd
  })
  names(predictedEstimates) <- forVarValues
  return(predictedEstimates)
}

getPredictedEstimates <- function(coefs, forVarName, forVarValues, selectVarName, selectVarValue, averageVarNames){
  predictedEstimates <- lapply(forVarValues, function(pvn){
    pvToAdd <- sapply(coefs, function(coef){
      if("(Intercept)" %in% names(dimnames(coef))){
        return(coef[1])
      }
      coef <- applyByDimname(coef, dimname = selectVarName, FUN = function(cf){
        cf[[selectVarValue]]
      })
      coef <- applyByDimname(coef, dimname = forVarName, FUN = function(cf){
        cf[[pvn]]
      })
      for(covar in averageVarNames){
        coef <- applyByDimname(coef, dimname = covar, FUN = function(cf){
          mean(cf)
        })
      }
      as.vector(coef)
    })
    pvToAdd
  })
  names(predictedEstimates) <- forVarValues
  return(predictedEstimates)
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
   if(verbose) message("---> STARTING estimateJunctionSeqDispersions: (v",packageVersion("JunctionSeq"),") (",date(),")")
   stopifnot( inherits( jscs, "JunctionSeqCountSet" ) )
   if( all( is.na( sizeFactors( jscs )) ) ){
     stop("Please calculate size factors before estimating dispersions\n")
   }
   method.GLM <- match.arg(method.GLM)
   fData(jscs)$status <- rep("OK",nrow(fData(jscs)))
   
   attr(jscs,"CallStack") <- c(attr(jscs,"CallStack"), list(deparse(match.call())))
   attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.GLM.DispEst = method.GLM)
   
   testable <- ! fData(jscs)$allZero
   fData(jscs)$status[! testable] <- "ALL_ZERO"
   if(meanCountTestableThreshold != "auto"){
     fData(jscs)$status[testable & fData(jscs)$baseMean <= meanCountTestableThreshold] <- "LOW_COUNT"
     testable <- testable & fData(jscs)$baseMean > meanCountTestableThreshold
     if(verbose) message("-----> ejsd: ",length(testable) - sum(testable)," counting bins dropped due to low coverage (mean normalized count < ",meanCountTestableThreshold,")")
   }
   
   fData(jscs)$testable <- testable
   if(verbose) message("-----> ejsd: ",sum(testable)," counting bins are marked 'testable'. across ",length(unique(fData(jscs)$geneID[fData(jscs)$testable]))," genes.")
   if(verbose) message("             (",sum(fData(jscs)$featureType[fData(jscs)$testable] == "exonic_part")," exonic regions, ",
                                        sum(fData(jscs)$featureType[fData(jscs)$testable] == "splice_site")," known junctions, ",
                                        sum(fData(jscs)$featureType[fData(jscs)$testable] == "novel_splice_site")," novel junctions)")
   
   jscs@formulas[["formulaDispersion"]] <- deparse(test.formula1)
   
   if(method.GLM %in% c("simpleML","DEXSeq-v1.8.0-style")){
     myApply <- getMyApply(nCores)
     rows <- seq_len(nrow(jscs))
     modelFrame <- constructModelFrame( jscs )
     mm <- rmDepCols( model.matrix( test.formula1, modelFrame ) )
     disps <- myApply( rows, function(i) {
        if( verbose & i %% 1000 == 0 ){
           message(paste0("-------> estimateJunctionSeqDispersions: (Calculating dispersion for feature ",i," of ",nrow(jscs),")","(",date(),")"))
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

     message("-----> ejsd: Calculating mu...")
     mu <- myApply( rows, function(i) {
        if( verbose && i %% 1000 == 0 ){
           message(paste0("-------> ejsd: (Calculating fitted values for feature ",i," of ",nrow(jscs),")","(",date(),")"))
        }

        if( fData(jscs)$testable[i] ) {
           disp <- fData(jscs)$dispBeforeSharing[i]
           count <- jscs@countVectors[i,]
           fit <- try( glmnb.fit( mm, count, dispersion = disp, offset = log( modelFrame$sizeFactor ) ), silent=TRUE)
           if( inherits( fit, "try-error" ) ) {
              warning( paste0( sprintf("Unable to estimate mu for %s:%s", as.character( geneIDs(jscs)[i] ), countbinIDs(jscs)[i]),
                           "\nReason: ", fit
                      ))
              rep(NA, nrow(mm))
           } else {
              fit[["fitted.values"]]
           }
        } else{
            rep(NA, nrow(mm)); 
        }
     })
     mu <- do.call(rbind, mu)
     colnames(mu) <- colnames(jscs@countVectors)
     rownames(mu) <- rownames(jscs@countVectors)
     jscs@fittedMu <- as.matrix(mu)
     
     mcols(jscs@DESeqDataSet)$dispGeneEst <- fData(jscs)$dispBeforeSharing
     mcols(jscs@DESeqDataSet)$baseMean <- fData(jscs)$baseMean
     mcols(jscs@DESeqDataSet)$baseVar <- fData(jscs)$baseVar
   } else {
     dds <- jscs@DESeqDataSet
     if(verbose) message("---------> Executing DESeq2 call: estimateUnsharedDispersions");
     dds <- estimateUnsharedDispersions( dds, formula=test.formula1, BPPARAM=getBPParam(nCores), quiet= !verbose)
     if(verbose) message("---------> Finished with DESeq2 call.");
     fData(jscs)$dispBeforeSharing <- mcols(dds)$dispGeneEst
     mcols(dds)$baseMean <- fData(jscs)$baseMean
     mcols(dds)$baseVar <- fData(jscs)$baseVar
     mcols(dds)$allZero <- fData(jscs)$allZero
     
     jscs@DESeqDataSet <- dds
     jscs@fittedMu <- assays(jscs@DESeqDataSet)[["mu"]]
   }
   if(verbose) message("-----> ejsd: Dispersion estimation failed for ",sum(is.na( fData(jscs)$dispBeforeSharing[fData(jscs)$testable] )) ," out of ",sum(fData(jscs)$testable)," 'testable' counting bins. Setting these features to be 'untestable'")
   fData(jscs)$status[is.na( fData(jscs)$dispBeforeSharing ) & fData(jscs)$testable] <- "DISPERSION_EST_FAILED"
   fData(jscs)$testable[which( is.na( fData(jscs)$dispBeforeSharing ) )] <- FALSE
   
   if(verbose) message("---> FINISHED estimateJunctionSeqDispersions (",date(),")")
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
  method.GLM <- match.arg(method.GLM)
  attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.GLM.HTEST = method.GLM)
  attr(jscs,"CallStack") <- c(attr(jscs,"CallStack"), list(deparse(match.call())))
  
  keep.debug.model.data <- TRUE
  stopifnot( inherits( jscs, "JunctionSeqCountSet" ) )
   if( all( is.na( sizeFactors( jscs )))) {
     stop("Please calculate size factors before estimating dispersions\n")
   } 
   if( all( is.na( fData(jscs)[,dispColumn] ) ) ){
     stop("Please estimate dispersions before calling this function\n")
   }
   myApply <- getMyApply(nCores)
   
   jscs@formulas[["test.formula0"]] <- deparse(test.formula0)
   jscs@formulas[["test.formula1"]] <- deparse(test.formula1)
   rows <- seq_len(nrow(jscs))
   modelFrame <- constructModelFrame( jscs )
   mm0 <- rmDepCols( model.matrix( test.formula0, modelFrame ) )
   mm1 <- rmDepCols( model.matrix( test.formula1, modelFrame ) )
   
#########################################
   fitExpToVar <- "condition"
   keepCoefs <- which(attr(mm1,"assign") == length( attr(terms(test.formula1),"term.labels"))  )
   keepCoefNames <- paste0("HtestCoef(",colnames(mm1)[keepCoefs],")")
   conditionLevels <- levels(modelFrame[[fitExpToVar]])
#########################################
   
   
   if(verbose && INTERNALDEBUGMODE) simpleReportMem()
   if(method.GLM %in% c("advanced","DESeq2-style")){
   
     if(keep.hypothesisTest.fit){
       message("Saving Hypothesis Test Model Fits is Deprecated with DESeq2-style hypothesis testing.")
       warning("Saving Hypothesis Test Model Fits is Deprecated with DESeq2-style hypothesis testing.")
     }
   
     BPPARAM <- getBPParam(nCores);
     reducedModelMatrix <- mm0
     fullModelMatrix <- mm1
     object <- jscs@DESeqDataSet
     mcols(object)$allZero <- ! fData(jscs)$testable
     
     if(verbose && INTERNALDEBUGMODE) simpleReportMem()
     splitParts <- sort(
         rep(seq_len(BPPARAM$workers), 
         length.out=nrow(object) ) )
     splitObject <- split( object, splitParts )
     
     if(verbose && INTERNALDEBUGMODE) simpleReportMem()
     message(paste0("-------> testJunctionsForDiffUsage: Starting hypothesis test iteration. ","(",date(),")"))
     splitObject <- bplapply( splitObject,
       function(x){
         x <- nbinomLRT( x, reduced = reducedModelMatrix, full=fullModelMatrix )
     }, BPPARAM=BPPARAM )
     if(verbose && INTERNALDEBUGMODE) simpleReportMem()
     message(paste0("-------> testJunctionsForDiffUsage: Finished hypothesis test iteration. ","(",date(),")"))
     
     mergeObject <- do.call(rbind, splitObject)
     if(verbose && INTERNALDEBUGMODE) simpleReportMem()
     matchedNames <- match( rownames(object), rownames(mergeObject)) 
     if(verbose && INTERNALDEBUGMODE) simpleReportMem()
     mcols(object) <- mcols( mergeObject )[matchedNames,]
     assays(object) <- assays(mergeObject[matchedNames,])
     extraAttributes <- setdiff( names( attributes(splitObject[[1]]) ),  names( attributes(object) ) )
     if(verbose && INTERNALDEBUGMODE) simpleReportMem()
     
     for( atr in extraAttributes ){
       attr( object, atr ) <- attr( splitObject[[1]], atr )
     }
     message(paste0("-------> testJunctionsForDiffUsage: Finished compiling hypothesis test results. ","(",date(),")"))
     jscs@DESeqDataSet <- object
     
     fData(jscs)$pvalue <- mcols(jscs@DESeqDataSet)$LRTPvalue
     fData(jscs)$testable <- fData(jscs)$testable & (! is.na(fData(jscs)$pvalue))
     fData(jscs)$padjust_noFilter <- rep(NA,nrow(fData(jscs)))
     fData(jscs)$padjust_noFilter[fData(jscs)$testable] <- p.adjust(fData(jscs)$pvalue[fData(jscs)$testable] , method = pAdjustMethod)
     if(! all(is.na(mcols(jscs@DESeqDataSet)$maxCooks))){
       fData(jscs)$maxCooks <- mcols(jscs@DESeqDataSet)$maxCooks
     } else {
       if(verbose) message("---> tJfDU(): No non-NA maxCooks values. Ignoring cooks.")
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
                        filterThreshold =meanCountTestableThreshold,
                        verbose = verbose)
        fData(jscs)$padjust <- cfp[["res"]]$padjust
        fData(jscs)$status <- cfp[["res"]]$status
        fData(jscs)$testable <- cfp[["res"]]$testable
        attr(jscs, "filterThreshold") <- cfp[["filterThreshold"]]
        attr(jscs, "filterNumRej") <- cfp[["filterNumRej"]]

     if(verbose & INTERNALDEBUGMODE) simpleReportMem()

     return(jscs)
   } else if(method.GLM %in% c("simpleML","DEXSeq-v1.8.0-style")){
     mdl.out <- myApply( rows,
       function(i) {
         if( verbose && i %% 1000 == 0 ){
           message(paste0("-------> testJunctionsForDiffUsage: (testing for DJU on feature ",i," of ",nrow(jscs),")","(",date(),")"))
         }

         if( fData(jscs)$testable[i] ) {
            out <- testFeatureForDJU.fromRow(test.formula1, jscs, i, modelFrame, mm0, mm1, fData(jscs)[i, dispColumn] , keepCoefs = keepCoefs)
            if(! keep.hypothesisTest.fit) out[["fit"]] <- list(fitH0 = "FIT_NOT_SAVED", fitH1 = "FIT_NOT_SAVED") 
            out
         } else {
           return(list(coefficient = rep(NA,length(keepCoefs)), 
                       logFC = rep(NA, length(conditionLevels) - 1),
                       pval = NA, 
                       disp = NA, 
                       countVector = jscs@countVectors[i,], 
                       fit = list(fitH0 = NA, fitH1 = NA)
                  ))
         }
       })

      pvals <- sapply(mdl.out, "[[", "pval")
      coefficient <- do.call(rbind.data.frame, lapply(mdl.out,"[[", "coefficient"))
      logFC <- do.call(rbind.data.frame, lapply(mdl.out,"[[", "logFC"))

      modelFits <- lapply(mdl.out, "[[", "fit"); #Currently nonfunctional
      names(modelFits) <- featureNames( jscs )
      jscs@modelFitForHypothesisTest <- modelFits

      message("dim(coefficient) = ",paste0(dim(coefficient),collapse=","))
      message("dim(logFC) = ",paste0(dim(logFC),collapse=","))

      colnames(coefficient) <- keepCoefNames
      conditionLevels <- levels(pData(jscs)[["condition"]])
      colnames(logFC) <- paste0("logFC(",conditionLevels[-1],"/",conditionLevels[1],")")

      for(CN in colnames(coefficient)){
        fData(jscs)[[CN]] <- coefficient[[CN]]
      }

      fData(jscs)$pvalue <- pvals
      fData(jscs)$simple_padjust <- p.adjust( fData(jscs)$pvalue, method=pAdjustMethod )
      
      fData(jscs)$padjust <- p.adjust( fData(jscs)$pvalue, method=pAdjustMethod )
      
      if(meanCountTestableThreshold == "auto"){
         warning("Automatic threshold selection is NOT compatible with simpleML mode!")
         meanCountTestableThreshold <- -1
      }
      attr(jscs, "filterThreshold") <- meanCountTestableThreshold
      attr(jscs, "filterNumRej") <- sum( fData(jscs)$padjust < optimizeFilteringForAlpha, na.rm = TRUE)

      return(jscs)
    }
}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################

#Work-in-progress, currently nonfunctional:
testForDiffUsage.simpleNormDist <- function( jscs,
                                test.formula0 = formula(~ sample + countbin), 
                                test.formula1 = formula(~ sample + countbin + condition : countbin),
                                #method.GLM = c(c("advanced","DESeq2-style"), c("simpleML","DEXSeq-v1.8.0-style")),
                                dispColumn="dispersion", nCores=1 , 
                                keep.hypothesisTest.fit = FALSE,
                                meanCountTestableThreshold = "auto",
                                optimizeFilteringForAlpha = 0.01,
                                method.cooksFilter = TRUE, cooksCutoff,
                                pAdjustMethod = "BH",
                                verbose = TRUE){
  #method.GLM <- match.arg(method.GLM)
  #attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.GLM.HTEST = method.GLM)
  #attr(jscs,"CallStack") <- c(attr(jscs,"CallStack"), list(deparse(match.call())))
  
  keep.debug.model.data <- TRUE
  stopifnot( inherits( jscs, "JunctionSeqCountSet" ) )
   if( all( is.na( sizeFactors( jscs )))) {
     stop("Please calculate size factors before estimating dispersions\n")
   } 
   #if( all( is.na( fData(jscs)[,dispColumn] ) ) ){
   #  stop("Please estimate dispersions before calling this function\n")
   #}
   myApply <- getMyApply(nCores)
   
   jscs@formulas[["test.formula0"]] <- deparse(test.formula0)
   jscs@formulas[["test.formula1"]] <- deparse(test.formula1)
   rows <- seq_len(nrow(jscs))
   modelFrame <- constructModelFrame( jscs )
   mm0 <- rmDepCols( model.matrix( test.formula0, modelFrame ) )
   mm1 <- rmDepCols( model.matrix( test.formula1, modelFrame ) )
   
#########################################
   fitExpToVar <- "condition"
   keepCoefs <- which(attr(mm1,"assign") == length( attr(terms(test.formula1),"term.labels"))  )
   keepCoefNames <- paste0("HtestCoef(",colnames(mm1)[keepCoefs],")")
   conditionLevels <- levels(modelFrame[[fitExpToVar]])
#########################################
   
   
   if(verbose && INTERNALDEBUGMODE) simpleReportMem()
     mdl.out <- myApply( rows,
       function(i) {
         if( verbose && i %% 1000 == 0 ){
           message(paste0("-------> testJunctionsForDiffUsage: (testing for DJU on feature ",i," of ",nrow(jscs),")","(",date(),")"))
         }

         if( fData(jscs)$testable[i] ) {
            out <- testFeatureForDJU.fromRow.simpleNormDist(test.formula1, jscs, i=i, modelFrame=modelFrame, mm0=mm0, mm1=mm1, keepCoefs = keepCoefs)
            if(! keep.hypothesisTest.fit) out[["fit"]] <- list(fitH0 = "FIT_NOT_SAVED", fitH1 = "FIT_NOT_SAVED") 
            out
         } else {
           return(list(coefficient = rep(NA,length(keepCoefs)), 
                       logFC = rep(NA, length(conditionLevels) - 1),
                       pval = NA, 
                       disp = NA, 
                       countVector = jscs@countVectors[i,], 
                       fit = list(fitH0 = NA, fitH1 = NA)
                  ))
         }
       })

      pvals <- sapply(mdl.out, "[[", "pval")
      coefficient <- do.call(rbind.data.frame, lapply(mdl.out,"[[", "coefficient"))
      logFC <- do.call(rbind.data.frame, lapply(mdl.out,"[[", "logFC"))

      modelFits <- lapply(mdl.out, "[[", "fit"); #Usually nonfunctional, unless keep.hypothesisTest.fit is TRUE
      names(modelFits) <- featureNames( jscs )
      jscs@modelFitForHypothesisTest <- modelFits

      message("dim(coefficient) = ",paste0(dim(coefficient),collapse=","))
      message("dim(logFC) = ",paste0(dim(logFC),collapse=","))

      colnames(coefficient) <- keepCoefNames
      conditionLevels <- levels(pData(jscs)[["condition"]])
      colnames(logFC) <- paste0("logFC(",conditionLevels[-1],"/",conditionLevels[1],")")

      for(CN in colnames(coefficient)){
        fData(jscs)[[CN]] <- coefficient[[CN]]
      }

      fData(jscs)$pvalue <- pvals
      fData(jscs)$simple_padjust <- p.adjust( fData(jscs)$pvalue, method=pAdjustMethod )
      
      fData(jscs)$padjust <- p.adjust( fData(jscs)$pvalue, method=pAdjustMethod )
      
      if(meanCountTestableThreshold == "auto"){
         warning("Automatic threshold selection is NOT compatible with simpleML mode!")
         meanCountTestableThreshold <- -1
      }
      attr(jscs, "filterThreshold") <- meanCountTestableThreshold
      attr(jscs, "filterNumRej") <- sum( fData(jscs)$padjust < optimizeFilteringForAlpha, na.rm = TRUE)

      return(jscs)
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
        mm <- rmDepCols( model.matrix( formula(jscs@formulas[["formulaDispersion"]]), modelFrame ) )
        
        attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.cooksFilter = method.cooksFilter)
        
        testable <- fData(jscs)$testable
        baseMean <- fData(jscs)$baseMean
        if(meanCountTestableThreshold != "auto"){
          if(verbose) message(">     (Explicit baseMean filter. Removing ",sum(testable & (! f.na(baseMean > meanCountTestableThreshold)))," features with baseMean <= ",meanCountTestableThreshold,")")
          testable <- testable & f.na(baseMean > meanCountTestableThreshold)
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
                        verbose = verbose)
        return(cfp)
        
}
