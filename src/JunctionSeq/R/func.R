



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
                                analysis.type = c("junctionsAndExons","junctionsOnly","exonsOnly"),
                                use.exons = NULL, use.junctions = NULL, 
                                use.novel.junctions = TRUE, 
                                meanCountTestableThreshold = 7.5,
                                nCores = 1,
                                use.covars = NULL, 
                                test.formula0 = formula(~ sample + countbin), 
                                test.formula1 = formula(~ sample + countbin + condition : countbin),
                                effect.formula = formula(~ condition + countbin + condition : countbin),
                                geneLevel.formula = formula(~ condition),
                                gzip.output = TRUE,
                                test.aggregated.genes = FALSE,
                                fitDispersionsForExonsAndJunctionsSeparately = TRUE,
                                use.alternate.method = TRUE,
                                verbose = TRUE){
  keep.debug.model.data <- TRUE;
  gtf.format <- TRUE;
  analysis.type <- match.arg(analysis.type);
  if(is.null(use.junctions) & is.null(use.exons)){
    if(analysis.type == "junctionsAndExons"){
      use.junctions <- TRUE;
      use.exons <- TRUE;
      #if(is.null(fitDispersionsForExonsAndJunctionsSeparately)){
        fitDispersionsForExonsAndJunctionsSeparately = TRUE;
      #}
    } else if(analysis.type == "junctionsOnly"){
      use.junctions <- TRUE;
      use.exons <- FALSE;
    } else if(analysis.type == "exonsOnly"){
      use.junctions <- FALSE;
      use.exons <- TRUE;
    }
  } else {
    if(is.null(use.junctions) | is.null(use.exons)){
      stop(paste0("Illegal syntax! If parameter use.junctions or use.exons are used, then BOTH must be set!\n use.junctions = '",use.junctions,"', use.exons = '",use.exons,"'"));
    }
    
    if(use.junctions & use.exons){
      analysis.type <- "junctionsAndExons";
    } else if(use.junctions & (! use.exons)){
      analysis.type <- "junctionsOnly";
    } else  if((! use.junctions) & use.exons){
      analysis.type <- "exonsOnly";
    } else {
      stop("Illegal syntax! Parameters use.exons and use.junctions cannot both be false!");
    }
  }
  
  if(verbose){
    message("> STARTING runJunctionSeqAnalyses (",date(),")");
    message(paste("> rJSA: sample.files: ", paste0(sample.files,collapse=", ")));
    message(paste("> rJSA: sample.names: ",  paste0(sample.names,collapse=", ")));
    message(paste("> rJSA: condition: ",    paste0(condition,collapse=", ")));
    message(paste("> rJSA: analysis.type: ",analysis.type));
    message(paste("> rJSA: use.junctions: ",use.junctions));
    message(paste("> rJSA: use.novel.junctions: ",use.novel.junctions));
    message(paste("> rJSA: use.exons: ",use.exons));
    message(paste("> rJSA: nCores: ",nCores));
    message(paste("> rJSA: use.covars: ",use.covars));
    message(paste("> rJSA: test.formula0: ",paste0(test.formula0,collapse=" ")));
    message(paste("> rJSA: test.formula1: ",paste0(test.formula1,collapse=" ")));
    message(paste("> rJSA: gzip.output: ",gzip.output));
    message(paste("> rJSA: use.alternate.method: ",use.alternate.method));
    message(paste("> rJSA: test.aggregated.genes: ",test.aggregated.genes ));
  }
  #require("DEXSeq");
  
  if(! is.factor(condition)){
     condition <- factor(condition, levels = sort(unique(condition)));
  }
  
  design <- data.frame(condition = condition);
  if(! is.null(use.covars)){
    message(paste0("> rJSA: using covars:"," ",date()));
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

    if(verbose) { message(paste0("> rJSA: Reading Count files..."," ",date(),".")) }
  jscs = readJunctionSeqCounts(countfiles = as.character(sample.files),
                             samplenames = sample.names,
                             design = design,
                             flat.gff.file = flat.gff.file, 
                             verbose = verbose,
                             use.junctions = use.junctions,
                             use.novel.junctions = use.novel.junctions,
                             use.exons = use.exons,
                             nCores = nCores
                             );
  #gct <- htsc[[1]];
  #jscs <- htsc[[2]];
  if(verbose) {  message(paste0("> rJSA: Count files read."," ",date(),".")) }
    
    if(verbose) { message(paste0("> rJSA: Estimating Size Factors..."," ",date(),".")) }
  jscs <- estimateSizeFactors(jscs)
    if(verbose) { message(paste0("> rJSA: Size Factors Done. Size Factors are:","."))
                  message(paste0("> rJSA: ",paste0(names(sizeFactors(jscs)),collapse=",")))
                  message(paste0("> rJSA: ",paste0(sizeFactors(jscs),collapse=",")) )
                  message(paste0("> rJSA: Estimating Dispersions..."," ",date(),"."))
                }
  #jscs <- estimateDispersions(jscs, nCores = nCores, formula=test.formula1);

  #if(use.gene.counts.variant){
    jscs <- estimateJunctionSeqDispersions(jscs, 
                                    nCores = nCores, 
                                    test.formula1=test.formula1, 
                                    meanCountTestableThreshold = meanCountTestableThreshold,
                                    use.alternate.method  =  use.alternate.method, 
                                    test.aggregated.genes = test.aggregated.genes, 
                                    verbose = verbose);
  #} else {
  #  jscs <- estimateDispersions(jscs, nCores = nCores, formula=test.formula1);
  #}

  #jscs <- estimateDispersions(jscs, nCores = nCores);
    if(verbose) { message(paste0("> rJSA: Dispersions estimated."," ",date(),".")) }
    if(verbose) { message(paste0("> rJSA: Fitting Dispersion Fcn..."," ",date(),".")) }
  jscs <- fitDispersionFunction(jscs, verbose = verbose, fitDispersionsForExonsAndJunctionsSeparately = fitDispersionsForExonsAndJunctionsSeparately)
  
    if(verbose) { message(paste0("> rJSA: Dispersions Fcn Fitted."," ",date(),".")) }
    if(verbose) { message(paste0("> rJSA: Testing for DEU..."," ",date(),".")) }
  
  #if(use.gene.counts.variant){
    jscs <- testForDiffUsage(jscs, nCores = nCores, test.formula0 = test.formula0, test.formula1 = test.formula1, use.alternate.method  =  use.alternate.method, verbose = verbose); 
  #} else {
  #  jscs <- testForDEU(jscs, nCores = nCores, formula0 = test.formula0, formula1 = test.formula1); 
  #}
  #jscs <- testForDEU(jscs, nCores = nCores, formula0 = test.formula0, formula1 = test.formula1); 
  #jscs <- testForDEU(jscs, nCores = nCores); 
  
  
    if(verbose) { message(paste0("> rJSA: DEU tests complete."," ",date(),".")) }
    if(verbose) { message(paste0("> rJSA: Estimating effect sizes using effects models..."," ",date(),".")) }
  jscs <- estimateEffectSizes( jscs , effect.formula = effect.formula, geneLevel.formula = geneLevel.formula, nCores = nCores)
    if(verbose) { message(paste0("> rJSA: Effect Sizes estimated.",".")) }
    if(verbose) { message(paste0("> rJSA: Generating results table..."," ",date(),".")) }
  
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
    if(verbose) { message(paste0("> rJSA: Saving state..."," ",date(),".")) }
    save(jscs,file=paste(outfile.prefix,"jscs.RData",sep=''));
    #save(res,file=paste(outfile.prefix,"res.RData",sep='')); 
    if(verbose) { message(paste0("> rJSA: State saved."," ",date(),".")) }
  } else {
    if(verbose) { message(paste0("> rJSA: Skipping save state because saveState = ", saveState, " and/or is.null(outfile.prefix) = ", is.null(outfile.prefix),". ",date(),".")) }
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
                        geneLevel.formula = formula(~ condition), calculate.geneLevel.expression = TRUE,
                        nCores=1, 
                        dispColumn="dispersion",
                        verbose = TRUE){
  stopifnot(is(jscs, "JunctionSeqCountSet"))
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
  
  
  #To implement later: progress function!
  #prog.fcn <- make.progress.report.fcn(length(runOnFeatures), 20, "-------> estimateEffectSizes: (Calculating effect size and predicted values for feature ");
  
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
        #exprEstimateVST = rep(NA,conditionCt),
        otherExprEstimate = rep(NA,conditionCt),
        #otherExprEstimateVST = rep(NA,conditionCt),
        relExprEstimate = rep(NA,conditionCt),
        #relExprEstimateVST = rep(NA,conditionCt),
        normCounts = countVector / modelFrame$sizeFactor#,
        #normCountsVST = vst(countVector / modelFrame$sizeFactor, jscs)
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
          #exprEstimateVST = rep(NA,conditionCt),
          otherExprEstimate = rep(NA,conditionCt),
          #otherExprEstimateVST = rep(NA,conditionCt),
          relExprEstimate = rep(NA,conditionCt),
          #relExprEstimateVST = rep(NA,conditionCt),
          normCounts = countVector / modelFrame$sizeFactor#,
          #normCountsVST = vst(countVector / modelFrame$sizeFactor, jscs)
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
        #exprEstimateVST = vst(exprEstimate, jscs),
        otherExprEstimate = otherExprEstimate,
        #otherExprEstimateVST = vst(otherExprEstimate, jscs),
        relExprEstimate = relExprEstimate,
        #relExprEstimateVST = vst(relExprEstimate, jscs),
        normCounts = countVector / modelFrame$sizeFactor#,
        #normCountsVST = vst(countVector / modelFrame$sizeFactor, jscs)
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
  #exprEstimateVST      <- extractPlottingEstimates_helper(jscs, mdl.out, "exprEstimateVST",      paste0("exprVST","_",conditionLevels));
  otherExprEstimate    <- extractPlottingEstimates_helper(jscs, mdl.out, "otherExprEstimate",    paste0("geneExpr","_",conditionLevels));
  #otherExprEstimateVST <- extractPlottingEstimates_helper(jscs, mdl.out, "otherExprEstimateVST", paste0("geneExprVST","_",conditionLevels));
  
  relExprEstimate      <- extractPlottingEstimates_helper(jscs, mdl.out, "relExprEstimate",      paste0("relExpr","_",conditionLevels));
  #relExprEstimateVST   <- extractPlottingEstimates_helper(jscs, mdl.out, "relExprEstimateVST",   paste0("relExprVST","_",conditionLevels));
  
  sampleNames <- colnames(counts(jscs));
  normCounts           <- extractPlottingEstimates_helper(jscs, mdl.out, "normCounts",   c(paste0("normCount_",sampleNames), paste0("normGeneCount_",sampleNames)));
  #normCountsVST        <- extractPlottingEstimates_helper(jscs, mdl.out, "normCountsVST",   c(paste0("normCountsVST_",sampleNames), paste0("normGeneCountsVST_",sampleNames)));
  
  jscs@plottingEstimates    <- list(exprEstimate = exprEstimate,
                                     geneExprEstimate = otherExprEstimate,
                                     relExprEstimate = relExprEstimate, 
                                     normCounts = normCounts);
  #jscs@plottingEstimatesVST <- list(exprEstimateVST = exprEstimateVST, 
  #                                   geneExprEstimateVST = otherExprEstimateVST, 
  #                                   relExprEstimateVST = relExprEstimateVST, 
  #                                   normCountsVST = normCountsVST);
  
  if(calculate.geneLevel.expression){
    ##NOTE: THIS SECTION PROBABLY NEEDS OPTIMIZATION. DO SOME BENCHMARKING.
    #
    ##geneLevelEstSimple <- do.call(cbind.data.frame, lapply(conditionLevels, function(cond){
    ##  #rowMeans(jscs@geneCountData[,pData(jscs)$condition == cond] * sizeFactors(jscs)[pData(jscs)$condition == cond]);
    ##}));
    ##colnames(geneLevelEstSimple) <- conditionLevels;
    ##rownames(geneLevelEstSimple) <- rownames(jscs@geneCountData);
    #
    ##jscs@geneLevelPlottingEstimates <- list(geneLevelEstSimple = geneLevelEstSimple);
    
    #geneLevelEstSimple <- do.call(rbind.data.frame,myApply(seq_len(nrow(jscs@geneCountData)), function(i){
    #  #if( verbose & i %% 100 == 0 ){
    #  #   message(paste0("-------> estimateEffectSizes: (Calculating gene-level effect size and predicted values for gene ",i," of ",nrow(jscs@geneCountData),")","(",date(),")"));
    #  #}
    #  geneLevelEstSimple <- sapply(conditionLevels, function(cond){
    #    mean(jscs@geneCountData[i,pData(jscs)$condition == cond] * sizeFactors(jscs)[pData(jscs)$condition == cond]);
    #  })
    #  geneLevelEstSimple;
    #}));
    #colnames(geneLevelEstSimple) <- paste0("geneLevelEst_",conditionLevels);
    #rownames(geneLevelEstSimple) <- rownames(jscs@geneCountData);
    #
    #geneLevelEstAveraged <- do.call(rbind.data.frame,myApply(seq_len(nrow(jscs@geneCountData)), function(i){
    #  #if( verbose & i %% 100 == 0 ){
    #  #   message(paste0("-------> estimateEffectSizes: (Calculating gene-level effect size and predicted values for gene ",i," of ",nrow(jscs@geneCountData),")","(",date(),")"));
    #  #}
    #  geneID <- rownames(jscs@geneCountData)[i];
    #  feature.idx <- which(fData(jscs)$geneID == geneID);
    #  X <- jscs@plottingEstimates[["exprEstimate"]][feature.idx,] + jscs@plottingEstimates[["geneExprEstimate"]][feature.idx,];
    #  geneLevelEstAveraged <- colMeans(X,na.rm=TRUE);
    #  geneLevelEstMin <- apply(X,2,min,na.rm=TRUE);
    #  geneLevelEstMax <- apply(X,2,max,na.rm=TRUE);
    #  num.est <- sum( ! is.na( X[,1] ) );
    #  
    #  #names(geneLevelEstAveraged) <- paste0("geneLevelEst_",conditionLevels);
    #  #names(geneLevelEstMin) <- paste0("geneLevelEstMin_",conditionLevels);
    #  #names(geneLevelEstMax) <- paste0("geneLevelEstMax_",conditionLevels);
    #  c(geneLevelEstAveraged,geneLevelEstMin,geneLevelEstMax, (geneLevelEstMax - geneLevelEstMin) / geneLevelEstMin, num.est);
    #}));
    #colnames(geneLevelEstAveraged) <- c(paste0("geneLevelEst_",conditionLevels), paste0("geneLevelEstMin_",conditionLevels), paste0("geneLevelEstMax_",conditionLevels),paste0("rangePct_",conditionLevels), "numAveraged");
    #rownames(geneLevelEstAveraged) <- rownames(jscs@geneCountData);
    
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
  
  return(jscs);
}

extractPlottingEstimates_helper <- function(jscs, mdl.out, estName, columnNames){
  x <- do.call( rbind.data.frame , lapply(mdl.out, "[[", estName) );
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
                            save.allGenes = TRUE, save.sigGenes = TRUE, save.fit = FALSE, save.VST = FALSE,
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
   
   if(save.VST){
     expression.data.vst <- cbind.data.frame(featureID = rownames(fData(jscs)), do.call(cbind.data.frame,lapply(jscs@plottingEstimates, function(pe){ vst(pe,jscs) }))  );
     colnames(expression.data.vst) <- c("featureID", paste0( unlist( lapply( jscs@plottingEstimates, colnames ) ), "VST" ) )
   }
   
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
     if(sum(sig.features) == 0){
       if(verbose) message("> wcr: Zero Significant Features! (at adjusted-p-value threshold ",FDR.threshold,")");
     } else {
       gene.list <- unique(as.character(fData(jscs)$geneID[sig.features]));
       if(verbose) message("> wcr: Writing results for ", length(gene.list), " genes with 1 or more significant junctions (at adjusted-p-value threshold ", FDR.threshold,")");
     
       sig.rows <- which( as.character(fData(jscs)$geneID) %in% gene.list );
       if(verbose) message("> wcr:     Found ", length(sig.rows), " counting bins belonging to those genes.");
     
     
                    write.simple.table.gz(expression.data[sig.rows,, drop=FALSE],     file=paste0(outfile.prefix,"sigGenes.expression.data.txt"),      use.gzip=gzip.output,row.names=F,col.names=T,quote=F, sep = '\t');
       if(save.VST) write.simple.table.gz(expression.data.vst[sig.rows,, drop=FALSE], file=paste0(outfile.prefix,"sigGenes.expression.data.VST.txt"),  use.gzip=gzip.output,row.names=F,col.names=T,quote=F, sep = '\t');
                    write.table.gz(fData(jscs)[sig.rows,, drop=FALSE],         file=paste0(outfile.prefix,"sigGenes.results.txt"), use.gzip=gzip.output,row.names="featureID",col.names=T,quote=F);
     }
   }
   
   if(save.bedTracks){
     if(save.allGenes){
       if(any(fData(jscs)$featureType == "splice_site" | fData(jscs)$featureType == "novel_splice_site")){
         writeExprBedTrack(paste0(outfile.prefix,"allGenes.junctionCoverage.bed.gz"), 
                         jscs = jscs, 
                         only.with.sig.gene = FALSE, plot.exons = FALSE, plot.junctions = TRUE,
                         trackLine = paste0("track name='JctExprAll' description='Splice Junction Coverage Estimates, by group' itemRgb='On' visibility=3"));
       }
       if(any(fData(jscs)$featureType == "exonic_part")){
         writeExprBedTrack(paste0(outfile.prefix,"allGenes.exonCoverage.bed.gz"), 
                         jscs = jscs, 
                         only.with.sig.gene = FALSE, plot.exons = TRUE, plot.junctions = FALSE,
                         trackLine = paste0("track name='ExonExprAll' description='Exonic Region Coverage Estimates, by group' itemRgb='On' visibility=3"));
       }
     }
     if(save.sigGenes){
       sig.features <- which(fData(jscs)$padjust < FDR.threshold);
       if(sum(sig.features) > 0){
         if(any(fData(jscs)$featureType == "splice_site" | fData(jscs)$featureType == "novel_splice_site")){
           writeExprBedTrack(paste0(outfile.prefix,"sigGenes.junctionCoverage.bed.gz"), 
                           jscs = jscs, 
                           only.with.sig.gene = TRUE, plot.exons = FALSE, plot.junctions = TRUE,
                           trackLine = paste0("track name='JctExprAll' description='Sig genes splice Junction Coverage Estimates, by group' itemRgb='On' visibility=3"));
         }
         if(any(fData(jscs)$featureType == "exonic_part")){
           writeExprBedTrack(paste0(outfile.prefix,"sigGenes.exonCoverage.bed.gz"), 
                           jscs = jscs, 
                           only.with.sig.gene = TRUE, plot.exons = TRUE, plot.junctions = FALSE,
                           trackLine = paste0("track name='ExonExprAll' description='Sig genes exonic Region Coverage Estimates, by group' itemRgb='On' visibility=3"));
         }
         writeSigBedTrack(paste0(outfile.prefix,"sigGenes.pvalues.bed.gz"), 
                          jscs = jscs,
                          trackLine = paste0("track name='JctPvals' description='Significant Splice Junctions' useScore=1 visibility=3")
                          );
       }
     }
   }
   
   if(verbose) message("> DONE writeCompleteResults (",date(),")");
   
   if(save.jscs) save( jscs, file = paste0(outfile.prefix, "jscs.RData") );
}

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
                          #flat.gff.data = NULL, flat.gff.file = NULL, 
                          gene.list = NULL, FDR.threshold = 0.01, 
                          use.plotting.device = "png", 
                          use.vst=FALSE,use.log = TRUE, truncateBelowOne = TRUE, 
                          exon.rescale.factor = 0.3,  
                          subdirectories.by.type = TRUE,
                          ma.plot=FALSE, variance.plot=FALSE,
                          with.TX=TRUE,without.TX=TRUE,
                          expr.plot=TRUE,normCounts.plot=TRUE,
                          rExpr.plot=FALSE,rawCounts.plot=FALSE,
                          colorRed.FDR.threshold = FDR.threshold, 
                          color=NULL,
                          plot.gene.level.expression = NULL,
                          plot.exon.results = NULL, plot.junction.results = NULL, plot.novel.junction.results = NULL,
                          plot.untestable.results = FALSE,
                          plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                          par.cex = 1, anno.cex.text = 1, anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                          drawCoordinates = TRUE, yAxisLabels.inExponentialForm = FALSE,
                          show.strand.arrows = 10, arrows.length = 0.125,
                          graph.margins = c(2,3,3,2),
                          base.plot.height = 12, base.plot.width = 12, 
                          base.plot.units = "in", 
                          GENE.annotation.relative.height = 0.2, TX.annotation.relative.height = 0.025, 
                          autoscale.height.to.fit.TX.annotation = TRUE,
                          autoscale.width.to.fit.bins = 35,
                          plotting.device.params = list(), 
                          number.plots = TRUE,
                          condition.legend.text = NULL, include.TX.names = TRUE, draw.start.end.sites = TRUE,
                          openPlottingDeviceFunc = NULL, closePlottingDeviceFunc = NULL,
                          verbose=TRUE, 
                          ...){
  condition <- jscs@phenoData$condition;
  flat.gff.data <- jscs@flatGffData;
  
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
    #if(is.null(flat.gff.data)){
    #  if(is.null(flat.gff.file)){
    #    stop("FATAL ERROR: buildAllPlots: either flat.gff.file or flat.gff.data must be specified!");
    #  }
    #  flat.gff.data <- readAnnotationData(flat.gff.file);
    #}

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
    
    if(subdirectories.by.type){
      if(! file.exists(outfile.prefix)){
        dir.create(outfile.prefix);
      }
      
      if(expr.plot & with.TX          & (! file.exists(paste0(outfile.prefix,"/exprTX"))))       dir.create(paste0(outfile.prefix,"/exprTX"));
      if(expr.plot & without.TX       & (! file.exists(paste0(outfile.prefix,"/expr"))))         dir.create(paste0(outfile.prefix,"/expr"));
      if(normCounts.plot & with.TX    & (! file.exists(paste0(outfile.prefix,"/normCountsTX")))) dir.create(paste0(outfile.prefix,"/normCountsTX"));
      if(normCounts.plot & without.TX & (! file.exists(paste0(outfile.prefix,"/normCounts"))))   dir.create(paste0(outfile.prefix,"/normCounts"));
      if(rExpr.plot & with.TX         & (! file.exists(paste0(outfile.prefix,"/rExprTX"))))      dir.create(paste0(outfile.prefix,"/rExprTX"));
      if(rExpr.plot & without.TX      & (! file.exists(paste0(outfile.prefix,"/rExpr"))))        dir.create(paste0(outfile.prefix,"/rExpr"));
      if(rawCounts.plot & with.TX     & (! file.exists(paste0(outfile.prefix,"/rawCountsTX"))))  dir.create(paste0(outfile.prefix,"/rawCountsTX"));
      if(rawCounts.plot & without.TX  & (! file.exists(paste0(outfile.prefix,"/rawCounts"))))    dir.create(paste0(outfile.prefix,"/rawCounts"));
    }
    
    geneNum <- 1;
    for(geneID in gene.list){
      message(paste("> buildAllPlots: starting geneID:",geneID));
      geneNum.string <- geneNum.strings[geneNum];
      
      if(subdirectories.by.type){
        outfile.prefixes <- paste0(outfile.prefix, c("/exprTX/","/expr/",
                                                    "/normCountsTX/","/normCounts/",
                                                    "/rExprTX/","/rExpr/",
                                                    "/rawCountsTX/","/rawCounts/"),geneNum.string);
      } else {
        outfile.prefixes <- paste0(outfile.prefix, geneNum.string);
      }
      
      buildAllPlotsForGene(geneID = geneID, jscs = jscs,
                            outfile.prefix = outfile.prefixes,
                            #flat.gff.data = flat.gff.data, 
                            use.plotting.device = use.plotting.device,
                            use.vst=use.vst, use.log = use.log, truncateBelowOne = truncateBelowOne,
                            exon.rescale.factor = exon.rescale.factor, plot.gene.level.expression = plot.gene.level.expression,
                            with.TX=with.TX,without.TX=without.TX,
                            expr.plot=expr.plot,normCounts.plot=normCounts.plot,
                            rExpr.plot=rExpr.plot,rawCounts.plot=rawCounts.plot,
                            colorRed.FDR.threshold = colorRed.FDR.threshold, 
                            color= color,
                            plot.exon.results = plot.exon.results, 
                            plot.junction.results = plot.junction.results,
                            plot.novel.junction.results = plot.novel.junction.results,
                            plot.untestable.results = plot.untestable.results,
                            plot.lwd = plot.lwd, axes.lwd = axes.lwd, anno.lwd = anno.lwd, 
                            drawCoordinates = drawCoordinates,
                            arrows.length = arrows.length, show.strand.arrows = show.strand.arrows,
                            graph.margins = graph.margins,
                            par.cex = par.cex,
                            anno.cex.text = anno.cex.text,
                            anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,
                            base.plot.height = base.plot.height, base.plot.width = base.plot.width,
                            base.plot.units = base.plot.units,
                            GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height,
                            autoscale.height.to.fit.TX.annotation = autoscale.height.to.fit.TX.annotation,
                            autoscale.width.to.fit.bins = autoscale.width.to.fit.bins,
                            plotting.device.params = plotting.device.params,
                            yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm,
                            condition.legend.text = condition.legend.text, include.TX.names = include.TX.names, draw.start.end.sites=draw.start.end.sites,
                            verbose=verbose, 
                            ...);

      geneNum <- geneNum + 1;
    }
    if(verbose) message("> buildAllPlots: Plotting complete.");
  }
  if(verbose) message("> buildAllPlots: Plotting and data writing complete.");
}


#pointsize = 18, res = 150
buildAllPlotsForGene <- function(geneID,jscs,
                          outfile.prefix = "./", 
                          #flat.gff.data = NULL, flat.gff.file = NULL,
                          use.plotting.device = "png", 
                          use.vst=FALSE, use.log = TRUE, truncateBelowOne = TRUE,
                          exon.rescale.factor = 0.3,   
                          with.TX=TRUE,without.TX=TRUE,
                          expr.plot=TRUE, normCounts.plot=TRUE,
                          rExpr.plot=FALSE, rawCounts.plot=FALSE,
                          colorRed.FDR.threshold = 0.01, 
                          color=NULL,
                          plot.gene.level.expression = NULL,
                          plot.exon.results = NULL, plot.junction.results = NULL, plot.novel.junction.results = NULL,
                          plot.untestable.results = FALSE,
                          plot.lwd=3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                          par.cex = 1, anno.cex.text = 1,
                          anno.cex.axis = anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                          drawCoordinates = TRUE, yAxisLabels.inExponentialForm = FALSE,
                          show.strand.arrows = 10, arrows.length = 0.125,
                          graph.margins = c(2,3,3,2),
                          base.plot.height = 12, base.plot.width = 12, 
                          base.plot.units = "in", 
                          GENE.annotation.relative.height = 0.2, TX.annotation.relative.height = 0.025,
                          autoscale.height.to.fit.TX.annotation = TRUE,
                          autoscale.width.to.fit.bins = 35,
                          plotting.device.params = list(),
                          condition.legend.text = NULL, include.TX.names = TRUE, draw.start.end.sites = TRUE,
                          openPlottingDeviceFunc = NULL, closePlottingDeviceFunc = NULL,
                          verbose=TRUE,  ...){
    message(paste("starting buildAllPlotsForGene() for geneID:",geneID));
    condition <- jscs@phenoData$condition
    flat.gff.data <- jscs@flatGffData;
    
    if(length(outfile.prefix) == 1){
      outfile.prefix <- rep(outfile.prefix, 8);
    }
    
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
    #if(is.null(flat.gff.data)){
    #   if(is.null(flat.gff.file)){
    #      stop("FATAL ERROR: buildAllPlotsForGene: either flat.gff.file or flat.gff.data must be specified!");
    #   }
    #   flat.gff.data <- readAnnotationData(flat.gff.file);
    #}
    
    transcripts <- sapply(sapply(flat.gff.data$transcripts[which(flat.gff.data$gene_id==geneID & flat.gff.data$featureType == "exonic_part")],toString), function(x){strsplit(x, "+",fixed=TRUE)})
    trans <- Reduce(union, transcripts)
    tx.ct <- length(trans);
    
    if(autoscale.height.to.fit.TX.annotation){
      GENE.annotation.height <- GENE.annotation.relative.height * 10;
      TX.annotation.height <- TX.annotation.relative.height * 10;
      withTxPlot.height.multiplier <- (10+1+GENE.annotation.height+TX.annotation.height*tx.ct) / (10+1+GENE.annotation.height);
    } else {
      withTxPlot.height.multiplier <- 1;
    }
    
    if(is.na(autoscale.width.to.fit.bins) | autoscale.width.to.fit.bins == 0){
      width.multiplier <- 1;
    } else {
      num.cols <- getColCt(geneID = geneID, merged.data = fData(jscs), 
                                plot.exon.results = plot.exon.results, plot.junction.results = plot.junction.results, plot.novel.junction.results=plot.novel.junction.results, 
                                plot.untestable.results=plot.untestable.results);
      if(num.cols > autoscale.width.to.fit.bins){
        width.multiplier <- num.cols / autoscale.width.to.fit.bins;
      } else {
        width.multiplier <- 1;
      }
    }
    
    
    if(expr.plot){
      plot.type <- "expr"
      if(with.TX){
        outfile <- paste(outfile.prefix[1],geneID,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs, colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd , par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results , plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[2],geneID,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows,...)
        closePlottingDeviceFunc();  
      }
    }
    if(normCounts.plot){
      plot.type <- "normCounts"
      if(with.TX){
        outfile <- paste(outfile.prefix[3],geneID,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[4],geneID,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows,...)
        closePlottingDeviceFunc();  
      }
    }
    if(rExpr.plot){
      plot.type <- "rExpr"
      if(with.TX){
        outfile <- paste(outfile.prefix[5],geneID,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[6],geneID,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=use.vst,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows,...)
        closePlottingDeviceFunc();  
      }
    }
    if(rawCounts.plot){
      plot.type <- "rawCounts"
      if(with.TX){
        outfile <- paste(outfile.prefix[7],geneID,"-",plot.type,"-withTx","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=withTxPlot.height.multiplier,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=FALSE,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,displayTranscripts=TRUE,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows,...)
        closePlottingDeviceFunc();  
      }
      if(without.TX){
        outfile <- paste(outfile.prefix[8],geneID,"-",plot.type,"","",sep="")
        openPlottingDeviceFunc(outfile,heightMult=1,widthMult=width.multiplier);
        plotJunctionSeqResultsForGene(geneID, jscs,  colorRed.FDR.threshold = FDR,color = color, plot.type = plot.type, use.vst=FALSE,use.log=use.log,truncateBelowOne=truncateBelowOne ,exon.rescale.factor=exon.rescale.factor,plot.lwd=plot.lwd,axes.lwd = axes.lwd, anno.lwd = anno.lwd, par.cex = par.cex, anno.cex.text = anno.cex.text, plot.exon.results = plot.exon.results , plot.junction.results = plot.junction.results , plot.novel.junction.results = plot.novel.junction.results, plot.untestable.results = plot.untestable.results,anno.cex.axis = anno.cex.axis, anno.cex.main = anno.cex.main,drawCoordinates = drawCoordinates,arrows.length = arrows.length, graph.margins = graph.margins, yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, plot.gene.level.expression = plot.gene.level.expression,condition.legend.text = condition.legend.text, include.TX.names = include.TX.names,GENE.annotation.relative.height = GENE.annotation.relative.height, TX.annotation.relative.height = TX.annotation.relative.height, draw.start.end.sites=draw.start.end.sites, show.strand.arrows = show.strand.arrows,...)
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

plotDispEsts <- function( jscs, ylim, xlim, linecol=c("#ff000080","#ff990080"),
  xlab = "mean of normalized counts", ylab = "dispersion",
  miniTicks = TRUE,
  cex = 0.45, pch = c(1,4),
  use.smoothScatter = FALSE, smooth.nbin = 512, nrpoints = 100,
  plot.exon.results = TRUE, 
  plot.junction.results = TRUE, 
  ... )
{
  px = rowMeans( counts( jscs, normalized=TRUE ) );
  pch <- ifelse(fData(jscs)$featureType == "exonic_part", pch[1], pch[2]);

  sel <- rep(TRUE, length(px));
  if(! plot.exon.results){
    sel <- sel & fData(jscs)$featureType != "exonic_part";
  }
  if(! plot.junction.results){
    sel <- sel & (fData(jscs)$featureType != "splice_site" & fData(jscs)$featureType != "novel_splice_site");
  }
  py <- fData(jscs)$dispBeforeSharing
  sel <- sel & f.na(px>0);
  sel <- sel & f.na(py>0);
  sel <- sel & (! is.na(px));
  sel <- sel & (! is.na(py));
  
  px <- px[sel]
  py <- py[sel]

  ymin <- ((min(py, na.rm=TRUE)));
  ymax <- ((max(py, na.rm=TRUE)));
  xmin <- ((min(px, na.rm=TRUE)));
  xmax <- ((max(px, na.rm=TRUE)));
  
  message("abundance ranges from ",xmin, " to ",xmax);
  message("dispersion ranges from ",ymin, " to ",ymax);
  
  
  if(missing(ylim)){
     ylim <- c(ymin,ymax);
  }
  if(missing(xlim)){
     xlim <- c(xmin,xmax);
  }
  ylim <- log10(ylim);
  xlim <- log10(xlim);
  
  decade.min <- min(floor(xlim[1]),floor(ylim[1])) - 1;
  decade.max <- max(ceiling(ylim[2]),ceiling(xlim[2])) + 1;
  decade.at <- decade.min:decade.max;
  decade.labels <- as.expression(sapply(decade.at,function(yda){
    substitute(10 ^ x, list(x = yda));
  }));
  ticks.at <- unlist(lapply(decade.at, function(D){
    log10( (10 ^ D) * 2:9 );
  }));
  
  message("ylim = c(",paste0(ylim,collapse=","),")");
    
    if(use.smoothScatter){
      smoothScatter(log10(px),log10(py),nbin = smooth.nbin, nrpoints = nrpoints,
                     xlim = xlim, ylim=ylim,xlab="",ylab="",axes=F, xaxs="i", yaxs="i", ...);
    } else {
      plot(log10(px),log10(py),
                     xlim = xlim, ylim=ylim,xlab="",ylab="",axes=F, xaxs="i", yaxs="i", pch = pch, cex = cex, ...);
    }
    #points(log10(X[is.over & (! is.sig)]),limited.log2fc[is.over & (! is.sig)], col= point.color[is.over & (! is.sig)],pch= point.pch[is.over & (! is.sig)],cex= points.cex, ...);
    box(...);
    axis(1, at = decade.at, labels = decade.labels, tcl = -0.5, las = 1, ...);
    if(miniTicks) axis(1, at = ticks.at, labels = FALSE, tcl = -0.25,  ...);
    axis(2, at = decade.at, labels = decade.labels, tcl = -0.5, las = 1, ...);
    if(miniTicks) axis(2, at = ticks.at, labels = FALSE, tcl = -0.25,  ...);
    
    log.xg = seq( decade.min, decade.max, length.out=200 );
    xg <- 10 ^ log.xg;
    
    
    if(! is.null(jscs@dispFunctionJct)){
      lines(log.xg, log10(jscs@dispFunctionExon(xg)) , col=linecol[1], lwd=1, ...);
      lines(log.xg, log10(jscs@dispFunctionJct(xg)) , col=linecol[2], lwd=1, ...);
    } else {
      lines(log.xg, log10(jscs@dispFunction(xg)) , col=linecol[1], lwd=1, ...);
    }
    #  dispFunction = "function",
    #  dispFunctionJct  = "function",
    #  dispFunctionExon = "function",
    
    #fun = function(x) { log10(jscs@dispFitCoefs[1] + jscs@dispFitCoefs[2] / x )}
    #lines( log.xg, fun(xg), col=linecol, lwd=4, ...);
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
                                  analysis.type = c("junctionsAndExons","junctionsOnly","exonsOnly"),
                                  nCores = 1,
                                  use.exons = NULL, use.junctions = NULL, 
                                  use.novel.junctions = TRUE, 
                                  verbose = TRUE)
{
   if(verbose) {
      message("-> STARTING readJunctionSeqCounts (",date(),")");
      #message("---> RJSC: countfiles: ",paste0(countfiles, collapse=","));
      #message("---> RJSC: samplenames: ",paste0(samplenames, collapse=","));
      #message("---> RJSC: flat.gff.file: ",flat.gff.file);
      #message("---> RJSC: use.junctions:",use.junctions);
      #message("---> RJSC: use.novel.junctions:",use.novel.junctions);
      #message("---> RJSC: use.exons:",use.exons);
   }
   
   analysis.type <- match.arg(analysis.type);
   if(is.null(use.junctions) & is.null(use.exons)){
    if(analysis.type == "junctionsAndExons"){
      use.junctions <- TRUE;
      use.exons <- TRUE;
    } else if(analysis.type == "junctionsOnly"){
      use.junctions <- TRUE;
      use.exons <- FALSE;
    } else if(analysis.type == "exonsOnly"){
      use.junctions <- FALSE;
      use.exons <- TRUE;
    }
   } else {
    if(is.null(use.junctions) | is.null(use.exons)){
      stop(paste0("Illegal syntax! If parameter use.junctions or use.exons are used, then BOTH must be set!\n use.junctions = '",use.junctions,"', use.exons = '",use.exons,"'"));
    }
    
    if(use.junctions & use.exons){
      analysis.type <- "junctionsAndExons";
    } else if(use.junctions & (! use.exons)){
      analysis.type <- "junctionsOnly";
    } else  if((! use.junctions) & use.exons){
      analysis.type <- "exonsOnly";
    } else {
      stop("Illegal syntax! Parameters use.exons and use.junctions cannot both be false!");
    }
   }
   if(verbose){
      #message("---> RJSC: countfiles: ",paste0(countfiles, collapse=","));
      message("---> RJSC: samplenames: ",paste0(samplenames, collapse=","));
      message("---> RJSC: flat.gff.file: ",flat.gff.file);
      message("---> RJSC: use.exons:",use.exons);
      message("---> RJSC: use.junctions:",use.junctions);
      message("---> RJSC: use.novel.junctions:",use.novel.junctions);
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
   if(verbose) message(paste0("---> Removed gene features. Found: ",sum(use.bins), " features to be included so far."))
   
   if(! use.junctions & ! use.novel.junctions & ! use.exons){
      stop("FATAL ERROR: At least one of: use.junctions, use.novel.junctions, or use.exons must be set to TRUE. Otherwise you've got no data to test!");
   }
   if(! use.junctions){
      use.bins <- use.bins & bin.type != "J" & bin.type != "N";
      if(verbose) message(paste0("---> Removed splice junction features. Found: ",sum(use.bins), " features to be included so far."))
      
   } else if(! use.novel.junctions){
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
      
      jscs@analysisType <- analysis.type;
      featureChar <- substr(fData(jscs)$countbinID,1,1);
      fData(jscs)[["featureType"]] <- ifelse(featureChar == "E","exonic_part",ifelse(featureChar == "J", "splice_site", "novel_splice_site"))
      varMetadata( featureData(jscs) )[ "featureType", "labelDescription" ] <- "The type of feature (exonic_part,, splice_site, or novel_splice_site).";

      
      jscs@annotationFile <- flat.gff.file
      jscs@flatGffData <- anno.data;
      pData(jscs)$countfiles <- countfiles
      
      if(verbose) message("-----> generating count vectors... (",date(),")");
      jscs@countVectors <- getAllJunctionSeqCountVectors(jscs, nCores = nCores);
      if(verbose) message("-----> count vectors generated (",date(),")");

      if(verbose) message("-> FINISHED readJunctionSeqCounts (",date(),")");  
      
      return(jscs);
      #return(list(geneCountTable, jscs))

   }else{
      #return(list(geneCountTable, newJunctionSeqCountSet(countData=dcounts, design=design, geneIDs=genesrle, countbinIDs=exons)))
      if(verbose) message("-> FINISHED readJunctionSeqCounts (",date(),")"); 
      message("Warning: flat gff annotation not set (via parameter flat.gff.file)! While technically optional, running without the annotation data may make interpretation of the data difficult. Much of the plotting functionality will not work!");
      warning("Warning: flat gff annotation not set (via parameter flat.gff.file)! While technically optional, running without the annotation data may make interpretation of the data difficult. Much of the plotting functionality will not work!");
      jscs <- newJunctionSeqCountSet(countData=dcounts, design=design, geneIDs=genesrle, countbinIDs=exons);
      jscs@analysisType <- analysis.type;
      
      if(verbose) message("-----> generating count vectors... (",date(),")");
      jscs@countVectors <- getAllJunctionSeqCountVectors(jscs, nCores = nCores);
      if(verbose) message("-----> count vectors generated (",date(),")");
      
      return(jscs);
   }
}


#setMethod("estimateDispersions", signature(object="JunctionSeqCountSet"),
estimateJunctionSeqDispersions <- function( jscs, 
                                            test.formula1 = formula(~ sample + countbin + condition : countbin),
                                            meanCountTestableThreshold=7.5, 
                                            nCores=1, test.aggregated.genes = FALSE, 
                                            use.alternate.method  = TRUE,
                                            verbose = TRUE){
   if(verbose) message("---> STARTING estimateJunctionSeqDispersions: (",date(),")");
   stopifnot( inherits( jscs, "JunctionSeqCountSet" ) )
   if( all( is.na( sizeFactors( jscs )) ) ){
     stop("Please calculate size factors before estimating dispersions\n")
   }
   fData(jscs)$meanBase <- rowMeans(counts(jscs, normalized=TRUE));
   fData(jscs)$status <- rep("OK",nrow(fData(jscs)));
   
   testable <- fData(jscs)$meanBase >= meanCountTestableThreshold;
   fData(jscs)$status[! testable] <- "LOW_COUNT";
   if(verbose) message("-----> ejsd: ",length(testable) - sum(testable)," counting bins dropped due to low coverage (mean normalized count < ",meanCountTestableThreshold,")");
   
   is.aggregate <- grepl("+", geneIDs(jscs), fixed=TRUE)
   if(verbose) message("-----> ejsd: ",sum(is.aggregate)," counting bins are from ",length(unique(grep("+", geneIDs(jscs), fixed=TRUE)))," multigene aggregates (ie. overlapping genes).");
   
   if(! test.aggregated.genes){
     if(verbose) {
       message("-----> ejsd: ",sum(is.aggregate & testable)," counting bins dropped because they were from multigene aggregates.");
       message("             (To include these aggregate genes, set test.aggregated.genes = TRUE)");
     }
     fData(jscs)$status[ testable & is.aggregate ] <- "IS_AGGREGATE_GENE";
     testable <- testable & (! is.aggregate);
   }

   onlyOneJunctionDropCt <- 0;
   for( r in split( seq_len(nrow(jscs)), geneIDs(jscs) ) ) {
     if( sum( testable[r] ) == 1 ){
       fData(jscs)$status[ r ] <- "ONLY_ONE_BIN_TESTABLE";
       testable[r] <- FALSE;
       onlyOneJunctionDropCt <- onlyOneJunctionDropCt + 1;
     }
   }
   message("-----> ejsd: ",sum(onlyOneJunctionDropCt)," counting bins dropped because they belong to a gene with only 1 testable bin.");
   fData(jscs)$testable <- testable;
   if(verbose) message("-----> ejsd: ",sum(testable)," counting bins are marked 'testable'. across ",length(unique(fData(jscs)$geneID[fData(jscs)$testable]))," genes.");
   if(verbose) message("             (",sum(fData(jscs)$featureType[fData(jscs)$testable] == "exonic_part")," exonic regions, ",
                                        sum(fData(jscs)$featureType[fData(jscs)$testable] == "splice_site")," known junctions, ",
                                        sum(fData(jscs)$featureType[fData(jscs)$testable] == "novel_splice_site")," novel junctions)");
   
   
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
   if(verbose) message("-----> ejsd: Dispersion estimation failed for ",sum(is.na( fData(jscs)$dispBeforeSharing[fData(jscs)$testable] )) ," out of ",sum(fData(jscs)$testable)," 'testable' counting bins. Setting these features to be 'untestable'");   
   fData(jscs)$status[is.na( fData(jscs)$dispBeforeSharing ) & fData(jscs)$testable] <- "DISPERSION_EST_FAILED";
   fData(jscs)$testable[which( is.na( fData(jscs)$dispBeforeSharing ) )] <- FALSE
   if(verbose) message("---> FINISHED estimateJunctionSeqDispersions (",date(),")");
   jscs
} #)



testForDiffUsage <- function( jscs,
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

getColCt <- function(geneID, merged.data, 
                          plot.exon.results, plot.junction.results, plot.novel.junction.results, 
                          plot.untestable.results){
  if(is.null(plot.exon.results)){
    plot.exon.results <- any( merged.data$featureType == "exonic_part" );
  }
  if(is.null(plot.junction.results)){
    plot.junction.results <- any( merged.data$featureType == "splice_site" );
  }
  if(is.null(plot.novel.junction.results)){
    plot.novel.junction.results <- any( merged.data$featureType == "novel_splice_site" );
  }
  
  rt <- merged.data$geneID == geneID;
  if(! plot.exon.results){
   rt <- rt & merged.data$featureType != "exonic_part" ;
  }
  if(! plot.junction.results){
   rt <- rt & merged.data$featureType != "splice_site" ;
  }
  if(! plot.novel.junction.results){
   rt <- rt & merged.data$featureType != "novel_splice_site" ;
  }
  if(! plot.untestable.results){
   rt <- rt & merged.data$testable;
  }
  return(sum(rt));
}


USE.MARGIN.MEX <- FALSE;

plotJunctionSeqResultsForGene <- function(geneID, jscs, 
                              colorRed.FDR.threshold=0.05,
                              plot.type = "expr", 
                              displayTranscripts = FALSE,
                              color = NULL, 
                              use.vst = FALSE, use.log = TRUE,truncateBelowOne = TRUE,
                              exon.rescale.factor = 0.3,
                              label.p.vals = TRUE, 
                              plot.lwd = 3, axes.lwd = plot.lwd, anno.lwd = plot.lwd, 
                              par.cex = 1, anno.cex.text = 1,
                              anno.cex.axis=anno.cex.text, anno.cex.main = anno.cex.text * 1.2,
                              fit.countbin.names = TRUE,
                              plot.gene.level.expression = NULL,
                              plot.exon.results = NULL, plot.junction.results = NULL, plot.novel.junction.results = NULL, 
                              plot.untestable.results = FALSE, draw.untestable.annotation = TRUE,
                              show.strand.arrows = 10, arrows.length = 0.125,
                              sort.features = TRUE,
                              drawCoordinates = TRUE,
                              yAxisLabels.inExponentialForm = FALSE,
                              title.main=NULL, title.ylab=NULL, graph.margins = c(2,3,3,2),
                              GENE.annotation.relative.height = 0.2, TX.annotation.relative.height = 0.025,
                              condition.legend.text = NULL, include.TX.names = TRUE, draw.start.end.sites = TRUE,
                              verbose=TRUE, debug.mode = FALSE, 
                              ...)
{
    tryCatch({
      GENE.annotation.height <- GENE.annotation.relative.height * 10;
      TX.annotation.height <- TX.annotation.relative.height * 10;
      flat.gff.data <- jscs@flatGffData;

      if(is.null(plot.exon.results)){
        plot.exon.results <- any( fData(jscs)$featureType == "exonic_part" );
      }
      if(is.null(plot.junction.results)){
        plot.junction.results <- any( fData(jscs)$featureType == "splice_site" );
      }
      if(is.null(plot.novel.junction.results)){
        plot.novel.junction.results <- any( fData(jscs)$featureType == "novel_splice_site" );
      }

     #if(! USE.MARGIN.MEX){
     #  graph.margins[2:4] <- graph.margins[2:4] * anno.cex.axis;
     #}

     if(is.null(condition.legend.text)){
       condition.legend.text <- levels(jscs@phenoData$condition);
       names(condition.legend.text) <- condition.legend.text;
     } else {
       if(is.null(names(condition.legend.text))){
         warning("names(condition.legend.text is NULL! condition.legend.text mis-formatted. Must be a list or character vector, with element names equal to the levels of pData(jscs)$condition. Falling back.");
         condition.legend.text <- levels(jscs@phenoData$condition);
         names(condition.legend.text) <- condition.legend.text;
       } else if(! all(levels(jscs@phenoData$condition) %in% names(condition.legend.text))){
         warning("Not all levels contained in names(condition.legend.text)! condition.legend.text mis-formatted. Must be a list or character vector, with element names equal to the levels of pData(jscs)$condition. Falling back.");
         condition.legend.text <- levels(jscs@phenoData$condition);
         names(condition.legend.text) <- condition.legend.text;
       }
     }

     #graph.margins <- c(2, 4, 4, 2);
     gene.level.buffer <- 0.5;

     if(is.null(plot.gene.level.expression)){
       if(plot.type == "rExpr"){
         plot.gene.level.expression <- FALSE;
       } else if(use.vst){
         plot.gene.level.expression <- FALSE;
       } else {
         plot.gene.level.expression <- TRUE;
       }
     }
     if(plot.gene.level.expression & plot.type == "rExpr"){
       warning("WARNING: plotting of gene-level expression is not supported for relative expression plots (it doesn't make sense to do so). Errors are likely to follow.");
     }
     if(plot.gene.level.expression & use.vst){
       warning("WARNING: plotting of gene-level expression is not supported for vst-transformed plots. Errors are likely to follow.");
     }


     FDR <- colorRed.FDR.threshold;
     merged.data <- fData(jscs);
     condition <- jscs@phenoData$condition;

     if(verbose){
        message("> plotJunctionSeqResultsForGene(): ", geneID, ", plot.type: ", plot.type);
     }

     rt <- merged.data$geneID == geneID;
     if(! plot.exon.results){
       if(debug.mode) message(">     Removing ",sum(rt & merged.data$featureType == "exonic_part")," exonic_part features. ",sum(rt & merged.data$featureType != "exonic_part")," features remaining");
       rt <- rt & merged.data$featureType != "exonic_part" ;
     }
     if(! plot.junction.results){
       if(debug.mode) message(">     Removing ",sum(rt & merged.data$featureType == "splice_site")," splice_site features. ",sum(rt & merged.data$featureType != "splice_site")," features remaining");
       rt <- rt & merged.data$featureType != "splice_site" ;
     }
     if(! plot.novel.junction.results){
       if(debug.mode) message(">     Removing ",sum(rt & merged.data$featureType == "novel_splice_site")," novel_splice_site features. ",sum(rt & merged.data$featureType != "novel_splice_site")," features remaining");
       rt <- rt & merged.data$featureType != "novel_splice_site" ;
     }
     untestable.rt <- which(rt & (! merged.data$testable));
     if(! plot.untestable.results){
       if(debug.mode) message(">     Removing ",sum(rt & ! merged.data$testable)," untestable features. ",sum(rt & merged.data$testable)," features remaining");
       rt <- rt & merged.data$testable;
     }
     if(debug.mode)   message(">     Plotting ",sum(rt), " features.");

     rt <- which(rt);
     if(length(rt) == 0){
       message("NO FEATURES TO PLOT!");
     } else {

       if(sort.features){
         #sort by start, then end (depreciated):
         #rt <- rt[order(merged.data$start[rt],merged.data$end[rt])];
         #New method: sort by midpoint. Looks better?
         rt <- rt[order(((merged.data$end[rt] - merged.data$start[rt])/2) + merged.data$start[rt] )];
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

       rt.allJunction <- which(flat.gff.data$gene_id==geneID & (flat.gff.data$featureType == "splice_site" | flat.gff.data$featureType == "novel_splice_site"));
       rango.allJunction <- 1:length(rt.allJunction);

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
                    vst( jscs@plottingEstimates[["relExprEstimate"]][rt,, drop=FALSE], jscs);
                  } else  if(use.log) {
                    apply(log10(jscs@plottingEstimates[["relExprEstimate"]][rt,, drop=FALSE]),c(1,2),FUN=convertY);
                  } else {
                    jscs@plottingEstimates[["relExprEstimate"]][rt,, drop=FALSE];
                  }
         plot.gene.level.expression <- FALSE;
         color.count <- rep(color[condition.names],each=nrow(count))
         y.axis.title <- "Relative Coverage";
         main.title <- paste0("Relative Coverage (",geneID,")");
       } else if(plot.type == "normCounts"){
         #count <- if(use.vst) merged.data[rt,paste("normCountVST_",sample.names,sep='')]
         #         else merged.data[rt,paste("normCount_",sample.names,sep='')]
         count <- if(use.vst){ 
                    vst(jscs@plottingEstimates[["normCounts"]][rt,1:length(sample.names), drop=FALSE], jscs) 
                  } else if(use.log) {
                    apply(log10(jscs@plottingEstimates[["normCounts"]][rt,1:length(sample.names), drop=FALSE]),c(1,2),FUN=convertY)
                  } else {
                    jscs@plottingEstimates[["normCounts"]][rt,1:length(sample.names), drop=FALSE]
                  }
         if(plot.gene.level.expression) { 
           geneCount <- if(use.vst){ 
                    plot.gene.level.expression <- FALSE;
                    NULL;
                  } else if(use.log){
                    sapply(log10(jscs@geneCountData[rownames(jscs@geneCountData) == geneID,] / sizeFactors(jscs)),FUN=convertY)
                  } else {
                    jscs@geneCountData[rownames(jscs@geneCountData) == geneID,] / sizeFactors(jscs);
           }
           color.geneCount <- color[as.character(condition)];
         }
         color.count <- rep(color[as.character(condition)], each=nrow(count));

         y.axis.title <- "Normalized Counts";
         main.title <- paste0("Normalized Counts (",geneID,")");
       } else if(plot.type == "rawCounts"){
         count <- if(use.vst){ 
                    #jscs@plottingEstimatesVST[["normCountsVST"]][rt,1:length(sample.names), drop=FALSE];
                    vst(jscs@countVectors[rt,1:length(sample.names), drop=FALSE], jscs);
                  } else if(use.log) {
                    apply(log10(jscs@countVectors[rt,1:length(sample.names), drop=FALSE]),c(1,2),FUN=function(y){max(INTERNAL.NINF.VALUE,y)})
                  } else {
                    jscs@countVectors[rt,1:length(sample.names), drop=FALSE]
                  }
         if(plot.gene.level.expression) {
           geneCount <- if(use.vst){ 
                    plot.gene.level.expression <- FALSE;
                    NULL;
                  } else if(use.log){
                    sapply(log10(jscs@geneCountData[rownames(jscs@geneCountData) == geneID,]),FUN=convertY)
                  } else {
                    jscs@geneCountData[rownames(jscs@geneCountData) == geneID,]
           }
           color.geneCount <- color[as.character(condition)];
         }
         color.count <- rep(color[as.character(condition)], each=nrow(count))
         y.axis.title <- "log10 Raw Counts";
         main.title <- paste0("Raw Counts (",geneID,")");
       } else if(plot.type == "expr"){
         count <- if(use.vst){ 
                    vst(jscs@plottingEstimates[["exprEstimate"]][rt,, drop=FALSE], jscs);
                  } else if(use.log) {
                    apply(log10(jscs@plottingEstimates[["exprEstimate"]][rt,, drop=FALSE]),c(1,2),FUN=convertY);
                  } else {
                    jscs@plottingEstimates[["exprEstimate"]][rt,, drop=FALSE];
                  }
         if(plot.gene.level.expression) {
           geneCount <- if(use.vst){ 
                    plot.gene.level.expression <- FALSE;
                    NULL;
                  } else if(use.log){
                    sapply(log10(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]][rownames(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]]) == geneID ,]),FUN=convertY);
                  } else {
                    jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]][rownames(jscs@geneLevelPlottingEstimates[["geneLevelEstModeled"]]) == geneID ,];
           }
           color.geneCount <- color[condition.names];
         }
         color.count <- rep(color[condition.names],each=nrow(count))
         y.axis.title <- "Normalized Coverage";
         main.title <- paste0("Mean Normalized Coverage (",geneID,")");
         #print(count);
       } else {
         stop(paste0("FATAL ERROR: Unknown plot type! plot.type = \"",plot.type,"\""));
       }
       if(! is.null(title.ylab)) y.axis.title <- title.ylab;
       if(! is.null(title.main)) main.title <- title.main;
       if(! plot.gene.level.expression) {
         geneCount <- NULL;
         color.geneCount <- NULL;
       }
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

       #SET COLORS FOR ANNOTATION:
       SIG.VERTLINE.COLOR <- "#F219ED60";
       NOSIG.VERTLINE.COLOR <- "#AAAAAA60";
       UNTESTABLE.VERTLINE.COLOR <- "#22222260";

       SIG.FEATURE.COLOR <- "#F219ED";
       NOSIG.FEATURE.COLOR <- "#111111";
       UNTESTABLE.FEATURE.COLOR <- "#DDDDDD";
       EXCLUDED.FEATURE.COLOR <- "#DDDDDD";

       SIG.FEATURE.BORDER.COLOR <- "#000000";
       NOSIG.FEATURE.BORDER.COLOR <- "#000000";
       UNTESTABLE.FEATURE.BORDER.COLOR <- "#AAAAAA";
       EXCLUDED.FEATURE.BORDER.COLOR <- "#000000";

       SIG.FEATURE.FILL.COLOR <- "#F219ED";
       NOSIG.FEATURE.FILL.COLOR <- "#CCCCCC";
       UNTESTABLE.FEATURE.FILL.COLOR <- "#EEEEEE";
       EXCLUDED.FEATURE.FILL.COLOR <- "#CCCCCC";

       intervals<-(0:nrow(count))/nrow(count)

       #numcond<-length(unique(design(jscs, drop=FALSE)[[fitExpToVar]]))
       numexons<-nrow(count)
       each <- merged.data$padjust[rt]
       #exoncol<-ifelse(each<=FDR, "#8B0000", "dark green")
       #exoncol[is.na(exoncol)]<-"black"
       #colorlines <- ifelse(each<=FDR, "#FF000060", "lightgrey")
       #exoncol <- ifelse(merged.data$testable[rt], ifelse(f.na(merged.data$padjust[rt] <= FDR), SIG.FEATURE.COLOR, NOSIG.FEATURE.COLOR), UNTESTABLE.FEATURE.COLOR);
       vertline.col <- ifelse(merged.data$testable[rt], ifelse(f.na(merged.data$padjust[rt] <= FDR), SIG.VERTLINE.COLOR, NOSIG.VERTLINE.COLOR), UNTESTABLE.VERTLINE.COLOR);
       annolink.col <- ifelse(merged.data$testable[rt], ifelse(f.na(merged.data$padjust[rt] <= FDR), SIG.FEATURE.COLOR,  NOSIG.FEATURE.COLOR),  UNTESTABLE.FEATURE.COLOR);
       exonlty <- rep(1,length(vertline.col));
       exonlty[as.character(merged.data$featureType[rt]) == "novel_splice_site"] <- 2 
       #exoncol[is.na(exoncol)] <- color2transparent("#CCCCCC",25)

       is.sig.feature <- f.na(each <= FDR);
       sig.feature <- which(is.sig.feature);
       #print(paste("Index of sig splices: ", sig.feature));
       #print(merged.data[rt,][sig.feature,]);

    if(debug.mode & verbose) message(">    pJSRforGene(): ", "Reached step 5.");

       #colorlines <- exoncol;  #ifelse(each<=FDR, "#F219ED60", "#B3B3B360")   # vertical dashed lines
       #colorlines[is.na(colorlines)] <- "#B3B3B360"
       #colorlinesB <- ifelse(each<=FDR, "#9E109B", "#666666")  # slanted solid lines
       #colorlinesB[is.na(colorlinesB)] <- "#666666"

       #exrt <- which(merged.data$geneID==geneID & merged.data$featureType=="exonic_part");

       ################## DETERMINE THE LAYOUT OF THE PLOT DEPENDING OF THE OPTIONS THE USER PROVIDES ###########
          sub <- data.frame(start=merged.data$start[rt], 
                            end=merged.data$end[rt], 
                            chr=merged.data$chr[rt], 
                            strand=merged.data$strand[rt], 
                            is.exon = (merged.data$featureType[rt] == "exonic_part"), 
                            is.testable = merged.data$testable[rt],
                            is.sig = is.sig.feature,
                            col = vertline.col,
                            featureType = merged.data$featureType[rt], 
                            featureID = rownames(merged.data)[rt],
                            stringsAsFactors=F);
          #if(draw.untestable.annotation & (! plot.untestable.results)){
          #  sub.untestable <- data.frame(start=merged.data$start[rt.untestable], end=merged.data$end[rt.untestable], is.exon = (merged.data$featureType[rt.untestable] == "exonic_part"), featureType = merged.data$featureType[rt.untestable], col = rep(UNTESTABLE.FEATURE.COLOR, length(rt.untestable)),stringsAsFactors=F);
          #} else {
          #  sub.untestable <- NULL;
          #}
          sub.allExon <- data.frame(start=flat.gff.data$start[rt.allExon], 
                                    end=flat.gff.data$end[rt.allExon], 
                                    chr=flat.gff.data$chrom[rt.allExon], 
                                    strand=flat.gff.data$strand[rt.allExon], 
                                    is.exon = (flat.gff.data$featureType[rt.allExon] == "exonic_part"),
                                    featureID = as.character(flat.gff.data$featureName[rt.allExon]), 
                                    stringsAsFactors=F);
          sub.allJunction <- data.frame(start=flat.gff.data$start[rt.allJunction], 
                                    end=flat.gff.data$end[rt.allJunction], 
                                    chr=flat.gff.data$chrom[rt.allJunction], 
                                    strand=flat.gff.data$strand[rt.allJunction], 
                                    is.exon = (flat.gff.data$featureType[rt.allJunction] == "exonic_part"), 
                                    feature.type = flat.gff.data$featureType[rt.allJunction] ,
                                    is.novel = (flat.gff.data$featureType[rt.allJunction] == "novel_splice_site"), 
                                    featureID = as.character(flat.gff.data$featureName[rt.allJunction]), 
                                    stringsAsFactors=F);

       testable.featureIDs <- sub$featureID[sub$is.testable];
       sig.featureIDs <- sub$featureID[sub$is.sig];
       untestable.featureIDs <- rownames(merged.data)[untestable.rt];

       if(debug.mode & verbose){ message("testable.featureIDs: ",paste0(testable.featureIDs,collapse=",")) }
       if(debug.mode & verbose){ message("sig.featureIDs: ",paste0(sig.featureIDs,collapse=",")) }
       if(debug.mode & verbose){ message("untestable.featureIDs: ",paste0(untestable.featureIDs,collapse=",")) }


       sub.allExon$is.testable   <- sub.allExon$featureID %in% testable.featureIDs;
       sub.allExon$is.sig        <- sub.allExon$featureID %in% sig.featureIDs;
       sub.allExon$is.untestable <- sub.allExon$featureID %in% untestable.featureIDs;
       sub.allJunction$is.testable <- sub.allJunction$featureID %in% testable.featureIDs;
       sub.allJunction$is.sig      <- sub.allJunction$featureID %in% sig.featureIDs;
       sub.allJunction$is.untestable <- sub.allJunction$featureID %in% untestable.featureIDs;

       sub.allExon$lineColor     <- ifelse(sub.allExon$is.testable,
                                           ifelse(sub.allExon$is.sig,
                                                  SIG.FEATURE.COLOR,
                                                  NOSIG.FEATURE.COLOR),
                                           ifelse(sub.allExon$is.untestable, 
                                                  UNTESTABLE.FEATURE.COLOR, 
                                                  EXCLUDED.FEATURE.COLOR));
       sub.allExon$fillColor     <- ifelse(sub.allExon$is.testable,ifelse(sub.allExon$is.sig,SIG.FEATURE.FILL.COLOR,NOSIG.FEATURE.FILL.COLOR),ifelse(sub.allExon$is.untestable, UNTESTABLE.FEATURE.FILL.COLOR, EXCLUDED.FEATURE.FILL.COLOR));
       sub.allExon$borderColor   <- ifelse(sub.allExon$is.testable,ifelse(sub.allExon$is.sig,SIG.FEATURE.BORDER.COLOR,NOSIG.FEATURE.COLOR),ifelse(sub.allExon$is.untestable, UNTESTABLE.FEATURE.BORDER.COLOR, EXCLUDED.FEATURE.BORDER.COLOR));
       sub.allJunction$lineColor <- ifelse(sub.allJunction$is.testable,ifelse(sub.allJunction$is.sig,SIG.FEATURE.COLOR,NOSIG.FEATURE.COLOR),ifelse(sub.allJunction$is.untestable, UNTESTABLE.FEATURE.COLOR, EXCLUDED.FEATURE.COLOR));
       sub.allJunction$lty <- ifelse(sub.allJunction$is.novel,2,1); 

       #print("##### sub.allJunction:");
       #print(sub.allJunction);

       sig.feature.names <- sub$featureID[is.sig.feature & sub$is.exon];
       allExon.isSig <- sub.allExon$featureID %in% sig.featureIDs;
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
          rel.calc.min <- min(sub.allJunction$start, sub.allExon$start)
          rel.calc.max <- max(sub.allJunction$end,   sub.allExon$end)

          if(is.na(exon.rescale.factor) | exon.rescale.factor <= 0 | exon.rescale.factor >= 1){
            rel <- (data.frame(sub.allExon$start, sub.allExon$end))-rel.calc.min;
            rel <- rel/(rel.calc.max - rel.calc.min);
            rescale.iv <- NULL;
          } else {
            #message("attempting rescale...");
            rescale.iv <- generate.interval.scale(
              data.frame(
                start = c(sub.allExon$start, sub.allJunction$start),
                end = c(sub.allExon$end, sub.allJunction$end),
                is.exon = c(sub.allExon$is.exon, sub.allJunction$is.exon)
              ),exon.rescale.factor);
            rel <- data.frame(start = rescale.coords(sub$start,rescale.iv), 
                              end   = rescale.coords(sub$end,  rescale.iv));
            #message("rescale.iv:");
            #print(rescale.iv);
            #message("rel:");
            #print(rel);
          }

         #print("rescale.iv");
         #print(rescale.iv);
          transcripts <- sapply(sapply(flat.gff.data$transcripts[rt.allExon],toString), function(x){strsplit(x, "+",fixed=TRUE)})
          trans <- Reduce(union, transcripts)
          if(displayTranscripts==TRUE){

             #if(length(trans) > 42){
             #   warning("This gene contains more than 42 transcripts annotated, only the first 42 will be plotted\n")
             #   trans <- trans[1:42];
             #} ## max support from transcripts is 45, which seems to be the max for the layout supported by graphics
             # NOTE - PROBLEM FIXED. All transcripts are now plotted on a single graphics frame, thereby avoiding this issue.
             mat <- 1:4
             hei<-c(10,1, GENE.annotation.height, TX.annotation.height * length(trans));
          }else{
             mat<-1:3
             hei<-c(10,1, GENE.annotation.height)
          }
          #Add another transcript height to the bottom, to serve as a lower margin.
          hei <- c(hei, TX.annotation.height)
          mat <- c(mat, length(mat)+1)
          layout(matrix(mat), heights=hei)
          ##par(mar=c(1, graph.margins[2:4]), cex = par.cex, mex = anno.cex.text)
          #if(USE.MARGIN.MEX){
          #  par(mar=graph.margins, cex = par.cex, mex = anno.cex.text)
          #} else {
          #  par(mar=graph.margins, cex = par.cex)
          #}


    if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached step 6.");

          ylimn <- c(min(min(count,na.rm=TRUE),0), max(count, na.rm=TRUE));
          #if(plot.type == "rawCounts" & use.log) ylimn[1] <- INTERNAL.NINF.VALUE;
          if((! use.vst) & use.log ) ylimn[1] <- INTERNAL.NINF.VALUE;

          p.values.labels <- ifelse(each<=FDR, format(each,digits=3), "");

          if(any(sub$is.exon)){
            italicize.label <- ! sub$is.exon;
          } else {
            italicize.label <- NULL;
          }

          intervals <- drawPlot(matr=count, ylimn,jscs, 
                   intervals, rango, textAxis=y.axis.title, 
                   rt=rt, color.count=color.count, 
                   colorlines=vertline.col, 
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
                   debug.mode = debug.mode, plot.gene.level.expression = plot.gene.level.expression, geneCount = geneCount, color.geneCount = color.geneCount,
                   yAxisLabels.inExponentialForm = yAxisLabels.inExponentialForm, italicize.label = italicize.label, condition.legend.text = condition.legend.text,
                   rel = rel, annolink.col = annolink.col, exonlty = exonlty, graph.margins = graph.margins,
                   ...);
    if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached end of step 6.");


    if(debug.mode & verbose) message(">    pJSRforGene(): ","Reached step 7.");

          #box()
          #########PLOT THE GENE MODEL:
          if(USE.MARGIN.MEX){
            par(mar=c(0, graph.margins[2], 0, graph.margins[4]), cex = par.cex, mex = anno.cex.text); 
          } else {
            par(mar=c(0, graph.margins[2], 0, graph.margins[4]), cex = par.cex); 
          }

          plot.new();

          # lines linking exons / splices to their column:
          segments(
                        apply((rbind(rel[rango,2], rel[rango, 1])), 2, median), 
                        0, #par("usr")[3], (old version: lines are connected.)
                        apply(rbind(intervals[rango], intervals[rango+1]-((intervals[rango+1]-intervals[rango])*0.2)), 2, median), 
                        1, col=annolink.col, lty = exonlty, lwd = plot.lwd, cex = anno.cex.text,cex.axis=anno.cex.axis, cex.main=anno.cex.main, xpd=NA, ...) #col=colorlinesB,...)


          #axes.lwd = axes.lwd, anno.lwd = anno.lwd,
          par(mar=c(2, graph.margins[2], 0, graph.margins[4]), cex = par.cex);

          #message("rel.calc.min = ", rel.calc.min);
          #message("rel.calc.max = ", rel.calc.max);

          #Get start/end sites:
          startSites <- c();
          endSites <- c();
          for(i in 1:length(trans)){
             logicexons <- sapply(transcripts, function(x){any(x == trans[i])})
             startSites <- c(startSites, sub.allExon$start[logicexons][1]);
             endSites <- c(endSites, sub.allExon$end[logicexons][sum(logicexons)]);
          }
          startSites <- unique(startSites);
          endSites <- unique(endSites);

          drawGene(rel.calc.min, 
                   rel.calc.max, 
                   tr=sub, 
                   tr.allExon=sub.allExon, 
                   tr.allJunction = sub.allJunction,
                   #tr.untestable=sub.untestable,
                   rango, 
                   rescale.iv = rescale.iv, 
                   exoncol=annolink.col,
                   allExon.exonCol=allExon.exonCol, 
                   names, 
                   trName="Gene model", 
                   anno.cex.text=anno.cex.text, par.cex = par.cex, 
                   exonlty = exonlty, 
                   plot.lwd = plot.lwd, anno.lwd=anno.lwd, 
                   show.strand.arrows = show.strand.arrows, 
                   cex.axis=anno.cex.axis, cex.main=anno.cex.main, 
                   arrows.length = arrows.length, 
                   draw.untestable.annotation = draw.untestable.annotation,
                   draw.start.end.sites = draw.start.end.sites, startSites = startSites, endSites = endSites,
                   ...)
     # graph.margins = graph.margins,
           #Maybe make these options external at some point?
           include.endpoints.on.coordinates <- FALSE;
           num.coord.miniticks.per.tick <- 10;

           if(drawCoordinates){
             if(! is.null(rescale.iv)){
               pretty.x <- pretty(c(rel.calc.min,rel.calc.max),n=5);
               pretty.interval <- pretty.x[2] - pretty.x[1];
               pretty.x <- pretty.x[pretty.x > rel.calc.min & pretty.x < rel.calc.max];
               rescaled.pretty.x <- rescale.coords(pretty.x,rescale.iv);
               #if(any(is.na(rescaled.pretty.x))){
               #  print(pretty.x);
               #  print(rescaled.pretty.x);
               #  print(rescale.iv);
               #}

               if(num.coord.miniticks.per.tick > 0){
                 rel.coord.miniticks <- pretty.interval * (1:(num.coord.miniticks.per.tick-1)) / num.coord.miniticks.per.tick;
                 unscaled.coord.miniticks <- unlist(lapply(pretty.x, function(a){ a + rel.coord.miniticks }));
                 unscaled.coord.miniticks <- c(pretty.x[1] - rel.coord.miniticks, unscaled.coord.miniticks);
                 #scaled.coord.miniticks <- unscaled.coord.miniticks;
                 rescaled.coord.miniticks <- rescale.coords(unscaled.coord.miniticks, rescale.iv) * (rel.calc.max-rel.calc.min) + rel.calc.min;
               } else {
                 rescaled.coord.miniticks <- FALSE;
               }

               if(include.endpoints.on.coordinates){
                 if(min(rescaled.pretty.x) > 0.05){
                   pretty.x <- c(rel.calc.min, pretty.x);
                   rescaled.pretty.x <- c(0,rescaled.pretty.x);
                 }
                 if(max(rescaled.pretty.x) < 0.95){
                   pretty.x <- c(pretty.x, rel.calc.max);
                   rescaled.pretty.x <- c(rescaled.pretty.x,1);
                 }
               }
               rescaled.pretty.x <- rescaled.pretty.x * (rel.calc.max-rel.calc.min) + rel.calc.min;
             } else {
               pretty.x <- pretty(c(rel.calc.min,rel.calc.max),n=5);
               coord.miniticks <- FALSE;
             }
             #axis(2,   at=logticks ,labels=FALSE, las=2, pos=0, tcl=0.25, cex.axis = anno.cex.text, ...)
             axis(1, at = rescaled.pretty.x,labels=pretty.x,cex.axis=anno.cex.text, tcl = -0.5, lwd = axes.lwd, ...);
             axis(1, at = c(rel.calc.min,rel.calc.max),labels=FALSE, cex.axis = anno.cex.text, tcl = 0, lwd = axes.lwd,...);
             axis(1, at = rescaled.coord.miniticks,labels=FALSE,cex.axis=anno.cex.text, tcl = -0.25, lwd = axes.lwd, ...);
       }

       if(displayTranscripts){
         if(USE.MARGIN.MEX){
           par(cex = par.cex, mar=c(0, graph.margins[2], 0, graph.margins[4]), mex = anno.cex.text);
         } else {
           par(cex = par.cex, mar=c(0, graph.margins[2], 0, graph.margins[4]));
         }

         plot.new();
         plot.window(xlim=c(rel.calc.min,rel.calc.max),ylim=c(0,length(trans)));

         for(i in 1:length(trans)){
            ymin <- length(trans) - i;
            if(include.TX.names){
              trName = trans[i];
            } else {
              trName = NULL;
            } 

            #logicexons <- sapply(transcripts, function(x){length(which(x==trans[i]))})
            logicexons <- sapply(transcripts, function(x){any(x == trans[i])})
            tr <- sub.allExon[logicexons,]  #   data.frame(start = sub.allExon$start[logicexons==1], end = sub.allExon$end[logicexons==1], featureType = sub.allExon$featureType[logicexons==1], stringsAsFactors = F);
            #tr <- tr[tr$featureType == "exonic_part",, drop=FALSE]
            curr.exoncol <- ifelse(allExon.isSig[logicexons],"#F219ED", "#CCCCCC");
            drawTranscript(rel.calc.min, 
                           rel.calc.max, 
                           ymin = ymin,
                           tr=tr, 
                           tr.allJunction = sub.allJunction,
                           rango=1:nrow(tr), 
                           rescale.iv=rescale.iv, 
                           #exoncol=curr.exoncol, 
                           names=c(), 
                           trName=trName, 
                           par.cex = par.cex, 
                           anno.cex.text = anno.cex.text,
                           sub.sig = sub.sig,
                           anno.lwd=anno.lwd, 
                           cex.axis=anno.cex.axis, 
                           cex.main=anno.cex.main, 
                           ...)
         }
         #box(lwd = axes.lwd, lty = 3);
      }
      par(mar = c(0,0,0,0));
      plot.new();
      par(mar = c(0,0,0,0));
      if(debug.mode & verbose) message("> plotJunctionSeqResultsForGene(): "," Done.");
    }
  }, error = function(e){
    message("Error caught while attempting plotJunctionSeqResultsForGene");
    message("---------------------");
    message("     Error text:");
    message("     ",e);
    message("");
    message("---------------------");
    message("     Input parameters:")
    message("     geneID = ",geneID);
    message("     colorRed.FDR.threshold = ", colorRed.FDR.threshold);
    message("     plot.type = ",plot.type);
    message("     displayTranscripts = ",displayTranscripts);
    message("     color = ",color);
    message("     use.vst = ",use.vst);
    message("     use.log = ",use.log);
    message("     truncateBelowOne = ",truncateBelowOne);
    message("     exon.rescale.factor = ",exon.rescale.factor);
    message("     label.p.vals = ",label.p.vals);
    message("     plot.lwd = ",plot.lwd);
    message("     axes.lwd = ",axes.lwd);
    message("     anno.lwd = ",anno.lwd);
    message("     fit.countbin.names = ",fit.countbin.names);
    message("     plot.gene.level.expression = ",plot.gene.level.expression);
    message("     plot.exon.results = ",plot.exon.results);
    message("     plot.junction.results = ",plot.junction.results);
    message("     plot.novel.junction.results = ",plot.novel.junction.results);
    message("     plot.untestable.results = ",plot.untestable.results);
    message("     draw.untestable.annotation = ",draw.untestable.annotation);
    message("     show.strand.arrows = ",show.strand.arrows);
    message("     sort.features = ",sort.features);
    message("     drawCoordinates = ",drawCoordinates);
    message("     yAxisLabels.inExponentialForm = ",yAxisLabels.inExponentialForm);
    message("---------------------");
    
    warning("Error caught while attempting plotJunctionSeqResultsForGene");
    warning("Error text:");
    warning(e);
  })
}
