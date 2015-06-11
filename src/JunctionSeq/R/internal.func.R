

generateSigExpressionEstimates <- function(ecs,nCores = 1,fitExpToVar="condition",verbose=TRUE, sig.gene.list){
   #require(statmod)
   #require(DEXSeq)
   ##require(multicore)
   
   myApply <- getMyApply(nCores);
   
   levelCt <- length(levels(ecs@phenoData$condition))
   sampleCt <- length(sampleNames(ecs@phenoData))
   
   out <- rbind( c(0.0,0.0,0.0,rep(0.0,levelCt),rep(0.0,levelCt),rep(0.0,sampleCt),rep(0.0,levelCt),rep(0.0,levelCt),rep(0.0,sampleCt)))
   #out <- rbind( c(0,0,0,rep(0,levelCt),rep(0,levelCt),rep(0,sampleCt),rep(0,levelCt),rep(0,levelCt),rep(0,sampleCt)))

   out <- data.frame(out)[0,]
      
   #names.out <- c("featureID","geneID","countbinID",
   #                paste("exprVST",sort(levels(ecs@phenoData$condition)),sep='_'),paste("expr",sort(levels(ecs@phenoData$condition)),sep='_'),
   #                paste("rExprVST",sort(levels(ecs@phenoData$condition)),sep='_'),paste("rExpr",sort(levels(ecs@phenoData$condition)),sep='_'),
   #                paste("normCountVST",sampleNames(ecs@phenoData),sep='_'),paste("normCount",sampleNames(ecs@phenoData),sep='_'),
   #                paste("rawCount",sampleNames(ecs@phenoData),sep='_')
   #                );
   names.out <- c("featureID","geneID","countbinID",
                   paste("exprVST",(levels(ecs@phenoData$condition)),sep='_'),
                   paste("expr",(levels(ecs@phenoData$condition)),sep='_'),
                   paste("rExprVST",(levels(ecs@phenoData$condition)),sep='_'),
                   paste("rExpr",(levels(ecs@phenoData$condition)),sep='_'),
                   paste("normCountVST",sampleNames(ecs@phenoData),sep='_'),
                   paste("normCount",sampleNames(ecs@phenoData),sep='_'),
                   paste("rawCount",sampleNames(ecs@phenoData),sep='_')
                   );
   #message(dim(out));
   gene.list <- unique(ecs@featureData$geneID);
   blank.line.length <- length(levels(ecs@phenoData$condition)) * 4;
   #rep(NA,length(levels(ecs@phenoData$condition)) * 4);
   
   if(verbose){
      message(paste("Starting expression estimate generation on ",length(sig.gene.list)," significant top-level features (eg genes)."));
   }
   get.one.genes.data <- function(geneIndex){
      geneID <- sig.gene.list[geneIndex];
      es <- fitAndArrangeCoefs( ecs, geneID, frm=as.formula(paste("count ~", "condition",  "* countbin")) )
      if(! is.null(es)){
          #If the model converged properly, then output the results:
          curr <- data.frame(getAllData(ecs,es,geneID))
          if(geneIndex - 1 %% 10000 == 0 & verbose){
            message("colnames:");
            message(paste0(colnames(curr),collapse=" "));
          }
          names(curr) <- names.out;
          return(curr);
      } else {
          #If the model does not converge, return NA for all fitted values.
          if(verbose){
            message(paste("glm fit failed for gene", geneID, " index: ",geneIndex));
          }
          curr <- featureData(ecs)@data[featureData(ecs)@data$geneID == geneID,c("geneID","countbinID")]
          rn <- row.names(curr);
          curr <- cbind.data.frame(rn,curr);
          names(curr)[1] <- "featureID";
          exonCt <- dim(curr)[1]
          NA.buffer <- matrix(rep(NA,blank.line.length * exonCt), nrow=exonCt, ncol=blank.line.length );
          #blankList <- lapply(1:exonCt, function(x){ blank.line });
          norCounts <- cbind.data.frame(get.norcounts.data(ecs,geneID),get.norcounts.data(ecs,geneID,vst.xform=FALSE), get.rawcounts.data(ecs,geneID))          
          curr <- cbind.data.frame(curr,NA.buffer,norCounts);
          names(curr) <- names.out;
          row.names(curr) <- rn;
          return( curr );
      }
   }
   
   genes.data.list <- myApply(as.list(1:length(sig.gene.list)), FUN=get.one.genes.data);
   
   if(verbose){
     message(paste("generated expression data for: ",length(genes.data.list),"top-level features (eg genes). Collapsing data..."));
     #message("head(genes.data.list):");
     #print(head(genes.data.list));
   }

   out <- do.call(rbind.data.frame,genes.data.list);
   names(out) <- names.out;
   
   if(verbose){
     message(paste("Collapsed expression data into a single data frame. generateAllExpressionEstimates() complete."));
   }
   #out <- out[! is.na(out$geneID) ,];
   
   #for(i in 4:length(names(out))){
   #  out[,i] <- as.numeric(as.character(out[,i]));
   #}
   
   out;
}



getAllData <- function(ecs,es,geneID,fitExpToVar="condition"){
   temp <- cbind.data.frame(get.expression.data(ecs,es,geneID),   get.expression.data(ecs,es,geneID,vst.xform=FALSE),
                 get.rexpr.data(ecs,es,geneID),   get.rexpr.data(ecs,es,geneID,vst.xform=FALSE),
                 get.norcounts.data(ecs,geneID),  get.norcounts.data(ecs,geneID,vst.xform=FALSE),
                 get.rawcounts.data(ecs,geneID))
   countbinID <- as.character(row.names(temp));
   #temp <- t(apply(temp,MAR=c(2),FUN=function(x){ as.numeric(as.character(x)) }));
   #temp <- apply(temp,MAR=c(2),FUN=function(x){ as.numeric(as.character(x)) });
   geneIDs <- rep(as.character(geneID),length(countbinID));
   featureID <- paste(geneID,countbinID,sep=':');
   out <- cbind.data.frame(featureID = featureID,geneID=geneIDs,countbinID = countbinID,temp);
   row.names(out) <- featureID;
   return(out);
}

##########################################################################
######### coefficient calculation functions from DEXSeq:
##########################################################################


doesArrayContainDim <- function(a, dimname){
  if(is.vector(a)){
    return(FALSE);
  } else {
    return(   any(  dimname == names(dimnames(a)))   );
  }
}

applyByDimnameOLD <- function(a, dimname, FUN){
  if(length(dimnames(a)) == 1){
    if(names(dimnames(a)) == dimname){
      return(FUN(a));
    } else {
      stop("Internal Error. Impossible state! (applyByDimname) Dimension not found: ",dimname);
    }
  } else {
    otherDim <- which(dimname != names(dimnames(a)));
    out <- array(apply(coef, MARGIN = otherDim, FUN = FUN));
    dimnames(out) <- dimnames(a)[otherDim];
    return(out);
  }
}

applyByDimname <- function(a, dimname, FUN){
  if(is.vector(a)){
    return(a);
  } else if(! is.array(a)){
    stop("Can only use applyByDimname on vectors or arrays.");
  } else {
    if(  any(  dimname == names(dimnames(a)))   ){
      if(length(dimnames(a)) == 1){
        if(names(dimnames(a)) == dimname){
          return(as.vector(FUN(a)));
        } else {
          stop("Internal Error. Impossible state! (applyByDimname) Dimension not found: ",dimname);
        }
      } else {
        otherDim <- which(dimname != names(dimnames(a)));
        out <- array(apply(a, MARGIN = otherDim, FUN = FUN));
        dimnames(out) <- dimnames(a)[otherDim];
        return(out);
      }
    } else {
      return(a);
    }
  }
}

extractFromArray <- function(a, dimname, val){
  if(length(dimnames(a)) == 1){
    return(a[[val]]);
  } else {
    otherDim <- which(dimname != names(dimnames(a)));
    out <- array(apply(coef, MARGIN = otherDim, FUN = function(x){
      x[[val]];
    }));
    dimnames(out) <- dimnames(a)[otherDim];
    return(out);
  }
}

#########################
# Developer's Note:
#       The function getAllData2 is a replacement for the
#       DEXSeq-based visualization modeling/estimation. Instead of modeling the
#       entire gene at once, it would model each splice junction individually
#       compared with the rest of the gene (similar to the methods used for 
#       the hypothesis tests, except using a main effect for the "condition" variable). 
#########################


#complex effect
#formula1 <- formula(~ condition + countbin + covar : countbin + condition : countbin);
getAllData2 <- function(jscs, runOnFeatures = seq_len(nrow(fData(jscs))),
                        fitExpToVar="condition", 
                        formula1 = formula(paste0("~ ",fitExpToVar," + countbin + ",fitExpToVar," : countbin")),
                        nCores=1, dispColumn="dispersion",verbose = TRUE){
  
  varlist <- rownames(attr(terms(formula1), "factors"));
  covarlist <- varlist[ ! varlist %in% c(fitExpToVar, "countbin") ]
  
  myApply <- getMyApply(nCores);
  
  modelFrame <- constructModelFrame( jscs )
  modelFrame$countbin <- factor(modelFrame$countbin, levels = c("others","this"));
  mm <- rmDepCols( model.matrix( formula1, modelFrame ) );
  conditionLevels <- levels(modelFrame[[fitExpToVar]]);
  conditionCt <- length(conditionLevels);
                  
  fitList <- myApply(runOnFeatures, function(i){
    if( verbose & i %% 1000 == 0 ){
       message(paste0("-------> generateCompleteResults: (Calculating effect size and predicted values for feature ",i," of ",nrow(jscs),")","(",date(),")"));
    }
    geneID <- fData(jscs)$geneID[i];
    countbinID <- fData(jscs)$countbinID[i];
    #countVector <- getJunctionSeqCountVector( jscs, geneID = geneID, countbinID = countbinID);
    countVector <- jscs@countVectors[i,];
    
    if(! fData(jscs)$testable[i]){
      return(c(rep(NA,conditionCt - 1), rep(NA,conditionCt), rep(NA,conditionCt), 
                                        rep(NA,conditionCt), rep(NA,conditionCt), 
                                        rep(NA,conditionCt), rep(NA,conditionCt), 
                                        countVector * modelFrame$sizeFactor, 
                                        vst(countVector * modelFrame$sizeFactor, jscs),
                                        countVector));
    } else {
      disp <- fData(jscs)[i, dispColumn];
      fit <- try( {
        glmnb.fit( mm,  countVector, dispersion = disp, offset = log( modelFrame$sizeFactor ) );
      });
      if( any(inherits( fit, "try-error" ) )) {
        warning( sprintf("glmnb.fit failed for %s:%s\n", as.character( geneID ), countbinID) );
        return(list(fitConverged = FALSE, fit = NULL));

        return(c(rep(NA,conditionCt - 1), rep(NA,conditionCt), rep(NA,conditionCt), 
                                          rep(NA,conditionCt), rep(NA,conditionCt), 
                                          rep(NA,conditionCt), rep(NA,conditionCt), 
                                          countVector * modelFrame$sizeFactor, 
                                          vst(countVector * modelFrame$sizeFactor, jscs),
                                          countVector));
      }
      
      coefs <- arrangeCoefs( formula1, modelFrame, mm, fit = fit, insertValues = TRUE );
      
      predictedEstimates <- lapply(conditionLevels, function(pvn){
        pvToAdd <- sapply(coefs, function(coef){
          if("(Intercept)" %in% names(dimnames(coef))){
            return(coef[1]);
          }
          if(doesArrayContainDim(coef, "countbin")){
            coef <- applyByDimname(coef, dimname = "countbin", FUN = function(cf){
              cf[["this"]];
            });
          }
          if(doesArrayContainDim(coef, fitExpToVar)){
            coef <- applyByDimname(coef, dimname = fitExpToVar, FUN = function(cf){
              cf[[pvn]];
            });
          }
          for(covar in covarlist){
            if(doesArrayContainDim(coef, covar)){
              coef <- applyByDimname(coef, dimname = covar, FUN = function(cf){
                mean(cf);
              });
            }
          }
          #message("class(coef) = ", class(coef));
          #message("length(coef) = ", length(coef));
          as.vector(coef);
        });
        pvToAdd;
      });
      names(predictedEstimates) <- conditionLevels;
      
      predictedEstimatesOther <- lapply(conditionLevels, function(pvn){
        pvToAdd <- sapply(coefs, function(coef){
          if("(Intercept)" %in% names(dimnames(coef))){
            return(coef[1]);
          }
          if(doesArrayContainDim(coef, "countbin")){
            coef <- applyByDimname(coef, dimname = "countbin", FUN = function(cf){
              cf[["others"]];
            });
          }
          if(doesArrayContainDim(coef, fitExpToVar)){
            coef <- applyByDimname(coef, dimname = fitExpToVar, FUN = function(cf){
              cf[[pvn]];
            });
          }
          for(covar in covarlist){
            if(doesArrayContainDim(coef, covar)){
              coef <- applyByDimname(coef, dimname = covar, FUN = function(cf){
                mean(cf);
              });
            }
          }
          #message("class(coef) = ", class(coef));
          #message("length(coef) = ", length(coef));
          as.vector(coef);
        });
        pvToAdd;
      });
      names(predictedEstimatesOther) <- conditionLevels;
      
      meanCondition <- mean(sapply(predictedEstimates, function(pe){
          pe[[fitExpToVar]];
      }));
      
      relativeEstimates <- predictedEstimates;
      for(j in seq_len(length(predictedEstimates))){
        relativeEstimates[[j]] <- predictedEstimates[[j]] - predictedEstimatesOther[[j]];
      }
      relativeLogExprEstimate <- sapply(relativeEstimates, sum);
      logFCs <- relativeLogExprEstimate[-1] - relativeLogExprEstimate[1];
      
      exprEstimate <- exp(sapply(predictedEstimates, sum));
      otherExprEstimate <- exp(sapply(predictedEstimatesOther, sum));
      relExprEstimate <- exp(relativeLogExprEstimate + predictedEstimates[[1]][["(Intercept).(Intercept)"]] + predictedEstimates[[1]][["countbin"]]);
      
      c(logFCs, 
        exprEstimate, vst(exprEstimate, jscs),
        otherExprEstimate, vst(otherExprEstimate, jscs),
        relExprEstimate, vst(relExprEstimate, jscs),
        countVector * modelFrame$sizeFactor, vst(countVector * modelFrame$sizeFactor, jscs),
        countVector);
      
      #return(list(fitConverged = TRUE, 
      #            fit = fit, 
      #            exprEstimate = exprEstimate, 
      #            otherExprEstimate = otherExprEstimate, 
      #            relativeExprEstimate = relativeExprEstimate, 
      #            logFCs = logFCs));
    }
  });
  
  sampleNames <- colnames(counts(jscs));
  
  
  
  out <- do.call(rbind.data.frame, fitList);
  names(out) <- c(paste0("logFC(",conditionLevels[-1],"/",conditionLevels[1],")"),
                  paste0("expr_",conditionLevels), paste0("exprVST_",conditionLevels),
                  paste0("exprGene_",conditionLevels),  paste0("exprGeneVST_",conditionLevels),
                  paste0("rexpr_",conditionLevels), paste0("rexprVST_",conditionLevels),
                  
                  paste0("normCount_",sampleNames), paste0("normGeneCount_",sampleNames), 
                    paste0("normCountVST_",sampleNames), paste0("normGeneCountVST_",sampleNames), 
                  paste0("rawCount_",sampleNames), paste0("rawGeneCounts_",sampleNames) 
                  );
  
  return(out);
}

getLogFoldChangeFromModel <- function(formula1, modelFrame, mm, fit){
  fitExpToVar <- "condition"
  conditionLevels <- levels(modelFrame[[fitExpToVar]])
  
  varlist <- rownames(attr(terms(formula1), "factors"));
  covarlist <- varlist[ ! varlist %in% c(fitExpToVar, "countbin") ]
  #conditionLevels <- levels(pData(jscs)[[fitExpToVar]]);
  
  coefs <- arrangeCoefs( formula1, modelFrame, mm, fit = fit, insertValues = TRUE );
  
  predictedEstimates <- lapply(conditionLevels, function(pvn){
    pvToAdd <- sapply(coefs, function(coef){
      if("(Intercept)" %in% names(dimnames(coef))){
        return(coef[1]);
      }
      if(doesArrayContainDim(coef, "countbin")){
        coef <- applyByDimname(coef, dimname = "countbin", FUN = function(cf){
          cf[["this"]];
        });
      }
      if(doesArrayContainDim(coef, fitExpToVar)){
        coef <- applyByDimname(coef, dimname = fitExpToVar, FUN = function(cf){
          cf[[pvn]];
        });
      }
      for(covar in covarlist){
        if(doesArrayContainDim(coef, covar)){
          coef <- applyByDimname(coef, dimname = covar, FUN = function(cf){
            mean(cf);
          });
        }
      }
      #message("class(coef) = ", class(coef));
      #message("length(coef) = ", length(coef));
      as.vector(coef);
    });
    pvToAdd;
  });
  names(predictedEstimates) <- conditionLevels;

  predictedEstimatesOther <- lapply(conditionLevels, function(pvn){
    pvToAdd <- sapply(coefs, function(coef){
      if("(Intercept)" %in% names(dimnames(coef))){
        return(coef[1]);
      }
      if(doesArrayContainDim(coef, "countbin")){
        coef <- applyByDimname(coef, dimname = "countbin", FUN = function(cf){
          cf[["others"]];
        });
      }
      if(doesArrayContainDim(coef, fitExpToVar)){
        coef <- applyByDimname(coef, dimname = fitExpToVar, FUN = function(cf){
          cf[[pvn]];
        });
      }
      for(covar in covarlist){
        if(doesArrayContainDim(coef, covar)){
          coef <- applyByDimname(coef, dimname = covar, FUN = function(cf){
            mean(cf);
          });
        }
      }
      #message("class(coef) = ", class(coef));
      #message("length(coef) = ", length(coef));
      as.vector(coef);
    });
    pvToAdd;
  });
  names(predictedEstimatesOther) <- conditionLevels;

  relativeEstimates <- predictedEstimates;
  for(j in seq_len(length(predictedEstimates))){
    relativeEstimates[[j]] <- predictedEstimates[[j]] - predictedEstimatesOther[[j]];
  }
  relativeLogExprEstimate <- sapply(relativeEstimates, sum);
  logFCs <- relativeLogExprEstimate[-1] - relativeLogExprEstimate[1];
  
  names(logFCs) <- paste0("logFC(",conditionLevels[-1],"/",conditionLevels[1],")")
  
  return(logFCs)
}




  # names.out <- c("featureID","geneID","countbinID",
   #                paste("exprVST",(levels(jscs@phenoData$condition)),sep='_'),
   #                paste("expr",(levels(jscs@phenoData$condition)),sep='_'),
   #                paste("rExprVST",(levels(jscs@phenoData$condition)),sep='_'),
   #                paste("rExpr",(levels(jscs@phenoData$condition)),sep='_'),
   #                paste("normCountVST",sampleNames(jscs@phenoData),sep='_'),
   #                paste("normCount",sampleNames(jscs@phenoData),sep='_'),
   #                paste("rawCount",sampleNames(jscs@phenoData),sep='_')
   #                );


get.expression.data <- function(ecs, es, geneID, fitExpToVar="condition", vst.xform = TRUE){
      #es <- fitAndArrangeCoefs( ecs, geneID, frm=as.formula(paste("count ~", fitExpToVar,  "* exon")) )
      if(is.null(es)){
          warning(sprintf("glm fit failed for gene %s", geneID))
          return()
      }
      coeff <- as.matrix( t( getEffectsForPlotting(es, averageOutExpression=FALSE, groupingVar=fitExpToVar) ) )
      coeff <- exp(coeff)
      #ylimn <- c(0, max(coeff, na.rm=TRUE))
      if(vst.xform) coeff <- vst( coeff, ecs )
      return(coeff)
}

get.rexpr.data <- function(ecs,es,geneID,fitExpToVar="condition", vst.xform = TRUE){
      #coeff <- as.matrix( t( getEffectsForPlotting( fitAndArrangeCoefs( ecs, geneID, frm=as.formula(paste("count ~", fitExpToVar,  "* exon")) ), averageOutExpression=TRUE, groupingVar=fitExpToVar) ) )
      coeff <- as.matrix( t( getEffectsForPlotting( es, averageOutExpression=TRUE, groupingVar=fitExpToVar) ) )
      
      coeff <- exp(coeff)
      #ylimn <- c(0, max(coeff, na.rm=TRUE))
      if(vst.xform) coeff <- vst( coeff, ecs )
      return(coeff)
}

get.norcounts.data <- function(ecs,geneID,fitExpToVar="condition", vst.xform = TRUE){
      count <- countTableForGene(ecs, geneID, normalized=TRUE)
      #ylimn <- c(0, max(count, na.rm=TRUE))
      if(vst.xform) count <- vst( count, ecs );
      return(count)
}
get.rawcounts.data <- function(ecs,geneID,fitExpToVar="condition", vst.xform = FALSE){
      count <- countTableForGene(ecs, geneID, normalized=FALSE)
      #ylimn <- c(0, max(count, na.rm=TRUE))
      if(vst.xform) count <- vst( count, ecs );
      return(count)
}

countTableForGene <- function( ecs, geneID, normalized=FALSE, withDispersion=FALSE ) {
   stopifnot( is( ecs, "JunctionSeqCountSet") )
   rows <- geneIDs(ecs) == geneID
   if( sum(rows) == 0 )
      stop( "geneID not found" )
	ans <- counts( ecs )[ rows, , drop=FALSE ]
	rownames( ans ) <- countbinIDs(ecs)[ rows ]
	attr( ans, "geneID" ) <- geneID
	if( normalized )
	   if(any(is.na(sizeFactors(ecs)))){
   	   stop("Please first call function estimateSizeFactors\n")
	   }else{
	   ans <- t( t(ans) / sizeFactors(ecs) )
	   }
	if( withDispersion )
      attr( ans, "dispersion" ) <- fData(ecs)$dispersion[rows]
   ans
}



getEffectsForPlotting <- function( coefs, groupingVar = "condition", averageOutExpression=FALSE )
{
   groupingExonInteraction <- which( sapply( coefs, function(x) 
      all( c( groupingVar, "countbin") %in% names(dimnames(x)) ) & length(dim(x)) == 2 ) ) 
   fittedValues <- coefs[[ groupingExonInteraction ]]
   if( names(dimnames(fittedValues))[1] == "countbin" )
      fittedValues <- t( fittedValues )
   stopifnot( identical( names(dimnames(fittedValues)), c( groupingVar, "countbin" ) ) )
   for( x in coefs[ -groupingExonInteraction ] ) {
      if( all( c( groupingVar, "countbin") %in% names(dimnames(x)) ) )
         stop( "Cannot yet deal with third-order terms." )
      if( !any( c( groupingVar, "countbin") %in% names(dimnames(x)) ) ) {
         fittedValues <- fittedValues + mean( x )
      } else if( averageOutExpression & identical( names(dimnames(x)), groupingVar ) ) {
         fittedValues <- fittedValues + mean( x )
      } else if( groupingVar %in% names(dimnames(x)) ) {
         groupMeans <- apply2( x, groupingVar, mean )
         stopifnot( identical( names(groupMeans), dimnames(fittedValues)[[1]] ) )
         fittedValues <- fittedValues + groupMeans
      } else if( "countbin" %in% names(dimnames(x)) ) {
         exonMeans <- apply2( x, "countbin", mean )
         fittedValues <- t( t(fittedValues) + exonMeans )
      } else {
         message( x )
         stop( "Unexpected term encountered." )
      }
   }
   fittedValues
}

fitAndArrangeCoefs <- function( ecs, geneID, frm = count ~ condition * countbin, balanceFeatures = TRUE )
{
   mf <- modelFrameForGene( ecs, geneID )
   if( length(levels(mf$countbin)) <= 1 )
      return( NULL )
   mm <- model.matrix( frm, mf )
   fit <- try( glmnb.fit(mm, mf$count, dispersion=mf$dispersion, offset=log(mf$sizeFactor)), silent=TRUE)
   if( is( fit, "try-error" ) )
      return( NULL )
   coefs <- arrangeCoefs( frm, mf, mm, fit )
   if( balanceFeatures ) 
      balanceFeatures( coefs, tapply( mf$dispersion, mf$countbin, `[`, 1 ) )
   else
      coefs
}


#################################
### Internal Calc Functions:


getAllJunctionSeqCountVectors <- function( ecs,  use.alternate.method = TRUE, nCores = 1 ) {
   stopifnot( inherits( ecs, "JunctionSeqCountSet" ) )
   gct <- ecs@geneCountData;
   
   myApply <- getMyApply(nCores);
   
   message("    getAllJunctionSeqCountVectors: dim(counts) = ",paste(dim(counts(ecs)),  sep=",",collapse=",")  , " (",date(),")" );
   message("    getAllJunctionSeqCountVectors: dim(gct) = ",paste(dim(gct),sep=",",collapse=",") );
   
   if(use.alternate.method){
      out.list <- myApply(1:nrow(fData(ecs)), function(i){
        geneID <- geneIDs(ecs)[i];
        #countbinID <- countbinIDs(ecs)[i];
        rowWithGeneID <- which(rownames(gct) == geneID);
        stopifnot( length(rowWithGeneID) == 1 );
        binCounts <- counts(ecs)[i,];
        geneCounts <- gct[rowWithGeneID,];
        return( c(binCounts, geneCounts - binCounts) );
      });
      message("    getAllJunctionSeqCountVectors: out.list generated. length = ",length(out.list), " (",date(),")");
      out <- as.matrix(do.call(rbind, out.list));
      
      #out <- as.matrix(t(sapply(1:nrow(fData(ecs)),function(i){
      #  geneID <- geneIDs(ecs)[i];
      #  countbinID <- countbinIDs(ecs)[i];
      #  rowWithGeneID <- which(rownames(gct) == geneID);
      #  stopifnot( length(rowWithGeneID) == 1 );
      #  binCounts <- counts(ecs)[i,];
      #  geneCounts <- gct[rowWithGeneID,];
      #  return( c(binCounts, geneCounts - binCounts) );
      #})));
      message("    getAllJunctionSeqCountVectors: out generated. dim = ",paste(dim(out),sep=",",collapse=",")," (",date(),")");
      
      rownames(out) <- rownames(fData(ecs));
      colnames(out) <- c( paste0(colnames(counts(ecs)),"_thisBin")  , paste0(colnames(counts(ecs)),"_gene")  );
      return(out);
   } else {
      stop("Fallback method is depreciated and no longer supported!");
   }
}

#getJunctionSeqCountVector <- function( ecs, geneID, countbinID, use.alternate.method = TRUE ) {
#   stopifnot( inherits( ecs, "JunctionSeqCountSet" ) )
#   gct <- ecs@geneCountData;
#   
#   if(! geneID %in% rownames(gct))                stop("FATAL ERROR: getJunctionSeqCountVector: geneID ",geneID," not found in gct");
#   if(! geneID %in% geneIDs(ecs))                 stop("FATAL ERROR: getJunctionSeqCountVector: geneID ",geneID," not found in ecs");
#   
#   if(use.alternate.method){
#      rowWithcountbinID <- which(geneIDs(ecs) == geneID & countbinIDs(ecs) == countbinID);
#      stopifnot(length(rowWithcountbinID) == 1);
#      
#      rowWithGeneID <- which(rownames(gct) == geneID);
#      stopifnot( length(rowWithGeneID) == 1 );
#      
#      binCounts <- counts(ecs)[rowWithcountbinID,];
#      geneCounts <- gct[rowWithGeneID,];
#      
#      return( c(binCounts, geneCounts - binCounts) );
#   } else {
#      rowsWithGeneID <- which(geneIDs(ecs) == geneID);
#      rowWithcountbinID <- which(countbinIDs(ecs)[rowsWithGeneID] == countbinID);
#      
#      stopifnot( length(rowWithcountbinID) == 1 );
#      
#      ctfg <- counts(ecs)[rowsWithGeneID,];
#      
#      return(c( ctfg[rowWithcountbinID,], colSums( ctfg[ -rowWithcountbinID, , drop=FALSE ] ) ));
#   }
#}

modelFrameForGene.v2 <- function( ecs, geneID, onlyTestable=FALSE) {
   stopifnot( is( ecs, "JunctionSeqCountSet") )
   if( onlyTestable & any(colnames(fData(ecs)) %in% "testable")){
      rows <- geneIDs(ecs) == geneID & fData(ecs)$testable
   }else{
      rows <- geneIDs(ecs) == geneID
   }

   numJunctions <- sum(rows)
   junctionCol <- rep( factor( countbinIDs(ecs)[rows], levels=countbinIDs(ecs)[rows] ), ncol( counts(ecs) ) )
   modelFrame <- data.frame(
      sample = rep( factor( colnames( counts(ecs) ) ), each = numJunctions ),
      countbin = junctionCol,
      sizeFactor = rep( sizeFactors(ecs), each = numJunctions ) )
   #for( cn in colnames( design(ecs,drop=FALSE) ) )
   #   modelFrame[[cn]] <- factor(rep( design(ecs,drop=FALSE)[[cn]], each=numJunctions ), levels=sort(levels(design(ecs,drop=FALSE)[[cn]] )))
   for( cn in colnames( design(ecs,drop=FALSE) ) )
         modelFrame[[cn]] <- factor(rep( design(ecs,drop=FALSE)[[cn]], each=numJunctions ), levels=levels(design(ecs,drop=FALSE)[[cn]] ))
   modelFrame$dispersion <- fData(ecs)$dispersion[ rows ][
      match( modelFrame$countbin, countbinIDs(ecs)[rows] ) ]
   modelFrame$count <- as.vector( counts(ecs)[rows,] )
   attr( modelFrame, "geneID" ) <- geneID
   modelFrame
}



modelFrameForGene <- function( ecs, geneID, onlyTestable=FALSE) {
   stopifnot( is( ecs, "JunctionSeqCountSet") )
   if( onlyTestable & any(colnames(fData(ecs)) %in% "testable")){
      rows <- geneIDs(ecs) == geneID & fData(ecs)$testable
   }else{
      rows <- geneIDs(ecs) == geneID
   }

   numJunctions <- sum(rows)
   junctionCol <- rep( factor( countbinIDs(ecs)[rows], levels=countbinIDs(ecs)[rows] ), ncol( counts(ecs) ) )
   modelFrame <- data.frame(
      sample = rep( factor( colnames( counts(ecs) ) ), each = numJunctions ),
      countbin = junctionCol,
      sizeFactor = rep( sizeFactors(ecs), each = numJunctions ) )
   #for( cn in colnames( design(ecs,drop=FALSE) ) )
   #   modelFrame[[cn]] <- factor(rep( design(ecs,drop=FALSE)[[cn]], each=numJunctions ), levels=sort(levels(design(ecs,drop=FALSE)[[cn]] )))
   for( cn in colnames( design(ecs,drop=FALSE) ) )
         modelFrame[[cn]] <- factor(rep( design(ecs,drop=FALSE)[[cn]], each=numJunctions ), levels=levels(design(ecs,drop=FALSE)[[cn]] ))
   modelFrame$dispersion <- fData(ecs)$dispersion[ rows ][
      match( modelFrame$countbin, countbinIDs(ecs)[rows] ) ]
   modelFrame$count <- as.vector( counts(ecs)[rows,] )
   attr( modelFrame, "geneID" ) <- geneID
   modelFrame
}

testFeatureForDJU.fromRow <- function(formula1, ecs, i, modelFrame, mm0, mm1, disp, keepCoefs = ncol(mm1), use.alternate.method = TRUE){
  stopifnot( inherits( ecs, "JunctionSeqCountSet" ) )
  geneID <- geneIDs(ecs)[i];
  countbinID <- countbinIDs(ecs)[i];
  #stopifnot( any(geneIDs(ecs) %in% geneID & countbinIDs(ecs) %in% countbinID ) )
  stopifnot( !is.na(disp) )
  if( all( is.na( sizeFactors( ecs )) ) ){
    stop("Please calculate size factors first\n")
  }
  countVector <- ecs@countVectors[i,];
  
  #countVector <- getJunctionSeqCountVector( ecs, geneID, countbinID , use.alternate.method = use.alternate.method);
  #countVector <- ecs@countVectors[i,];
  conditionLevels <- levels(modelFrame[["condition"]]);
  
  if(is.na(disp)){
     warning("Dispersion is NA!");
     return(list(coefficient = rep(NA,length(keepCoefs)), logFC = rep(NA, length(conditionLevels) - 1),pval = NA, disp = NA, countVector = countVector, fit = list(fitH0 = NA, fitH1 = NA) ));
  }
  
  return(testFeatureForDJU.fromCountVector(formula1 = formula1, mm0 = mm0,mm1 = mm1,disp = disp, keepCoefs = keepCoefs, modelFrame = modelFrame, countVector = countVector, 
                                           geneID = geneID, countbinID = countbinID));
}


#testFeatureForDJU <- function(formula1, ecs, geneID, countbinID, modelFrame, mm0, mm1, disp, keepCoefs = ncol(mm1), use.alternate.method = TRUE){
#  stopifnot( inherits( ecs, "JunctionSeqCountSet" ) )
#  stopifnot( any(geneIDs(ecs) %in% geneID & countbinIDs(ecs) %in% countbinID ) )
#  stopifnot( !is.na(disp) )
#  if( all( is.na( sizeFactors( ecs )) ) ){
#    stop("Please calculate size factors first\n")
#  }
#  countVector <- getJunctionSeqCountVector( ecs, geneID, countbinID , use.alternate.method = use.alternate.method);
#  #countVector <- ecs@countVectors[i,];
#  conditionLevels <- levels(modelFrame[["condition"]]);
#  
#  if(is.na(disp)){
#     warning("Dispersion is NA!");
#     return(list(coefficient = rep(NA,length(keepCoefs)), logFC = rep(NA, length(conditionLevels) - 1),pval = NA, disp = NA, countVector = countVector, fit = list(fitH0 = NA, fitH1 = NA) ));
#  }
#  
#  return(testFeatureForDJU.fromCountVector(formula1 = formula1, mm0 = mm0,mm1 = mm1,disp = disp, keepCoefs = keepCoefs, modelFrame = modelFrame, countVector = countVector, 
#                                           geneID = geneID, countbinID = countbinID));
#}



testFeatureForDJU.fromCountVector <- function(formula1, mm0,mm1,disp, keepCoefs = ncol(mm1), modelFrame,countVector, 
                                              # The remaining parameters have no effect on functionality, only impact error reporting:
                                              geneID = "", countbinID = ""){
  fitExpToVar <- "condition"
  conditionLevels <- levels(modelFrame[[fitExpToVar]])
  fitB <- try( {
      fit0 <- glmnb.fit( mm0,  countVector, dispersion = disp, offset = log( modelFrame$sizeFactor ) );
      fit1 <- glmnb.fit( mm1,  countVector, dispersion = disp, offset = log( modelFrame$sizeFactor ) );
      list(fit0,fit1);
  });
  if( any(inherits( fitB, "try-error" ) )) {
     warning( sprintf("glmnb.fit failed for %s:%s\n", as.character( geneID ), countbinID) );
     return(list(coefficient = rep(NA,length(keepCoefs)), logFC = rep(NA, length(conditionLevels) - 1),pval = NA, disp = NA, countVector = countVector, fit = list(fitH0 = NA, fitH1 = NA) ));
  }
  fit0 <- fitB[[1]];
  fit1 <- fitB[[2]];

  coefficient <- fit1$coefficients[keepCoefs];
  
  pval <- pchisq( deviance( fit0 ) - deviance( fit1 ), ncol( mm1 ) - ncol( mm0 ), lower.tail=FALSE );
  
  #logFC <- rep(NA, length(conditionLevels) - 1) # getLogFoldChangeFromModel(conditionLevels, formula1, modelFrame, mm1, fit1);
  logFC <- getLogFoldChangeFromModel(formula1, modelFrame, mm1, fit1);
  
  return(list(coefficient = coefficient, logFC = logFC, pval = pval, disp = disp, countVector = countVector, fit = list(fitH0 = fit0, fitH1 = fit1) ));
}


#testFeatureForDJU_OLDVERSION <- function(ecs, gct, geneID, countbinID, modelFrame, mm0, mm1, disp, use.alternate.method = TRUE){
#  stopifnot( inherits( ecs, "JunctionSeqCountSet" ) )
#  stopifnot( any(geneIDs(ecs) %in% geneID & countbinIDs(ecs) %in% countbinID ) )
#  stopifnot( !is.na(disp) )
#  if( all( is.na( sizeFactors( ecs )) ) ){
#    stop("Please calculate size factors first\n")
#  }
#  count <- getJunctionSeqCountVector( ecs, gcs, geneID, countbinID , use.alternate.method = use.alternate.method);
#  
#  fit0 <- glmnb.fit( mm0,  count, dispersion = disp, offset = log( modelFrame$sizeFactor ) );
#  fit1 <- glmnb.fit( mm1, count, dispersion = disp, offset = log( modelFrame$sizeFactor ) );
#
#  coefficient <- fit1$coefficients[length(fit1$coefficients)]
#  
#  pval <- pchisq( deviance( fit0 ) - deviance( fit1 ), ncol( mm1 ) - ncol( mm0 ), lower.tail=FALSE )
#  
#  return(list(coefficient, pval));
#}