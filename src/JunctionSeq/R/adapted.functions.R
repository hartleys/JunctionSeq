#These functions were extracted directly from the DEXSeq source-code.
#
#Most of these functions relate to setting up, running, and/or interpreting
#   the output from negative-binomial generalized linear models via the
#   

#setMethod("estimateSizeFactors", signature(object="JunctionSeqCountSet"),
#   function( object ){
#      cds <- object
#      stopifnot( is( cds, "JunctionSeqCountSet") )
#      geomeans <- exp( rowMeans( log( counts(cds) ) ) )
#      sizeFactors(cds) <-
#         apply( counts(cds), 2, function(cnts)
#            median( ( cnts / geomeans )[ geomeans>0 ] ) )
#      cds
#   }
#)


vst <- function(x, ecs)
{
  if(is.null(ecs@dispFitCoefs)){
    stop("Cannot perform variance stabilizing transform without parametrically-fitted dispersion.");
  }

suppressWarnings(
  ( 2 / ( sqrt(ecs@dispFitCoefs[1]) ) ) * 
     log( 2 * ecs@dispFitCoefs[1] * sqrt(x) + 
          2 * sqrt( ecs@dispFitCoefs[1] * ( ecs@dispFitCoefs[2] + 1 + ecs@dispFitCoefs[1] * x ) ) ) -
  ( 2 / ( sqrt(ecs@dispFitCoefs[1]) ) ) * 
     log( 2 * sqrt( ecs@dispFitCoefs[1] * ( ecs@dispFitCoefs[2] + 1 ) ) ) )
}

rmDepCols <- function(m) {
   q <- qr( m )
   if( q$rank < ncol(m) ){
      #Added functionality: pass on the "assign" attribute:
      out <- m[ , -q$pivot[ (q$rank+1) : ncol(m) ] ]
      attr(out,"assign") <- attr(m,"assign")[ -q$pivot[ (q$rank+1) : ncol(m) ] ]
      out;
   }else{
      m
   }
}


logConditionalLikelihood <- function( disp, mm, y, muhat )
{
   # calculate the log likelihood:
   if(length(disp) != length(y)){
      disp <- rep(disp, length(y))
   }

   ll <- sum( sapply( seq(along=y), function(i)
      dnbinom( y[i], mu=muhat[i], size=1/disp[i], log=TRUE ) ) )

   # transform the residuals, i.e., y - muhat, to the linear
   # predictor scale by multiplying them with the derivative
   # of the link function, i.e., by 1/muhat, and add this to the
   # linear predictors, log(muhat), to get the predictors that
   # are used in the IWLS regression
   z <- log(muhat) + ( y - muhat ) / muhat

   # the variance function of the NB is as follows
   v0 <- muhat + disp * muhat^2

   # transform the variance vector to linear predictor scale by
   # multiplying with the squared derivative of the link to
   # get the (reciprocal) weights for the IWLS
   w <- 1 / ( ( 1 / muhat )^2 * v0 )

   # All we need from the IRLS run is the QR decomposition of
   # its matrix
   qrres <- qr( mm*sqrt(w) )

   # from it, we extract we leverages and calculate the Cox-Reid
   # term:
   cr <- sum( log( abs( diag( qrres$qr )[ seq_len(qrres$rank) ] ) ) )

   # return the profile log likelihood:
   ll - cr
}

#estimateFeatureDispersion <- function( ecs, geneID, countbinID, modelFrame, mm , use.alternate.method  = TRUE ){
#   stopifnot( inherits( ecs, "JunctionSeqCountSet" ) )
#   if( all( is.na( sizeFactors( ecs ) )) ){
#     stop("Please calculate size factors before estimating dispersions\n")
#   }
#   
#   count <- getJunctionSeqCountVector( ecs, geneID, countbinID , use.alternate.method)
#   disp <- .1
#   for( i in 1:10 ) {
#     fit <- glmnb.fit( mm, count, dispersion = disp, offset = log( modelFrame$sizeFactor ) )
#     olddisp <- disp
#     disp <- exp( optimize( function(logalpha)
#        logConditionalLikelihood( exp(logalpha), mm, count, fitted.values(fit) ), 
#        log( c( 1e-11, 1e5 ) ), maximum=TRUE, tol=.01 )$maximum )
#     if( abs( log(disp) - log(olddisp) ) < .03 )
#        break
#     }
#  disp
#}

estimateFeatureDispersionFromRow <- function( ecs, i, modelFrame, mm , use.alternate.method  = TRUE ){
   stopifnot( inherits( ecs, "JunctionSeqCountSet" ) )
   if( all( is.na( sizeFactors( ecs ) )) ){
     stop("Please calculate size factors before estimating dispersions\n")
   }
   
   count <- ecs@countVectors[i,];
   disp <- .1
   for( i in 1:10 ) {
     fit <- glmnb.fit( mm, count, dispersion = disp, offset = log( modelFrame$sizeFactor ) )
     olddisp <- disp
     disp <- exp( optimize( function(logalpha)
        logConditionalLikelihood( exp(logalpha), mm, count, fitted.values(fit) ), 
        log( c( 1e-11, 1e5 ) ), maximum=TRUE, tol=.01 )$maximum )
     if( abs( log(disp) - log(olddisp) ) < .03 )
        break
     }
  disp
}

constructModelFrame <- function( ecs ){
  stopifnot( inherits( ecs, "JunctionSeqCountSet" ) )
  modelFrame <- cbind(
     sample = sampleNames(ecs),
     design(ecs, drop=FALSE),
     sizeFactor = sizeFactors(ecs) )
  modelFrame <- rbind(
     cbind( modelFrame, countbin="this" ),
     cbind( modelFrame, countbin="others" ) )
  rownames(modelFrame) <- NULL
  return( modelFrame )
}


arrangeCoefs <- function( frm, mf, mm = model.matrix( frm, mf ), fit = NULL, insertValues = TRUE ) {

   #if( any( attr( mm, "contrasts" ) != "contr.treatment" ) )
   #   stop( "Can only deal with standard 'treatment' contrasts." )   # Do I need this? # No.
   if( is.null(fit) & insertValues )
      stop( "If fit==NULL, returnCoefValues must be FALSE" )
   if( !is.null(fit) )
      stopifnot( all( colnames(mm) == names(coefficients(fit)) ) )

   fctTbl <- attr( terms(frm), "factors" )

   coefIndicesList <- 
   lapply( seq_len(ncol(fctTbl)), function( fctTblCol ) {
      termName <- colnames(fctTbl)[ fctTblCol ]
      varsInTerm <- stringr::str_split( termName, stringr::fixed(":") )[[1]] 
      #This check returns an error when a term appears as a interaction effect but not as a main effect (for example, for the hypothesis test model!)
      #stopifnot( all( fctTbl[ varsInTerm, fctTblCol ] == 1 ) )
      #stopifnot( sum( fctTbl[ , fctTblCol ] ) == length( varsInTerm ) )
      coefNames <- colnames(mm)[ attr( mm, "assign" ) == fctTblCol ]
      lvlTbl <- stringr::str_match( coefNames, 
         stringr::str_c( "^", stringr::str_c( varsInTerm, "([^:]*)", collapse=":" ), "$" ) )[ , -1, drop=FALSE ]
      stopifnot( ncol(lvlTbl) == length( varsInTerm ) )
      stopifnot( nrow(lvlTbl) == length( coefNames ) )
      if( !all( sapply( varsInTerm, function(v) is.factor(mf[[v]]) | is.character(mf[[v]]) ) ) )
         stop( "Non-factor in model frame" )

      varLevels <- lapply( varsInTerm, function(v) levels( factor( mf[[v]] ) ) ) 
      coefIndices <- array( NA_character_, dim = sapply( varLevels, length ), dimnames = varLevels )
      names( dimnames( coefIndices ) ) <- varsInTerm

      for( i in seq_len( nrow(lvlTbl) ) )
         coefIndices <- do.call( `[[<-`, c( quote(coefIndices), as.list( lvlTbl[ i, ] ), coefNames[i] ) )

      coefIndices
   } )
   names( coefIndicesList ) <- colnames( fctTbl )

   if( attr( terms(frm), "intercept" ) ) {
      a <- array( c( `(Intercept)` = "(Intercept)" ) )
      dimnames(a) <- list( `(Intercept)` = c( "(Intercept)" ) )
      coefIndicesList <- c( list( `(Intercept)` = a ), coefIndicesList )
   }

   if( !insertValues )
      ans <- coefIndicesList
   else
      ans <- lapply( coefIndicesList, function(coefIndices) {
         a <- ifelse( is.na(coefIndices), 0, coefficients(fit)[ coefIndices ] )
         attr( a, "variables" ) <- attr( coefIndices, "variables" )
         a } )
      
   lapply( ans, function(x) 
      if( is.array(x) ) 
         x 
      else { 
         y <- array( x, dim=length(x) )
         attr( y, "variables" ) <- attr( x, "variables" )
         dimnames(y) <- list( names(x) )
         y } )
}

balanceFeatures <- function( coefs, dispersions ) {
   stopifnot( any( sapply( coefs, function(x) 
      identical( names(dimnames(x)), "(Intercept)" ) ) ) )
   termsWithFeature <- sapply( coefs, function(x) "countbin" %in% names(dimnames(x)) )
   meanMainEffect <- sum( sapply( coefs[!termsWithFeature], mean, na.rm=TRUE ) )
   meanExonEffects <- rowSums( sapply( coefs[termsWithFeature], function(x) 
      apply2( x, "countbin", mean, na.rm=TRUE ) ) )

   meanExonFittedValue <- exp( meanMainEffect + meanExonEffects )

   exonWeights <-  1 / ( dispersions + 1 / meanExonFittedValue )

   shifts <- lapply( coefs[termsWithFeature], function(x) { 
      nonExonDims <- which(  names(dimnames(x)) != "countbin" )
      list(
         vars = names(dimnames(x))[ nonExonDims ],
         wmeans = apply2( x, nonExonDims, weighted.mean, exonWeights) ) } )

   lapply( coefs, function(x) {
      nonExonVars <- names(dimnames(x))[ names(dimnames(x)) != "countbin" ]
      if( identical( nonExonVars, "(Intercept)" ) )
         whichShift <- which( sapply( shifts, function(xx) length( xx$vars ) == 0 ) )
      else
         whichShift <- which( sapply( shifts, function(xx) identical( xx$vars, nonExonVars ) ) )
      if( length( whichShift ) == 0 )
         return( x )
      if( length( whichShift ) > 1 )
         stop( "Confused about selecting shift." )
      if( "countbin" %in% names(dimnames(x)) )
         x - shifts[[ whichShift ]]$wmeans
      else
         x + shifts[[ whichShift ]]$wmeans
    } )
}         

##################################################################

estimateSizeFactors <- function( jscs ){
   #Updated version, uses gene-level counts rather than sub-feature counts:
   #  This is preferred, since using the sub-feature counts in the normalization
   #  effectively over-weights genes with numerous annotated sub-features.
   #In practice, the results are almost always functionally-identical in most datasets.
   stopifnot( is( jscs, "JunctionSeqCountSet") );
   
   geomeans <- exp( rowMeans(log( jscs@geneCountData )) );
   sizeFactors(jscs) <- apply( jscs@geneCountData, 2, function(cnts){
     median( (cnts / geomeans)[ geomeans > 0] );
   });
   return(jscs);
   
   #Old version, uses feature counts.
   #stopifnot( is( jscs, "JunctionSeqCountSet") );
   #geomeans <- exp( rowMeans( log( counts(jscs) ) ) );
   #sizeFactors(jscs) <- apply( counts(jscs), 2, function(cnts){
   #   median( ( cnts / geomeans )[ geomeans>0 ] )
   #});
   #return(jscs);                                                      
}

fitDispersionFunction <- function( jscs , fitType = c("parametric", "local", "mean"), fitDispersionsForExonsAndJunctionsSeparately = TRUE, advancedMode = TRUE, verbose = TRUE){
   if(advancedMode){
     if(verbose) message("> fitDispersionFunction(): Advanced mode.");
     return(fitDispersionFunction_advancedMode(jscs, fitType = fitType, verbose = verbose, fitDispersionsForExonsAndJunctionsSeparately = fitDispersionsForExonsAndJunctionsSeparately));
   } else {
     #Options(warn = 1);
     if(verbose) message("> fitDispersionFunction(): Fallback mode.");
     if(verbose) message("> fitDispersionFunction() Starting (",date(),")");
     stopifnot(is(jscs, "JunctionSeqCountSet"))
     if(all(is.na(fData(jscs)$dispBeforeSharing))){
        stop("no CR dispersion estimations found, please first call estimateDispersions function")
     }
     #means <- colMeans( t(counts(jscs))/sizeFactors(jscs) )
     means <- fData(jscs)$meanBase
     disps <- fData(jscs)$dispBeforeSharing
     coefs <- c( .1, 1 )
     iter <- 0
     while(TRUE) {
        residuals <- disps / ( coefs[1] + coefs[2] / means )
        good <- which((residuals > 1e-4) & (residuals < 15))
        mm <- model.matrix(disps[good] ~ I(1/means[good]))

        #FIX ME! TEST ME! (Note: Done.)
        fit <- tryCatch({
            #testVar <- "ATEST!";
            glmgam.fit(mm, disps[good], coef.start=coefs, maxit=250);
          }, warning = function(w){
            message("> fitDispersionFunction(): warning encountered in glmgam.fit (iteration ",iter,")\n    ",w);
            #message(testVar);
            glmgam.fit(mm, disps[good], coef.start=coefs, maxit=250);
          }, error = function(e){
            message("> fitDispersionFunction(): Fatal Error encountered in glmgam.fit (iteration ",iter,")");
            message("> fitDispersionFunction(): Failed to fit the dispersion function!");
            stop(e);
          }
        );

        oldcoefs <- coefs
        coefs <- coefficients(fit)

        if(verbose) message("> (Iteration ",iter,") Dispersion Coefs: [",coefs[1], ",", coefs[2],"]");
        if(coefs[1] < 0){
           coefs[1] <- 0
           message("> fitDispersionFunction(): warning encountered on iteration ",iter,":\n    Negative intercept value in the dispersion function, it will be set to 0. Check fit diagnostics plot section from the vignette.");

           warning("Negative intercept value in the dispersion function, it will be set to 0. Check fit diagnostics plot section from the vignette.")
           break
        }
        if( sum( log( coefs / oldcoefs )^2 ) < .005 )
           break
        iter <- iter + 1
        if( iter > 25 ) {
           warning( "Dispersion fit did not converge." )
           break }
      }
      if(verbose) message("> fitDispersionFunction(): Finished on iteration ",iter,". Dispersion Coefs: [",coefs[1], ",", coefs[2],"]");
      jscs@dispFitCoefs <- coefs
      fData(jscs)$dispFitted <- jscs@dispFitCoefs[1] + jscs@dispFitCoefs[2] / colMeans( t(counts(jscs))/sizeFactors(jscs) )
      fData(jscs)$dispersion <- pmin(
         pmax(
            fData(jscs)$dispBeforeSharing,
            fData(jscs)$dispFitted,
            na.rm = TRUE ),
            1e8 )   # 1e8 as an arbitrary way-too-large value to capture infinities
      if(verbose) message("> fitDispersionFunction() Done. (",date(),")");
      return(jscs)
    }
}

fitDispersionFunction_advancedMode <- function( jscs, fitType = c("parametric", "local", "mean"), fitDispersionsForExonsAndJunctionsSeparately = TRUE, verbose = TRUE){
   #Options(warn = 1);
   if(verbose) message("> fitDispersionFunction() Starting (",date(),")");
   stopifnot(is(jscs, "JunctionSeqCountSet"))
   if(all(is.na(fData(jscs)$dispBeforeSharing))){
      stop("no CR dispersion estimations found, please first call estimateDispersions function")
   }
   #means <- colMeans( t(counts(jscs))/sizeFactors(jscs) )
   means <- fData(jscs)$meanBase
   disps <- fData(jscs)$dispBeforeSharing
   isExon <- fData(jscs)$featureType == "exonic_part";
   
   if((! any(isExon)) | (! any(! isExon))){
     fitDispersionsForExonsAndJunctionsSeparately <- FALSE;
   }
   
   if(verbose) message(">    fdf: Fitting dispersions:");
   dispFunction <- fitDispersionFunctionHelper_DESeq2(means = means,disps = disps, fitType = fitType, quiet = ! verbose);
   jscs@dispFunction <- dispFunction;
   if(fitType == "parametric") jscs@dispFitCoefs <- attr(jscs@dispFunction,"coefficients");
   
   if(fitDispersionsForExonsAndJunctionsSeparately){
     if(verbose) message(">    fdf: Fitting dispersions of exons and junctions to separate fitted trends.");
     if(verbose) message(">    fdf: Fitting exon dispersions:");
     dispFunctionExon <- fitDispersionFunctionHelper_DESeq2(means = means[ isExon],disps = disps[ isExon], fitType = fitType, quiet = ! verbose);
     if(verbose) message(">    fdf: Fitting splice-junction dispersions:");
     dispFunctionJct  <- fitDispersionFunctionHelper_DESeq2(means = means[!isExon],disps = disps[!isExon], fitType = fitType, quiet = ! verbose);
     
     jscs@dispFunctionExon <- dispFunctionExon;
     jscs@dispFunctionJct  <- dispFunctionJct;

     fData(jscs)$dispFitted <- ifelse(isExon,
                                      jscs@dispFunctionExon(means),
                                      jscs@dispFunctionJct(means)
                                      );
   } else {
     if(verbose) message(">    fdf: Fitting exons and junctions dispersions together to a single fitted trend.");
     fData(jscs)$dispFitted <- jscs@dispFunction(means);
   }
   fData(jscs)$dispersion <- pmin(
        pmax(
           fData(jscs)$dispBeforeSharing,
           fData(jscs)$dispFitted,
           na.rm = TRUE 
           ),
        1e8 # 1e8 as an arbitrary way-too-large value to capture infinities
   )   
   if(verbose) message("> fitDispersionFunction() Done. (",date(),")");
   return(jscs)
}

#This code is taken from DESeq2 (which is licensed under the Lesser-GPL v3), but adapted to work with our data structures:
fitDispersionFunctionHelper_DESeq2 <- function(means, disps, fitType = c("parametric", "local", "mean"),  minDisp = 1e-08, quiet = FALSE) {
    keep <- (! is.na(means)) & (! is.na(disps));
    means <- means[keep];
    disps <- disps[keep];
    
    #objectNZ <- object[!mcols(object)$allZero, ]
    #useForFit <- mcols(objectNZ)$dispGeneEstConv
    fitType <- match.arg(fitType);
    stopifnot(length(fitType) == 1)
    stopifnot(length(minDisp) == 1)
    
    if (fitType == "parametric") {
        trial <- try(dispFunction <- parametricDispersionFit(means, disps, quiet))
        if (inherits(trial, "try-error")) {
            warning("the parametric fit of dispersion estimates over the mean of counts\nfailed, which occurs when the trend is not well captured by the\nfunction y = a/x + b. A local regression fit is automatically performed,\nand the analysis can continue. You can specify fitType='local' or 'mean'\nto avoid this message if re-running the same data.\nWhen using local regression fit, the user should examine plotDispEsts(dds)\nto make sure the fitted line is not sharply curving up or down based on\nthe position of individual points.")
            fitType <- "local"
        }
    }
    if (fitType == "local") {
        dispFunction <- localDispersionFit(means = means, 
                                           disps = disps, minDisp = minDisp)
    }
    if (fitType == "mean") {
        useForMean <- disps > 10 * minDisp
        meanDisp <- mean(means, 
                         na.rm = TRUE, trim = 0.05)
        dispFunction <- function(means) meanDisp
    }
    if (!(fitType %in% c("parametric", "local", "mean"))) {
        stop("unknown fitType")
    }
    return(dispFunction);
    #attr(dispFunction, "fitType") <- fitType
    #dispersionFunction(object) <- dispFunction
    #dispDataFrame <- buildDataFrameWithNARows(list(dispFit = dispFit), 
    #    mcols(object)$allZero)
    #mcols(dispDataFrame) <- DataFrame(type = "intermediate", 
    #    description = "fitted values of dispersion")
    #mcols(object) <- cbind(mcols(object), dispDataFrame)
    #return(object)
}
# Estimate a parametric fit of dispersion to the mean intensity
parametricDispersionFit <- function( means, disps, quiet = FALSE ) {
   coefs <- c( .1, 1 )
   iter <- 0;
   while(TRUE) {
      residuals <- disps / ( coefs[1] + coefs[2] / means )
      good <- which( (residuals > 1e-4) & (residuals < 15) )
      # check for glm convergence below to exit while-loop
      suppressWarnings({fit <- glm( disps[good] ~ I(1/means[good]),
         family=Gamma(link="identity"), start=coefs )})
      oldcoefs <- coefs
      coefs <- coefficients(fit)
      if ( !all( coefs > 0 ) )
         stop(simpleError("parametric dispersion fit failed"))
      if ( ( sum( log( coefs / oldcoefs )^2 ) < 1e-6 )  & fit$converged )
         break
      iter <- iter + 1
      if(! quiet) message(">       (Iteration ",iter,") Parametric Dispersion Coefs: [",coefs[1], ",", coefs[2],"]");
      if ( iter > 25 ) 
        stop(simpleError("dispersion fit did not converge"))
    }
   if(! quiet)    message(">       (FINAL) Parametric Dispersion Coefs: [",coefs[1], ",", coefs[2],"]");
   names( coefs ) <- c( "asymptDisp", "extraPois" )
   ans <- function(q) coefs[1] + coefs[2] / q
   attr( ans, "coefficients" ) <- coefs
   ans
}
# Local fit of dispersion to the mean intensity
# fitting is done on log dispersion, log mean scale
localDispersionFit <- function( means, disps, minDisp ) {
  if (all(disps < minDisp*10)) {
    return(rep(minDisp,length(disps)))
  }
  d <- data.frame(logDisps = log(disps), logMeans = log(means))
  fit <- locfit(logDisps ~ logMeans, data=d[disps >= minDisp*10,,drop=FALSE],
                weights = means[disps >= minDisp*10])
  dispFunction <- function(means) exp(predict(fit, data.frame(logMeans=log(means))))
  return(dispFunction)
}


fitDispersionFunctionHelper_SIMPLE <- function(means, disps,  quiet = FALSE){
   verbose <- ! quiet;
   coefs <- c( .1, 1 );
   iter <- 0;
   while(TRUE) {
      residuals <- disps / ( coefs[1] + coefs[2] / means )
      good <- which((residuals > 1e-4) & (residuals < 15))
      mm <- model.matrix(disps[good] ~ I(1/means[good]))
      
      fit <- tryCatch({
          #testVar <- "ATEST!";
          glmgam.fit(mm, disps[good], coef.start=coefs, maxit=250);
        }, warning = function(w){
          message(">      fitDispersionFunction(): warning encountered in glmgam.fit (iteration ",iter,")\n    ",w);
          #message(testVar);
          glmgam.fit(mm, disps[good], coef.start=coefs, maxit=250);
        }, error = function(e){
          message(">     fitDispersionFunction(): Fatal Error encountered in glmgam.fit (iteration ",iter,")");
          message(">     fitDispersionFunction(): Failed to fit the dispersion function!");
          stop(e);
        }
      );
      
      oldcoefs <- coefs
      coefs <- coefficients(fit)
      
      if(verbose) message(">     (Iteration ",iter,") Dispersion Coefs: [",coefs[1], ",", coefs[2],"]");
      if(coefs[1] < 0){
        coefs[1] <- 0
        message(">     fitDispersionFunction(): warning encountered on iteration ",iter,":\n    Negative intercept value in the dispersion function, it will be set to 0. Check fit diagnostics plot section from the vignette.");
        warning("Negative intercept value in the dispersion function, it will be set to 0. Check fit diagnostics plot section from the vignette.")
        #break
      }
      if( sum( log( coefs / oldcoefs )^2 ) < .005 )
        break
      iter <- iter + 1
      if( iter > 25 ) {
        warning( "Dispersion fit did not converge." )
        break 
      }
    }
    if(verbose) message("> fitDispersionFunction(): Finished on iteration ",iter,". Dispersion Coefs: [",coefs[1], ",", coefs[2],"]");
    
    return( function(m){
      coefs[1] + coefs[2] / m;
    });
    
    #return(coefs);
}


estimatelog2FoldChanges <- function(ecs, fitExpToVar="condition", denominator="", getOnlyEffects=FALSE, averageOutExpression=TRUE, nCores=1, quiet=FALSE, file="", estimate.fc.using.genewide.model = TRUE)
{
   myApply <- getMyApply(nCores);
   
   stopifnot(is(ecs, "JunctionSeqCountSet"))
   if(any(is.na(sizeFactors(ecs)))){
      stop("Please estimate sizeFactors first\n")}
   if(!fitExpToVar %in% ecs@designColumns){
      stop("fitExpToVar parameter is not in the design columns, double check ecs@designColumns")}

   if(! estimate.fc.using.genewide.model){
     
     #OLD VERSION:
     #if( denominator == "" ){
     #  denominator <- as.character(levels(design(ecs,drop=FALSE)[[fitExpToVar]])[1]);
     #}
     #
     #varCoefCols <- which(substr(names(fData(ecs)), 1, 5 + nchar(fitExpToVar)) == paste0("Coef(",fitExpToVar)  );
     #if(length(varCoefCols) == 0){
     #  stop("No coefficients found! Cannot calculate log2foldchanges from junctionwise model! Run testJunctionsForDJU()!");
     #}
     #
     #varCoefNames <- names(fData(ecs))[varCoefCols];
     #numerator <- substr(varCoefNames, 6 + nchar(fitExpToVar), sapply(varCoefNames, nchar) - 1 );
     #
     #fc.names <- paste0("log2fold(",numerator,"/",denominator,")");
     #
     #log2fc.coef <- fData(ecs)[,varCoefCols, drop=F] / log(2);
     #names(log2fc.coef) <- fc.names;
     #
     #fData(ecs) <- cbind.data.frame(fData(ecs), log2fc.coef)
     #
     #return(ecs);
   } else {
     if(sum(is.na(featureData(ecs)$dispersion))==nrow(counts(ecs))){
        stop("No dispersion parameters found, first call function estimateDispersions...\n")
     }
     
     frm <- as.formula(paste("count ~", fitExpToVar,  "* countbin"))
     testablegenes <- as.character(unique(fData(ecs)[which(fData(ecs)$testable),]$geneID))

     geteffects <- function(geneID){
       coefficients <- fitAndArrangeCoefs(ecs, geneID=geneID, frm, balanceFeatures=TRUE)
       if( is.null( coefficients ) ){
          return(coefficients)
       }
       ret <- t(getEffectsForPlotting(coefficients, averageOutExpression=averageOutExpression, groupingVar=fitExpToVar))
       rownames(ret) <- paste(geneID, rownames(ret), sep=":")
       return(ret)
     }

     #if( nCores > 1 ){
     #   if(!is.loaded("mc_fork", PACKAGE="parallel")){
     #   stop("Please load first parallel package or set parameter nCores to 1...")}
     #   alleffects <- parallel:::mclapply( testablegenes, function(x){ geteffects(x) }, mc.cores=nCores )
     # }else{
     #   alleffects <- lapply( testablegenes, function(x){geteffects(x)})
     # }
     alleffects <- myApply( testablegenes, function(x){geteffects(x)});

      names(alleffects) <- testablegenes
      alleffects <- do.call(rbind, alleffects)
      alleffects <- vst(exp( alleffects ), ecs)
      toadd <- matrix(NA, nrow=nrow(ecs), ncol=ncol(alleffects))
      rownames(toadd) <- featureNames(ecs)

      if( getOnlyEffects ){
         colnames(toadd) <- colnames(alleffects)
         toadd[rownames(alleffects), colnames(alleffects)] <- alleffects
       }else{
         if( denominator == "" ){
            #denominator <- as.character(design(ecs, drop=FALSE)[[fitExpToVar]][1])
            denominator <- as.character(levels(design(ecs,drop=FALSE)[[fitExpToVar]])[1]);
         }
         stopifnot( any( colnames(alleffects) %in% denominator ) )
         denoCol <- which(colnames(alleffects) == denominator)
         alleffects <- log2(alleffects / alleffects[,denoCol])
         colnames(alleffects) <- sprintf("log2fold(%s/%s)", colnames(alleffects), denominator)
         colnames(toadd) <- colnames(alleffects)
         alleffects <- alleffects[,-denoCol, drop=FALSE]
         toadd <- toadd[,-denoCol, drop=FALSE]
         toadd[rownames(alleffects), colnames(alleffects)] <- alleffects
       }

      fData(ecs) <- cbind(fData(ecs), toadd)
      return(ecs);
    }
}
