#These functions were extracted directly from the DEXSeq and DESeq2 source-code.
#
# We unfortunetely could not use the functions directly because they relied on the internal structure of the
#   deseq data objects, which are different in JunctionSeq.
# Additionally, many of the DESeq2 and DEXSeq commands have been changed several times in the 
#   past few releases, breaking any code that calls these functions internally. Thus, for consistency
#   they are copied over here in static form.
#
# Note that DEXSeq is licensed under the GPL v3, and DESeq2 is licensed under the LGPL v3. Therefore this
#   code packaged together is licensed under the GPL v3, as noted in the LICENSE file.
# All additions to the base DEXSeq and DESeq2 code are "united states government work" and thus cannot be
#   copyrighted. See the LICENSE file for more information.
#
#Most of these functions relate to setting up, running, and/or interpreting
#   the output from negative-binomial generalized linear models.
#

#From DEXSeq:
JS.perGeneQValue = function(pvals, wTest, geneID, method = JS.perGeneQValueExact) {
  
  ## use only those exons that were testable
  pvals     = pvals[wTest]
  ## 'factor' removes ununsed levels
  geneID    = factor(geneID[wTest])
  geneSplit = split(seq(along=geneID), geneID)

  ## summarise p-values of exons for one gene: take the minimum
  pGene = sapply(geneSplit, function(i) min(pvals[i]))
  stopifnot(all(is.finite(pGene)))

  ## Determine the thetas to be used
  theta = unique(sort(pGene))

  ## compute q-values associated with each theta
  q = method(pGene, theta, geneSplit)

  ## return a named vector of q-values per gene
  res        = rep(NA_real_, length(pGene))
  res        = q[match(pGene, theta)]
  res = pmin(1, res)
  names(res) = names(geneSplit)
  #stopifnot(!any(is.na(res)))
  return(res)
}

##----------------------------------------------------------------------
## For each value of theta, determine how many minima of random per-exon
## p-values are smaller (using the number of exons per gene)
##----------------------------------------------------------------------
JS.perGeneQValueBySimulation = function(pGene, theta, geneSplit, nperm = 24) {
  nr = sum(listLen(geneSplit))
  pRand = apply(matrix(runif(nr*nperm), nrow=nr, ncol=nperm), 2,
    function(p) sapply(geneSplit, function(i) min(p[i])))

  ## check that the apply/sapply stuff worked as intended
  stopifnot(nrow(pRand)==length(pGene), ncol(pRand)==nperm)

  hTest   = hist(pGene, breaks=c(theta,+Inf), plot=FALSE)
  hRand   = hist(pRand, breaks=c(theta,+Inf), plot=FALSE)
  stopifnot(sum(hTest$counts)==length(pGene),
            sum(hRand$counts)==length(pRand))

  numPos       = cumsum(hTest$counts)
  numFalsePos  = cumsum(hRand$counts)/nperm

  return(numFalsePos/numPos)
}

##--------------------------------------------------
## Exact computation - see methods part of the paper
##---------------------------------------------------
JS.perGeneQValueExact = function(pGene, theta, geneSplit) {
  stopifnot(length(pGene)==length(geneSplit))

  ## Compute the numerator \sum_{i=1}^M 1-(1-theta)^{n_i}
  ## Below we first identify the summands which are the same
  ## (because they have the same n_i), then do the sum via the
  ## mapply
  numExons     = listLen(geneSplit)
  tab          = tabulate(numExons)
  notZero      = (tab>0)
  numerator    = mapply(function(m, n) m * (1 - (1-theta)^n),
                            m = tab[notZero],
                            n = which(notZero))
  numerator    = rowSums(numerator)

  ## Compute the denominator: for each value of theta, the number
  ## of genes with pGene <= theta[i].
  ## Note that in cut(..., right=TRUE), the intervals are
  ## right-closed (left open) intervals.
  bins   = cut(pGene, breaks=c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)
  counts = tabulate(bins, nbins = nlevels(bins))
  denom  = cumsum(counts)
  stopifnot(denom[length(denom)]==length(pGene))

  return(numerator/denom)
}

#######################################################################################################

calc.filtered.adjusted.p <- function(
    testable, status,
    independentFiltering = TRUE,
    filter, pvalue, maxCooks,
    theta, alpha = 0.1,
    dispModelMatrix,
    cooksFilter = TRUE, cooksCutoff,
    pAdjustMethod = "BH",
    verbose = TRUE
){
   #independentFiltering = TRUE;
   if(verbose) message("> Performing final p.adjust filtering.");
   pvalue[!testable] <- NA;
   status <- as.character(status);
   
   m <- nrow(dispModelMatrix);
   p <- ncol(dispModelMatrix);
  
  if(cooksFilter){
    if(is.null(maxCooks)){
      if(verbose) message(">      No cook's cutoffs found.");
    } else {
      if (m > p) {
        defaultCutoff <- qf(.99, p, m - p)
        if (missing(cooksCutoff)) {
          cooksCutoff <- defaultCutoff
        }
        stopifnot(length(cooksCutoff)==1)
        if (is.logical(cooksCutoff) & cooksCutoff) {
          cooksCutoff <- defaultCutoff
        }
      } else {
        cooksCutoff <- FALSE
      }
      if(verbose) message(">      Applying cooks cutoff outlier filtering.");
      # apply cutoff based on maximum Cook's distance
      performCooksCutoff <- (is.numeric(cooksCutoff) | cooksCutoff) 
      if ((m > p) & performCooksCutoff) {
        cooksOutlier <- f.na(maxCooks > cooksCutoff);

        if(verbose) message(">      Filtering out ",sum(cooksOutlier & testable)," features because maxCooks > ",cooksCutoff);
        pvalue[cooksOutlier] <- NA
        status[testable & cooksOutlier] <- "COOKS_OUTLIER";
        testable[testable & cooksOutlier] <- FALSE;
      }
    }
  }
  
  # perform independent filtering
  if(independentFiltering) {
    if(verbose) message(">      automatically selecting a filtering threshold to optimize results at the alpha < ",alpha," significance level.");
    
    if (missing(filter)) {
      stop("Must set a filter parameter!");
    }
    if (missing(theta)) {
      lowerQuantile <- mean(filter == 0)
      if (lowerQuantile < .95) upperQuantile <- .95 else upperQuantile <- 1
      theta <- seq(lowerQuantile, upperQuantile, length=20)
    }
    stopifnot(length(theta) > 1)
    #stopifnot(length(filter) == nrow(object))
    filtPadj <- genefilter::filtered_p(filter=filter, test=pvalue,
                                       theta=theta, method=pAdjustMethod) 
    numRej  <- colSums(filtPadj < alpha, na.rm = TRUE)
    j <- which.max(numRej)
    padj <- filtPadj[, j, drop=TRUE]
    cutoffs <- quantile(filter, theta)
    filterThreshold <- cutoffs[j];
    use <- filter >= filterThreshold;
    filterNumRej <- data.frame(theta=theta, numRej=numRej)
    
    
    if(verbose) message(">      Automatically selecting a filtering threshold of ",filterThreshold," to optimize results at the alpha < ",alpha," significance level.");
    if(verbose) message(">         (Automatic independent filtering: ",sum(testable & (! use))," out of ",sum(testable)," features filtered out, using baseMean < ",filterThreshold,")");
    if(verbose) message(">         (Rejected null hypothesis for ",sum(numRej[j]), " features at alpha < ",alpha,")");
    status[testable & (! use)] <- "LOW_COUNTS_INDEP_FILTER";
    testable[testable & (! use)] <- FALSE;
    
  } else {
    # regular p-value adjustment
    # which does not include those rows which were removed
    # by maximum Cook's distance
    padj <- p.adjust(pvalue,method=pAdjustMethod)
  }
  
  if(verbose) message("> Final p.adjust filtering complete.");
  
  return(list(
    res = data.frame(
      status = status,
      testable = testable,
      padjust = padj
    ),
    filterThreshold = filterThreshold,
    filterNumRej = filterNumRej
  ));
  
}





#FROM DESeq2:
#' Low-level function to estimate size factors with robust regression.
#' 
#' Given a matrix or data frame of count data, this function estimates the size
#' factors as follows: Each column is divided by the geometric means of the
#' rows. The median (or, if requested, another location estimator) of these
#' ratios (skipping the genes with a geometric mean of zero) is used as the size
#' factor for this column. Typically, one will not call this function directly, but use
#' \code{\link{estimateSizeFactors}}.
#' 
#' @param counts a matrix or data frame of counts, i.e., non-negative integer
#' values
#' @param locfunc a function to compute a location for a sample. By default, the
#' median is used. However, especially for low counts, the
#' \code{\link[genefilter]{shorth}} function from genefilter may give better results.
#' @param geoMeans by default this is not provided, and the
#' geometric means of the counts are calculated within the function.
#' A vector of geometric means from another count matrix can be provided
#' for a "frozen" size factor calculation
#' @param controlGenes optional, numeric or logical index vector specifying those genes to
#' use for size factor estimation (e.g. housekeeping or spike-in genes)
#' @return a vector with the estimates size factors, one element per column
#' @author Simon Anders
#' @seealso \code{\link{estimateSizeFactors}}
#' @examples
#' 
#' dds <- makeExampleDESeqDataSet()
#' estimateSizeFactorsForMatrix(counts(dds))
#' geoMeans <- exp(rowMeans(log(counts(dds))))
#' estimateSizeFactorsForMatrix(counts(dds),geoMeans=geoMeans)
#'
estimateSizeFactorsForMatrix <- function( counts, locfunc = stats::median, geoMeans, controlGenes )
{
  if (missing(geoMeans)) {
    loggeomeans <- rowMeans(log(counts))
  } else {
    if (length(geoMeans) != nrow(counts)) {
      stop("geoMeans should be as long as the number of rows of counts")
    }
    loggeomeans <- log(geoMeans)
  }
  if (all(is.infinite(loggeomeans))) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }
  sf <- if (missing(controlGenes)) {
    apply(counts, 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
    })
  } else {
    if ( !( is.numeric(controlGenes) | is.logical(controlGenes) ) ) {
      stop("controlGenes should be either a numeric or logical vector")
    }
    loggeomeansSub <- loggeomeans[controlGenes]
    apply(counts[controlGenes,,drop=FALSE], 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeansSub)[is.finite(loggeomeansSub) & cnts > 0]))
    })
  }
  sf
}



estimateUnsharedDispersions <- function(object,
                                        maxit=100, quiet=FALSE, formula=design(object), 
                                        BPPARAM=MulticoreParam(workers=1)) {
  
  allVars <- all.vars(formula)
  if( any(!allVars %in% colnames( colData(object) )) ){
     notPresent <- allVars[!allVars %in% colnames( colData(object) )]
     notPresent <- paste(notPresent, collapse=",")
     stop(sprintf("the variables '%s' of the parameter 'formula' are not specified in the columns of the colData", notPresent ) )
  }
  splitParts <- sort(
    rep(seq_len(BPPARAM$workers), 
    length.out=nrow(object) ) )
  splitObject <- split( object, splitParts )
  modelMatrix <- rmDepCols(
    model.matrix(formula, as.data.frame(colData(object))))
  splitObject <- bplapply( splitObject, 
      function(x){
        estimateDispersionsGeneEst(x, 
          maxit=maxit, quiet=quiet, 
          modelMatrix = modelMatrix, 
          niter = 10)}, 
    BPPARAM=BPPARAM )
  
  mergeObject <- do.call( rbind, splitObject )
  matchedNames <- match( rownames(object), rownames(mergeObject))  
  mcols(object) <- mcols( mergeObject )[matchedNames,]
  assays(object) <- assays(mergeObject[matchedNames,])
  
  return(object);
}

adapted.estimateDispersionsMAP <- function( jscs, useRows, #dispFn, 
                                    test.formula1 = formula(jscs@formulas[["formulaDispersion"]]),
                                    outlierSD = 2, dispPriorVar, minDisp = 1e-08, 
                                    kappa_0 = 1, dispTol = 1e-06, maxit = 100, modelMatrix, verbose = TRUE) {
    stopifnot(length(outlierSD) == 1)
    stopifnot(length(minDisp) == 1)
    stopifnot(length(kappa_0) == 1)
    stopifnot(length(dispTol) == 1)
    stopifnot(length(maxit) == 1)
    
    dispBeforeSharing <- fData(jscs)$dispBeforeSharing[useRows];
    dispFitted        <- fData(jscs)$dispFitted[useRows];
    #dispFn <- jscs@dispFunction;
    means <- fData(jscs)$baseMean[useRows];
    countVectors <- jscs@countVectors[useRows,];
    mu <- jscs@fittedMu[useRows,];
    
    #if (!is.null(mcols(object)$dispersion)) {
    #    if (!quiet) 
    #        message("found already estimated dispersions, removing these")
    #    removeCols <- c("dispersion", "dispOutlier", "dispMAP", 
    #        "dispIter", "dispConv")
    #    mcols(object) <- mcols(object)[, !names(mcols(object)) %in% 
    #        removeCols, drop = FALSE]
    #}
    #if (missing(modelMatrix)) {
    #    modelMatrix <- model.matrix(design(object), data = as.data.frame(colData(object)))
    #}
    if (missing(modelMatrix)) {
        #modelMatrix <- model.matrix(test.formula1, data = as.data.frame(pData(jscs)))
        modelFrame <- constructModelFrame( jscs )
        modelMatrix <- rmDepCols( model.matrix( test.formula1, modelFrame ) )
    }
    aboveMinDisp <- f.na(dispBeforeSharing > minDisp * 100);
    varLogDispEsts <- mad(log(dispBeforeSharing[aboveMinDisp]) - log(dispFitted[aboveMinDisp]), na.rm = TRUE)^2;
    #else {
    #    message("using supplied model matrix")
    #}
    if (missing(dispPriorVar)) {
    #   if (sum(mcols(object)$dispGeneEst >= minDisp * 100, na.rm = TRUE) == 0) {
        if (sum(dispBeforeSharing >= minDisp * 100, na.rm = TRUE) == 0) {
            stop(paste0("all genes have dispersion estimates < ", minDisp * 10, "!"));
            #warning(paste0("all genes have dispersion estimates < ", minDisp * 10, ", returning disp = ", minDisp *  10))
            #resultsList <- list(dispersion = minDisp * 10)
            #dispDataFrame <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
            #mcols(dispDataFrame) <- DataFrame(type = "intermediate",  description = "final estimates of dispersion")
            #mcols(object) <- cbind(mcols(object), dispDataFrame)
            #dispFn <- dispersionFunction(object)
            #attr(dispFn, "dispPriorVar") <- 0.25
            #dispersionFunction(object, estimateVar = FALSE) <- dispFn
            #return(object)
        }
        dispPriorVar <- adapted.estimateDispersionsPriorVar(dispBeforeSharing = dispBeforeSharing, dispFitted = dispFitted, varLogDispEsts = varLogDispEsts, modelMatrix = modelMatrix, minDisp = minDisp, verbose = verbose)
        #dispFn <- dispersionFunction(object)
        #attr(dispFn, "dispPriorVar") <- dispPriorVar
        #dispersionFunction(object, estimateVar = FALSE) <- dispFn
        
    } #else {
    #    dispFn <- dispersionFunction(object)
    #    attr(dispFn, "dispPriorVar") <- dispPriorVar
    #    dispersionFunction(object, estimateVar = FALSE) <- dispFn
    #}
    stopifnot(length(dispPriorVar) == 1)
    #objectNZ <- object[!mcols(object)$allZero, , drop = FALSE]
    #varLogDispEsts <- attr(dispFn, "varLogDispEsts")
    
    
    log_alpha_prior_sigmasq <- dispPriorVar
    
    #mu <- 
    dispInit <- ifelse(dispBeforeSharing > 0.1 * dispFitted, 
                       dispBeforeSharing, 
                       dispFitted)
    dispInit[is.na(dispInit)] <- dispFitted[is.na(dispInit)]
    dispResMAP <- fitDispWrapper(ySEXP = countVectors, xSEXP = modelMatrix, 
        mu_hatSEXP = mu, log_alphaSEXP = log(dispInit), log_alpha_prior_meanSEXP = log(dispFitted), 
        log_alpha_prior_sigmasqSEXP = log_alpha_prior_sigmasq, 
        min_log_alphaSEXP = log(minDisp/10), kappa_0SEXP = kappa_0, 
        tolSEXP = dispTol, maxitSEXP = maxit, use_priorSEXP = TRUE)
    dispMAP <- exp(dispResMAP$log_alpha)
    dispConv <- dispResMAP$iter < maxit
    refitDisp <- ! dispConv
    if (sum(refitDisp) > 0) {
        dispInR <- fitDispGridWrapper(y = countVectors[refitDisp, , drop = FALSE], 
                                      x = modelMatrix, mu = mu[refitDisp, , drop = FALSE], 
                                      logAlphaPriorMean = log(dispFitted)[refitDisp], 
                                      logAlphaPriorSigmaSq = log_alpha_prior_sigmasq, usePrior = TRUE)
        dispMAP[refitDisp] <- dispInR
    }
    maxDisp <- max(10, nrow(modelMatrix))
    dispMAP <- pmin(pmax(dispMAP, minDisp), maxDisp)
    dispersionFinal <- dispMAP
    dispOutlier <- f.na(log(dispBeforeSharing) > log(dispFitted) + outlierSD * sqrt(varLogDispEsts))
    #dispOutlier[is.na(dispOutlier)] <- FALSE
    dispersionFinal[dispOutlier] <- dispBeforeSharing[dispOutlier]
    resultsList <- list(dispersion = dispersionFinal, 
                        dispIter = dispResMAP$iter, 
                        dispOutlier = dispOutlier, 
                        dispMAP = dispMAP)
    #dispDataFrame <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
    #mcols(dispDataFrame) <- DataFrame(type = rep("intermediate", 
    #    ncol(dispDataFrame)), description = c("final estimate of dispersion", 
    #    "number of iterations", "dispersion flagged as outlier", 
    #    "maximum a posteriori estimate"))
    #mcols(object) <- cbind(mcols(object), dispDataFrame)
    #return(object)
    return(resultsList);
}



adapted.estimateDispersionsPriorVar <- function (dispBeforeSharing, dispFitted, varLogDispEsts, minDisp = 1e-08, modelMatrix, verbose = TRUE){
    #objectNZ <- object[!mcols(object)$allZero, , drop = FALSE]
    
    aboveMinDisp <- dispBeforeSharing >= minDisp * 100
    #if (missing(modelMatrix)) {
    #    modelMatrix <- model.matrix(design(object), data = as.data.frame(colData(object)))
    #}
    dispResiduals <- log(dispBeforeSharing) - log(dispFitted);
    if (sum(aboveMinDisp, na.rm = TRUE) == 0) {
        stop("no data found which is greater than minDisp")
    }
    #varLogDispEsts <- attr(dispersionFunction(object), "varLogDispEsts")
    m <- nrow(modelMatrix)
    p <- ncol(modelMatrix)
    if(((m - p) <= 3) & (m > p)) {
        
        if (exists(".Random.seed")) {
            oldRandomSeed <- .Random.seed
        }
        set.seed(2)
        obsDist <- dispResiduals[aboveMinDisp]
        brks <- -20:20/2
        obsDist <- obsDist[obsDist > min(brks) & obsDist < max(brks)]
        obsVarGrid <- seq(from = 0, to = 8, length = 200)
        obsDistHist <- hist(obsDist, breaks = brks, plot = FALSE)
        klDivs <- sapply(obsVarGrid, function(x) {
            randDist <- log(rchisq(10000, df = (m - p))) + rnorm(10000, 0, sqrt(x)) - log(m - p)
            randDist <- randDist[randDist > min(brks) & randDist < max(brks)]
            randDistHist <- hist(randDist, breaks = brks, plot = FALSE)
            z <- c(obsDistHist$density, randDistHist$density)
            small <- min(z[z > 0])
            kl <- sum(obsDistHist$density * (log(obsDistHist$density + small) - log(randDistHist$density + small)))
            kl
        })
        lofit <- loess(klDivs ~ obsVarGrid, span = 0.2)
        obsVarFineGrid <- seq(from = 0, to = 8, length = 1000)
        lofitFitted <- predict(lofit, obsVarFineGrid)
        argminKL <- obsVarFineGrid[which.min(lofitFitted)]
        expVarLogDisp <- trigamma((m - p)/2)
        dispPriorVar <- pmax(argminKL, 0.25)
        if (exists("oldRandomSeed")) {
            .Random.seed <<- oldRandomSeed
        }
        return(dispPriorVar)
    }
    if (m > p) {
        expVarLogDisp <- trigamma((m - p)/2)
        dispPriorVar <- pmax((varLogDispEsts - expVarLogDisp), 0.25)
    } else {
        dispPriorVar <- varLogDispEsts
        expVarLogDisp <- 0
    }
    return(dispPriorVar);
}

##############################


vst <- function(x, ecs)
{
  if(is.null(ecs@dispFunctionType[["fitType"]])){
    stop("vst: Cannot find dispersion fit data! Fit dispersions first.");
  } else if(ecs@dispFunctionType[["fitType"]] != "parametric"){
    warning("vst cannot be performed without parametric dispersion fit. Using un-transformed data instead.");
    return(x);
  } else {
    suppressWarnings(
      ( 2 / ( sqrt(ecs@dispFitCoefs[1]) ) ) * 
         log( 2 * ecs@dispFitCoefs[1] * sqrt(x) + 
              2 * sqrt( ecs@dispFitCoefs[1] * ( ecs@dispFitCoefs[2] + 1 + ecs@dispFitCoefs[1] * x ) ) ) -
      ( 2 / ( sqrt(ecs@dispFitCoefs[1]) ) ) * 
         log( 2 * sqrt( ecs@dispFitCoefs[1] * ( ecs@dispFitCoefs[2] + 1 ) ) ) 
    )
  }
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
  #modelFrame <- rbind(modelFrame, modelFrame);
  #modelFrame$countbin <- 
  modelFrame <- rbind(
     cbind( modelFrame, countbin="this" ),
     cbind( modelFrame, countbin="others" ) )
  modelFrame$countbin <- factor(modelFrame$countbin, levels = c("this","others"));
  rownames(modelFrame) <- NULL
  return( modelFrame )
}

arrangeCoefs <- function( frm, mf, mm = model.matrix( frm, mf ), fit = NULL, insertValues = TRUE ) {

   #if( any( attr( mm, "contrasts" ) != "contr.treatment" ) )
   #   stop( "Can only deal with standard 'treatment' contrasts." )   # Do I need this? # No.
   if( is.null(fit) & insertValues )
      stop( "If fit==NULL, returnCoefValues must be FALSE" )
   if( !is.null(fit) )
      stopifnot( all( colnames(mm) == names(coefficients(fit)) ) );

   fctTbl <- attr( terms(frm), "factors" );

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

fitDispersionFunction_simpleMode <- function(jscs, verbose = TRUE){
     #Options(warn = 1);
     if(verbose) message("> fitDispersionFunction(): Fallback mode.");
     if(verbose) message("> fitDispersionFunction() Starting (",date(),")");
     stopifnot(is(jscs, "JunctionSeqCountSet"))
     if(all(is.na(fData(jscs)$dispBeforeSharing))){
        stop("no CR dispersion estimations found, please first call estimateDispersions function")
     }
     #means <- colMeans( t(counts(jscs))/sizeFactors(jscs) )
     means <- fData(jscs)$baseMean
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
            
      jscs@dispFunctionType <- list(fitType = "parametric", 
                                    attemptedFitType = "parametric", 
                                    finalDispersionMethod = "max", 
                                    fitDispersionsForExonsAndJunctionsSeparately = FALSE);
                                 
      if(verbose) message("> fitDispersionFunction() Done. (",date(),")");
      return(jscs);
}

fitDispersionFunction_advancedMode <- function( jscs, 
                                                fitType = c("parametric", "local", "mean"), 
                                                finalDispersionMethod = c("shrink","max","fitted","noShare"),
                                                fitDispersionsForExonsAndJunctionsSeparately = TRUE, 
                                                nCores = 1,
                                                verbose = TRUE){
   #Options(warn = 1);
   minDisp <- 1e-8; # 1e8 as an arbitrary way-too-small value to capture infinities
   minFitDisp <- minDisp * 100;
   
   if(verbose) message("> fitDispersionFunction() Starting (",date(),")");
   fitType <- match.arg(fitType);
   finalDispersionMethod <- match.arg(finalDispersionMethod);
   if(verbose) message(">   (fitType = ",fitType,")");
   if(verbose) message(">   (finalDispersionMethod = ",finalDispersionMethod,")");
   if(verbose) message(">   (fitDispersionsForExonsAndJunctionsSeparately = ",fitDispersionsForExonsAndJunctionsSeparately,")");
   
   stopifnot(is(jscs, "JunctionSeqCountSet"))
   if(all(is.na(fData(jscs)$dispBeforeSharing))){
      stop("no CR dispersion estimations found, please first call estimateDispersions function")
   }
   #means <- colMeans( t(counts(jscs))/sizeFactors(jscs) )
   means <- fData(jscs)$baseMean;
   disps <- fData(jscs)$dispBeforeSharing;
   useForFit <- f.na(disps > minFitDisp) & (! fData(jscs)$allZero);
   isExon <- fData(jscs)$featureType == "exonic_part";
   
   if((! any(isExon)) | (! any(! isExon))){
     fitDispersionsForExonsAndJunctionsSeparately <- FALSE;
   }
   if(sum(useForFit) == 0){
     stop("No loci found with dispersion > ",minFitDisp, ". dispersion fit failed!");
   }
   if(sum(useForFit & isExon) == 0 & fitDispersionsForExonsAndJunctionsSeparately){
     stop("No exonic regions found with dispersion > ",minFitDisp, ". dispersion fit failed!");
   }
   if(sum(useForFit & (! isExon)) == 0 & fitDispersionsForExonsAndJunctionsSeparately){
     stop("No splice junctions found with dispersion > ",minFitDisp, ". dispersion fit failed!");
   }
   
   if(verbose) message(">    fdf: Fitting dispersions:");
   dispFunction <- adapted.estimateDispersionsFit(means = means[useForFit], disps = disps[useForFit], fitType = fitType, quiet = ! verbose);
   #dispFunction <- dispersionFunction( estimateDispersionsFit(jscs@DESeqDataSet,fitType=fitType, quiet=!verbose) );
   jscs@dispFunction <- dispFunction;
   if(attr(dispFunction,"fitType") == "parametric") jscs@dispFitCoefs <- attr(jscs@dispFunction,"coefficients");
   
   if(fitDispersionsForExonsAndJunctionsSeparately){
     if(verbose) message(">    fdf: Fitting dispersions of exons and junctions to separate fitted trends.");
     if(verbose) message(">    fdf: Fitting exon dispersions:");
     useCols <- useForFit & isExon;
     dispFunctionExon <- adapted.estimateDispersionsFit(means = means[useCols],disps = disps[useCols], fitType = fitType, quiet = ! verbose);
     #dispFunctionExon <- dispersionFunction( estimateDispersionsFit(jscs@DESeqDataSet[isExon],fitType=fitType, quiet=!verbose) );
     
     if(verbose) message(">    fdf: Fitting splice-junction dispersions:");
     useCols <- useForFit & (! isExon);
     dispFunctionJct  <- adapted.estimateDispersionsFit(means = means[useCols],disps = disps[useCols], fitType = fitType, quiet = ! verbose);
     #dispFunctionJct <- dispersionFunction( estimateDispersionsFit(jscs@DESeqDataSet[! isExon],fitType=fitType, quiet=!verbose) );

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
   
   if(finalDispersionMethod == "max"){
      if(verbose) message("> fdf(): Using the higher of the fitted or feature-specific dispersion estimates.");
      fData(jscs)$dispersion <- pmin(
           pmax(
              fData(jscs)$dispBeforeSharing,
              fData(jscs)$dispFitted,
              na.rm = TRUE 
              ),1e8)
   } else if(finalDispersionMethod == "shrink"){
      if(verbose) message("> fdf(): 'Shrinking' fitted and feature-specific dispersion estimates.");
      fData(jscs)$dispersion <- NA;
      if(fitDispersionsForExonsAndJunctionsSeparately){
        resultsListExon <- adapted.estimateDispersionsMAP( jscs = jscs, 
                                               useRows = (isExon) & fData(jscs)$testable,
                                               verbose = verbose)
        resultsListJct <- adapted.estimateDispersionsMAP( jscs = jscs, 
                                               useRows = (! isExon) & fData(jscs)$testable,
                                               verbose = verbose)
        fData(jscs)$dispersion[(isExon) & fData(jscs)$testable]   <- resultsListExon[["dispersion"]];
        fData(jscs)$dispersion[(! isExon) & fData(jscs)$testable] <- resultsListJct[["dispersion"]];
      } else {
        resultsList <- adapted.estimateDispersionsMAP( jscs = jscs, 
                                               useRows = fData(jscs)$testable,
                                               verbose = verbose)
        fData(jscs)$dispersion[fData(jscs)$testable] <- resultsList[["dispersion"]];
      }
   } else if(finalDispersionMethod == "fitted"){
      if(verbose) message("> fdf(): Using fitted dispersion estimates. (NOTE: NOT RECOMMENDED!)");
      fData(jscs)$dispersion <- pmin(fData(jscs)$dispFitted,1e8)
      #fData(jscs)$dispersion[is.na(fData(jscs)$dispersion)] <- 1e8;
   } else if(finalDispersionMethod == "noShare"){
      if(verbose) message("> fdf(): Using feature-specific dispersion estimates. (NOTE: NOT RECOMMENDED!)");
      fData(jscs)$dispersion <- pmin(fData(jscs)$dispBeforeSharing,1e8);
      #fData(jscs)$dispersion <- pmin(fData(jscs)$dispFitted,1e8);
      #fData(jscs)$dispersion[is.na(fData(jscs)$dispersion)] <- 1e8;
   }

   jscs@dispFunctionType <- list(fitType = attr(dispFunction,"fitType"), 
                                 attemptedFitType = fitType,
                                 finalDispersionMethod = finalDispersionMethod,
                                 fitDispersionsForExonsAndJunctionsSeparately = fitDispersionsForExonsAndJunctionsSeparately);
   
   failedDisp <- is.na(fData(jscs)$dispersion) & fData(jscs)$testable;
   
   if(verbose) message("> fdf() Dispersion estimate failed for ",sum(failedDisp)," out of ",sum(fData(jscs)$testable)," features.");
   
   
   if(verbose) message("> fitDispersionFunction() Done. (",date(),")");
   return(jscs)
}


#This code is taken from DESeq2 (which is licensed under the Lesser-GPL v3), but adapted to work with our data structures:
adapted.estimateDispersionsFit <- function(means, disps, fitType = c("parametric", "local", "mean"),  minDisp = 1e-08, quiet = FALSE) {
    #keep <- (! is.na(means)) & (! is.na(disps));
    #keep <- keep & f.na(means > 0);
    #keep <- keep & (disps > 100 * minDisp);
    #means <- means[keep];
    #disps <- disps[keep];
    
    #objectNZ <- object[!mcols(object)$allZero, ]
    #useForFit <- mcols(objectNZ)$dispGeneEstConv
    fitType <- match.arg(fitType);
    stopifnot(length(fitType) == 1)
    stopifnot(length(minDisp) == 1)
    
    if (fitType == "parametric") {
        trial <- try(dispFunction <- adapted.parametricDispersionFit(means, disps, quiet))
        if (inherits(trial, "try-error")) {
            warning("the parametric fit of dispersion estimates over the mean of counts\nfailed, which occurs when the trend is not well captured by the\nfunction y = a/x + b. A local regression fit is automatically performed,\nand the analysis can continue. You can specify fitType='local' or 'mean'\nto avoid this message if re-running the same data.\nWhen using local regression fit, the user should examine plotDispEsts(dds)\nto make sure the fitted line is not sharply curving up or down based on\nthe position of individual points.")
            fitType <- "local"
        }
    }
    if (fitType == "local") {
        dispFunction <- adapted.localDispersionFit(means = means, 
                                           disps = disps, minDisp = minDisp)
    }
    if (fitType == "mean") {
        #useForMean <- disps > 10 * minDisp
        meanDisp <- mean(means, 
                         na.rm = TRUE, trim = 0.05)
        dispFunction <- function(means) meanDisp
    }
    if (!(fitType %in% c("parametric", "local", "mean"))) {
        stop("unknown fitType")
    }
    attr(dispFunction, "fitType") <- fitType;
    
    varLogDispEsts <- mad(log(disps) - log(dispFunction(means)), na.rm = TRUE)^2;
    attr( dispFunction, "varLogDispEsts" ) <- varLogDispEsts;
    
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
adapted.parametricDispersionFit <- function( means, disps, quiet = FALSE ) {
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
adapted.localDispersionFit <- function( means, disps, minDisp ) {
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
