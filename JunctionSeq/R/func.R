#Some of these functions were based on similar functions created for the DEXSeq package.
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




##########################################################################
##########################################################################
##########################################################################
######### Front-end Interface functions:
##########################################################################
##########################################################################
##########################################################################

INTERNALDEBUGMODE = FALSE

##########################################################################
######### Running JunctionSeq Analyses:
##########################################################################

runJunctionSeqAnalyses <- function(sample.files, sample.names, condition, 
                                flat.gff.file = NULL,
                                analysis.type = c("junctionsAndExons","junctionsOnly","exonsOnly"),
                                meanCountTestableThreshold = "auto",
                                nCores = 1,
                                use.covars = NULL, 
                                test.formula0 = formula(~ sample + countbin), 
                                test.formula1 = formula(~ sample + countbin + condition : countbin),
                                effect.formula = formula(~ condition + countbin + condition : countbin),
                                geneLevel.formula = formula(~ condition),
                                use.exons = NULL, use.junctions = NULL, 
                                use.known.junctions = TRUE, use.novel.junctions = TRUE, 
                                use.multigene.aggregates = FALSE, 
                                gene.names = NULL,
                                method.GLM = c(c("advanced","DESeq2-style"), c("simpleML","DEXSeq-v1.8.0-style")),
                                method.dispFit = c("parametric", "local", "mean"), 
                                method.dispFinal = c("shrink","max","fitted","noShare"),
                                method.sizeFactors = c("byGenes","byCountbins"),
                                method.countVectors = c("geneLevelCounts","sumOfAllBinsForGene","sumOfAllBinsOfSameTypeForGene"),
                                method.expressionEstimation = c("feature-vs-gene","feature-vs-otherFeatures"),
                                method.cooksFilter = TRUE,
                                optimizeFilteringForAlpha = 0.01,
                                fitDispersionsForExonsAndJunctionsSeparately = TRUE,
                                keep.hypothesisTest.fit = FALSE,
                                keep.estimation.fit = FALSE,
                                replicateDEXSeqBehavior.useRawBaseMean = FALSE,
                                verbose = TRUE, debug.mode = FALSE
                                ){
  keep.debug.model.data <- TRUE
  gtf.format <- TRUE
  analysis.type <- match.arg(analysis.type)
  method.GLM <- match.arg(method.GLM)
  method.dispFit <- match.arg(method.dispFit)
  method.dispFinal <- match.arg(method.dispFinal)
  method.sizeFactors <- match.arg(method.sizeFactors)
  method.countVectors <- match.arg(method.countVectors)
  method.expressionEstimation <- match.arg(method.expressionEstimation)
  
  callStack <- list(deparse(match.call(), nlines=1))
  
  if(is.null(use.junctions) && is.null(use.exons)){
    if(analysis.type == "junctionsAndExons"){
      use.junctions <- TRUE
      use.exons <- TRUE
    } else if(analysis.type == "junctionsOnly"){
      use.junctions <- TRUE
      use.exons <- FALSE
    } else if(analysis.type == "exonsOnly"){
      use.junctions <- FALSE
      use.exons <- TRUE
    }
  } else {
    if(is.null(use.junctions) || is.null(use.exons)){
      stop(paste0("Illegal syntax! If parameter use.junctions or use.exons are used, then BOTH must be set!\n use.junctions = '",use.junctions,"', use.exons = '",use.exons,"'"))
    }
    
    if(isTRUE(use.junctions) && isTRUE(use.exons)){
      analysis.type <- "junctionsAndExons"
    } else if(isTRUE(use.junctions) && isTRUE(! use.exons)){
      analysis.type <- "junctionsOnly"
    } else  if(isTRUE(! use.junctions) && isTRUE(use.exons)){
      analysis.type <- "exonsOnly"
    } else {
      stop("Illegal syntax! Parameters use.exons and use.junctions cannot both be false!")
    }
  }
  
  if(isTRUE(verbose)){
    message("> STARTING runJunctionSeqAnalyses (v",packageVersion("JunctionSeq"),") (",date(),")")
    message(paste("> rJSA: sample.files: ", paste0(sample.files,collapse=", ")))
    message(paste("> rJSA: sample.names: ",  paste0(sample.names,collapse=", ")))
    message(paste("> rJSA: condition: ",    paste0(condition,collapse=", ")))
    message(paste("> rJSA: analysis.type: ",analysis.type))
    message(paste("> rJSA: use.junctions: ",use.junctions))
    message(paste("> rJSA: use.novel.junctions: ",use.novel.junctions))
    message(paste("> rJSA: use.exons: ",use.exons))
    message(paste("> rJSA: nCores: ",nCores))
    message(paste("> rJSA: use.covars: ",use.covars))
    message(paste("> rJSA: test.formula0: ",paste0(test.formula0,collapse=" ")))
    message(paste("> rJSA: test.formula1: ",paste0(test.formula1,collapse=" ")))
    message(paste("> rJSA: use.multigene.aggregates: ", use.multigene.aggregates))
    
  }
  #require("DEXSeq")
  
  if(! is.factor(condition)){
     condition <- factor(condition, levels = sort(unique(condition)))
  }
  
  design <- data.frame(condition = condition)
  if(! is.null(use.covars)){
    message(paste0("> rJSA: using covars:"," ",date()))
    if(class(use.covars) != "data.frame"){
       stop(paste0("FATAL ERROR: use.covars must be a data.frame! Instead it appears to be: ",class(use.covars)))
    }
    for(i in 1:ncol(use.covars)){
      if(! is.factor(use.covars[[i]]) ){
        use.covars[[i]] <- factor(use.covars[[i]], levels = sort(unique(use.covars[[i]]))  )
      }
    }
    
    design <- data.frame(cbind(design,use.covars))
    for(i in 1:length(names(use.covars))){
      message(paste0("      covar: ",names(use.covars)[i]))
      message(paste0(c("      ",paste0(use.covars[,i],collapse=", "))))
    }
    names(design) <- c("condition",names(use.covars))
  }
  row.names(design) <- sample.names

    if(isTRUE(verbose)) { message(paste0("> rJSA: Reading Count files..."," ",date(),".")) }
  jscs = readJunctionSeqCounts(countfiles = as.character(sample.files),
                             samplenames = sample.names,
                             design = design,
                             flat.gff.file = flat.gff.file, 
                             verbose = verbose,
                             use.junctions = use.junctions,
                             use.novel.junctions = use.novel.junctions,
                             use.known.junctions = use.known.junctions,
                             use.exons = use.exons,
                             use.multigene.aggregates = use.multigene.aggregates,
                             nCores = nCores,
                             method.countVectors = method.countVectors,
                             test.formula1 = test.formula1,
                             gene.names = gene.names
                             )
  attr(jscs,"CallStack") <- c(list(deparse(match.call())), attr(jscs,"CallStack"))
  
  if(isTRUE(verbose)) {  message(paste0("> rJSA: Count files read."," ",date(),".")) }
  if(isTRUE(debug.mode)) reportMem(jscs)
  if(isTRUE(verbose)) { message(paste0("> rJSA: Estimating Size Factors..."," ",date(),".")) }
  
  jscs <- estimateJunctionSeqSizeFactors(jscs, method.sizeFactors = method.sizeFactors, replicateDEXSeqBehavior.useRawBaseMean = replicateDEXSeqBehavior.useRawBaseMean, verbose = verbose)
    if(isTRUE(verbose)) { message(paste0("> rJSA: Size Factors Done. Size Factors are:","."))
                  message(paste0("> rJSA: ",paste0(names(sizeFactors(jscs)),collapse=",")))
                  message(paste0("> rJSA: ",paste0(sizeFactors(jscs),collapse=",")) )
                  if(isTRUE(debug.mode)) reportMem(jscs)
                  message(paste0("> rJSA: Estimating Dispersions..."," ",date(),"."))
                }

    jscs <- estimateJunctionSeqDispersions(jscs, 
                                    method.GLM = method.GLM,
                                    nCores = nCores, 
                                    test.formula1=test.formula1, 
                                    meanCountTestableThreshold = meanCountTestableThreshold, 
                                    verbose = verbose)

    if(isTRUE(verbose)) { message(paste0("> rJSA: Dispersions estimated."," ",date(),".")) }
    if(isTRUE(debug.mode)) reportMem(jscs)
    if(isTRUE(verbose)) { message(paste0("> rJSA: Fitting Dispersion Fcn..."," ",date(),".")) }
  jscs <- fitJunctionSeqDispersionFunction(jscs, 
                                method.GLM = method.GLM,
                                method.dispFit = method.dispFit,
                                method.dispFinal = method.dispFinal,
                                verbose = verbose, 
                                fitDispersionsForExonsAndJunctionsSeparately = fitDispersionsForExonsAndJunctionsSeparately)
  
    if(isTRUE(verbose)) { message(paste0("> rJSA: Dispersions Fcn Fitted."," ",date(),".")) }
    if(isTRUE(debug.mode)) reportMem(jscs)
    if(isTRUE(verbose)) { message(paste0("> rJSA: Testing for DEU..."," ",date(),".")) }
    jscs <- testForDiffUsage(jscs, 
                  nCores = nCores, 
                  test.formula0 = test.formula0, 
                  test.formula1 = test.formula1, 
                  method.GLM = method.GLM,
                  keep.hypothesisTest.fit = keep.hypothesisTest.fit,
                  meanCountTestableThreshold = meanCountTestableThreshold,
                  method.cooksFilter = method.cooksFilter,
                  optimizeFilteringForAlpha = optimizeFilteringForAlpha,
                  verbose = verbose); 
  
  
    if(isTRUE(verbose)) { message(paste0("> rJSA: DEU tests complete."," ",date(),".")) }
    if(isTRUE(debug.mode)) reportMem(jscs)
    if(isTRUE(verbose)) { message(paste0("> rJSA: Estimating effect sizes using effects models..."," ",date(),".")) }
  jscs <- estimateEffectSizes( jscs , 
                               effect.formula = effect.formula, 
                               geneLevel.formula = geneLevel.formula, 
                               method.expressionEstimation = method.expressionEstimation,
                               keep.estimation.fit = keep.estimation.fit,
                               nCores = nCores)
    if(isTRUE(verbose)) { message(paste0("> rJSA: Effect Sizes estimated.",".")) }
    if(isTRUE(debug.mode)) reportMem(jscs)
    if(isTRUE(verbose)) { message(paste0("> rJSA: Generating results table..."," ",date(),".")) }
  
  sampleNames(jscs) = as.character(sample.names)

  if(isTRUE(debug.mode)) reportMem(jscs)
  if(isTRUE(verbose)) message("> FINISHED runJunctionSeqAnalyses (",date(),")")
  return(jscs)
}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


writeCompleteResults <- function(jscs, outfile.prefix, 
                            gzip.output = TRUE, FDR.threshold = 0.05,
                            save.allGenes = TRUE, save.sigGenes = TRUE, save.fit = FALSE, save.VST = FALSE,
                            save.bedTracks = TRUE,
                            save.jscs = FALSE,
                            bedtrack.format = c("BED","GTF","GFF3"),
                            verbose = TRUE
                            ){
   bedtrack.format <- match.arg(bedtrack.format)
   if(isTRUE(verbose)){
       message("> STARTING writeCompleteResults (v",packageVersion("JunctionSeq"),") (",date(),")")
       message(paste("> wcr: outfile.prefix: ",outfile.prefix))
       message(paste("> wcr: FDR.threshold: ",FDR.threshold))
       message(paste("> wcr: save.allGenes: ",save.allGenes))
       message(paste("> wcr: save.sigGenes: ",save.sigGenes))
       message(paste("> wcr: save.fit: ",save.fit))
       message(paste("> wcr: save.VST: ",save.VST))
       message(paste("> wcr: bedtrack.format: ",bedtrack.format))
   }
   keep.debug.model.data <- TRUE
   
   if(isTRUE(save.jscs)) save( jscs, file = paste0(outfile.prefix, "jscs.RData") )
   
   if(isTRUE(verbose)) message("> wcr: Writing sizeFactors.")
   writeSizeFactors(jscs, file = paste0(outfile.prefix,"sizeFactors.txt"))
   
   countVectors <- jscs@countVectors
   sampleNames <- sampleNames(jscs@phenoData)
   colnames(countVectors) <- c(paste0("rawCounts_",sampleNames), paste0("rawGeneCounts_",sampleNames))
   
   expression.data <- cbind.data.frame(featureID = rownames(fData(jscs)), do.call(cbind.data.frame,jscs@plottingEstimates) , countVectors )
   colnames(expression.data) <- c("featureID", unlist( lapply( jscs@plottingEstimates, colnames ) ), colnames(countVectors) )
   
   if(isTRUE(save.VST)){
     expression.data.vst <- cbind.data.frame(featureID = rownames(fData(jscs)), do.call(cbind.data.frame,lapply(jscs@plottingEstimates, function(pe){ vst(pe,jscs) }))  )
     colnames(expression.data.vst) <- c("featureID", paste0( unlist( lapply( jscs@plottingEstimates, colnames ) ), "VST" ) )
   }
   
   if(isTRUE(save.allGenes) && isTRUE(verbose)) message("> wcr: Writing results for ",length(unique(as.character(fData(jscs)$geneID)))," genes.")
   if(isTRUE(save.allGenes) && isTRUE(verbose)) message("> wcr:     Found ",nrow(fData(jscs))," counting bins belonging to these genes.")
   if(isTRUE(save.allGenes)) write.simple.table.gz(expression.data,     file=paste0(outfile.prefix,"allGenes.expression.data.txt"),      use.gzip=gzip.output,row.names=FALSE,col.names=TRUE,quote=FALSE, sep = '\t')
   if(isTRUE(save.allGenes) && isTRUE(save.VST)) write.simple.table.gz(expression.data.vst, file=paste0(outfile.prefix,"allGenes.expression.data.VST.txt"),  use.gzip=gzip.output,row.names=FALSE,col.names=TRUE,quote=FALSE, sep = '\t')
   
   if(isTRUE(save.allGenes)) write.table.gz(fData(jscs),         file=paste0(outfile.prefix,"allGenes.results.txt"), use.gzip=gzip.output,row.names="featureID",col.names=TRUE,quote=FALSE)
   
   modelFitForHypothesisTest <- jscs@modelFitForHypothesisTest
   modelFitForEffectSize <- jscs@modelFitForEffectSize
   if(isTRUE(save.fit)) save( modelFitForHypothesisTest, file = paste0(outfile.prefix, "modelFitForHypothesisTest.RData") )
   if(isTRUE(save.fit)) save( modelFitForEffectSize,     file = paste0(outfile.prefix, "modelFitForEffectSize.RData") )
   
   sig.features <- which( f.na(fData(jscs)$padjust < FDR.threshold) )
   gene.list <- unique(as.character(fData(jscs)$geneID[sig.features]))
   

   
   if(isTRUE(save.sigGenes)){
     if(length(sig.features) == 0){
       if(isTRUE(verbose)) message("> wcr: Zero Significant Features! (at adjusted-p-value threshold ",FDR.threshold,")")
     } #else {
       
       try({
         genewiseTable <- makeGeneWiseTable(jscs, gene.list, FDR.threshold = FDR.threshold, verbose = verbose)
         write.simple.table.gz(pData(genewiseTable),file=paste0(outfile.prefix,"sigGenes.genewiseResults.txt"),      use.gzip=gzip.output,row.names=FALSE,col.names=TRUE,quote=FALSE, sep = '\t')
       })
       
       if(isTRUE(verbose)) message("> wcr: Writing results for ", length(gene.list), " genes with 1 or more significant junctions (at adjusted-p-value threshold ", FDR.threshold,")")
     
       sig.rows <- which( as.character(fData(jscs)$geneID) %in% gene.list )
       if(isTRUE(verbose)) message("> wcr:     Found ", length(sig.rows), " counting bins belonging to those genes.")
       
       
                    write.simple.table.gz(expression.data[sig.rows,, drop=FALSE],     file=paste0(outfile.prefix,"sigGenes.expression.data.txt"),      use.gzip=gzip.output,row.names=FALSE,col.names=TRUE,quote=FALSE, sep = '\t')
       if(isTRUE(save.VST)) write.simple.table.gz(expression.data.vst[sig.rows,, drop=FALSE], file=paste0(outfile.prefix,"sigGenes.expression.data.VST.txt"),  use.gzip=gzip.output,row.names=FALSE,col.names=TRUE,quote=FALSE, sep = '\t')
                    write.table.gz(fData(jscs)[sig.rows,, drop=FALSE],         file=paste0(outfile.prefix,"sigGenes.results.txt"), use.gzip=gzip.output,row.names="featureID",col.names=TRUE,quote=FALSE)
     #}
   }
   
   if(isTRUE(gzip.output)){
     bedFileExt <- ".bed.gz"
   } else {
     bedFileExt <- ".bed"
   }
   
   if(isTRUE(save.bedTracks)){
     if(isTRUE(save.allGenes)){
       if(any(fData(jscs)$featureType == "splice_site" | fData(jscs)$featureType == "novel_splice_site")){
         writeExprBedTrack(paste0(outfile.prefix,"allGenes.junctionCoverage",bedFileExt), 
                         jscs = jscs, only.with.sig.gene = FALSE,
                         plot.exons = FALSE, plot.junctions = TRUE,
                         output.format = bedtrack.format, verbose = verbose, use.gzip = gzip.output,
                         trackLine = paste0("track name='JctExprAll' description='Splice Junction Coverage Estimates, by group' itemRgb='On' visibility=3"))
       }
       if(any(fData(jscs)$featureType == "exonic_part")){
         writeExprBedTrack(paste0(outfile.prefix,"allGenes.exonCoverage",bedFileExt), 
                         jscs = jscs, only.with.sig.gene = FALSE,
                         plot.exons = TRUE, plot.junctions = FALSE,
                         output.format = bedtrack.format, verbose = verbose, use.gzip = gzip.output,
                         trackLine = paste0("track name='ExonExprAll' description='Exonic Region Coverage Estimates, by group' itemRgb='On' visibility=3"))
       }
     }
     if(isTRUE(save.sigGenes)){
       sig.features <- which( f.na(fData(jscs)$padjust < FDR.threshold) )
       #if(length(sig.features) > 0){
         if(any(fData(jscs)$featureType == "splice_site" | fData(jscs)$featureType == "novel_splice_site")){
           writeExprBedTrack(paste0(outfile.prefix,"sigGenes.junctionCoverage",bedFileExt), 
                           jscs = jscs, only.with.sig.gene = TRUE,
                           plot.exons = FALSE, plot.junctions = TRUE,
                           output.format = bedtrack.format, verbose = verbose, FDR.threshold= FDR.threshold, use.gzip = gzip.output,
                           trackLine = paste0("track name='JctExprAll' description='Sig genes splice Junction Coverage Estimates, by group' itemRgb='On' visibility=3"))
         }
         if(any(fData(jscs)$featureType == "exonic_part")){
           writeExprBedTrack(paste0(outfile.prefix,"sigGenes.exonCoverage",bedFileExt), 
                           jscs = jscs, only.with.sig.gene = TRUE, 
                           plot.exons = TRUE, plot.junctions = FALSE,
                           output.format = bedtrack.format, verbose = verbose, FDR.threshold= FDR.threshold, use.gzip = gzip.output,
                           trackLine = paste0("track name='ExonExprAll' description='Sig genes exonic Region Coverage Estimates, by group' itemRgb='On' visibility=3"))
         }
         writeSigBedTrack(paste0(outfile.prefix,"sigGenes.pvalues",bedFileExt), 
                          jscs = jscs,
                          output.format = bedtrack.format, verbose = verbose, FDR.threshold= FDR.threshold, use.gzip = gzip.output,
                          trackLine = paste0("track name='JctPvals' description='Significant Splice Junctions' useScore=1 visibility=3")
                          )
       #}
     }
   }
   
   if(isTRUE(verbose)) message("> DONE writeCompleteResults (",date(),")")
   

}

writeSizeFactors <- function(jscs, file = NULL){
   if( all( is.na( sizeFactors( jscs )) ) ){
     stop("No size factors to write.\n")
   }
   sf <- data.frame(sample.ID = names(sizeFactors(jscs)), size.factor = sizeFactors(jscs))
   if(ncol(jscs@altSizeFactors) > 0){
     sf <- cbind.data.frame(sf, jscs@altSizeFactors)
   }
   if(! is.null(file)){
     write.table(sf,file = file, quote=FALSE,row.names=FALSE,col.names=TRUE,sep='\t')
   }
   return(sf)
}


##########################################################################
##########################################################################
##########################################################################
######### Top-Level Coefficient calculation functions
##########################################################################
##########################################################################
##########################################################################


generateAllExpressionEstimates.v2 <- function(jscs,nCores = 1,fitExpToVar="condition",verbose=TRUE){

   myApply <- getMyApply(nCores)

   levelCt <- length(levels(jscs@phenoData$condition))
   sampleCt <- length(sampleNames(jscs@phenoData))
   
   out <- getAllData2(jscs, runOnFeatures = seq_len(nrow(fData(jscs))),
                            fitExpToVar=fitExpToVar, 
                            formula1 = formula(paste0("~ ",fitExpToVar," + countbin + ",fitExpToVar," : countbin")),
                            nCores=nCores, dispColumn="dispersion",verbose = verbose)
   
   out <- cbind.data.frame(featureID = rownames(fData(jscs)),
                           geneID = fData(jscs)$geneID,
                           countbinID = fData(jscs)$countbinID,
                           out)
   
   out
}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


generateAllExpressionEstimates <- function(jscs,nCores = 1,fitExpToVar="condition",verbose=TRUE){
   
   myApply <- getMyApply(nCores)

   levelCt <- length(levels(jscs@phenoData$condition))
   sampleCt <- length(sampleNames(jscs@phenoData))
   
   out <- rbind( c(0.0,0.0,0.0,rep(0.0,levelCt),rep(0.0,levelCt),rep(0.0,sampleCt),rep(0.0,levelCt),rep(0.0,levelCt),rep(0.0,sampleCt)))
   #out <- rbind( c(0,0,0,rep(0,levelCt),rep(0,levelCt),rep(0,sampleCt),rep(0,levelCt),rep(0,levelCt),rep(0,sampleCt)))

   out <- data.frame(out)[0,]
   names.out <- c("featureID","geneID","countbinID",
                   paste("exprVST",(levels(jscs@phenoData$condition)),sep='_'),
                   paste("expr",(levels(jscs@phenoData$condition)),sep='_'),
                   paste("rExprVST",(levels(jscs@phenoData$condition)),sep='_'),
                   paste("rExpr",(levels(jscs@phenoData$condition)),sep='_'),
                   paste("normCountVST",sampleNames(jscs@phenoData),sep='_'),
                   paste("normCount",sampleNames(jscs@phenoData),sep='_'),
                   paste("rawCount",sampleNames(jscs@phenoData),sep='_')
                   )
   gene.list <- unique(jscs@featureData$geneID)
   blank.line.length <- length(levels(jscs@phenoData$condition)) * 4
   
   if(isTRUE(verbose)){
      message(paste("Starting expression estimate generation on",length(gene.list),"top-level features (eg genes)."))
   }
   get.one.genes.data <- function(geneIndex){
      geneID <- gene.list[geneIndex]
      es <- fitAndArrangeCoefs( jscs, geneID, frm=as.formula(paste("count ~", fitExpToVar,  "* countbin")) )
      if(! is.null(es)){
          #If the model converged properly, then output the results:
          curr <- data.frame(getAllData(jscs,es,geneID))
          if(geneIndex - 1 %% 10000 == 0 & verbose){
            message("colnames:")
            message(paste0(colnames(curr),collapse=" "))
          }
          names(curr) <- names.out
          return(  list(TRUE , curr)  )
      } else {
          #If the model does not converge, return NA for all fitted values.
          if(isTRUE(verbose)){
            message(paste("glm fit failed for gene", geneID, " index: ",geneIndex))
          }
          curr <- featureData(jscs)@data[featureData(jscs)@data$geneID == geneID,c("geneID","countbinID"), drop=FALSE]
          rn <- row.names(curr)
          curr <- cbind.data.frame(rn,curr)
          names(curr)[1] <- "featureID"
          exonCt <- dim(curr)[1]
          NA.buffer <- matrix(rep(NA,blank.line.length * exonCt), nrow=exonCt, ncol=blank.line.length )
          norCounts <- cbind.data.frame(get.norcounts.data(jscs,geneID),get.norcounts.data(jscs,geneID,vst.xform=FALSE), get.rawcounts.data(jscs,geneID))          
          curr <- cbind.data.frame(curr,NA.buffer,norCounts)
          names(curr) <- names.out
          row.names(curr) <- rn
          return( list(FALSE, curr) )
      }
   }
   
   pair.buffer <- myApply(as.list(1:length(gene.list)), FUN=get.one.genes.data)
   
   genes.data.list <- lapply(pair.buffer, "[[", 2)
   did.gene.converge <- sapply(pair.buffer, "[[", 1)
   
   if(isTRUE(verbose)){
     message(paste("generated expression data for: ",length(genes.data.list),"top-level features (ie genes). ",sum(did.gene.converge)," converged. ",sum(! did.gene.converge)," failed."))
   }

   out <- do.call(rbind.data.frame,genes.data.list)
   names(out) <- names.out
   
   #############################
   #NEW CODE: save data to jscs!
   message("Saving data to JunctionSeqCountSet. (",date(),")")
   conditionLevels <- levels(jscs@phenoData$condition)
   names(conditionLevels) <- conditionLevels
   sampleNames <- colnames(counts(jscs))
   
   exprEstimate          <- out[,  strStartsWith(names(out), "expr_") ]
   exprEstimateVST       <- out[,  strStartsWith(names(out), "exprVST_") ]
   otherExprEstimate     <- do.call(cbind.data.frame, lapply( conditionLevels, function(x){ rep(NA,nrow(jscs))} )  )
   otherExprEstimateVST  <- do.call(cbind.data.frame, lapply( conditionLevels, function(x){ rep(NA,nrow(jscs))} )  )
   relExprEstimate       <- out[,  strStartsWith(names(out), "rExpr_") ]
   relExprEstimateVST    <- out[,  strStartsWith(names(out), "rExprVST_") ]
   
   normCounts      <- cbind( out[,  strStartsWith(names(out),"normCount_") ],     do.call(cbind.data.frame, lapply( 1:sampleCt, function(x){ rep(NA,nrow(jscs))} )  )   )
   normCountsVST    <- cbind( out[,  strStartsWith(names(out),"normCountVST") ],   do.call(cbind.data.frame, lapply( 1:sampleCt, function(x){ rep(NA,nrow(jscs))} )  )   )
   
   colnames(otherExprEstimate) <- paste0("geneExpr","_",conditionLevels)
   colnames(otherExprEstimateVST) <- paste0("geneExprVST","_",conditionLevels)
   colnames(relExprEstimate) <- paste0("relExpr","_",conditionLevels)
   colnames(relExprEstimateVST) <- paste0("relExprVST","_",conditionLevels)
   
   colnames(normCounts) <- c(paste0("normCount_",sampleNames), paste0("normGeneCount_",sampleNames))
   colnames(normCountsVST) <- c(paste0("normCountsVST_",sampleNames), paste0("normGeneCountsVST_",sampleNames))
   
   jscs@plottingEstimates    <- list(exprEstimate = exprEstimate,
                                     geneExprEstimate = otherExprEstimate,
                                     relExprEstimate = relExprEstimate, 
                                     normCounts = normCounts)
   jscs@plottingEstimatesVST <- list(exprEstimateVST = exprEstimateVST, 
                                     geneExprEstimateVST = otherExprEstimateVST, 
                                     relExprEstimateVST = relExprEstimateVST, 
                                     normCountsVST = normCountsVST)
   message("Done saving data to JunctionSeqCountSet. (",date(),")")
   #############################
   
   if(isTRUE(verbose)){
     message(paste("Collapsed expression data into a single data frame. generateAllExpressionEstimates() complete."))
   }
   
   list(jscs, out)
}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


readJunctionSeqCounts <- function(countfiles = NULL, countdata = NULL,
                                  samplenames,  design,
                                  flat.gff.file=NULL, 
                                  test.formula1 = formula(~ sample + countbin + condition : countbin),
                                  analysis.type = c("junctionsAndExons","junctionsOnly","exonsOnly"),
                                  nCores = 1,
                                  use.exons = NULL, use.junctions = NULL, 
                                  use.known.junctions = TRUE,
                                  use.novel.junctions = TRUE, 
                                  use.multigene.aggregates = FALSE,
                                  gene.names = NULL,
                                  verbose = TRUE,
                                  method.countVectors = c("geneLevelCounts","sumOfAllBinsForGene","sumOfAllBinsOfSameTypeForGene"))
{
   method.countVectors <- match.arg(method.countVectors)
   
   if(isTRUE(verbose)) {
      message("-> STARTING readJunctionSeqCounts (",date(),")")
   }
   
   
   analysis.type <- match.arg(analysis.type)
   if(is.null(use.junctions) && is.null(use.exons)){
    if(analysis.type == "junctionsAndExons"){
      use.junctions <- TRUE
      use.exons <- TRUE
    } else if(analysis.type == "junctionsOnly"){
      use.junctions <- TRUE
      use.exons <- FALSE
    } else if(analysis.type == "exonsOnly"){
      use.junctions <- FALSE
      use.exons <- TRUE
    }
   } else {
    if(is.null(use.junctions) || is.null(use.exons)){
      stop(paste0("Illegal syntax! If parameter use.junctions or use.exons are used, then BOTH must be set!\n use.junctions = '",use.junctions,"', use.exons = '",use.exons,"'"))
    }
    
    if(use.junctions && use.exons){
      analysis.type <- "junctionsAndExons"
    } else if(use.junctions && (! use.exons)){
      analysis.type <- "junctionsOnly"
    } else  if((! use.junctions) && use.exons){
      analysis.type <- "exonsOnly"
    } else {
      stop("Illegal syntax! Parameters use.exons and use.junctions cannot both be false!")
    }
   }
   if(isTRUE(verbose)){
      message("---> RJSC; (v",packageVersion("JunctionSeq"),")")
      message("---> RJSC: samplenames: ",paste0(samplenames, collapse=","))
      message("---> RJSC: flat.gff.file: ",flat.gff.file)
      message("---> RJSC: use.exons:",use.exons)
      message("---> RJSC: use.junctions:",use.junctions)
      message("---> RJSC: use.novel.junctions:",use.novel.junctions)
   }
   
   
   if((is.null(countfiles) && is.null(countdata))){
     stop("Fatal error: Either countfiles OR countdata must be set! Both are null!")
   }
   if(  (!is.null(countfiles)) && (!is.null(countdata))   ){
     stop("Fatal error: Either countfiles OR countdata must be set! Both are non-null!")
   }
   
   stopifnot( class(design) == "data.frame" )
   
   for(i in 1:ncol(design)){
     if( ! is.factor(design[[i]])){
        stop("ERROR: design must be a data.frame composed entirely of factors!")
     }
   }
   
   if(! is.null(countfiles)){
     lf <- lapply( countfiles, function(x)
        read.table( x, header=FALSE,stringsAsFactors=FALSE ) )
   } else {
     lf <- countdata
   }
   
   if( !all( sapply( lf[-1], function(x) all( x$V1 == lf[1]$V1 ) ) ) )
      stop( "Count files have differing gene ID column." )
   if(isTRUE(verbose)) message("---> File read complete.");  

   dcounts <- sapply( lf, `[[`, "V2" )
   rownames(dcounts) <- lf[[1]][,1]
   dcounts <- dcounts[ substr(rownames(dcounts),1,1)!="_", ]

   bin.type <- sapply( rownames(dcounts), 
                     function(x){
                        substr(strsplit(x, ":",fixed=TRUE)[[1]][2],1,1)
                     })
   raw.geneID <- sapply( rownames(dcounts), 
                     function(x){
                        strsplit(x, ":",fixed=TRUE)[[1]][1]
                     })
   
   if(isTRUE(verbose)) message(paste0("---> Extracted counts. Found ",dim(dcounts)[1]," features so far."));  
   
   geneCountTable <- dcounts[bin.type == "A",, drop=FALSE]
   rownames(geneCountTable) <- sapply(strsplit(rownames(geneCountTable), ":"),"[[",1)
   colnames(geneCountTable) <- as.character(samplenames)
   use.bins <- bin.type != "A"
   
   if(isTRUE(verbose)) message(paste0("---> Extracted gene-level counts. Found: ",dim(geneCountTable)[1], " genes and aggregate-genes."))
   if(isTRUE(verbose)) message(paste0("---> Removed gene features. Found: ",sum(use.bins), " features to be included so far."))
   
   if(isTRUE(! use.junctions) && isTRUE(! use.exons)){
      stop("FATAL ERROR: At least one of: use.junctions or use.exons must be set to TRUE. Otherwise you've got no data to test!")
   }
   
   if(isTRUE(! use.junctions)){
      use.bins <- use.bins & bin.type != "J" & bin.type != "N"
      if(isTRUE(verbose)) message(paste0("---> Removed splice junction features. Found: ",sum(use.bins), " features to be included so far."))
   } 
   if(isTRUE(! use.novel.junctions)){
      use.bins <- use.bins & bin.type != "N"
      if(isTRUE(verbose)) message(paste0("---> Removed novel splice junction features. Found: ",sum(use.bins), " features to be included so far."))
   }
   if(isTRUE(! use.known.junctions)){
     use.bins <- use.bins & bin.type != "J"
     if(isTRUE(verbose)) message(paste0("---> Removed known splice junction features. Found: ",sum(use.bins), " features to be included so far."))
   }
   if(isTRUE(! use.exons)){
      use.bins <- use.bins & bin.type != "E"
      if(isTRUE(verbose)) message(paste0("---> Removed exon features. Found: ",sum(use.bins), " features to be included so far."))
   }
   
   is.multiGene.aggregate <- grepl("+", raw.geneID, fixed=TRUE)
   if(isTRUE(verbose)) message("---> Note: ",sum(is.multiGene.aggregate[use.bins])," counting bins are from ",length(unique(grep("+", raw.geneID[use.bins & is.multiGene.aggregate], fixed=TRUE)))," multigene aggregates (ie. overlapping genes).")
   if(isTRUE(! use.multigene.aggregates)){
      use.bins <- use.bins & (! is.multiGene.aggregate)
      if(isTRUE(verbose)) message(paste0("---> Removed multigene-aggregate features. Found: ",sum(use.bins), " features to be included so far."))
   }
   
   
   if(isTRUE(verbose)) message(paste0("---> Final feature count: ",sum(use.bins), " features to be included in the analysis."))
   dcounts <- dcounts[use.bins,, drop=FALSE]
   bin.type <- bin.type[use.bins]
   
   if(isTRUE(verbose)) message("---> Extracted feature counts.");  
   
   colnames(dcounts) <- as.character(samplenames)
   splitted <- strsplit(rownames(dcounts), ":")
   exons <- sapply(splitted, "[[", 2)
   genesrle <- sapply( splitted, "[[", 1)
   
   if(isTRUE(verbose)) message("---> counts complete.");  
   
   

   if(! is.null(flat.gff.file)){
      if(isTRUE(verbose)) message("-----> reading annotation...");  
      anno.data <- readAnnotationData(flat.gff.file)
      if(isTRUE(verbose)) message("-----> formatting annotation...");  
      #featureName     featureType     chrom   start   end     strand  gene_id part_number     transcripts
      
      exoninfo <- data.frame(chr = anno.data$chrom, start = anno.data$start, end = anno.data$end, strand = anno.data$strand)
      if(isTRUE(verbose)) message("-----> initial generation...");  
      rownames(exoninfo) <- anno.data$featureName
      
      transcripts <- anno.data$transcripts
      transcripts <- gsub("\\+", ";", transcripts)
      names(transcripts) <- rownames(exoninfo)
      
      matching <- match(rownames(dcounts), rownames(exoninfo))
      if(any(is.na(matching))){
         stop("FATAL ERROR! Annotation file appears to be missing information! Are you sure you are using the correct flattened annotation file, created by prepare_annotation_with_splices.py?")
      }
      if(isTRUE(verbose)) message("-----> creating jscs...");  
      jscs <- newJunctionSeqCountSet(countData=dcounts, geneCountData = geneCountTable, design=design, geneIDs=genesrle, countbinIDs=exons, featureIntervals=exoninfo[matching,], transcripts=transcripts[matching])
      jscs@annotationFile <- flat.gff.file
      jscs@flatGffData <- anno.data
      
      jscs@flatGffGeneData <- readGeneInfo(flat.gff.file)
      jscs <- mapGeneNames(jscs, gene.names)
   } else {
      if(isTRUE(verbose)) message("-> FINISHED readJunctionSeqCounts (",date(),")"); 
      message("Warning: flat gff annotation not set (via parameter flat.gff.file)! While technically optional, running without the annotation data may make interpretation of the data difficult. Much of the plotting functionality will not work!")
      warning("Warning: flat gff annotation not set (via parameter flat.gff.file)! While technically optional, running without the annotation data may make interpretation of the data difficult. Much of the plotting functionality will not work!")
      jscs <- newJunctionSeqCountSet(countData=dcounts, design=design, geneIDs=genesrle, countbinIDs=exons)
   }
   attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.countVectors = method.countVectors)
   attr(jscs,"CallStack") <- list(deparse(match.call()))

   jscs@analysisType <- analysis.type
   featureChar <- substr(fData(jscs)$countbinID,1,1)
   fData(jscs)[["featureType"]] <- ifelse(featureChar == "E","exonic_part",ifelse(featureChar == "J", "splice_site", "novel_splice_site"))
   varMetadata( featureData(jscs) )[ "featureType", "labelDescription" ] <- "The type of feature (exonic_part,, splice_site, or novel_splice_site)."
   pData(jscs)$countfiles <- countfiles
   if(isTRUE(verbose)) message("-----> generating count vectors... (",date(),")")
   jscs@countVectors <- getAllJunctionSeqCountVectors(jscs, nCores = nCores, method.countVectors); #use.alternate.method = use.alternate.method)
   if(isTRUE(verbose)) message("-----> count vectors generated (",date(),")");   
   
   if(isTRUE(verbose)) message("-----> generating DESeqDataSet... (",date(),")")
   jscs <- makeDESeqDataSetFromJSCS(jscs, test.formula1 = test.formula1)
   if(isTRUE(verbose)) message("-----> DESeqDataSet generated (",date(),")")
   
   fData(jscs)[["allZero"]] <- (rowSums(counts(jscs)) == 0) | 
                               (rowSums(counts(jscs@DESeqDataSet)[, colData(jscs@DESeqDataSet)$countbin == "others"]) ==0)
   mcols(jscs@DESeqDataSet)$allZero <- fData(jscs)[["allZero"]]
   
   return(jscs)
   if(isTRUE(verbose)) message("-> FINISHED readJunctionSeqCounts (",date(),")");  
}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################

mapGeneNames <- function(jscs, gene.names = NULL, gene.name.separator = "+", gene.multimap.separator = ","){
  jscs@flatGffGeneData$gene_name <- mapGeneNamesToList(jscs@flatGffGeneData$geneID, 
                                                       gene.names = gene.names,
                                                       gene.name.separator = gene.name.separator,
                                                       gene.multimap.separator = gene.multimap.separator)
  return(jscs)
}

mapGeneNamesToList <- function(geneIDs, gene.names = NULL, gene.name.separator = "+", gene.multimap.separator = ","){
  if(is.null(gene.names)){
    out = geneIDs
    names(out) = geneIDs
    return(out)
  }
  if(class(gene.names) != "data.frame"){
    stop("Error: gene.names must be a data frame!")
  }
  oldIDs <- as.character(geneIDs)
  oldIDs.map <- as.character(gene.names[,1])
  newIDs <- as.character(gene.names[,2])
  
  oldID.list <- strsplit(geneIDs, "+", fixed=TRUE)
  
  out <- sapply(oldID.list, function(gs){
     gs <- as.character(gs)
     paste0(sapply(gs, function(g){
       idx <- which(as.character(oldIDs.map) == g)
       if(length(idx) == 0){
         g
       } else if(length(idx) == 1){
         newIDs[idx]
       } else {
         paste0(newIDs[idx], collapse= gene.multimap.separator)
       }
     }), collapse= gene.name.separator)
  })
  return(out)
}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


readAnnotationData <- function(flat.gff.file){
     aggregates<-read.delim(flat.gff.file, stringsAsFactors=FALSE, header=FALSE)
     colnames(aggregates)<-c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
     aggregates<-aggregates[which(aggregates$class != "aggregate_gene"),]
     aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
     gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
     transcripts <- gsub(".*tx_set\\s(\\S+).*", "\\1", aggregates$attr)
     part_number <- gsub(".*num\\s(\\S+).*", "\\1", aggregates$attr)
     
     featureCode <- ifelse(aggregates$class == "exonic_part","E", ifelse(aggregates$class == "splice_site","J", ifelse(aggregates$class == "novel_splice_site","N","?")))
     
     featureName <- paste0(gene_id,":",featureCode,part_number)
     
     out <- data.frame(featureName = featureName,
                       featureType = aggregates$class,
                       chrom = aggregates$chr,
                       start = aggregates$start,
                       end = aggregates$end,
                       strand = aggregates$strand,
                       gene_id = gene_id,
                       part_number = part_number,
                       transcripts = transcripts
                       )
                       
     out$start <- as.integer(out$start - 1)
     
     return(out)
}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


readGeneInfo <- function(flat.gff.file){
     aggregates<-read.delim(flat.gff.file, stringsAsFactors=FALSE, header=FALSE)
     colnames(aggregates)<-c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
     aggregates <- aggregates[which(aggregates$class == "aggregate_gene"),]
     
     attrSimple <- gsub("\"|=|;", "", aggregates$attr)
     attrCells <- strsplit(attrSimple,"\\s")
     #Backwards compatibility with older versions of QoRTs:
     if( ! grepl(".*tx_set\\s\\S+.*",attrSimple[[1]])){
       gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", attrSimple)
       part_number <- gsub(".*num\\s(\\S+).*", "\\1", attrSimple)
       out <- data.frame(geneID = gene_id,num = part_number, stringsAsFactors=FALSE)
     } else {
       gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", attrSimple)
       transcripts <- gsub(".*tx_set\\s(\\S+).*", "\\1", attrSimple)
       part_number <- gsub(".*num\\s(\\S+).*", "\\1", attrSimple)
       tx_strands <- gsub(".*tx_strands\\s(\\S+).*", "\\1", attrSimple)
       aggregateGeneStrand <- gsub(".*aggregateGeneStrand\\s(\\S+).*", "\\1", attrSimple)

       out <- data.frame(geneID = gene_id,aggregateGeneStrand = aggregateGeneStrand, tx_set = transcripts, num = part_number, tx_strands = tx_strands, stringsAsFactors=FALSE)
     }
     
     return(out)

}

##############################################################################################################################################################################################################################
######### 
##############################################################################################################################################################################################################################


getColCt <- function(geneID, merged.data, 
                          plot.exon.results, plot.junction.results, plot.novel.junction.results, 
                          plot.untestable.results){
  if(is.null(plot.exon.results)){
    plot.exon.results <- any( merged.data$featureType == "exonic_part" )
  }
  if(is.null(plot.junction.results)){
    plot.junction.results <- any( merged.data$featureType == "splice_site" | merged.data$featureType == "novel_splice_site" )
  }
  if(is.null(plot.novel.junction.results)){
        if(plot.junction.results){
          plot.novel.junction.results <- any( merged.data$featureType == "novel_splice_site" )
        } else {
          plot.novel.junction.results <- FALSE
        }
  }
  
  rt <- merged.data$geneID == geneID
  if(isTRUE(! plot.exon.results)){
   rt <- rt & merged.data$featureType != "exonic_part" 
  }
  if(isTRUE(! plot.junction.results)){
   rt <- rt & merged.data$featureType != "splice_site" 
  }
  if(isTRUE(! plot.novel.junction.results)){
   rt <- rt & merged.data$featureType != "novel_splice_site" 
  }
  if(isTRUE(! plot.untestable.results)){
   rt <- rt & merged.data$testable
  }
  return(sum(rt))
}


