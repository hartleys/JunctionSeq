#These functions were extracted directly from the DEXSeq and DESeq2 source-code.
#
# We unfortunetely could not use the functions directly because they relied on the internal structure of the
#   deseq data objects, which are different in JunctionSeq.
# Additionally, many of the DESeq2 and DEXSeq internal data structures have changed several times in the 
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

#' @rdname dispersionFunction
#' @export
setGeneric("dispersionFunction", function(object,...) standardGeneric("dispersionFunction"))

#' @rdname dispersionFunction
#' @export
setGeneric("dispersionFunction<-", function(object,...,value) standardGeneric("dispersionFunction<-"))

#' @rdname dispersions
#' @export
setGeneric("dispersions", function(object,...) standardGeneric("dispersions"))

#' @rdname dispersions
#' @export
setGeneric("dispersions<-", function(object,...,value) standardGeneric("dispersions<-"))

#' @rdname normalizationFactors
#' @export
setGeneric("normalizationFactors", function(object,...) standardGeneric("normalizationFactors"))

#' @rdname normalizationFactors
#' @export
setGeneric("normalizationFactors<-", function(object,...,value) standardGeneric("normalizationFactors<-"))
