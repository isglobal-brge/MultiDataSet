#' MethylationSet instances
#'
#' Container with the data needed to perform methylation analysis. \code{MethylationSet}
#' inherits from \code{eSet} and contains \code{meth} matrix as assay data member.
#'
#' FeatureData, which contains annotation data, is required to perform any of the
#' analysis.
#'
#' @export
#' @rdname MethylationSet-class
#' @name MethylationSet
#' @aliases MethylationSet-class MethylationSet-methods
#' @slot assayData Contains matrices with equal dimensions, and with column number
#' equal to nrow(phenoData). assayData must contain a matrix meth with rows representing
#' features (e.g., methylation probes sets) and columns representing samples.
#' @slot phenoData See \linkS4class{eSet}
#' @slot annotation See \linkS4class{eSet}
#' @slot featureData See \linkS4class{eSet}. fData should contain at least chromosome
#' and positions columns.
setClass (
    Class = "MethylationSet",
    contains = "eSet",
    prototype = prototype(new("VersionedBiobase",
        versions = c(classVersion("eSet"), MethylationSet = "1.0.0")))
)

## ------------------------------------------------------------------------- ##

#' @export
setGeneric("betas", function(object){
    standardGeneric("betas")
})

#' Transforms beta values to M-values
#'
#' Given a \code{MethylationSet} or a \code{AnalysisResults} returns
#' the matrix of M values using a logit2 transformation. Betas equal to 0 will
#' be transformed to threshold and betas equal to 1, to 1 - threshold.
#'
#' @name getMs
#' @rdname getMs-methods
#' @export
#' @param object \code{MethylationSet} or \code{AnalysisResults}
#' @param threshold Numeric with the threshold to avoid 0s and 1s.
#' @return Matrix with the M values.
#' @examples
#' if(require(MEALData)){
#'  data(mset)
#'  Ms <- getMs(mset)
#' }
setGeneric("getMs", function(object, threshold = 0.0001) standardGeneric("getMs"))