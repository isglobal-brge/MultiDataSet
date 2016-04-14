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

#' MultiDataSet instances
#'
#' The class \code{MultiDataSet} is a superior class to store multiple
#' datasets in form of triplets (assayData-phenoData-featureData). The datasets must be \code{eSet} or
#' \code{SummarizedExperiment}.
#'
#' The names of the three lists (\code{assayData}, \code{phenoData} and
#' \code{featureData}) must be the same.
#'
#' @name MultiDataSet-class
#' @rdname MultiDataSet-class
#' @exportClass MultiDataSet
#' @slot assayData List of \code{assayData} elements.
#' @slot phenoData List of \code{AnnotatedDataFrame} containing the phenoData
#' of each \code{assayData}.
#' @slot featureData List of \code{AnnotatedDataFrame} containing the featureData
#' of each \code{assayData}.
#' @slot return_method List of functions used to create the original object.
#' @seealso \code{\link{add_eset}}, \code{\link{add_rse}}
setClass(
  Class = "MultiDataSet",
  representation = representation(
    assayData = "list",
    phenoData = "list",
    featureData = "list",
    return_method = "list"
  )
)