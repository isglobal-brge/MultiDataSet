#' @describeIn MultiDataSet Method to add a \code{RangedSummarizedExperiment} to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param object \code{MultiDataSet}
#' @param dataset.type Character with the type of data of the omic set (e.g. expression, methylation...)
#' @param overwrite Logical to indicate if the set stored in the slot will be overwritten.
#' @return \code{MultiDataSet}
setMethod(
  f = "add_rse",
  signature = c("MultiDataSet", "RangedSummarizedExperiment"),
  definition = function(object, set, dataset.type, dataset.name, warnings = TRUE, overwrite = FALSE) {
    validObject(set)
    dataset.name <- paste(c(dataset.type, dataset.name), collapse = "+")

#     if (!"id" %in% colnames(colData(set))){
#       stop("colData of set must contain a column called id.")
#     }
    if(dataset.name %in% names(object)){
      if (!overwrite){
        stop("There is already an object in this slot. Set overwrite = TRUE to overwrite the previous set.")
      }
      if (warnings) {
      warning("Slot '", dataset.name, "' is already set in 'MultiDataSet'. Previous content will be overwritten.")
      }
    }
    
    env <- new("environment")
    lapply(names(SummarizedExperiment::assays(set)), function(x) assign(x, SummarizedExperiment::assays(set)[[x]], env))
    object@assayData[[dataset.name]] <- env
    
    pheno <- Biobase::AnnotatedDataFrame(as.data.frame(SummarizedExperiment::colData(set)))
    if (!"id" %in% colnames(pheno)){
      warning("No id column found in rowRanges. The id will be equal to the sampleNames")
      pheno$id <- rownames(pheno)
    }
    object@phenoData[[dataset.name]] <- pheno
    object@featureData[[dataset.name]] <- Biobase::AnnotatedDataFrame(as.data.frame(SummarizedExperiment::rowRanges(set)))
    
    returnfunc <- function(env, phe, fet) {
      assays <- SummarizedExperiment::Assays(as.list(env))
      new(class(set), assays = assays, colData = S4Vectors::DataFrame(as(phe, "data.frame")), 
          rowRanges = GenomicRanges::makeGRangesFromDataFrame(as(fet, "data.frame"), keep.extra.columns=TRUE), 
          elementMetadata = S4Vectors::DataFrame(matrix(nrow = nrow(assays), ncol = 0 )))
    }
    
    object@return_method[[dataset.name]] <- returnfunc
    return(object)
  }
)
