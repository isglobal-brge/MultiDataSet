#' @describeIn MultiDataSet Method to add a slot of methylation to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param methySet \code{MethylationSet} to be used to fill the slot.
#' @param ... Further arguments passed to \code{add_eset}.
setMethod(
  f = "add_methy",
  signature = c("MultiDataSet", "MethylationSet"),
  definition = function(object, methySet, ...) {
    
    if (!all(c("chromosome", "position") %in% fvarLabels(methySet))){
      stop("methyset must contain a fData with columns chromosome and position.")
    }
    
    object <- add_eset(object, methySet, dataset.type = "methylation", ...)
    
    return(object)
  }
)

#' @describeIn MultiDataSet Method to add a slot of methylation to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
setMethod(
  f = "add_methy",
  signature = c("MultiDataSet", "RatioSet"),
  definition = function(object, methySet, ...) {
    
    if (!all(c("chromosome", "position") %in% fvarLabels(methySet))){
      stop("methyset must contain a fData with columns chromosome and position.")
    }
    
    object <- add_eset(object, methySet, dataset.type = "methylation", ...)
    
    return(object)
  }
)
