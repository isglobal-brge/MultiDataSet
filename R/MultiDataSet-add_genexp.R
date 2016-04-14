#' @describeIn MultiDataSet Method to add a slot of expression to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param gexpSet \code{ExpressionSet} to be used to fill the slot.
setMethod(
  f = "add_genexp",
  signature = c("MultiDataSet", "ExpressionSet"),
  definition = function(object, gexpSet, ...) {
    object <- add_eset(object, gexpSet, dataset.type = "expression", ...)
    
    return(object)
  }
)

#' @describeIn MultiDataSet Method to add a slot of (RNASeq) expression to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
setMethod(
  f = "add_rnaseq",
  signature = c("MultiDataSet", "ExpressionSet"),
  definition = function(object, gexpSet, ...) {
    object <- add_eset(object, gexpSet, dataset.type = "rnaseq", ...)
    
    return(object)
  }
)

